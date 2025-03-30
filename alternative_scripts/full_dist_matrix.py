import numpy as np

def full_dist_matrix(atom1, box_dim, atom2=None):
    """
    Computes the NxN distance matrix between two sets of atoms (or within one set),
    applying periodic boundary conditions (PBC) for orthogonal or triclinic boxes.

    Parameters
    ----------
    atom1 : list of dict
        Each dict must have 'x', 'y', 'z' coordinates.
    box_dim : list or array
        Box dimensions. Can be:
          - length 3: [lx, ly, lz] (orthogonal box)
          - length 6: [xlo, ylo, zlo, xhi, yhi, zhi] (commonly used in some data formats)
          - length 9: [lx, ly, lz, 0, 0, xy, 0, xz, yz] (triclinic tilt factors)
    atom2 : list of dict or None, optional
        If provided, computes distances between atom1 and atom2.
        If None, distances are computed within atom1 (default).

    Returns
    -------
    dist_matrix : np.ndarray (float32)
        An (N2 x N1) array of distances, transposed at the end to (N1 x N2).
    X_dist : np.ndarray (float32)
        X-displacement components for all pairs (N2 x N1, then transposed).
    Y_dist : np.ndarray (float32)
        Y-displacement components for all pairs.
    Z_dist : np.ndarray (float32)
        Z-displacement components for all pairs.

    Notes
    -----
    - This is a Python translation of a similar MATLAB function. 
    - For orthogonal boxes (3-length `box_dim`), a simple minimum-image
      approach is used by shifting coordinate differences that exceed half
      the box length.
    - For 6-length boxes, it assumes [xlo, ylo, zlo, xhi, yhi, zhi] and
      converts them to (lx, ly, lz).
    - For 9-length boxes (triclinic), the tilt factors (xy, xz, yz) are
      used in a straightforward adaptation of the minimum-image approach.
      This portion is not extensively tested here.
    """

    # If atom2 not given, compute distances within atom1
    if atom2 is None:
        atom2 = atom1

    # Convert to numpy float32 arrays of coordinates
    coords1 = np.array([[atm['x'], atm['y'], atm['z']] for atm in atom1], dtype=np.float32)
    coords2 = np.array([[atm['x'], atm['y'], atm['z']] for atm in atom2], dtype=np.float32)

    n1 = coords1.shape[0]
    n2 = coords2.shape[0]

    # Default to orthonormal box
    # We'll interpret box_dim length to set lx, ly, lz, xy, xz, yz
    box_dim = np.array(box_dim, dtype=np.float32)
    if box_dim.size == 3:
        # [lx, ly, lz]
        lx, ly, lz = box_dim
        xy, xz, yz = 0.0, 0.0, 0.0
    elif box_dim.size == 6:
        # [xlo, ylo, zlo, xhi, yhi, zhi]
        xlo, ylo, zlo, xhi, yhi, zhi = box_dim
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        xy, xz, yz = 0.0, 0.0, 0.0
    elif box_dim.size == 9:
        # [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        lx, ly, lz = box_dim[0], box_dim[1], box_dim[2]
        xy, xz, yz = box_dim[5], box_dim[7], box_dim[8]
    else:
        raise ValueError("box_dim must have length 3, 6, or 9.")

    # Prepare output arrays: shape (n2 x n1), single precision
    dist_matrix = np.zeros((n2, n1), dtype=np.float32)
    X_dist = np.zeros((n2, n1), dtype=np.float32)
    Y_dist = np.zeros((n2, n1), dtype=np.float32)
    Z_dist = np.zeros((n2, n1), dtype=np.float32)

    # For each atom i in atom1, compute displacement to all atoms in atom2
    for i in range(n1):
        rx = coords1[i,0] - coords2[:,0]
        ry = coords1[i,1] - coords2[:,1]
        rz = coords1[i,2] - coords2[:,2]

        # Orthogonal minimum-image or partial triclinic approach
        if box_dim.size == 3 or (xy == 0 and xz == 0 and yz == 0):
            # Simple orthogonal PBC
            # Shift rx, ry, rz to be within [-L/2, L/2]
            rx[rx >  lx/2] -= lx
            rx[rx < -lx/2] += lx
            ry[ry >  ly/2] -= ly
            ry[ry < -ly/2] += ly
            rz[rz >  lz/2] -= lz
            rz[rz < -lz/2] += lz

        else:
            # Triclinic approach (not heavily tested)
            # First handle z displacement
            zmask_up = (rz >  lz/2)
            zmask_dn = (rz < -lz/2)
            rz[zmask_up] -= lz
            rz[zmask_dn] += lz
            rx[zmask_up] -= xz
            rx[zmask_dn] += xz
            ry[zmask_up] -= yz
            ry[zmask_dn] += yz

            # Then y displacement
            ymask_up = (ry > ly/2)
            ymask_dn = (ry < -ly/2)
            ry[ymask_up] -= ly
            ry[ymask_dn] += ly
            rx[ymask_up] -= xy
            rx[ymask_dn] += xy

            # Finally x displacement
            xmask_up = (rx > lx/2)
            xmask_dn = (rx < -lx/2)
            rx[xmask_up] -= lx
            rx[xmask_dn] += lx

        # Distances
        r = np.sqrt(rx*rx + ry*ry + rz*rz)
        dist_matrix[:,i] = r
        X_dist[:,i] = rx
        Y_dist[:,i] = ry
        Z_dist[:,i] = rz

    # Transpose everything to match the MATLAB output: (n1 x n2)
    dist_matrix = dist_matrix.T
    X_dist = X_dist.T
    Y_dist = Y_dist.T
    Z_dist = Z_dist.T

    return dist_matrix, X_dist, Y_dist, Z_dist
