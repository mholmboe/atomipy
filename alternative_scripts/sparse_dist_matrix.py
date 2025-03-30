import numpy as np

def sparse_dist_matrix(atom, box_dim, rmaxshort=None, rmaxlong=12.0, r_atom_type='H'):
    """
    Computes the NxN distance matrix for all atom pairs using a cell-list approach
    under periodic boundary conditions (PBC) for a triclinic (or orthogonal) box.
    
    Parameters
    ----------
    atom : list of dict
        Each dict must contain at least:
         - 'x', 'y', 'z' (atom coordinates in Angstrom)
         - 'type' (atom type, e.g. 'H', 'O', etc.)
    box_dim : list or array
        The box dimensions, which can be:
         - length 3: [lx, ly, lz] (orthogonal box)
         - length 6: [a, b, c, alpha, beta, gamma] (cell parameters)
         - length 9: [lx, ly, lz, 0,0, xy, 0, xz, yz] (triclinic)
    rmaxshort : float or None, optional
        If set, a shorter cutoff distance for specified atom types (e.g. 'H').
        If None, only one cutoff (rmaxlong) is used.
    rmaxlong : float, optional
        Long-range cutoff distance (default 12.0 Å).
    r_atom_type : str, optional
        The atom 'type' to which rmaxshort applies (default 'H').

    Returns
    -------
    dist_matrix : np.ndarray
        NxN distance matrix (0 where distance > cutoff or i=j).
    bond_list : list of tuple
        Pairs of (i_atom, j_atom) that lie within the relevant cutoff.
    dist_list : list of float
        Distances corresponding to bond_list.
    X_dist, Y_dist, Z_dist : np.ndarray
        NxN displacement matrices (component-wise) for pairs within cutoff.
    
    Notes
    -----
    - This function uses a cell-list algorithm for efficiency, subdividing the box
      into smaller cells of size ~rmaxlong.
    - If rmaxshort is provided, any pair involving `r_atom_type` uses that shorter cutoff.
    - The final result includes full NxN symmetric distance and displacement matrices.
    - This is a Pythonic rewrite of the original MATLAB function 
      `cell_list_dist_matrix_atom.m`.
    
    Example
    -------
    >>> # Suppose 'atom' is a list of dicts, each with 'x','y','z','type'
    >>> # and box_dim is [lx, ly, lz]
    >>> dm, bonds, dists, Xd, Yd, Zd = cell_list_dist_matrix_atom(atom, box_dim, 1.25, 2.25, 'H')
    >>> print(dm.shape)  # NxN
    >>> print(len(bonds), len(dists))
    """
    # ----------------------------------------------------------------------
    # 1) Helper function to parse box_dim into [a,b,c,alpha,beta,gamma].
    # ----------------------------------------------------------------------
    def box_dim2cell(bdim):
        """Minimal function to convert various box_dim forms into a 6-length cell array."""
        bdim = np.asarray(bdim, dtype=float)
        if bdim.size == 9:
            # [lx, ly, lz, 0, 0, xy, 0, xz, yz]
            lx, ly, lz = bdim[0], bdim[1], bdim[2]
            xy, xz, yz = bdim[5], bdim[7], bdim[8]
            a = lx
            b = np.sqrt(ly**2 + xy**2)
            c = np.sqrt(lz**2 + xz**2 + yz**2)
            # angles
            alpha = np.degrees(np.arccos((ly*yz + xy*xz)/(b*c))) if (b*c) != 0 else 90.0
            beta  = np.degrees(np.arccos(xz/c)) if c != 0 else 90.0
            gamma = np.degrees(np.arccos(xy/b)) if b != 0 else 90.0
            return [a, b, c, alpha, beta, gamma]
        elif bdim.size == 6:
            # [a, b, c, alpha, beta, gamma]
            return bdim.tolist()
        elif bdim.size == 3:
            # [lx, ly, lz] => orthogonal => alpha=beta=gamma=90
            return [bdim[0], bdim[1], bdim[2], 90.0, 90.0, 90.0]
        else:
            raise ValueError("box_dim must have length 3, 6, or 9.")

    # ----------------------------------------------------------------------
    # 2) Parse input
    # ----------------------------------------------------------------------
    # If only two inputs in MATLAB, they set rmaxlong=2.45 by default. 
    # We'll keep the default 12.0 here, but you can override it below:
    # (Just a note: in the original code, they used 2.45 or 12. 
    #  We'll keep 12 for a more general case, but you can change.)
    # If you want to emulate exactly: rmaxlong = 2.45 if rmaxshort is None else ...
    # For clarity, we won't exactly replicate that. We'll trust the user to pass what they want.

    # Convert box_dim => cell
    cell = box_dim2cell(box_dim)
    a, b, c, alpha, beta, gamma = cell

    # Angles in radians
    alpha_r = np.radians(alpha)
    beta_r  = np.radians(beta)
    gamma_r = np.radians(gamma)

    # Number of atoms
    N = len(atom)
    if N == 0:
        # Quick return if no atoms
        return (np.zeros((0,0)), [], [], np.zeros((0,0)), np.zeros((0,0)), np.zeros((0,0)))

    # Coordinates in Nx3
    coords = np.array([[atm['x'], atm['y'], atm['z']] for atm in atom], dtype=float)

    # ----------------------------------------------------------------------
    # 3) Construct the triclinic box matrix H, and its inverse
    #    (For orthogonal box, gamma=beta=alpha=90, so it reduces neatly.)
    # ----------------------------------------------------------------------
    ax = a
    bx = b * np.cos(gamma_r)
    by = b * np.sin(gamma_r)
    cx = c * np.cos(beta_r)
    cy = c*(np.cos(alpha_r) - np.cos(beta_r)*np.cos(gamma_r)) / (np.sin(gamma_r) + 1e-15)
    cz = np.sqrt(c**2 - cx**2 - cy**2)

    H = np.array([[ax, bx, cx],
                  [0.0, by, cy],
                  [0.0, 0.0, cz]], dtype=float)
    Hinv = np.linalg.inv(H)

    # ----------------------------------------------------------------------
    # 4) Cell-list construction
    # ----------------------------------------------------------------------
    cell_size = rmaxlong
    bounding_box = H @ np.array([1,1,1], dtype=float)
    # get box extents in each direction
    # (avoid overhead, just approximate)
    bounding_box_size = np.abs(bounding_box)
    n_cells = np.maximum(np.floor(bounding_box_size / cell_size).astype(int), 1)
    num_cells = np.prod(n_cells)

    # Convert coords -> fractional in [0,1)
    frac_coords = (Hinv @ coords.T).T
    frac_coords = frac_coords - np.floor(frac_coords)

    # Assign each atom to a cell
    cell_indices = np.floor(frac_coords * n_cells).astype(int)
    # Clip boundaries
    cell_indices = np.clip(cell_indices, [0,0,0], n_cells - 1)

    # Flatten cell indices
    lin_idx = (cell_indices[:,0] +
               cell_indices[:,1]*n_cells[0] +
               cell_indices[:,2]*n_cells[0]*n_cells[1])
    cell_list = [[] for _ in range(num_cells)]
    for i_atom, c_id in enumerate(lin_idx):
        cell_list[c_id].append(i_atom)

    # Neighbor offsets
    neighbor_offsets = [-1, 0, 1]
    offsets = []
    for ix in neighbor_offsets:
        for iy in neighbor_offsets:
            for iz in neighbor_offsets:
                offsets.append((ix, iy, iz))

    # ----------------------------------------------------------------------
    # 5) Prepare outputs
    # ----------------------------------------------------------------------
    dist_matrix = np.zeros((N, N), dtype=float)
    X_dist = np.zeros((N, N), dtype=float)
    Y_dist = np.zeros((N, N), dtype=float)
    Z_dist = np.zeros((N, N), dtype=float)
    bond_list = []
    dist_list = []

    # ----------------------------------------------------------------------
    # 6) Loop over each cell & its neighbors
    # ----------------------------------------------------------------------
    for c_id in range(num_cells):
        cell_atoms = cell_list[c_id]
        if not cell_atoms:
            continue
        # Convert c_id -> (cx, cy, cz)
        cx = c_id % n_cells[0]
        cy = (c_id // n_cells[0]) % n_cells[1]
        cz = c_id // (n_cells[0]*n_cells[1])

        for offset in offsets:
            nx = (cx + offset[0]) % n_cells[0]
            ny = (cy + offset[1]) % n_cells[1]
            nz = (cz + offset[2]) % n_cells[2]
            neighbor_lin = nx + ny*n_cells[0] + nz*n_cells[0]*n_cells[1]
            neighbor_atoms = cell_list[neighbor_lin]
            if not neighbor_atoms:
                continue

            # avoid double counting
            if neighbor_lin < c_id:
                continue

            # pairwise distance
            for i_atom in cell_atoms:
                for j_atom in neighbor_atoms:
                    if neighbor_lin == c_id and j_atom <= i_atom:
                        continue

                    # Choose local cutoff
                    if rmaxshort is not None:
                        # If either is r_atom_type, use shorter cutoff
                        if (atom[i_atom]['type'] == r_atom_type or
                            atom[j_atom]['type'] == r_atom_type):
                            cutoff = rmaxshort
                        else:
                            cutoff = rmaxlong
                    else:
                        cutoff = rmaxlong

                    # fractional difference => min image => back to cart
                    diff_frac = (Hinv @ (coords[j_atom] - coords[i_atom])).T
                    # shift into [-0.5, 0.5)
                    diff_frac = diff_frac - np.round(diff_frac)
                    d_vec = (H @ diff_frac).T
                    dist_val = np.linalg.norm(d_vec)

                    if dist_val <= cutoff:
                        bond_list.append((i_atom, j_atom))
                        dist_list.append(dist_val)
                        dist_matrix[i_atom, j_atom] = dist_val
                        dist_matrix[j_atom, i_atom] = dist_val
                        X_dist[i_atom, j_atom] = d_vec[0]
                        X_dist[j_atom, i_atom] = -d_vec[0]
                        Y_dist[i_atom, j_atom] = d_vec[1]
                        Y_dist[j_atom, i_atom] = -d_vec[1]
                        Z_dist[i_atom, j_atom] = d_vec[2]
                        Z_dist[j_atom, i_atom] = -d_vec[2]

    return dist_matrix, bond_list, dist_list, X_dist, Y_dist, Z_dist
