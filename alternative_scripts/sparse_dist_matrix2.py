import numpy as np
from atomipy import box_dim2cell

def cell_list_dist_matrix2(atom1, atom2, box_dim, rmaxshort=None, rmaxlong=12.0, r_atom_type='H'):
    """
    Computes the NxN distance matrix between two sets of atoms (atom1 and atom2)
    using a cell-list approach under periodic boundary conditions (PBC) for a
    triclinic (or orthogonal) box. Intra-set distances (atom1-atom1 or atom2-atom2)
    are skipped, so only cross-set distances are recorded.

    Parameters
    ----------
    atom1, atom2 : list of dict
        Each dict must contain at least:
          - 'x', 'y', 'z' for coordinates
          - 'type' (e.g., 'H', 'O', etc.)
    box_dim : list or array
        Box dimensions. Can be:
          - length 3: [lx, ly, lz] (orthogonal box)
          - length 6: [a, b, c, alpha, beta, gamma] (cell parameters)
          - length 9: [lx, ly, lz, 0, 0, xy, 0, xz, yz] (triclinic)
    rmaxshort : float or None, optional
        If set, a shorter cutoff distance specifically for pairs involving
        atoms of type `r_atom_type`. If None, only `rmaxlong` is used.
    rmaxlong : float, optional
        The default (long) cutoff distance (default = 12.0 Å).
    r_atom_type : str, optional
        The atom type to which `rmaxshort` applies (default = 'H').

    Returns
    -------
    dist_matrix : np.ndarray
        An (N1+N2)x(N1+N2) matrix of distances. Only cross-set pairs (atom1-atom2)
        within the cutoff are non-zero. Distances for i=j or same-set pairs are 0.
    bond_list : list of tuple
        Pairs (i_atom, j_atom) that lie within the respective cutoff, where
        i_atom is from 0..N1-1 and j_atom from N1..(N1+N2-1).
    dist_list : list of float
        Distances corresponding to each pair in bond_list.
    X_dist, Y_dist, Z_dist : np.ndarray
        The displacement matrices (component-wise) for all valid pairs.

    Notes
    -----
    - Uses a cell-list algorithm for efficiency, subdividing the box by ~rmaxlong.
    - If both atoms in a pair are type `r_atom_type`, or if either is, that pair
      will use the shorter cutoff `rmaxshort`.
    - Intra-set distances are excluded; only cross-set pairs are returned.
    - This function is a Pythonic port of the MATLAB function
      `cell_list_dist_matrix_atom1atom2.m`.

    Example
    -------
    >>> # Suppose atom1, atom2 each is a list of dicts with 'x','y','z','type'
    >>> # and box_dim is e.g. [lx, ly, lz]
    >>> dm, bonds, dists, Xd, Yd, Zd = cell_list_dist_matrix_atom1atom2(atom1, atom2, box_dim, 1.25, 2.25, 'H')
    >>> print(dm.shape)  # => (N1+N2, N1+N2)
    >>> print(len(bonds), len(dists))
    """

    # ----------------------------------------------------------------------
    # 1) Helper function to parse box_dim into [a,b,c,alpha,beta,gamma].
    # ----------------------------------------------------------------------
    # def box_dim2cell(bdim):
    #     """Minimal function to convert various box_dim forms into a 6-length cell array."""
    #     bdim = np.asarray(bdim, dtype=float)
    #     if bdim.size == 9:
    #         # [lx, ly, lz, 0, 0, xy, 0, xz, yz]
    #         lx, ly, lz = bdim[0], bdim[1], bdim[2]
    #         xy, xz, yz = bdim[5], bdim[7], bdim[8]
    #         a = lx
    #         b = np.sqrt(ly**2 + xy**2)
    #         c = np.sqrt(lz**2 + xz**2 + yz**2)
    #         # angles
    #         alpha = np.degrees(np.arccos((ly*yz + xy*xz)/(b*c))) if (b*c) != 0 else 90.0
    #         beta  = np.degrees(np.arccos(xz/c)) if c != 0 else 90.0
    #         gamma = np.degrees(np.arccos(xy/b)) if b != 0 else 90.0
    #         return [a, b, c, alpha, beta, gamma]
    #     elif bdim.size == 6:
    #         # [a, b, c, alpha, beta, gamma]
    #         return bdim.tolist()
    #     elif bdim.size == 3:
    #         # [lx, ly, lz] => orthogonal => alpha=beta=gamma=90
    #         return [bdim[0], bdim[1], bdim[2], 90.0, 90.0, 90.0]
    #     else:
    #         raise ValueError("box_dim must have length 3, 6, or 9.")

    # ----------------------------------------------------------------------
    # 2) Interpret box_dim => cell parameters
    # ----------------------------------------------------------------------
    cell = box_dim2cell(box_dim)
    a, b, c, alpha, beta, gamma = cell
    alpha_r = np.radians(alpha)
    beta_r  = np.radians(beta)
    gamma_r = np.radians(gamma)

    # ----------------------------------------------------------------------
    # 3) Merge the two atom lists so we can build a single cell list
    # ----------------------------------------------------------------------
    atom1 = list(atom1)  # ensure it's a list
    atom2 = list(atom2)
    N1 = len(atom1)
    N2 = len(atom2)
    all_atoms = atom1 + atom2  # combined
    N = N1 + N2

    if N == 0:
        # If no atoms, return empty results
        return (np.zeros((0,0)), [], [], np.zeros((0,0)), np.zeros((0,0)), np.zeros((0,0)))

    # ----------------------------------------------------------------------
    # 4) Build coordinate array
    # ----------------------------------------------------------------------
    coords = np.array([[atm['x'], atm['y'], atm['z']] for atm in all_atoms], dtype=float)

    # ----------------------------------------------------------------------
    # 5) Construct the triclinic box matrix H
    # ----------------------------------------------------------------------
    ax = a
    bx = b * np.cos(gamma_r)
    by = b * np.sin(gamma_r)
    cx = c * np.cos(beta_r)
    # to avoid divide-by-zero in extremely orthonormal cases, add a small epsilon
    denominator = np.sin(gamma_r) + 1e-15
    cy = c * (np.cos(alpha_r) - np.cos(beta_r)*np.cos(gamma_r)) / denominator
    cz = np.sqrt(c**2 - cx**2 - cy**2)

    H = np.array([[ax, bx, cx],
                  [0.0, by, cy],
                  [0.0, 0.0, cz]])
    Hinv = np.linalg.inv(H)

    # ----------------------------------------------------------------------
    # 6) Build the cell list
    #    We pick a cellSize ~ rmaxlong*1.5 to ensure we capture neighbor pairs.
    # ----------------------------------------------------------------------
    if rmaxlong is None:
        cell_size = 3.0
    else:
        cell_size = 1.5*rmaxlong

    # bounding box estimate
    bounding_box = H @ np.array([1.0, 1.0, 1.0])
    bounding_box_size = np.abs(bounding_box)
    n_cells = np.maximum(np.floor(bounding_box_size / cell_size).astype(int), 1)
    num_cells = np.prod(n_cells)

    # Convert to fractional coords in [0,1)
    frac_coords = (Hinv @ coords.T).T
    frac_coords = frac_coords - np.floor(frac_coords)

    # Assign each atom to a cell
    cell_idx = np.floor(frac_coords * n_cells).astype(int)
    # Clip if they land exactly on the boundary
    cell_idx = np.clip(cell_idx, [0,0,0], n_cells - 1)

    # Flatten
    lin_idx = (cell_idx[:,0] +
               cell_idx[:,1]*n_cells[0] +
               cell_idx[:,2]*n_cells[0]*n_cells[1])
    cell_list = [[] for _ in range(num_cells)]
    for i_atom, c_id in enumerate(lin_idx):
        cell_list[c_id].append(i_atom)

    # ----------------------------------------------------------------------
    # 7) Prepare output arrays
    # ----------------------------------------------------------------------
    dist_matrix = np.zeros((N, N), dtype=float)
    X_dist = np.zeros((N, N), dtype=float)
    Y_dist = np.zeros((N, N), dtype=float)
    Z_dist = np.zeros((N, N), dtype=float)
    bond_list = []
    dist_list = []

    # neighbor offsets => 27 neighbors for each cell (PBC)
    offsets = [(ix, iy, iz) for ix in (-1,0,1) for iy in (-1,0,1) for iz in (-1,0,1)]

    # ----------------------------------------------------------------------
    # 8) Loop over cells & neighbors, compute cross-set distances
    # ----------------------------------------------------------------------
    for c_id in range(num_cells):
        atom_list_c = cell_list[c_id]
        if not atom_list_c:
            continue

        # convert c_id -> (cx, cy, cz) in cell space
        cx_idx = c_id % n_cells[0]
        cy_idx = (c_id // n_cells[0]) % n_cells[1]
        cz_idx = c_id // (n_cells[0]*n_cells[1])

        for (dx, dy, dz) in offsets:
            nx = (cx_idx + dx) % n_cells[0]
            ny = (cy_idx + dy) % n_cells[1]
            nz = (cz_idx + dz) % n_cells[2]
            neighbor_lin = nx + ny*n_cells[0] + nz*(n_cells[0]*n_cells[1])
            atom_list_n = cell_list[neighbor_lin]
            if not atom_list_n:
                continue

            # avoid double-counting
            if neighbor_lin < c_id:
                continue

            for i_atom in atom_list_c:
                for j_atom in atom_list_n:
                    # skip duplicates in same cell
                    if neighbor_lin == c_id and j_atom <= i_atom:
                        continue

                    # skip pairs that are within the same set
                    # i.e., both in atom1 or both in atom2
                    i_in_first = (i_atom < N1)
                    j_in_first = (j_atom < N1)
                    if (i_in_first and j_in_first) or (not i_in_first and not j_in_first):
                        continue

                    # Choose cutoff based on whether either atom is r_atom_type
                    if rmaxshort is not None:
                        if (all_atoms[i_atom]['type'] == r_atom_type or
                            all_atoms[j_atom]['type'] == r_atom_type):
                            cutoff = rmaxshort
                        else:
                            cutoff = rmaxlong
                    else:
                        cutoff = rmaxlong

                    # minimum-image
                    diff_frac = (Hinv @ (coords[j_atom] - coords[i_atom])).T
                    diff_frac -= np.round(diff_frac)
                    d_vec = (H @ diff_frac).T
                    d_val = np.linalg.norm(d_vec)

                    if d_val <= cutoff:
                        bond_list.append((i_atom, j_atom))
                        dist_list.append(d_val)
                        dist_matrix[i_atom, j_atom] = d_val
                        dist_matrix[j_atom, i_atom] = d_val
                        # displacement
                        X_dist[i_atom, j_atom] = d_vec[0]
                        X_dist[j_atom, i_atom] = -d_vec[0]
                        Y_dist[i_atom, j_atom] = d_vec[1]
                        Y_dist[j_atom, i_atom] = -d_vec[1]
                        Z_dist[i_atom, j_atom] = d_vec[2]
                        Z_dist[j_atom, i_atom] = -d_vec[2]

    return dist_matrix, bond_list, dist_list, X_dist, Y_dist, Z_dist
