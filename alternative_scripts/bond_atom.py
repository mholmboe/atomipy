import math
import numpy as np

def bond_atom(atom, box_dim, rmaxlong=2.25, distance_factor=0.65, rmaxshort=1.18):
    """
    Assigns bonds and angles to each atom in 'atom' by computing pairwise
    distances under PBC and comparing against radii-based cutoffs.

    Parameters
    ----------
    atom : list of dict
        Each dict represents an atom, including:
          - 'x', 'y', 'z' (coordinates)
          - 'type' (initial guess for element or atom type)
          - 'molid' (molecule ID)
          - 'index' (atom index)
        Optionally 'element' can be present; otherwise guessed from 'type'.
    box_dim : list or array
        Box dimensions (orthogonal or triclinic). Used by
        `sparse_dist_matrix` for PBC.
    rmaxlong : float, optional
        Default max cutoff for bonds (default=2.25 Å).
    distance_factor : float, optional
        Factor to scale the sum of atomic radii in deciding bonds (default=0.65).
    rmaxshort : float, optional
        Shorter cutoff for special cases (e.g. 'H') (default=1.18 Å).

    Returns
    -------
    atom : list of dict
        Updated with fields:
          - atom[i]['neigh']: neighbor distances, indices, etc.
          - atom[i]['bond']:  bond distances, indices, etc.
          - atom[i]['angle']: angle info if enough neighbors
    bond_index : np.ndarray
        An Nx3 array: [i_atom, j_atom, distance].
    angle_index : np.ndarray
        An Mx10 array describing angles at the central atom i.
    dist_matrix, X_dist, Y_dist, Z_dist : np.ndarray
        NxN matrices for distances and displacement components.
    neigh_index : np.ndarray
        [j, i, distance] for neighbor pairs within rmaxlong in the same molecule.

    Notes
    -----
    - Relies on a function `cell_list_dist_matrix_atom(atom, box_dim, rmaxshort, rmaxlong)`
      that returns (dist_matrix, bond_list, dist_list, X_dist, Y_dist, Z_dist).
    - We guess 'element' from 'type' if missing, then set atom[i]['type'] = element.
    - The code assigns approximate vdw radii, then zeros out pairs where distance
      exceeds (r1 + r2)*distance_factor (capped by rmaxlong).
    - We replicate the neighbor/bond angle logic from the MATLAB code.
    """

    # Helper to guess element if not provided
    def guess_element(atom_type):
        t = atom_type.strip().capitalize()
        if t.startswith('Si') or t.startswith('Sy') or t.startswith('Sc'):
            return 'Si'
        elif t.startswith('S'):
            return 'S'
        elif t.startswith('Al'):
            return 'Al'
        elif t.startswith('Fe'):
            return 'Fe'
        elif t.startswith('Mg'):
            return 'Mg'
        elif t.startswith('O'):
            return 'O'
        elif t.startswith('H'):
            return 'H'
        elif t.startswith('C'):
            return 'C'
        elif t.startswith('Na'):
            return 'Na'
        elif t.startswith('K'):
            return 'K'
        return t[:2]

    # 1) Ensure each atom has 'element'
    for a in atom:
        if 'element' not in a:
            a['element'] = guess_element(a['type'])
        # Now set 'type' = 'element'
        a['type'] = a['element']
        # Also add 'fftype' if missing
        if 'fftype' not in a:
            a['fftype'] = a['type']

    # 2) Compute distance matrix (and displacements) with the cell-list approach
    #    You need a function cell_list_dist_matrix_atom(atom, box_dim, rmaxshort, rmaxlong)
    from sparse_dist_matrix import sparse_dist_matrix

    dist_matrix, bond_list, dist_list, X_dist, Y_dist, Z_dist = \
        sparse_dist_matrix(atom, box_dim, rmaxshort, rmaxlong)
    n_atoms = len(atom)

    # 3) Approximate radii (vdw or ionic):
    vdw_radii = {
        'H': 0.30,   # special for hydrogen
        'O': 0.66,
        'C': 0.77,
        'Si':1.46,
        'Na':1.02,
        'K': 1.38,
        'Al':1.25,
        'S': 1.70,
        'Fe':1.20,
    }
    fallback_radius = 1.00

    # Build array of radii
    radii = np.zeros(n_atoms, dtype=np.float32)
    for i, a in enumerate(atom):
        el = a['element']
        radii[i] = vdw_radii.get(el, fallback_radius)

    # Combine radii in outer sum, scale by distance_factor
    radius_limit = np.add.outer(radii, radii) * distance_factor
    # Cap by rmaxlong
    radius_limit[radius_limit > rmaxlong] = rmaxlong
    # Zero out distance_matrix if dist>radius_limit
    mask = dist_matrix > radius_limit
    dist_matrix[mask] = 0.0

    # 4) Clear existing 'neigh', 'bond', 'angle' if any
    for a in atom:
        for field in ('neigh','bond','angle'):
            if field in a:
                del a[field]

    # 5) Build neighbor/bond info
    coords = np.array([[a['x'], a['y'], a['z']] for a in atom], dtype=np.float32)
    bond_index_list = []  # for storing [i_min, i_max, dist]
    neigh_indices = [[] for _ in range(n_atoms)]   # store neighbor indices for angles
    neigh_vectors = [[] for _ in range(n_atoms)]   # store vectors for angles

    for i in range(n_atoms):
        dcol = dist_matrix[:, i]
        j_list = np.where(dcol>0)[0]
        neigh_data = {'dist':[], 'index':[], 'type':[], 'coords':[], 'r_vec':[]}
        bond_data = {'dist':[], 'index':[], 'type':[]}

        for j in j_list:
            if atom[i]['molid'] == atom[j]['molid']:  # same molecule => bond
                dist_ij = dcol[j]
                neigh_data['dist'].append(dist_ij)
                neigh_data['index'].append(j)
                neigh_data['type'].append(atom[j]['type'])
                neigh_data['coords'].append([atom[j]['x'], atom[j]['y'], atom[j]['z']])
                # i->j vector is -(X_dist[j,i], Y_dist[j,i], Z_dist[j,i]) since X_dist[j,i] is j->i
                vec_ij = [-X_dist[j,i], -Y_dist[j,i], -Z_dist[j,i]]
                neigh_data['r_vec'].append(vec_ij)

                bond_data['dist'].append(dist_ij)
                bond_data['index'].append([i,j])
                bond_data['type'].append(1)  # arbitrary

                i_min, i_max = min(i,j), max(i,j)
                bond_index_list.append((i_min, i_max, dist_ij))

                neigh_indices[i].append(j)
                neigh_vectors[i].append(vec_ij)

        atom[i]['neigh'] = neigh_data
        atom[i]['bond']  = bond_data

    # 6) Angle detection
    angle_index_list = []
    for i in range(n_atoms):
        n_inds = neigh_indices[i]
        n_vecs = neigh_vectors[i]
        if len(n_inds) < 2:
            continue

        for v in range(len(n_inds)):
            for w in range(v+1, len(n_inds)):
                vec_v = np.array(n_vecs[v], dtype=np.float32)
                vec_w = np.array(n_vecs[w], dtype=np.float32)
                cross_mag = np.linalg.norm(np.cross(vec_v, vec_w))
                dot_vw = float(np.dot(vec_v, vec_w))
                denom = float(np.linalg.norm(vec_v)*np.linalg.norm(vec_w))
                if denom<1e-12:
                    continue
                angle_deg = math.degrees(math.atan2(cross_mag, dot_vw))
                if 0.0 < angle_deg <= 180.0:
                    i_v = n_inds[v]
                    i_w = n_inds[w]
                    angle_index_list.append([
                        i_v, i, i_w, angle_deg,
                        vec_v[0], vec_v[1], vec_v[2],
                        vec_w[0], vec_w[1], vec_w[2]
                    ])

        if angle_index_list:
            # gather angle info for this atom i
            # we filter those that have center i
            # but let's do it outside the main loop for clarity
            pass

    # 7) Clean up bond_index => unique
    bond_index = np.array(bond_index_list, dtype=np.float32)
    if bond_index.size == 0:
        bond_index = bond_index.reshape((0,3))
    else:
        # sort by first col, then second
        order = np.lexsort((bond_index[:,1], bond_index[:,0]))
        bond_index = bond_index[order,:]
        # unique => round dist for grouping
        dist_rounded = np.round(bond_index[:,2], 4)
        combined = np.column_stack((bond_index[:,:2], dist_rounded))
        _, unique_idx = np.unique(combined, axis=0, return_index=True)
        bond_index = bond_index[sorted(unique_idx), :]

    # 8) Clean up angle_index => unique
    angle_index = np.array(angle_index_list, dtype=np.float32)
    if angle_index.size == 0:
        angle_index = angle_index.reshape((0,10))
    else:
        order = np.lexsort((angle_index[:,0], angle_index[:,1]))  # sort by center -> j1
        angle_index = angle_index[order,:]
        angle_rounded = np.round(angle_index[:,3], 3)
        combined = np.column_stack((angle_index[:,:3], angle_rounded))
        _, unique_idx = np.unique(combined, axis=0, return_index=True)
        angle_index = angle_index[sorted(unique_idx), :]

    # 9) Build neigh_index from dist_matrix => [j, i, dist]
    #    filtering out dist>rmaxlong or different molid
    tri_mask = dist_matrix>0
    rows, cols = np.where(np.tril(tri_mask))
    dvals = dist_matrix[rows, cols]
    neigh_data = np.column_stack((cols, rows, dvals))  # => [j, i, dist]

    # sort by j => col=0
    order = np.argsort(neigh_data[:,0])
    neigh_data = neigh_data[order,:]
    # unique => columns 0:2
    _, uniq_idx = np.unique(neigh_data[:,:2], axis=0, return_index=True)
    neigh_data = neigh_data[sorted(uniq_idx), :]

    # remove if dist>rmaxlong or different molid
    to_remove = []
    for idx, row in enumerate(neigh_data):
        j_, i_, dist_ = row
        j_, i_ = int(j_), int(i_)
        if dist_>rmaxlong:
            to_remove.append(idx)
        elif atom[j_]['molid'] != atom[i_]['molid']:
            to_remove.append(idx)
    neigh_index = np.delete(neigh_data, to_remove, axis=0)

    # 10) Store angle info in each central atom 
    #     (Angle_index is global, but we can store a subset in atom[i]['angle'] if desired)
    angle_map = {}
    for row in angle_index:
        j1, i_ctr, j2, ang_deg = row[:4]
        i_ctr = int(i_ctr)
        if i_ctr not in angle_map:
            angle_map[i_ctr] = []
        angle_map[i_ctr].append(row)

    for i_ctr, rows in angle_map.items():
        atom[i_ctr]['angle'] = {
            'type': 1,
            'index': [r[:3].tolist() for r in rows],
            'angle': [r[3] for r in rows],
            'vec1':  [r[4:7].tolist() for r in rows],
            'vec2':  [r[7:10].tolist() for r in rows]
        }

    # Return
    return atom, bond_index, angle_index, dist_matrix, X_dist, Y_dist, Z_dist, neigh_index
