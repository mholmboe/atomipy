"""
Analysis module for atomipy - provides functions for structural analysis like RDF,
coordination numbers, and unwrapping coordinates.
"""

import numpy as np
from .distances import dist_matrix, get_neighbor_list
from .transform import wrap_coordinates

def unwrap_coordinates(atoms, Box, molid=None):
    """
    Unwrap atom coordinates to fix molecules split across periodic boundaries.
    
    This function shifts atoms within the same molecule such that they are
    clustered together according to the minimum image convention relative
    to the first atom of the molecule.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
    Box : list of float
        Box dimensions (1x3, 1x6, or 1x9).
    molid : int or list of int, optional
        If provided, only unwrap molecules with these IDs.
        
    Returns
    -------
    list of dict
        Atoms with unwrapped coordinates.
    """
    if not atoms:
        return []
    
    # Identify molecules
    molids = np.array([a['molid'] for a in atoms])
    unique_molids = np.unique(molids)
    
    if molid is not None:
        if isinstance(molid, int):
            unique_molids = [molid]
        else:
            unique_molids = molid

    # Create a copy to avoid modifying originals
    unwrapped_atoms = [a.copy() for a in atoms]
    
    # We need to handle this molecule by molecule
    for mid in unique_molids:
        indices = np.where(molids == mid)[0]
        if len(indices) <= 1:
            continue
            
        # Use the first atom as reference
        ref_idx = indices[0]
        ref_pos = np.array([unwrapped_atoms[ref_idx]['x'], 
                           unwrapped_atoms[ref_idx]['y'], 
                           unwrapped_atoms[ref_idx]['z']])
        
        for idx in indices[1:]:
            pos = np.array([unwrapped_atoms[idx]['x'], 
                            unwrapped_atoms[idx]['y'], 
                            unwrapped_atoms[idx]['z']])
            
            # Find the shift needed for minimum image relative to ref_pos
            diff = pos - ref_pos
            # Use wrap_coordinates logic or similar minimum image shift
            # For simplicity, we assume we want the nearest image
            # We can use the logic from distances.py if we want code reuse
            # But here we just need to shift current atom's coordinates
            
            # TODO: Full triclinic support for unwrapping if Box is 1x9
            # For now, let's start with orthogonal support and expand
            if len(Box) == 3:
                lx, ly, lz = Box
                shift_x = -lx * np.round(diff[0] / lx)
                shift_y = -ly * np.round(diff[1] / ly)
                shift_z = -lz * np.round(diff[2] / lz)
                
                unwrapped_atoms[idx]['x'] += shift_x
                unwrapped_atoms[idx]['y'] += shift_y
                unwrapped_atoms[idx]['z'] += shift_z
            elif len(Box) == 9:
                # Triclinic unwrapping (Gromacs Box_dim format: lx ly lz 0 0 xy 0 xz yz)
                lx, ly, lz = Box[0], Box[1], Box[2]
                xy, xz, yz = Box[5], Box[7], Box[8]
                
                # Z direction
                sz = np.round(diff[2] / lz)
                unwrapped_atoms[idx]['z'] -= sz * lz
                unwrapped_atoms[idx]['y'] -= sz * yz
                unwrapped_atoms[idx]['x'] -= sz * xz
                
                # Recalculate diff for Y
                diff = np.array([unwrapped_atoms[idx]['x'], unwrapped_atoms[idx]['y'], unwrapped_atoms[idx]['z']]) - ref_pos
                sy = np.round(diff[1] / ly)
                unwrapped_atoms[idx]['y'] -= sy * ly
                unwrapped_atoms[idx]['x'] -= sy * xy
                
                # Recalculate diff for X
                diff = np.array([unwrapped_atoms[idx]['x'], unwrapped_atoms[idx]['y'], unwrapped_atoms[idx]['z']]) - ref_pos
                sx = np.round(diff[0] / lx)
                unwrapped_atoms[idx]['x'] -= sx * lx
                
    return unwrapped_atoms

def calculate_rdf(atoms, Box, rmax=15.0, dr=0.1, atom_types=None, pair_types=None, typeA=None, typeB=None):
    """
    Calculate the Radial Distribution Function g(r) for specified atom pairs.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
    Box : list of float
        Box dimensions.
    rmax : float
        Maximum distance for RDF.
    dr : float
        Bin width for RDF.
    atom_types : list of str, optional
        Filter atoms by these types for both particles in the pair.
    pair_types : tuple of lists, optional
        (types_A, types_B) to compute RDF between elements of A and B.
        
    Returns
    -------
    tuple
        (r, g_r) - Bin centers and RDF values.
    """
    if not atoms:
        return np.array([]), np.array([])

    # Map typeA/typeB to pair_types if provided
    if typeA and typeB:
        pair_types = ([typeA] if isinstance(typeA, str) else typeA, 
                      [typeB] if isinstance(typeB, str) else typeB)
    elif typeA:
        atom_types = [typeA] if isinstance(typeA, str) else typeA

    if pair_types:
        types_a, types_b = pair_types
        indices_a = [i for i, a in enumerate(atoms) if a.get('type') in types_a]
        indices_b = [i for i, a in enumerate(atoms) if a.get('type') in types_b]
        is_cross = True
    elif atom_types:
        indices_a = [i for i, a in enumerate(atoms) if a.get('type') in atom_types]
        indices_b = indices_a
        is_cross = False
    else:
        indices_a = list(range(len(atoms)))
        indices_b = indices_a
        is_cross = False

    if not indices_a or not indices_b:
        return np.array([]), np.array([])

    # Use dist_matrix for all-to-all distances
    # For large systems, we should use cell_list for performance, but dist_matrix is easier for RDF
    # actually distances.py has get_neighbor_list which is better for large systems
    
    # Construct histogram
    bins = np.arange(0, rmax + dr, dr)
    hist = np.zeros(len(bins) - 1)
    
    # Average density rho = N/V
    # Get box volume
    if len(Box) == 3:
        vol = np.prod(Box)
    elif len(Box) == 9:
        vol = Box[0] * Box[1] * Box[2]
    else:
        # Fallback for 1x6 Cell
        from .cell_utils import Cell2Box_dim
        bdim = Cell2Box_dim(Box)
        vol = bdim[0] * bdim[1] * bdim[2]
        
    rho = len(indices_b) / vol
    
    # Calculate distances
    # If system is small, use full matrix
    if len(atoms) < 5000:
        dists, _, _, _ = dist_matrix(atoms, Box)
        # Subset to indices
        subset_dists = dists[np.ix_(indices_a, indices_b)]
        
        # Flatten and filter
        valid_dists = subset_dists[(subset_dists > 1e-7) & (subset_dists <= rmax)]
        h, _ = np.histogram(valid_dists, bins=bins)
        hist += h
    else:
        # For large systems, use neighbor list to avoid memory overflow
        # get_neighbor_list returns (i, j, d, dx, dy, dz)
        # We might need to iterate over indices_a if it's too big
        # For now, let's assume dist_matrix is used as a baseline and we optimizes later
        dists, _, _, _ = dist_matrix(atoms, Box)
        subset_dists = dists[np.ix_(indices_a, indices_b)]
        valid_dists = subset_dists[(subset_dists > 1e-7) & (subset_dists <= rmax)]
        h, _ = np.histogram(valid_dists, bins=bins)
        hist += h

    # Normalize RDF
    r = bins[:-1] + dr/2
    # shell volume V(r) = 4/3 * pi * ((r+dr)^3 - r^3)
    shell_v = 4/3 * np.pi * (bins[1:]**3 - bins[:-1]**3)
    
    # g(r) = hist(r) / (N_a * rho * V(r))
    # where N_a is number of atoms in type A
    g_r = hist / (len(indices_a) * rho * shell_v)
    
    return r, g_r

def coordination_number(atoms, Box, cutoff=3.0, atom_types=None, neighbor_types=None, typeA=None, typeB=None):
    """
    Calculate the coordination number for each atom.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
    Box : list of float
        Box dimensions.
    cutoff : float, optional
        Distance cutoff for neighbors.
    atom_types : list of str, optional
        Filter target atoms.
    neighbor_types : list of str, optional
        Filter neighbor atoms.
        
    Returns
    -------
    list of int
        Coordination numbers for each atom in `atoms` (filtered or all).
    """
    # Map typeA/typeB to atom_types/neighbor_types if provided
    if typeA: atom_types = [typeA] if isinstance(typeA, str) else typeA
    if typeB: neighbor_types = [typeB] if isinstance(typeB, str) else typeB

    # Use neighbor list for performance
    i_idx, j_idx, d, _, _, _ = get_neighbor_list(atoms, Box, cutoff)
    
    cn = np.zeros(len(atoms), dtype=int)
    
    # If we have filters, we need to handle them
    if neighbor_types:
        neighbor_mask = np.array([a.get('type') in neighbor_types for a in atoms])
        # j_idx refers to indices in 'atoms'
        valid_mask = neighbor_mask[j_idx]
        i_idx = i_idx[valid_mask]
        j_idx = j_idx[valid_mask]
        
    # Count occurrences in i_idx
    unique_i, counts = np.unique(i_idx, return_counts=True)
    cn[unique_i] = counts
    
    # Also count in j_idx (since neighbor list might be half-matrix if implemented that way)
    # Actually get_neighbor_list in distances.py uses upper triangle mask
    # so we MUST count both i and j.
    unique_j, counts_j = np.unique(j_idx, return_counts=True)
    cn[unique_j] += counts_j
    
    if atom_types:
        target_indices = [i for i, a in enumerate(atoms) if a.get('type') in atom_types]
        return cn[target_indices].tolist()
        
    return cn.tolist()

def closest_atom(atoms, reference, Box=None):
    """
    Return the atom closest to a reference coordinate or atom.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries to search.
    reference : dict or list of float
        Either an atom dictionary or [x, y, z] coordinates.
    Box : list of float, optional
        Box dimensions for PBC distance.
        
    Returns
    -------
    dict
        The closest atom dictionary.
    """
    if not atoms:
        return None
        
    if isinstance(reference, dict):
        ref_pos = np.array([reference['x'], reference['y'], reference['z']])
    else:
        ref_pos = np.array(reference)
        
    positions = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    diff = positions - ref_pos
    
    if Box is not None:
        # PBC correction
        if len(Box) == 3:
            diff -= np.round(diff / np.array(Box)) * np.array(Box)
        elif len(Box) == 9:
            # Simple orthogonal approximation for Box[0..2] if not handling full triclinic here
            # But let's use the logic from distances if possible or just implement orthogonal-lite
            lx, ly, lz = Box[0], Box[1], Box[2]
            diff[:, 0] -= lx * np.round(diff[:, 0] / lx)
            diff[:, 1] -= ly * np.round(diff[:, 1] / ly)
            diff[:, 2] -= lz * np.round(diff[:, 2] / lz)

    dists_sq = np.sum(diff**2, axis=1)
    min_idx = np.argmin(dists_sq)
    
    return atoms[min_idx]


def min_distances(atoms, Box, group_by='molid', n_pairs=10, cutoff=None):
    """
    Calculate the minimum distances between different molecules/residues.

    For each pair of distinct groups (molecules or residues), finds the
    single closest atom–atom contact, then returns the `n_pairs` shortest
    contacts sorted by distance.  Useful for detecting clashes after
    solvation or verifying inter-molecular separation.

    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries, each with at least 'x', 'y', 'z',
        'molid' (int), and 'resname' (str) keys.
    Box : list of float
        Box dimensions (1×3, 1×6, or 1×9).
    group_by : str
        'molid' (default) — group atoms by their molid integer.
        'resname' — group atoms by their residue name.
    n_pairs : int
        Maximum number of closest inter-group pairs to return (default 10).
    cutoff : float or None
        If given, only consider pairs whose minimum distance is ≤ cutoff Å.

    Returns
    -------
    list of dict
        Sorted list (ascending distance) of inter-group contacts.  Each
        entry has keys:
          'group_a'   – label of the first group
          'group_b'   – label of the second group
          'atom_a'    – index of the closest atom in group A
          'atom_b'    – index of the closest atom in group B
          'type_a'    – atom type of the closest atom in group A
          'type_b'    – atom type of the closest atom in group B
          'distance'  – minimum distance in Å (float, 2 d.p.)
    """
    if not atoms:
        return []

    # Determine box vectors for PBC
    if len(Box) == 3:
        lx, ly, lz = float(Box[0]), float(Box[1]), float(Box[2])
        box_vec = np.array([lx, ly, lz])
    elif len(Box) >= 9:
        lx, ly, lz = float(Box[0]), float(Box[1]), float(Box[2])
        box_vec = np.array([lx, ly, lz])
    elif len(Box) == 6:
        from .cell_utils import Cell2Box_dim
        bdim = Cell2Box_dim(Box)
        lx, ly, lz = float(bdim[0]), float(bdim[1]), float(bdim[2])
        box_vec = np.array([lx, ly, lz])
    else:
        lx, ly, lz = 1e9, 1e9, 1e9
        box_vec = None

    # Build groups
    key_func = (lambda a: a.get('molid', 0)) if group_by == 'molid' else (lambda a: a.get('resname', '?'))
    groups = {}
    for idx, atom in enumerate(atoms):
        key = key_func(atom)
        if key not in groups:
            groups[key] = []
        groups[key].append(idx)

    group_keys = list(groups.keys())
    n_groups = len(group_keys)
    if n_groups < 2:
        return []

    # Pre-build position array
    pos = np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=float)

    results = []
    for gi in range(n_groups):
        for gj in range(gi + 1, n_groups):
            key_a = group_keys[gi]
            key_b = group_keys[gj]
            idx_a = groups[key_a]
            idx_b = groups[key_b]

            pos_a = pos[idx_a]   # shape (Na, 3)
            pos_b = pos[idx_b]   # shape (Nb, 3)

            # Broadcast pairwise difference
            diff = pos_a[:, None, :] - pos_b[None, :, :]   # (Na, Nb, 3)

            # Minimum image convention (orthogonal approximation)
            if box_vec is not None:
                diff -= np.round(diff / box_vec) * box_vec

            dists = np.sqrt(np.sum(diff ** 2, axis=2))  # (Na, Nb)
            flat_min = np.argmin(dists)
            ia_local, ib_local = np.unravel_index(flat_min, dists.shape)
            d_min = float(dists[ia_local, ib_local])

            if cutoff is not None and d_min > cutoff:
                continue

            atom_a_idx = idx_a[ia_local]
            atom_b_idx = idx_b[ib_local]
            results.append({
                'group_a': key_a,
                'group_b': key_b,
                'atom_a': atom_a_idx,
                'atom_b': atom_b_idx,
                'type_a': atoms[atom_a_idx].get('type', '?'),
                'type_b': atoms[atom_b_idx].get('type', '?'),
                'distance': round(d_min, 2),
            })

    results.sort(key=lambda r: r['distance'])
    return results[:n_pairs]

