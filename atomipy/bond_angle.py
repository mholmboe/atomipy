import numpy as np
from .dist_matrix import dist_matrix
from .cell_list_dist_matrix import cell_list_dist_matrix
from .cell_list_dist_matrix_fast import cell_list_dist_matrix_fast

def bond_angle(atoms, Box, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True, calculate_coordination=True, neighbor_element=None, dm_method=None):
    """Compute bonds and angles for a given atomic structure.
    
    ... [docs] ...
    
        dm_method: Optional manual selection of distance matrix method: 
                   'direct', 'old_cl', or 'fast_cl'. If None, auto-select based on size.
    """
    # Get the number of atoms
    N = len(atoms)

    # Check the size of the system to determine which method to use
    # For large systems (>2000 atoms), cell_list_dist_matrix_fast is more memory efficient and now highly vectorized
    if dm_method == 'direct' or (dm_method is None and len(atoms) <= 2000):
        # Smaller system (or forced direct) - use standard distance matrix
        print(f"Using Direct distance matrix method (size: {len(atoms)})")
        dmat, dx, dy, dz = dist_matrix(atoms, Box)
    elif dm_method == 'old_cl':
        # Forced old cell list
        print(f"Using Old Cell List method (size: {len(atoms)})")
        dmat, dx, dy, dz, _, _ = cell_list_dist_matrix(atoms, Box, cutoff=max(rmaxH, rmaxM))
    else:
        # Default or forced fast_cl
        print(f"Using Memory-Efficient Fast Cell List method (size: {len(atoms)})")
        dmat, dx, dy, dz, _, _ = cell_list_dist_matrix_fast(atoms, Box, cutoff=max(rmaxH, rmaxM), rmaxH=rmaxH)
    
    # Since dist_matrix doesn't provide precalculated bond lists, we'll create them based on cutoffs
    precalc_bond_list = []
    dist_list = []
    
    # Try to import tqdm for progress bar
    try:
        from tqdm import tqdm
        has_tqdm = True
    except ImportError:
        print("Note: Install tqdm package for progress bars (pip install tqdm)")
        has_tqdm = False

    # Create bond lists based on the distance matrix and appropriate cutoffs
    total_pairs = N * (N - 1) // 2  # Total number of pairs to check
    pair_count = 0
    
    # Vectorized bond identification
    # Determine which atoms are hydrogen
    types = np.array([atom.get('type', atom.get('name', '')) for atom in atoms])
    is_h = np.array([bool(t and t[0].upper() == 'H') for t in types])
    
    # Create a cutoff matrix: rmaxH if either atom is H, otherwise rmaxM
    # is_h[:, None] | is_h[None, :] is True if i is H OR j is H
    cutoff_matrix = np.where(is_h[:, np.newaxis] | is_h[np.newaxis, :], rmaxH, rmaxM)
    
    # Find indices where distance is within cutoff and i < j (upper triangle)
    # Using np.triu to only get each pair once
    mask = (dmat > 0) & (dmat <= cutoff_matrix)
    ii, jj = np.where(np.triu(mask, k=1))
    
    # Create the bond list and distance list
    precalc_bond_list = np.column_stack((ii, jj))
    dist_list = dmat[ii, jj]
    
    print(f"  Vectorized bond finding: identified {len(precalc_bond_list)} potential bonds")

    # Initialize lists for all atoms
    for i in range(N):
        atoms[i]['neigh'] = []
        atoms[i]['bonds'] = []
        atoms[i]['angles'] = []
    
    # Filter bonds based on element types and molecule IDs
    bond_pairs = []  # Store bonds as (atom1_idx, atom2_idx, distance)
    
    # Process the precalculated bonds from cell_list_dist_matrix
    if len(precalc_bond_list) > 0:
        for k in range(len(precalc_bond_list)):
            i, j = precalc_bond_list[k]
            distance = dist_list[k]
            
            # Ensure i < j for consistency - smaller index always in first column
            if i > j:
                i, j = j, i
                
            el_i = atoms[i].get('element','X')
            el_j = atoms[j].get('element','X')
            
            # Get molecule IDs if available, otherwise use None
            molid_i = atoms[i].get('molid', None)
            molid_j = atoms[j].get('molid', None)
            
            # Apply element check and molecule check if needed
            molecule_condition = True if not same_molecule_only else (molid_i == molid_j)
            element_condition = same_element_bonds or el_i != el_j
            
            if element_condition and molecule_condition:
                # Add to both atoms' neighbor and bond lists
                atoms[i]['neigh'].append(j)
                atoms[i]['bonds'].append((j, distance))
                
                atoms[j]['neigh'].append(i)
                atoms[j]['bonds'].append((i, distance))
                
                # Store bond information as tuple (low_idx, high_idx, distance)
                bond_pairs.append((i, j, distance))
    
    # Calculate angles for atoms with bonds
    angle_data = []  # Store angle data
    
    for i in range(N):
        # Skip if atom has less than 2 bonds
        if len(atoms[i]['neigh']) < 2:
            continue
            
        # Compute angles for each pair of bonded neighbors
        for m in range(len(atoms[i]['neigh'])):
            for n in range(m+1, len(atoms[i]['neigh'])):
                j = atoms[i]['neigh'][m]
                k = atoms[i]['neigh'][n]
                
                # Get vectors from atom i to atoms j and k with PBC correction
                rij = np.array([dx[i, j], dy[i, j], dz[i, j]])
                rik = np.array([dx[i, k], dy[i, k], dz[i, k]])
                
                # Normalize vectors
                rij_norm = np.linalg.norm(rij)
                rik_norm = np.linalg.norm(rik)
                
                # Calculate angle using dot product
                cos_angle = np.dot(rij, rik) / (rij_norm * rik_norm)
                
                # Clamp to valid range to prevent numerical errors
                cos_angle = max(min(cos_angle, 1.0), -1.0)
                angle = np.degrees(np.arccos(cos_angle))
                
                # Add angle to atom's data
                atoms[i]['angles'].append(((j, k), angle))
                
                # Store angle data with proper ordering for Angle_index
                # Ensure first atom has lower index than third atom
                if j < k:
                    atom1, atom3 = j, k
                    # Vector from middle atom (i) to lowest index atom (j)
                    dx12, dy12, dz12 = dx[i, j], dy[i, j], dz[i, j]
                    # Vector from middle atom (i) to highest index atom (k)
                    dx23, dy23, dz23 = dx[i, k], dy[i, k], dz[i, k]
                else:
                    atom1, atom3 = k, j
                    # Vector from middle atom (i) to lowest index atom (k)
                    dx12, dy12, dz12 = dx[i, k], dy[i, k], dz[i, k]
                    # Vector from middle atom (i) to highest index atom (j)
                    dx23, dy23, dz23 = dx[i, j], dy[i, j], dz[i, j]
                
                # Store the angle data in consistent format
                angle_data.append((atom1, i, atom3, angle, 
                                 dx12, dy12, dz12, 
                                 dx23, dy23, dz23))
    
    # Convert bond_pairs list to Nx3 numpy array
    Bond_index = np.array(bond_pairs)
    
    # For each bond, ensure the smaller atom index is in the first column
    if len(Bond_index) > 0:
        # This should already be taken care of when creating bond_pairs,
        # but let's make sure by doing a final check
        for i in range(len(Bond_index)):
            if Bond_index[i, 0] > Bond_index[i, 1]:
                # Swap indices to put smaller first
                Bond_index[i, 0], Bond_index[i, 1] = Bond_index[i, 1], Bond_index[i, 0]
        
        # Now sort rows based on first column (atom1_idx) and then second column (atom2_idx)
        sorted_indices = np.lexsort((Bond_index[:, 1], Bond_index[:, 0]))
        Bond_index = Bond_index[sorted_indices]
    
    # Convert angle_data list to Mx10 numpy array
    Angle_index = np.array(angle_data)
    
    # Sort Angle_index row-wise by atomic indices
    if len(Angle_index) > 0:
        # Sort the entire array: first by middle atom (column 1), 
        # then by lowest bonded atom (column 0), then by highest bonded atom (column 2)
        sorted_indices = np.lexsort((Angle_index[:, 2], Angle_index[:, 0], Angle_index[:, 1]))
        Angle_index = Angle_index[sorted_indices]
    
    # Calculate coordination numbers if requested
    if calculate_coordination:
        for i, atom in enumerate(atoms):
            neighbors = atom.get('neigh', [])
            
            # Filter neighbors by element if requested
            if neighbor_element and neighbors:
                neighbors = [idx for idx in neighbors if atoms[idx].get('element') == neighbor_element]
                
            atom['cn'] = len(neighbors)
    
    return atoms, Bond_index, Angle_index


def bond_angle_dihedral(atoms, Box, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False,
                        same_molecule_only=True, calculate_coordination=True, neighbor_element=None):
    """
    Calculate bonds, angles, dihedrals, and 1-4 pair interactions.

    This mirrors the MATLAB bond_angle_dihedral_atom.m flow: bonds/angles are built first,
    then dihedrals are assembled from adjacency and refined using the angle vectors
    (including reversed orientations) to compute dihedral angles. The pairlist is derived
    from terminal atoms of the dihedrals with existing bonds removed.

    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries containing Cartesian coordinates (keys: ``x``, ``y``, ``z``).
    Box : list or array-like
        Simulation cell dimensions in one of the supported formats:
        - 1x3: [lx, ly, lz] for orthogonal boxes
        - 1x6: [a, b, c, alpha, beta, gamma] for cell parameters
        - 1x9: [lx, ly, lz, 0, 0, xy, 0, xz, yz] for GROMACS triclinic format
    rmaxH : float, optional
        Cutoff distance (Å) for bonds involving hydrogen (default: 1.2).
    rmaxM : float, optional
        Cutoff distance (Å) for bonds between non-hydrogen atoms (default: 2.45).
    same_element_bonds : bool, optional
        If False, bonds between atoms of the same element are ignored (default: False).
    same_molecule_only : bool, optional
        If True, restrict bonds to atoms sharing the same ``molid`` (default: True).
    calculate_coordination : bool, optional
        If True, add coordination numbers to each atom dictionary (default: True).
    neighbor_element : str, optional
        If provided, coordination counts only neighbors of this element.

    Returns
    -------
    atoms : list of dict
        Atom list updated with neighbor, bond, angle, and coordination information.
    Bond_index : numpy.ndarray
        Array with shape (n_bonds, 3) containing bonded atom indices (0-based) and distances.
    Angle_index : numpy.ndarray
        Array with shape (n_angles, 10) describing angles as returned by :func:`bond_angle`.
    Dihedral_index : numpy.ndarray
        Array with shape (n_dihedrals, 5): [atom1, atom2, atom3, atom4, dihedral_angle_deg].
    Pairlist : numpy.ndarray
        Array with shape (n_pairs, 2) listing unique 1-4 pairs not directly bonded.

    Raises
    ------
    ValueError
        If ``atoms`` is empty or ``Box`` is missing/invalid.
    """
    if atoms is None or len(atoms) == 0:
        raise ValueError("atoms list cannot be empty")

    if Box is None:
        raise ValueError("Box parameter must be provided")

    atoms, Bond_index, Angle_index = bond_angle(
        atoms,
        Box,
        rmaxH=rmaxH,
        rmaxM=rmaxM,
        same_element_bonds=same_element_bonds,
        same_molecule_only=same_molecule_only,
        calculate_coordination=calculate_coordination,
        neighbor_element=neighbor_element,
    )

    Bond_index = np.atleast_2d(np.array(Bond_index))
    Angle_index = np.atleast_2d(np.array(Angle_index))

    if Bond_index.size == 0 or Angle_index.size == 0:
        empty_dihedral = np.empty((0, 5))
        empty_pairlist = np.empty((0, 2), dtype=int)
        return atoms, Bond_index, Angle_index, empty_dihedral, empty_pairlist

    # Build adjacency from bonds
    nbonds = Bond_index.shape[0]
    adjacency = [[] for _ in range(len(atoms))]
    for i in range(nbonds):
        a = int(Bond_index[i, 0])
        b = int(Bond_index[i, 1])
        adjacency[a].append(b)
        adjacency[b].append(a)

    # Oriented angle list using full neighbor permutations (not sorted) for dihedral construction
    from .dist_matrix import dist_matrix  # local import to avoid cycles
    _, dx, dy, dz = dist_matrix(atoms, Box)
    oriented_angles = []
    for center, neighs in enumerate(adjacency):
        if len(neighs) < 2:
            continue
        for n1 in neighs:
            for n2 in neighs:
                if n1 == n2:
                    continue
                v1 = np.array([dx[center, n1], dy[center, n1], dz[center, n1]])
                v2 = np.array([dx[center, n2], dy[center, n2], dz[center, n2]])
                norm1 = np.linalg.norm(v1)
                norm2 = np.linalg.norm(v2)
                if norm1 == 0 or norm2 == 0:
                    continue
                cos_angle = np.dot(v1, v2) / (norm1 * norm2)
                cos_angle = max(min(cos_angle, 1.0), -1.0)
                angle = np.degrees(np.arccos(cos_angle))
                oriented_angles.append(
                    (
                        n1,
                        center,
                        n2,
                        angle,
                        v1[0],
                        v1[1],
                        v1[2],
                        v2[0],
                        v2[1],
                        v2[2],
                    )
                )
    oriented_angles = np.array(oriented_angles)
    if oriented_angles.size == 0:
        empty_dihedral = np.empty((0, 5))
        empty_pairlist = np.empty((0, 2), dtype=int)
        return atoms, Bond_index, Angle_index, empty_dihedral, empty_pairlist
    angle_indices = oriented_angles[:, :3].astype(int)

    # Prepare angle stack (forward + reversed) to compute dihedral angles
    angles_rev = np.zeros_like(oriented_angles)
    angles_rev[:, 0] = angle_indices[:, 2]
    angles_rev[:, 1] = angle_indices[:, 1]
    angles_rev[:, 2] = angle_indices[:, 0]
    angles_rev[:, 3] = oriented_angles[:, 3]
    angles_rev[:, 4:7] = oriented_angles[:, 7:10]
    angles_rev[:, 7:10] = oriented_angles[:, 4:7]

    angle_stack = np.vstack((angles_rev, oriented_angles))
    angle_stack_indices = angle_stack[:, :3].astype(int)

    dihedral_list = []
    for i in range(angle_stack.shape[0]):
        for j in range(i, angle_stack.shape[0]):
            if (angle_stack_indices[i, 1] == angle_stack_indices[j, 0] and
                    angle_stack_indices[i, 2] == angle_stack_indices[j, 1]):
                A = np.cross(angle_stack[i, 4:7], angle_stack[i, 7:10])
                B = np.cross(angle_stack[j, 4:7], angle_stack[j, 7:10])
                normA = np.linalg.norm(A)
                normB = np.linalg.norm(B)
                if normA < 1e-8 or normB < 1e-8:
                    continue
                cos_phi = np.dot(A, B) / (normA * normB)
                cos_phi = max(min(cos_phi, 1.0), -1.0)
                theta = float(np.round(np.degrees(np.arccos(cos_phi)), 2))

                if angle_stack_indices[i, 1] < angle_stack_indices[i, 2]:
                    dihedral = (
                        angle_stack_indices[i, 0],
                        angle_stack_indices[i, 1],
                        angle_stack_indices[i, 2],
                        angle_stack_indices[j, 2],
                        theta,
                    )
                else:
                    dihedral = (
                        angle_stack_indices[j, 2],
                        angle_stack_indices[i, 2],
                        angle_stack_indices[i, 1],
                        angle_stack_indices[i, 0],
                        theta,
                    )
                dihedral_list.append(dihedral)

    if dihedral_list:
        Dihedral_index = np.unique(np.array(dihedral_list, dtype=float), axis=0)
        sort_idx = np.lexsort((Dihedral_index[:, 2], Dihedral_index[:, 1]))
        Dihedral_index = Dihedral_index[sort_idx]
    else:
        Dihedral_index = np.empty((0, 5))

    # Pairlist from dihedral terminals, removing bonded pairs
    if len(Dihedral_index) > 0:
        Pairlist = np.sort(Dihedral_index[:, [0, 3]].astype(int), axis=1)
        Pairlist = np.unique(Pairlist, axis=0)
        if Bond_index.size > 0:
            bond_pairs = np.sort(Bond_index[:, :2].astype(int), axis=1)
            keep_mask = []
            for pair in Pairlist:
                if np.any(np.all(bond_pairs == np.sort(pair), axis=1)):
                    keep_mask.append(False)
                else:
                    keep_mask.append(True)
            Pairlist = Pairlist[np.array(keep_mask, dtype=bool)]
    else:
        Pairlist = np.empty((0, 2), dtype=int)

    return atoms, Bond_index, Angle_index, Dihedral_index, Pairlist
