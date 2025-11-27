import numpy as np
from .dist_matrix import dist_matrix
from .cell_list_dist_matrix import cell_list_dist_matrix

def bond_angle(atoms, Box, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True, calculate_coordination=True, neighbor_element=None):
    """Compute bonds and angles for a given atomic structure.
    
    
    For each atom, bonds are determined based on a distance threshold:
      - rmaxH (default 1.2 Å) if either atom is hydrogen
      - rmaxM (default 2.45 Å) for bonds between non-hydrogen atoms.
    
    Angles are then computed for each pair of bonds at the central atom using the periodic boundary condition (PBC)
    corrected vectors. The function updates each atom's 'neigh', 'bonds', and 'angles' fields in-place.
    The same cutoffs are applied to both the neighbor list and bond list.

    Args:
       atoms: list of atom dictionaries (coordinates in Angstroms).
       Box: a 1x3, 1x6 or 1x9 list representing Cell dimensions (in Angstroms):
           - For orthogonal boxes, a 1x3 list [lx, ly, lz] where Box = Box_dim, and Cell would be [lx, ly, lz, 90, 90, 90]
           - For Cell parameters, a 1x6 list [a, b, c, alpha, beta, gamma] (Cell format)
           - For triclinic boxes, a 1x9 list [lx, ly, lz, 0, 0, xy, 0, xz, yz] (GROMACS Box_dim format)
       rmaxH: cutoff distance for bonds involving hydrogen (default 1.2 Å).
       rmaxM: cutoff distance for bonds between all other atoms (default 2.45 Å).
       same_element_bonds: if False, bonds between atoms of the same element are ignored (default True).
       same_molecule_only: if True, bonds are only formed between atoms with the same 'molid' (default False).

    Returns:
       tuple: (atoms, Bond_index, Angle_index)
           - atoms: Updated atom list with 'neigh', 'bonds', and 'angles'.
           - Bond_index: Nx3 numpy array where each row contains [atom1_idx, atom2_idx, distance]
             with atom indices sorted from low to high.
           - Angle_index: Mx10 numpy array with the following columns:
             [atom1_idx, atom2_idx, atom3_idx, angle, dx12, dy12, dz12, dx23, dy23, dz23]
             where atom2_idx is the middle/center atom of the angle, atom1_idx is the bonded atom
             with the lowest index, and atom3_idx is the bonded atom with the highest index.
             The dx,dy,dz values represent the distance vector components between the respective atoms.
        calculate_coordination: If True, calculate coordination numbers for each atom and store in 'cn' field.
        neighbor_element: Optional filter to only count neighbors of a specific element when calculating
                         coordination numbers.
    """
    # Get the number of atoms
    N = len(atoms)

    # Check the size of the system to determine which method to use
    # For large systems (>20000 atoms), cell_list_dist_matrix is more memory efficient
    # For smaller systems, dist_matrix is faster
    if len(atoms) > 15000:
        # Large system - use Cell list method which is more memory efficient
        print(f"Large system - using Cell list method for the distance matrix")
        dmat, dx, dy, dz, _, _ = cell_list_dist_matrix(atoms, Box, cutoff=max(rmaxH, rmaxM))
    else:
        # Smaller system - use standard distance matrix which is faster
        print(f"Small system - calculating the full distance matrix")
        dmat, dx, dy, dz = dist_matrix(atoms, Box)
    
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
    
    # Use tqdm if available, otherwise use a basic counter with percentage updates
    if has_tqdm:
        iterator = tqdm(range(N), desc="Finding bonds", unit="atom")
    else:
        print("Finding bonds...")
        iterator = range(N)
        last_percent = -1
    
    for i in iterator:
        for j in range(i+1, N):  # Only consider each pair once
            pair_count += 1
            
            # If no tqdm, show percentage updates
            if not has_tqdm and N > 1000:
                percent = int(100 * pair_count / total_pairs)
                if percent > last_percent and percent % 10 == 0:
                    print(f"  {percent}% complete...")
                    last_percent = percent
                    
            if dmat[i, j] > 0:  # Skip diagonal and zero distances
                # Determine if either atom is hydrogen by checking if type starts with H/h
                type_i = atoms[i].get('type', atoms[i].get('name', ''))
                type_j = atoms[j].get('type', atoms[j].get('name', ''))
                
                # Get first character of type and check if it's 'H' or 'h'
                isH_i = type_i and type_i[0].upper() == 'H'
                isH_j = type_j and type_j[0].upper() == 'H'
                
                # Apply appropriate cutoff based on atom types
                if isH_i or isH_j:
                    cutoff = rmaxH  # Use hydrogen cutoff
                else:
                    cutoff = rmaxM  # Use non-hydrogen cutoff
                    
                # If within cutoff, add to bond list
                if dmat[i, j] <= cutoff:
                    precalc_bond_list.append([i, j])
                    dist_list.append(dmat[i, j])

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
