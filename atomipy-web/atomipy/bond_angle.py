import numpy as np
from .dist_matrix import dist_matrix
from .cell_list_dist_matrix import cell_list_dist_matrix


def bond_angle(atoms, Box_dim, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True, calculate_coordination=True, neighbor_element=None):
    """Compute bonds and angles for a given atomic structure.
    
    For each atom, bonds are determined based on a distance threshold:
      - rmaxH (default 1.2 Å) if either atom is hydrogen
      - rmaxM (default 2.45 Å) for bonds between non-hydrogen atoms.
    
    Angles are then computed for each pair of bonds at the central atom using the periodic boundary condition (PBC)
    corrected vectors. The function updates each atom's 'neigh', 'bonds', and 'angles' fields in-place.
    The same cutoffs are applied to both the neighbor list and bond list.

    Args:
       atoms: list of atom dictionaries (coordinates in Angstroms).
       Box_dim: 1x9 list representing triclinic cell dimensions (in Angstroms).
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
        # Large system - use cell list method which is more memory efficient
        print(f"Large system - using cell list method for the distance matrix")
        dmat, dx, dy, dz, _, _ = cell_list_dist_matrix(atoms, cutoff=max(rmaxH, rmaxM), Box_dim=Box_dim)
    else:
        # Smaller system - use standard distance matrix which is faster
        print(f"Small system - calculating the full distance matrix")
        dmat, dx, dy, dz = dist_matrix(atoms, Box_dim=Box_dim)
    
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

