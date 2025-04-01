import numpy as np
from .bond_angle import bond_angle
from .charge_minff import charge_minff
from .charge_formal import assign_formal_charges
from .element import element  # Correct function name is 'element' not 'set_element'
from .mass import set_atomic_masses

def minff(atoms, Box_dim, ffname='minff', rmaxlong=2.45, rmaxH=1.2):
    """Assign MINFF forcefield specific atom types to atoms based on their coordination environment.
    
    This function updates the 'fftype' field based on the atom's element and its bonding environment,
    using a two-pass approach to first determine coordination numbers and then assign types based on
    structural environment.
    
    For details, see the MINFF forcefield documentation at github.com/mholmboe/minff.
    
    Args:
        atoms: A list of atom dictionaries, each atom is expected to have position coordinates
              and element/type information.
        Box_dim: Box dimensions for periodic boundary conditions.
        ffname: The forcefield name, default is 'minff'.
        rmaxlong: Maximum bond distance for non-hydrogen bonds, default is 2.45 Å.
        rmaxH: Maximum bond distance for hydrogen bonds, default is 1.2 Å.
    
    Returns:
        The updated atoms list with 'fftype' fields assigned.
    """
    # Set the atoms chemical element names
    atoms = element(atoms)  # Use correct function name 'element'

   # First assign formal charges to all atoms (especially for ions and water)
    # This sets appropriate charges based on atom types and residue names
    atoms = assign_formal_charges(atoms)

    # Set atom masses using the mass.py module
    atoms = set_atomic_masses(atoms)
    
    # Run the entire process twice to ensure all atoms have proper typing
    # This is especially important for oxygen atoms which need to know
    # whether their metal neighbors are tetrahedral or octahedral
    for _ in range(2):  # Run the typing process twice
        # First, ensure all atoms have element types defined
        for atom in atoms:
            #if 'element' not in atom:
            # Try to extract element from atom type
            atom_type = atom.get('type', 'X')
            
            # Convert to lowercase for case-insensitive comparison
            atom_type_lower = atom_type.lower()
            
            # Map atom types to elements based on first 1-3 characters
            if atom_type_lower.startswith('si'):  
                atom['element'] = 'Si'
            elif atom_type_lower.startswith('sc'):  
                atom['element'] = 'Si'
            elif atom_type_lower.startswith('ale'): 
                atom['element'] = 'Ale'
            elif atom_type_lower.startswith('alt'): 
                atom['element'] = 'Alt'
            elif atom_type_lower.startswith('al'):  
                atom['element'] = 'Al'
            elif atom_type_lower.startswith('mg'):  
                atom['element'] = 'Mg'
            elif atom_type_lower.startswith('fee'): 
                atom['element'] = 'Fee'
            elif atom_type_lower.startswith('fet'): 
                atom['element'] = 'Fet'
            elif atom_type_lower.startswith('fe'):  
                atom['element'] = 'Fe'
            elif atom_type_lower.startswith('f'):   
                atom['element'] = 'F'
            elif atom_type_lower.startswith('li'):  
                atom['element'] = 'Li'
            elif atom_type_lower.startswith('ow'):  
                atom['element'] = 'Ow'
            elif atom_type_lower.startswith('hw'):  
                atom['element'] = 'Hw'
            elif atom_type_lower.startswith('o'):   
                atom['element'] = 'O'
            elif atom_type_lower.startswith('h'):   
                atom['element'] = 'H'
            elif atom_type_lower.startswith('ti'):   
                atom['element'] = 'Ti'
            elif atom_type_lower.startswith('ca'):   
                atom['element'] = 'Ca'
            else:
                atom['element'] = atom_type
        
        # Initialize atom types and fftypes to match element type
        for atom in atoms:
            atom['type'] = atom['element']
            atom['fftype'] = atom['element']
        
        # Only calculate bonds in the first pass
        if _ == 0:
            # Get bonds and angles using bond_angle function (this also calculates coordination numbers)
            atoms, bond_index, angle_index = bond_angle(atoms, Box_dim, rmaxH=rmaxH, rmaxM=rmaxlong)
            
            # Store bond information and prepare for atom typing
            for i, atom in enumerate(atoms):
                # Skip water and ion residues if present
                if atom.get('resname') in ['SOL', 'ION']:
                    continue
                    
                # Get neighbors from the bonds
                neighbors = atom.get('neigh', [])
                if not neighbors:
                    continue
                    
                # Use the coordination number already calculated by bond_angle
                atom['coord_num'] = atom.get('cn', 0)
                
                # For Fe atoms, calculate average Fe-O bond distance to determine oxidation state
                if atom['element'] == 'Fe':
                    # Get bonds to this atom
                    bonds = atom.get('bonds', [])
                    if bonds:
                        # Extract distances and calculate average
                        bond_distances = [dist for _, dist in bonds]
                        avg_bond_distance = sum(bond_distances) / len(bond_distances)
                        atom['avg_bond_dist'] = avg_bond_distance
        
        # Assign atom types based on coordination and bond information
        for i, atom in enumerate(atoms):
            # Skip water and ion residues
            if atom.get('resname') in ['SOL', 'ION']:
                continue
                
            # Get neighbors from the bonds
            neighbors = atom.get('neigh', [])
            if not neighbors:
                continue
                
            # Get neighbor types for pattern matching
            neighbor_types = [atoms[neigh_idx]['element'] for neigh_idx in neighbors]
            neighbor_types.sort()
            neighbors_str = ''.join(neighbor_types)
            
            # Number of neighbors (coordination number)
            coord_num = atom.get('coord_num', 0)
            
            # Determine fftype based on element and coordination environment
            el = atom['element']
            
            # Lithium assignments
            if el == 'Li':
                if coord_num == 6:
                    atom['fftype'] = 'Lio'
                elif coord_num == 4:
                    atom['fftype'] = 'Lio'
                elif coord_num > 6:
                    atom['fftype'] = 'Lio_ov'  # Over-coordinated
                elif 4 < coord_num < 6:
                    atom['fftype'] = 'Lio_un'  # Under-coordinated
            
            # Silicon assignments
            elif el == 'Si':
                o_neighbors = neighbor_types.count('O')
                if o_neighbors == 4:
                    atom['fftype'] = 'Sit'
                elif o_neighbors == 3:
                    atom['fftype'] = 'Site'  # As in Stishovite
                elif o_neighbors == 6:
                    atom['fftype'] = 'Sio'   # As in Stishovite
                elif coord_num > 4:
                    atom['fftype'] = 'Si_ov'  # Over-coordinated
                elif coord_num < 4:
                    atom['fftype'] = 'Si_un'  # Under-coordinated
            
            # Aluminum assignments
            elif el == 'Al':
                o_neighbors = neighbor_types.count('O')
                if o_neighbors == 6:
                    atom['fftype'] = 'Al'     # Octahedral Al
                elif o_neighbors == 5:
                    atom['fftype'] = 'Ale'    # 5-coordinated Al
                elif o_neighbors == 4:
                    atom['fftype'] = 'Alt'    # Tetrahedral Al
                elif coord_num > 6:
                    atom['fftype'] = 'Al_ov'  # Over-coordinated
                elif coord_num < 4:
                    atom['fftype'] = 'Al_un'  # Under-coordinated
            
            # Magnesium assignments
            elif el == 'Mg':
                if coord_num == 6:
                    # Check if there are more Mg than Si (e.g. in forsterite)
                    mg_count = sum(1 for a in atoms if a.get('element') == 'Mg')
                    si_count = sum(1 for a in atoms if a.get('element') == 'Si')
                    if mg_count > si_count:
                        atom['fftype'] = 'Mgo'  # E.g. in forsterite
                    else:
                        atom['fftype'] = 'Mg'
                elif coord_num > 6:
                    atom['fftype'] = 'Mg_ov'  # Over-coordinated
                elif coord_num < 6:
                    atom['fftype'] = 'Mg_un'  # Under-coordinated

            # Titanium assignments
            elif el == 'Ti':
                o_neighbors = neighbor_types.count('O')
                if coord_num == 6:
                    atom['fftype'] = 'Tio'  # Rutile/anatase type (TiO2)
                elif coord_num == 4:
                    atom['fftype'] = 'Tit'  # Tetrahedral Ti
                elif coord_num > 6:
                    atom['fftype'] = 'Ti_ov'  # Over-coordinated
                elif coord_num < 4:
                    atom['fftype'] = 'Ti_un'  # Under-coordinated
                    
            # Calcium assignments
            elif el == 'Ca':
                o_neighbors = neighbor_types.count('O')
                f_neighbors = neighbor_types.count('F')
                
                if o_neighbors == 6:
                    atom['fftype'] = 'Cao'  # Octahedral Ca
                elif o_neighbors == 4:
                    atom['fftype'] = 'Cah'  # 4-coordinated Ca
                elif f_neighbors == 8:
                    # Likely in fluorite (CaF2) structure
                    print(f"Ca in CaF2 Fluorite? (atom index: {atom.get('index', '?')})")
                    atom['fftype'] = 'Cah'
                elif coord_num > 6:
                    print(f"Ca atom over coordinated (atom index: {atom.get('index', '?')})")
                    print(f"Neighbors: {neighbor_types}")
                    atom['fftype'] = 'Cao_ov'  # Over-coordinated
                else:
                    # Fall back for other cases
                    print(f"Ca with unusual coordination (atom index: {atom.get('index', '?')})")
                    print(f"Neighbors: {neighbor_types}")
            
            # Iron assignments with Fe3+/Fe2+ distinction based on bond distance
            elif el == 'Fe':
                o_neighbors = neighbor_types.count('O')
                avg_bond_dist = atom.get('avg_bond_dist', 0)
                
                if o_neighbors == 6:  # Octahedral Fe
                    if avg_bond_dist < 2.07:  # Fe3+ site
                        atom['fftype'] = 'Feo'
                    else:  # Fe2+ site
                        atom['fftype'] = 'Fe2'
                elif o_neighbors == 4:  # Tetrahedral Fe
                    if avg_bond_dist < 2.0:  # Fe3+ site (typical distance cutoff for tetrahedral)
                        atom['fftype'] = 'Fet'
                    else:  # Fe2+ site
                        atom['fftype'] = 'Fe2'
                        input()
                elif coord_num > 6:
                    atom['fftype'] = 'Fe_ov'  # Over-coordinated
                elif coord_num < 4:
                    atom['fftype'] = 'Fe_un'  # Under-coordinated
            
            # Fluoride assignments
            elif el == 'F':
                if coord_num == 3:
                    atom['fftype'] = 'Fs'
                elif coord_num == 4:
                    print(f"Fs atom as in CaF2 - Fluorite? (atom index: {atom.get('index', '?')})")
                    atom['fftype'] = 'Fs'
                elif coord_num > 4:
                    print(f"Fs atom over coordinated (atom index: {atom.get('index', '?')})")
                    print(f"Neighbors: {neighbor_types}")
                    atom['fftype'] = 'Fs_ov'  # Over-coordinated
                elif coord_num > 0 and coord_num < 3:
                    print(f"Fs atom under coordinated (atom index: {atom.get('index', '?')})")
                    print(f"Neighbors: {neighbor_types}")
                    atom['fftype'] = 'Fs_un'  # Under-coordinated
                else:
                    print(f"F with unusual coordination (atom index: {atom.get('index', '?')})")
                    print(f"Neighbors: {neighbor_types}")
            
            # Hydrogen assignments
            elif el == 'H':
                if coord_num == 1:
                    atom['fftype'] = 'H'
                elif coord_num > 1:
                    atom['fftype'] = 'H_ov'  # Over-coordinated
            
            # Oxygen assignments - based on neighbor pattern
            elif el == 'O':
                # Begin with basic cases based on key neighbor patterns
                if neighbors_str == 'AlAlAl' or neighbors_str == 'AlAlAlAl':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'AlAlAlAlt':
                    atom['fftype'] = 'Obt'
                elif neighbors_str == 'AlAlAlH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str == 'AlAlAlt' or neighbors_str == 'AlAlt':
                    atom['fftype'] = 'Ops'
                elif neighbors_str == 'AlAlFe' or neighbors_str == 'AlAlFeo':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'AlAlH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str == 'AlAlSi':
                    atom['fftype'] = 'Op'
                elif neighbors_str == 'AlAlSiSi':
                    atom['fftype'] = 'Oz'
                elif neighbors_str == 'AlAleH' or neighbors_str == 'AleH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str == 'AlAleSi' or neighbors_str == 'AleSi':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'AlAltH':
                    atom['fftype'] = 'Ops'
                elif neighbors_str == 'AlFeFe' or neighbors_str == 'AlFeoFeo' or neighbors_str == 'AltFeFe' or neighbors_str == 'AltFeoFeo':
                    atom['fftype'] = 'Ops'
                elif neighbors_str == 'AlFeH' or neighbors_str == 'AlFeoH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str == 'AlFeSi' or neighbors_str == 'AlFeoSi':
                    atom['fftype'] = 'Op'
                elif neighbors_str == 'AlFet' or neighbors_str == 'AlAlFet':
                    atom['fftype'] = 'Ops'
                elif neighbors_str == 'AlH':
                    atom['fftype'] = 'Oalh'  # Al-O-H or Al-O-Si
                elif neighbors_str == 'AlHH':
                    atom['fftype'] = 'Oalhh'
                elif neighbors_str == 'AlHSi':
                    atom['fftype'] = 'Oahs'  # Al-OH-Si for acidic edge
                elif neighbors_str.startswith('AlHMg'):
                    atom['fftype'] = 'Ohmg'
                elif neighbors_str == 'AlMgSi' or neighbors_str == 'AlMgoSi':
                    atom['fftype'] = 'Omg'
                elif neighbors_str == 'AlOmg':
                    atom['fftype'] = 'Odsub'
                elif neighbors_str == 'AltH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str.startswith('AltMgh') or neighbors_str == 'AltMgoMgoMgo':
                    atom['fftype'] = 'Ops'
                elif neighbors_str == 'AltSi':
                    atom['fftype'] = 'Obs'
                # Calcium environments
                elif neighbors_str == 'CaCaCaCaCaCa' or neighbors_str == 'CaoCaoCaoCaoCaoCao':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'CaCaCaH' or neighbors_str == 'CahCahCahH' or neighbors_str == 'CaoCaoCaoH':
                    atom['fftype'] = 'Oh'
                # Iron environments
                elif neighbors_str == 'Fe2Fe2Fe2Fe2Fe2Fe2' or neighbors_str == 'FeFeFeFeFeFe':
                    atom['fftype'] = 'Op'
                elif neighbors_str == 'FeFe' or neighbors_str == 'FeoFeo':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'FeFeFe' or neighbors_str == 'FeoFeoFeo':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'FeFeFet' or neighbors_str == 'FeoFeoFet':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'FeFeFeFet' or neighbors_str == 'FeoFeoFeoFet':
                    atom['fftype'] = 'Obt'
                elif neighbors_str == 'FeFeFeFe' or neighbors_str == 'FeoFeoFeoFeo':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'FeFeFeH' or neighbors_str == 'FeoFeoFeoH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str == 'FeFeH' or neighbors_str == 'FeoFeoH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str == 'FeFeSi' or neighbors_str == 'FeoFeoSi':
                    atom['fftype'] = 'Op'
                elif neighbors_str == 'FeH' or neighbors_str == 'FeoH':
                    atom['fftype'] = 'Oh'
                elif neighbors_str.startswith('FeHMg') or neighbors_str.startswith('FeoHMg'):
                    atom['fftype'] = 'Ohmg'
                elif neighbors_str == 'FeMgSi' or neighbors_str.startswith('FeMgoSi') or neighbors_str.startswith('FeoMgSi'):
                    atom['fftype'] = 'Omg'
                elif neighbors_str == 'FeSi' or neighbors_str == 'FeoSi' or neighbors_str == 'FetSi':
                    atom['fftype'] = 'Oalt'
                elif neighbors_str == 'FetFet' or neighbors_str == 'FetFetH' or neighbors_str == 'FetH':
                    atom['fftype'] = 'Oh'
                # Lithium/Magnesium environments
                elif neighbors_str == 'HLiMgMg' or neighbors_str == 'HLiMgoMgo' or neighbors_str == 'HLioMghMgh':
                    atom['fftype'] = 'Ohli'
                elif neighbors_str == 'HHMg' or neighbors_str == 'HHMgh' or neighbors_str.startswith('HHMgo'):
                    atom['fftype'] = 'Omhh'
                elif neighbors_str == 'HMg' or neighbors_str == 'HMgh':
                    atom['fftype'] = 'Ome'
                elif neighbors_str == 'HMgMg' or neighbors_str == 'HMghMgh':
                    atom['fftype'] = 'Ohmg'
                elif neighbors_str == 'HMgMgMg' or neighbors_str == 'HMghMghMgh':
                    atom['fftype'] = 'Oh'
                elif neighbors_str == 'HSi':
                    atom['fftype'] = 'Osih'
                elif neighbors_str.startswith('LiLiLiLi'):
                    atom['fftype'] = 'Oli'
                elif neighbors_str == 'LiMgMgSi' or neighbors_str == 'LioMgMgSi' or neighbors_str == 'LioMgoMgoSi':
                    atom['fftype'] = 'Oli'
                # More magnesium environments
                elif neighbors_str == 'MgMgMgMgMgMg' or neighbors_str == 'MghMghMghMghMghMgh' or neighbors_str == 'MgoMgoMgoMgoMgoMgo':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'MgMgMgSi' or neighbors_str == 'MghMghMghSi' or neighbors_str == 'MgoMgoMgoSi':
                    atom['fftype'] = 'Op'
                elif neighbors_str == 'MgMgSi' or neighbors_str == 'MgoMgoSi':
                    atom['fftype'] = 'Odsub'
                elif neighbors_str == 'MgSi':
                    atom['fftype'] = 'Omg'
                # Silicon and Titanium environments
                elif neighbors_str == 'SiSi' or neighbors_str == 'SiSiSi':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'TiTiTi' or neighbors_str == 'TioTioTio':
                    atom['fftype'] = 'Ob'
                # Special case for AlSi with nested condition
                elif neighbors_str == 'AlSi':
                    # Check if any neighbor is Alt type
                    has_alt_neighbor = False
                    for neigh_idx in atom.get('neigh', []):
                        if atoms[neigh_idx].get('fftype') == 'Alt':
                            has_alt_neighbor = True
                            break
                    
                    if has_alt_neighbor:
                        atom['fftype'] = 'Oalt'
                    else:
                        atom['fftype'] = 'Oalsi'  # Special zeolite case
                # Water molecule 
                elif neighbors_str == 'HH':
                    atom['fftype'] = 'Ow'
                    print(f"Water molecule detected (atom index: {atom.get('index', '?')})")
                # Over/under coordination cases
                elif coord_num > 2:
                    print(f"O atom overcoordinated (atom index: {atom.get('index', '?')})")
                    print(f"Neighbors: {neighbor_types}")
                    atom['fftype'] = 'O_ov'
                elif coord_num == 1 or neighbors_str == 'AlAl' or neighbors_str == 'AlMg':
                    if neighbors_str == 'Si':
                        atom['fftype'] = 'Osi'
                    elif neighbors_str.startswith('AlAl'):
                        atom['fftype'] = 'Oal'
                    else:
                        print(f"O atom undercoordinated (atom index: {atom.get('index', '?')})")
                        print(f"Neighbors: {neighbor_types}")
                        atom['fftype'] = 'O_un'
                # Basic silicon case
                elif neighbors_str == 'Si':
                    atom['fftype'] = 'Osi'
                elif neighbors_str == 'SiH':
                    atom['fftype'] = 'Osih'
                # If nothing matched, assign a basic O type
                else:
                    atom['fftype'] = 'Osi'
        
        # Update atom types to match their new fftype
        for atom in atoms:
            atom['type'] = atom['fftype']
    
    # Apply charges based on the MINFF forcefield after atom typing is complete
    atom_labels = ['Al', 'Alt', 'Ale', 'Tio', 'Feo', 'Fet', 'Fee', 'Fe3e', 'Fe2', 'Fe2e', 
                  'Na', 'K', 'Cs', 'Mgo', 'Mgh', 'Mge', 'Cao', 'Cah', 'Sit', 'Si', 
                  'Sio', 'Site', 'Lio', 'H']
    charges = [1.782, 1.782, 1.985, 2.48, 1.5, 1.5, 1.75, 1.75, 1.184, 1.32, 
               1.0, 1.0, 1.0, 1.562, 1.74, 1.635, 1.66, 1.52, 1.884, 1.884, 
               1.884, 2.413, 0.86, 0.4]
    # By default, apply charges to all atoms (set resname=None)
    # To limit charge assignment to specific residues, provide a resname (e.g., 'MIN')
    atoms = charge_minff(atoms, Box_dim, atom_labels, charges, resname=None)
    
    # Find unique types of atomtypes and their neighbors
    # This is equivalent to the MATLAB code for finding unique types after atom type assignment
    all_neighbors = []
    
    # Gather atom types and neighbor information
    for atom in atoms:
        if atom.get('neigh', []) and 'fftype' in atom:
            atom_type = atom['fftype']
            
            # Get all neighbor atom types (not just elements)
            neighbor_types = [atoms[neigh_idx].get('fftype', atoms[neigh_idx].get('element', '')) for neigh_idx in atom.get('neigh', [])]
            
            # Create a concatenated string of all neighbor atom types in sorted order
            all_neighbor_str = ''.join(sorted(neighbor_types))
            
            # For compatibility with existing code, still compute unique elements and counts
            # (used for consolidation logic)
            neighbor_elements = [atoms[neigh_idx]['element'] for neigh_idx in atom.get('neigh', [])]
            neighbor_counts = {}
            for n_type in neighbor_elements:
                if n_type in neighbor_counts:
                    neighbor_counts[n_type] += 1
                else:
                    neighbor_counts[n_type] = 1
                    
            # Sort neighbors by element name for consistent ordering (for consolidation)
            sorted_neighbors = sorted(neighbor_counts.items())
            neighbor_str_with_counts = ''.join([f"{n_type}{count}" for n_type, count in sorted_neighbors])
            
            # Get charge if available
            charge = atom.get('charge', 0)
            
            # Combine atom type, neighbors, and charge - store both versions of neighbor info
            all_neighbors.append([atom_type, neighbor_str_with_counts, all_neighbor_str, 1, None, charge])
    
    # Consolidate duplicate entries (same atom type with same neighbor pattern)
    i = 0
    while i < len(all_neighbors) - 1:
        # If current and next row have same atom type and neighbor pattern
        if (all_neighbors[i][0] == all_neighbors[i+1][0] and 
            all_neighbors[i][1] == all_neighbors[i+1][1]):
            # Increment count in the current row
            all_neighbors[i][3] += 1
            # Remove the duplicate row
            all_neighbors.pop(i+1)
        else:
            i += 1
    
    # Note: Removed suffix numbers for duplicate atom types since we removed the Label column
    
    # Convert charges to unique values per atom type
    # Pre-compute a dictionary of charges by atom type
    charges_by_type = {}
    for row in all_neighbors:
        atom_type = row[0]
        if atom_type not in charges_by_type:
            charges_by_type[atom_type] = []
        charges_by_type[atom_type].append(row[5])  # Updated index for charges

    # Then update each row with unique charges
    for i, row in enumerate(all_neighbors):
        atom_type = row[0]
        try:
            # Try direct set conversion
            unique_charges = list(set(charges_by_type[atom_type]))
        except TypeError:
            # Fall back to string conversion if needed
            unique_charges = list(set(str(c) if isinstance(c, (list, dict, set)) else c 
                                  for c in charges_by_type[atom_type]))
        all_neighbors[i][5] = unique_charges
    
    # Create a dictionary to consolidate truly unique atom types, neighbor patterns, and their counts
    unique_patterns = {}
    for row in all_neighbors:
        atom_type = row[0]
        count = row[3]
        neighbor_pattern = row[2]  # Use neighbor atom types for display
        charges = row[5]
        
        # Create a key using atom type and neighbor pattern
        key = (atom_type, neighbor_pattern)
        
        if key in unique_patterns:
            # If this pattern already exists, add to its count
            unique_patterns[key]['count'] += count
        else:
            # First time seeing this pattern
            unique_patterns[key] = {
                'count': count,
                'charges': charges
            }
    
    # Print a compact table of truly unique atom types and their neighbors
    print("\nUnique Atom Types and Their Coordination Environment")
    print("-" * 70)
    print(f"{'Type':<10} {'Count':<6} {'Neighbors':<25} {'Charge':<15}")
    print("-" * 70)
    
    # Sort by atom type for a more organized display
    for key in sorted(unique_patterns.keys()):
        atom_type, neighbor_pattern = key
        count = unique_patterns[key]['count']
        charges = unique_patterns[key]['charges']
        charge_str = ', '.join([f"{c:.3f}" if isinstance(c, float) else str(c) for c in charges])
        print(f"{atom_type:<10} {count:<6} {neighbor_pattern:<25} {charge_str:<15}")
    print("-" * 70)
    
    return atoms, all_neighbors
