import numpy as np
from .bond_angle import bond_angle
from .charge_minff import charge_minff
from .charge_formal import assign_formal_charges

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
    
    # Run the entire process twice to ensure all atoms have proper typing
    # This is especially important for oxygen atoms which need to know
    # whether their metal neighbors are tetrahedral or octahedral
    for _ in range(2):  # Run the typing process twice
        # First, ensure all atoms have element types defined
        for atom in atoms:
            if 'element' not in atom:
                # Try to extract element from atom type
                atom_type = atom.get('type', 'X')
                
                # Map atom types to elements based on first 1-3 characters
                if atom_type.startswith('Si'):  
                    atom['element'] = 'Si'
                elif atom_type.startswith('SC'):  
                    atom['element'] = 'Si'
                elif atom_type.startswith('Ale'): 
                    atom['element'] = 'Ale'
                elif atom_type.startswith('Alt'): 
                    atom['element'] = 'Alt'
                elif atom_type.startswith('Al'):  
                    atom['element'] = 'Al'
                elif atom_type.startswith('Mg'):  
                    atom['element'] = 'Mg'
                elif atom_type.startswith('Fee'): 
                    atom['element'] = 'Fee'
                elif atom_type.startswith('Fet'): 
                    atom['element'] = 'Fet'
                elif atom_type.startswith('Fe'):  
                    atom['element'] = 'Fe'
                elif atom_type.startswith('F'):   
                    atom['element'] = 'F'
                elif atom_type.startswith('Li'):  
                    atom['element'] = 'Li'
                elif atom_type.startswith('Ow'):  
                    atom['element'] = 'Ow'
                elif atom_type.startswith('Hw'):  
                    atom['element'] = 'Hw'
                elif atom_type.startswith('O'):   
                    atom['element'] = 'O'
                elif atom_type.startswith('H'):   
                    atom['element'] = 'H'
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
            
            # Iron assignments with Fe3+/Fe2+ distinction based on bond distance
            elif el == 'Fe':
                o_neighbors = neighbor_types.count('O')
                avg_bond_dist = atom.get('avg_bond_dist', 0)
                
                if o_neighbors == 6:  # Octahedral Fe
                    if avg_bond_dist < 2.07:  # Fe3+ site
                        atom['fftype'] = 'Feo3'
                    else:  # Fe2+ site
                        atom['fftype'] = 'Feo2'
                elif o_neighbors == 4:  # Tetrahedral Fe
                    if avg_bond_dist < 2.0:  # Fe3+ site (typical distance cutoff for tetrahedral)
                        atom['fftype'] = 'Fet3'
                    else:  # Fe2+ site
                        atom['fftype'] = 'Fet2'
                elif coord_num > 6:
                    atom['fftype'] = 'Fe_ov'  # Over-coordinated
                elif coord_num < 4:
                    atom['fftype'] = 'Fe_un'  # Under-coordinated
            
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
                    atom['fftype'] = 'Obh'
                elif neighbors_str == 'AlAlH':
                    atom['fftype'] = 'Oalh'
                elif neighbors_str == 'AlAlHH':
                    atom['fftype'] = 'Oalhh'
                elif neighbors_str == 'AlH':
                    atom['fftype'] = 'Oalh'
                elif neighbors_str == 'AlSi':
                    atom['fftype'] = 'Oas'
                elif neighbors_str == 'AlAlSi':
                    atom['fftype'] = 'Oas'
                elif neighbors_str == 'SiSi':
                    atom['fftype'] = 'Ob'
                elif neighbors_str == 'SiSiSi':
                    atom['fftype'] = 'Obt'
                elif neighbors_str == 'Si':
                    atom['fftype'] = 'Osi'
                elif neighbors_str == 'SiH':
                    atom['fftype'] = 'Osih'
                # If nothing matched, assign a basic O type
                elif atom['fftype'] == 'O':
                    atom['fftype'] = 'Osi'
        
        # Update atom types to match their new fftype
        for atom in atoms:
            atom['type'] = atom['fftype']

    # First assign formal charges to all atoms (especially for ions and water)
    # This sets appropriate charges based on atom types and residue names
    atoms = assign_formal_charges(atoms)
    
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
    
    return atoms
