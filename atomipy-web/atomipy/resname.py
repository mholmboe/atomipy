"""
This module provides functions for assigning residue names to atoms.

The main function guesses and assigns residue names (resname) based on atom types,
with special handling for water and ions.
"""

def assign_resname(atoms, default_resname='MIN'):
    """
    Assign residue names to atoms based on their types.
    
    This function assigns 'SOL' to water atoms, 'ION' to ion atoms,
    and the specified default residue name to all other atoms.
    
    Args:
        atoms: List of atom dictionaries
        default_resname: The residue name to assign to non-water, non-ion atoms (default: 'MIN')
    
    Returns:
        The updated list of atom dictionaries with 'resname' field assigned
    
    Example:
        # Assign residue names to all atoms
        atoms = ap.resname.assign_resname(atoms)
        
        # Assign 'CLAY' as the default residue name
        atoms = ap.resname.assign_resname(atoms, default_resname='CLAY')
    """
    # Define atom type patterns for water and ions
    water_types = ['Hw', 'Ow', 'OW', 'HW']
    
    ion_types = [
        # Monovalent cations
        'Li', 'Li+', 'LI+', 'Na', 'NA', 'Na+', 'NA+', 'K', 'K+', 
        'Rb', 'RB', 'Rb+', 'RB+', 'Cs', 'CS', 'Cs+', 'CS+',
        
        # Divalent cations
        'Ca', 'CA', 'Ca2+', 'CA2+', 'Cu', 'CU', 'Cu2+', 'CU2+',
        'Ni', 'NI', 'Ni2+', 'NI2+', 'Zn', 'ZN', 'Zn2+', 'ZN2+',
        'Sr', 'SR', 'Sr2+', 'SR2+', 'Ba', 'BA', 'Ba2+', 'BA2+',
        
        # Anions
        'F-', 'Cl', 'CL', 'Cl-', 'CL-', 'Br', 'BR', 'Br-', 'BR-', 'I', 'I-'
    ]
    
    # Create list of atom indices for each category
    sol_indices = []
    ion_indices = []
    
    # Process each atom to identify water and ions
    for i, atom in enumerate(atoms):
        atom_type = atom.get('type', '')
        
        # Check for water atoms (case-insensitive prefix match)
        is_water = any(atom_type.lower().startswith(water_type.lower()) for water_type in water_types)
        
        # Check for ion atoms (case-insensitive exact match)
        is_ion = atom_type in ion_types or atom_type.upper() in ion_types
        
        if is_water:
            sol_indices.append(i)
        elif is_ion:
            ion_indices.append(i)
    
    # Sort indices
    sol_indices.sort()
    ion_indices.sort()
    
    # Assign 'SOL' to water atoms
    for i in sol_indices:
        atoms[i]['resname'] = 'SOL'
    
    # Assign 'ION' to ion atoms
    for i in ion_indices:
        atoms[i]['resname'] = 'ION'
    
    # Assign default resname to remaining atoms
    other_indices = [i for i in range(len(atoms)) if i not in sol_indices and i not in ion_indices]
    for i in other_indices:
        # Only assign if resname is not already set
        if 'resname' not in atoms[i] or not atoms[i]['resname']:
            atoms[i]['resname'] = default_resname
    
    # Print summary
    print(f"Assigned resnames: {len(sol_indices)} water (SOL), {len(ion_indices)} ions (ION), "
          f"{len(other_indices)} other atoms ({default_resname})")
    
    return atoms


def change_default_resname(atoms, new_resname, current_resname='MIN'):
    """
    Change the residue name for atoms with a specific current residue name.
    
    This is useful to reclassify groups of atoms without affecting water and ions.
    
    Args:
        atoms: List of atom dictionaries
        new_resname: The new residue name to assign
        current_resname: The current residue name to change (default: 'MIN')
    
    Returns:
        The updated list of atom dictionaries
    
    Example:
        # Change all 'MIN' residues to 'CLAY'
        atoms = ap.resname.change_default_resname(atoms, 'CLAY')
    """
    count = 0
    for atom in atoms:
        if atom.get('resname') == current_resname:
            atom['resname'] = new_resname
            count += 1
    
    print(f"Changed {count} atoms from resname '{current_resname}' to '{new_resname}'")
    return atoms
