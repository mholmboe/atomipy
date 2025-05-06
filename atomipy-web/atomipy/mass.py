def mass():
    """Return dictionary mapping elements to their atomic masses in atomic mass units (amu).
    Example values provided.
    """
    masses = {
        'Li': 6.941,
        'Na': 22.9898,
        'Mg': 24.305,
        'Si': 28.085,
        'Al': 26.982,
        'K': 39.102,
        'Ca': 40.078,
        'Fe': 55.845,
        'Ti': 47.867,
        'F': 18.9984,
        'O': 15.9994,
        'H': 1.00784
    }
    return masses


def set_atomic_masses(atoms):
    """Set the mass attribute for each atom in the atoms list based on its element.
    
    Parameters
    ----------
    atoms : list of dictionaries
        List of atom dictionaries. If 'element' key is missing, it will be determined using element.py.
        
    Returns
    -------
    atoms : list of dictionaries
        The same list of atoms with updated 'mass' values.
    """
    # Import element module for determining element types if needed
    from . import element as element_module
    
    # Get the mass dictionary
    mass_dict = mass()
    
    # Count how many masses were set
    mass_set_count = 0
    
    # First, ensure all atoms have element information
    atoms_with_missing_elements = [atom for atom in atoms if 'element' not in atom or atom['element'] is None]
    if atoms_with_missing_elements:
        print(f"Determining element types for {len(atoms_with_missing_elements)} atoms without element information...")
        element_module.element(atoms)
    
    # Iterate through atoms and set masses
    for atom in atoms:
        # Get the element 
        element = atom.get('element')
        if element is None:
            print(f"Warning: Atom {atom.get('id', '?')} has no element information even after element assignment")
            continue
            
        # Handle special element cases from element.py
        if element == 'Ale' or element == 'Alt':
            element = 'Al'
        elif element == 'Fee' or element == 'Fet':
            element = 'Fe'
        elif element == 'Ow':
            element = 'O'
        elif element == 'Hw':
            element = 'H'
            
        # Get the first 1-2 characters which represent the element symbol
        element_symbol = element[:1] if len(element) == 1 else element[:2]
        
        # Try to find the element in the mass dictionary (case insensitive)
        element_upper = element_symbol.upper()
        element_title = element_symbol.title()
        
        # Try different ways to match the element
        if element_upper in mass_dict:
            atom['mass'] = mass_dict[element_upper]
            mass_set_count += 1
        elif element_title in mass_dict:
            atom['mass'] = mass_dict[element_title]
            mass_set_count += 1
        else:
            # If we couldn't find a mass for this element, print a warning
            print(f"Warning: No mass found for element '{element}' in atom {atom.get('id', '?')}")
    
    # Summary statistics
    if mass_set_count == len(atoms):
        print(f"✓ Successfully set atomic masses for all {mass_set_count} atoms")
    else:
        print(f"Set atomic masses for {mass_set_count} out of {len(atoms)} atoms")
        
    return atoms
