def element_atom(atom):
    """Guess the chemical element for an atom entry.
    
    It uses the 'resname' field and known element symbols to assign the 'element' key.
    Expected elements: Si, Al, Fe, Mg, Ti, Li, F, O, H.
    
    Args:
       atom: A dictionary representing an atom, with at least a 'resname' key.

    Returns:
       The element symbol as a string, and also updates the atom's 'element' field.
    """
    valid_elements = ['Si', 'Al', 'Fe', 'Mg', 'Ti', 'Li', 'F', 'O', 'H']
    resname = atom.get('resname', '').strip()
    # Check if resname itself matches any expected element
    if resname in valid_elements:
        atom['element'] = resname
    else:
        # Otherwise, try to use the first two characters (capitalized) if they form a known element
        candidate = resname[:2].title()
        if candidate in valid_elements:
            atom['element'] = candidate
        else:
            # Fallback: use first character capitalized
            candidate = resname[:1].upper()
            if candidate in valid_elements:
                atom['element'] = candidate
            else:
                # As last resort, leave as unknown
                atom['element'] = 'X'
    return atom
