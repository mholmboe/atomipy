def vdw_radius():
    """Return dictionary mapping elements to their van der Waals radii.
    Values in angstroms.
    """
    radii = {
        'Si': 2.10,
        'Al': 1.84,
        'Fe': 2.10,
        'Mg': 1.73,
        'Ti': 2.00,
        'Li': 1.82,
        'F': 1.47,
        'O': 1.52,
        'H': 1.20
    }
    return radii


def ionic_radius():
    """Return dictionary mapping elements to their ionic radii.
    Values in angstroms for the common oxidation states.
    """
    radii = {
        'Si': 0.40,  # Si4+
        'Al': 0.53,  # Al3+
        'Fe': 0.78,  # Fe2+ (high spin)
        'Mg': 0.72,  # Mg2+
        'Ti': 0.61,  # Ti4+
        'Li': 0.76,  # Li+
        'F': 1.33,   # F-
        'O': 1.40,   # O2-
        'H': 0.25    # Estimate for H+ (often used as 0 since it's a bare proton)
    }
    return radii


def radius(radius_type='vdw'):
    """Return dictionary mapping elements to their radii based on the specified type.
    
    Args:
        radius_type: String, either 'vdw' for van der Waals radii or 'ionic' for ionic radii.
        
    Returns:
        Dictionary mapping element symbols to their radii in angstroms.
    """
    if radius_type.lower() == 'vdw':
        return vdw_radius()
    elif radius_type.lower() == 'ionic':
        return ionic_radius()
    else:
        raise ValueError("radius_type must be either 'vdw' or 'ionic'")
