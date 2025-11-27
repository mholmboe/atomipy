def vdw_radius():
    """
    Return a dictionary mapping elements to their van der Waals radii (Å).

    Returns
    -------
    dict
        Element symbol -> van der Waals radius in Angstroms.

    Examples
    --------
    # Access oxygen VdW radius
    from atomipy.radius import vdw_radius
    radii = vdw_radius()
    radii['O']  # 1.52
    """
    radii = {
        # Alkali metals
        'H': 1.20,
        'Li': 1.82,
        'Na': 2.27,
        'K': 2.75,
        'Rb': 3.03,
        'Cs': 3.43,
        'Fr': 3.48,  # Estimated
        
        # Alkaline earth metals
        'Be': 1.53,
        'Mg': 1.73,
        'Ca': 2.31,
        'Sr': 2.49,
        'Ba': 2.68,
        'Ra': 2.83,  # Estimated
        
        # Transition metals
        'Sc': 2.15,
        'Ti': 2.00,
        'V': 2.05,
        'Cr': 2.05,
        'Mn': 2.05,
        'Fe': 2.10,
        'Co': 2.00,
        'Ni': 1.97,
        'Cu': 1.96,
        'Zn': 2.01,
        'Y': 2.32,
        'Zr': 2.23,
        'Nb': 2.18,
        'Mo': 2.17,
        'Tc': 2.16,  # Estimated
        'Ru': 2.13,
        'Rh': 2.10,
        'Pd': 2.10,
        'Ag': 2.11,
        'Cd': 2.18,
        'Hf': 2.23,
        'Ta': 2.22,
        'W': 2.18,
        'Re': 2.16,
        'Os': 2.16,
        'Ir': 2.13,
        'Pt': 2.13,
        'Au': 2.14,
        'Hg': 2.23,
        
        # Post-transition metals
        'Al': 1.84,
        'Ga': 1.87,
        'In': 1.93,
        'Sn': 2.17,
        'Tl': 1.96,
        'Pb': 2.02,
        'Bi': 2.07,
        
        # Metalloids
        'B': 1.92,
        'Si': 2.10,
        'Ge': 2.11,
        'As': 1.85,
        'Sb': 2.06,
        'Te': 2.06,
        'Po': 1.97,  # Estimated
        
        # Non-metals
        'C': 1.70,
        'N': 1.55,
        'O': 1.52,
        'F': 1.47,
        'P': 1.80,
        'S': 1.80,
        'Cl': 1.75,
        'Se': 1.90,
        'Br': 1.85,
        'I': 1.98,
        'At': 2.02,  # Estimated
        
        # Noble gases
        'He': 1.40,
        'Ne': 1.54,
        'Ar': 1.88,
        'Kr': 2.02,
        'Xe': 2.16,
        'Rn': 2.20   # Estimated
    }
    return radii


def ionic_radius():
    """Return dictionary mapping elements to their ionic radii.
    Values in angstroms for the common oxidation states.
    
    Returns
    -------
    dict
        Dictionary mapping element symbols to their ionic radii in Angstroms.
        For elements with multiple oxidation states, the most common one is used.
    """
    radii = {
        # Alkali metals (usually +1)
        'H': 0.25,  # H+
        'Li': 0.76, # Li+
        'Na': 1.02, # Na+
        'K': 1.38,  # K+
        'Rb': 1.52, # Rb+
        'Cs': 1.67, # Cs+
        'Fr': 1.80, # Fr+ (estimated)
        
        # Alkaline earth metals (usually +2)
        'Be': 0.45, # Be2+
        'Mg': 0.72, # Mg2+
        'Ca': 1.00, # Ca2+
        'Sr': 1.18, # Sr2+
        'Ba': 1.35, # Ba2+
        'Ra': 1.43, # Ra2+ (estimated)
        
        # Transition metals (common oxidation states)
        'Sc': 0.75, # Sc3+
        'Ti': 0.61, # Ti4+
        'V': 0.54,  # V5+
        'Cr': 0.52, # Cr6+
        'Mn': 0.67, # Mn2+
        'Fe': 0.65, # Fe3+
        'Co': 0.65, # Co2+
        'Ni': 0.69, # Ni2+
        'Cu': 0.73, # Cu2+
        'Zn': 0.74, # Zn2+
        'Y': 0.90,  # Y3+
        'Zr': 0.72, # Zr4+
        'Nb': 0.64, # Nb5+
        'Mo': 0.59, # Mo6+
        'Tc': 0.70, # Tc4+
        'Ru': 0.68, # Ru3+
        'Rh': 0.67, # Rh3+
        'Pd': 0.76, # Pd2+
        'Ag': 1.15, # Ag+
        'Cd': 0.95, # Cd2+
        'Hf': 0.71, # Hf4+
        'Ta': 0.64, # Ta5+
        'W': 0.60,  # W6+
        'Re': 0.72, # Re4+
        'Os': 0.69, # Os4+
        'Ir': 0.63, # Ir4+
        'Pt': 0.63, # Pt4+
        'Au': 0.85, # Au3+
        'Hg': 1.02, # Hg2+
        
        # Post-transition metals
        'Al': 0.54, # Al3+
        'Ga': 0.62, # Ga3+
        'In': 0.80, # In3+
        'Sn': 0.69, # Sn4+
        'Tl': 0.88, # Tl3+
        'Pb': 1.19, # Pb2+
        'Bi': 1.03, # Bi3+
        
        # Metalloids
        'B': 0.23,  # B3+
        'Si': 0.40, # Si4+
        'Ge': 0.53, # Ge4+
        'As': 0.46, # As5+
        'Sb': 0.76, # Sb3+
        'Te': 0.56, # Te6+
        'Po': 0.94, # Po4+ (estimated)
        
        # Non-metals
        'C': 0.16,  # C4+
        'N': 0.16,  # N5+
        'O': 1.40,  # O2-
        'F': 1.33,  # F-
        'P': 0.38,  # P5+
        'S': 0.29,  # S6+
        'Cl': 1.81, # Cl-
        'Se': 0.42, # Se6+
        'Br': 1.96, # Br-
        'I': 2.20,  # I-
        'At': 2.27  # At- (estimated)
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
