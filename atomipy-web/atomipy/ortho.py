"""
This module provides functions for handling orthogonal cell transformations.

The main function converts atom coordinates from triclinic to orthogonal
coordinate systems using either box dimensions or angle parameters.
"""

import numpy as np

def triclinic_to_orthogonal(atoms, box_dim, angleparam=None, angletype=None, preserve_fractional=True):
    """
    Convert atoms in a triclinic box to an orthogonal box.

    Parameters
    ----------
    atoms : list of dict
        Each dict should have at least the keys 'x', 'y', 'z'.
    box_dim : list or array
        Box dimensions that can be:
          - length 3: [lx, ly, lz] (orthogonal)
          - length 6: [a, b, c, alpha, beta, gamma] (cell parameters)
          - length 9: [lx, ly, lz, 0, 0, xy, 0, xz, yz] (triclinic)
    angleparam : list or array, optional
        If provided along with `angletype='angle'`, it should be [alpha, beta, gamma] in degrees.
        If provided along with `angletype='tilt'`, it should be [xy, xz, yz].
    angletype : str, optional
        - 'angle' => interpret angleparam as angles [alpha, beta, gamma].
        - 'tilt'  => interpret angleparam as tilt factors [xy, xz, yz].
    preserve_fractional : bool, optional
        If True, also adds fractional coordinates as 'xfrac', 'yfrac', 'zfrac' to atoms.
        Default is True.

    Returns
    -------
    atoms : list of dict
        The same list of dicts, with updated 'x', 'y', 'z' for orthogonal coordinates.
        If preserve_fractional is True, 'xfrac', 'yfrac', 'zfrac' are also assigned.
    
    Examples
    --------
    # Convert a triclinic system to orthogonal:
    atoms = ap.ortho.triclinic_to_orthogonal(atoms, [10, 10, 10, 0, 0, 0.5, 0, 0.2, 0.1])
    
    # Use explicit angles:
    atoms = ap.ortho.triclinic_to_orthogonal(
        atoms, [10, 10, 10], angleparam=[90, 110, 120], angletype='angle'
    )
    """
    # Convert input to numpy array for calculations
    box_dim = np.array(box_dim, dtype=float)
    
    # Set default values
    # a, b, c are lengths, alpha, beta, gamma are angles in degrees
    # lx, ly, lz, xy, xz, yz are triclinic tilt factors
    a = b = c = 0.0
    alpha = beta = gamma = 90.0
    lx = ly = lz = 0.0
    xy = xz = yz = 0.0

    # If angleparam is provided, use that plus box_dim[:3]
    if angleparam is not None and angletype is not None:
        angleparam = np.array(angleparam, dtype=float)
        a, b, c = box_dim[0], box_dim[1], box_dim[2]

        if angletype.lower().startswith('angle'):
            # Interpret angleparam as [alpha, beta, gamma]
            alpha, beta, gamma = angleparam
            lx = a
            xy = b * np.cos(np.radians(gamma))
            ly = np.sqrt(b**2 - xy**2)
            xz = c * np.cos(np.radians(beta))
            yz = (b*c*np.cos(np.radians(alpha)) - xy*xz)/ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

        elif angletype.lower().startswith('tilt'):
            # Interpret angleparam as [xy, xz, yz]
            xy, xz, yz = angleparam
            lx = a
            ly = np.sqrt(b**2 - xy**2)
            lz = np.sqrt(c**2 - xz**2 - yz**2)
            # Compute alpha, beta, gamma for reference
            alpha = np.degrees(np.arccos((ly*yz + xy*xz) / (b*c))) if b*c else 90
            beta  = np.degrees(np.arccos(xz / c)) if c else 90
            gamma = np.degrees(np.arccos(xy / b)) if b else 90

        else:
            print("Warning: Unrecognized angletype. Use 'angle' or 'tilt'.")
            return atoms

    else:
        # No angleparam given, interpret box_dim by length.
        if box_dim.size == 3:
            # Orthogonal box
            lx, ly, lz = box_dim
            a = lx
            b = ly
            c = lz
            alpha = beta = gamma = 90.0
            xy = xz = yz = 0.0

        elif box_dim.size == 6:
            # We assume these are cell parameters: [a, b, c, alpha, beta, gamma]
            a, b, c, alpha, beta, gamma = box_dim
            lx = a
            xy = b * np.cos(np.radians(gamma))
            ly = np.sqrt(b**2 - xy**2)
            xz = c * np.cos(np.radians(beta))
            yz = (b*c*np.cos(np.radians(alpha)) - xy*xz)/ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

        elif box_dim.size == 9:
            # Triclinic box: [lx, ly, lz, 0, 0, xy, 0, xz, yz]
            lx, ly, lz = box_dim[0], box_dim[1], box_dim[2]
            xy, xz, yz = box_dim[5], box_dim[7], box_dim[8]
            a = lx
            b = np.sqrt(ly**2 + xy**2)
            c = np.sqrt(lz**2 + xz**2 + yz**2)
            # angles for reference
            alpha = np.degrees(np.arccos((ly*yz + xy*xz)/(b*c))) if b*c else 90
            beta  = np.degrees(np.arccos(xz/c)) if c else 90
            gamma = np.degrees(np.arccos(xy/b)) if b else 90

        else:
            print("Warning: Invalid box_dim format. Expected length 3, 6, or 9.")
            return atoms

    # Zero out near-zero tilt values for numerical stability
    if abs(xy) < 1e-5: xy = 0.0
    if abs(xz) < 1e-5: xz = 0.0
    if abs(yz) < 1e-5: yz = 0.0
    
    # Compute cell volume using standard formula
    # volume = a*b*c * sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma)
    #                       + 2*cos(alpha)*cos(beta)*cos(gamma))
    alpha_rad = np.radians(alpha)
    beta_rad  = np.radians(beta)
    gamma_rad = np.radians(gamma)
    cos_a, cos_b, cos_g = np.cos(alpha_rad), np.cos(beta_rad), np.cos(gamma_rad)
    volume_factor = np.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g)
    box_volume = a * b * c * volume_factor

    # Define transformation matrices
    sin_g = np.sin(gamma_rad)
    fromFrac = np.array([
        [a, b*cos_g, c*cos_b],
        [0, b*sin_g, c*(cos_a - cos_b*cos_g)/sin_g if sin_g != 0 else 0],
        [0, 0, c*volume_factor/sin_g if sin_g != 0 else c]
    ])

    toFrac = np.linalg.inv(fromFrac)  # Inverse of fromFrac
    
    # For large atom lists, vectorized approach is more efficient
    if len(atoms) > 1000 and all('x' in atom and 'y' in atom and 'z' in atom for atom in atoms):
        # Extract coordinates as arrays
        coords = np.array([[atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)] for atom in atoms])
        
        # Convert to fractional coordinates (matrix multiplication for entire array)
        frac_coords = np.dot(coords, toFrac.T)
        
        # Convert fractional to orthogonal
        ortho_coords = np.zeros_like(frac_coords)
        ortho_coords[:, 0] = frac_coords[:, 0] * lx
        ortho_coords[:, 1] = frac_coords[:, 1] * ly
        ortho_coords[:, 2] = frac_coords[:, 2] * lz
        
        # Update atoms
        for i, atom in enumerate(atoms):
            if preserve_fractional:
                atom['xfrac'] = float(round(frac_coords[i, 0], 4))
                atom['yfrac'] = float(round(frac_coords[i, 1], 4))
                atom['zfrac'] = float(round(frac_coords[i, 2], 4))
            atom['x'] = float(round(ortho_coords[i, 0], 4))
            atom['y'] = float(round(ortho_coords[i, 1], 4))
            atom['z'] = float(round(ortho_coords[i, 2], 4))
    else:
        # Standard approach for smaller atom lists
        for atom in atoms:
            x, y, z = atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)
            
            # Convert to fractional
            frac = np.dot(toFrac, [x, y, z])
            
            # Convert back to orthogonal using only lx, ly, lz
            x_ortho = frac[0] * lx
            y_ortho = frac[1] * ly
            z_ortho = frac[2] * lz
            
            if preserve_fractional:
                atom['xfrac'] = float(round(frac[0], 4))
                atom['yfrac'] = float(round(frac[1], 4))
                atom['zfrac'] = float(round(frac[2], 4))
            
            atom['x'] = float(round(x_ortho, 4))
            atom['y'] = float(round(y_ortho, 4))
            atom['z'] = float(round(z_ortho, 4))

    # Return the modified atoms and optionally store orthogonal box info
    ortho_box_dim = [lx, ly, lz]
    
    # For compatibility, also calculate and set a box_volume variable
    # (This can be accessed via locals() if needed externally)
    _ = box_volume
    
    return atoms


def get_orthogonal_box(box_dim, angleparam=None, angletype=None):
    """
    Calculate orthogonal box dimensions from triclinic box parameters.
    
    Parameters
    ----------
    box_dim : list or array
        Box dimensions that can be:
          - length 3: [lx, ly, lz] (orthogonal) 
          - length 6: [a, b, c, alpha, beta, gamma] (cell parameters)
          - length 9: [lx, ly, lz, 0, 0, xy, 0, xz, yz] (triclinic)
    angleparam : list or array, optional
        If provided along with `angletype='angle'`, it should be [alpha, beta, gamma] in degrees.
        If provided along with `angletype='tilt'`, it should be [xy, xz, yz].
    angletype : str, optional
        - 'angle' => interpret angleparam as angles [alpha, beta, gamma].
        - 'tilt'  => interpret angleparam as tilt factors [xy, xz, yz].
        
    Returns
    -------
    ortho_box : list
        Orthogonal box dimensions [lx, ly, lz]
    """
    # Convert input to numpy array
    box_dim = np.array(box_dim, dtype=float)
    
    # Set default values
    a = b = c = 0.0
    alpha = beta = gamma = 90.0
    lx = ly = lz = 0.0
    xy = xz = yz = 0.0

    # If angleparam is provided, use that plus box_dim[:3]
    if angleparam is not None and angletype is not None:
        angleparam = np.array(angleparam, dtype=float)
        a, b, c = box_dim[0], box_dim[1], box_dim[2]

        if angletype.lower().startswith('angle'):
            # Interpret angleparam as [alpha, beta, gamma]
            alpha, beta, gamma = angleparam
            lx = a
            xy = b * np.cos(np.radians(gamma))
            ly = np.sqrt(b**2 - xy**2)
            xz = c * np.cos(np.radians(beta))
            yz = (b*c*np.cos(np.radians(alpha)) - xy*xz)/ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

        elif angletype.lower().startswith('tilt'):
            # Interpret angleparam as [xy, xz, yz]
            xy, xz, yz = angleparam
            lx = a
            ly = np.sqrt(b**2 - xy**2)
            lz = np.sqrt(c**2 - xz**2 - yz**2)

        else:
            print("Warning: Unrecognized angletype. Use 'angle' or 'tilt'.")
            return [0, 0, 0]

    else:
        # No angleparam given, interpret box_dim by length.
        if box_dim.size == 3:
            # Orthogonal box
            lx, ly, lz = box_dim
            
        elif box_dim.size == 6:
            # We assume these are cell parameters: [a, b, c, alpha, beta, gamma]
            a, b, c, alpha, beta, gamma = box_dim
            lx = a
            xy = b * np.cos(np.radians(gamma))
            ly = np.sqrt(b**2 - xy**2)
            xz = c * np.cos(np.radians(beta))
            yz = (b*c*np.cos(np.radians(alpha)) - xy*xz)/ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

        elif box_dim.size == 9:
            # Triclinic box: [lx, ly, lz, 0, 0, xy, 0, xz, yz]
            lx, ly, lz = box_dim[0], box_dim[1], box_dim[2]
            
        else:
            print("Warning: Invalid box_dim format. Expected length 3, 6, or 9.")
            return [0, 0, 0]

    return [lx, ly, lz]


# For backward compatibility with the original orto_atom function
def orto_atom(atom, box_dim, angleparam=None, angletype=None):
    """
    Legacy function for compatibility with the original orto_atom.
    
    See triclinic_to_orthogonal for documentation on parameters.
    """
    return triclinic_to_orthogonal(atom, box_dim, angleparam, angletype)
