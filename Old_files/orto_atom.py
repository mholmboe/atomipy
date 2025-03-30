import numpy as np

def orto_atom(atoms, box_dim, angleparam=None, angletype=None):
    """
    Convert atoms in a triclinic box to an orthonormal (orthogonal) box.

    Parameters
    ----------
    atoms : list of dict
        Each dict should have at least the keys 'x', 'y', 'z', 'type'.
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
    atoms : list of dict
        The same list of dicts, with updated 'x', 'y', 'z' for orthonormal coordinates.
        Additional keys 'xfrac', 'yfrac', 'zfrac' are also assigned (fractional coords).
    Notes
    -----
    - The final orthogonal box dimensions [lx, ly, lz] are stored in 'orto_Box_dim' (not returned, 
      but can be extracted if desired).
    - 'box_volume' is computed and stored (as a float) in the local variable scope (not returned).

    Example
    -------
    >>> # Suppose 'atoms' is a list of dicts with x, y, z, type, ...
    >>> # box_dim has 9 elements for a triclinic cell
    >>> atoms_ortho = orto_atom(atoms, box_dim)
    >>> # Or define angles explicitly:
    >>> atoms_ortho = orto_atom(atoms, [lx, ly, lz], angleparam=[90, 110, 120], angletype='angle')
    """

    def deg2rad(d):
        return np.radians(d)

    # Set default values
    # a, b, c are lengths, alpha, beta, gamma are angles in degrees
    # lx, ly, lz, xy, xz, yz are triclinic tilt factors
    a = b = c = 0.0
    alpha = beta = gamma = 90.0
    lx = ly = lz = 0.0
    xy = xz = yz = 0.0

    box_dim = np.array(box_dim, dtype=float)

    # If angleparam is provided, use that plus box_dim[:3]
    if angleparam is not None and angletype is not None:
        angleparam = np.array(angleparam, dtype=float)
        a, b, c = box_dim[0], box_dim[1], box_dim[2]

        if angletype.lower().startswith('angle'):
            # Interpret angleparam as [alpha, beta, gamma]
            alpha, beta, gamma = angleparam
            lx = a
            xy = b * np.cos(deg2rad(gamma))
            ly = np.sqrt(b**2 - xy**2)
            xz = c * np.cos(deg2rad(beta))
            yz = (b*c*np.cos(deg2rad(alpha)) - xy*xz)/ly
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
            print("Unrecognized angletype. Use 'angle' or 'tilt'.")
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
            xy = b * np.cos(deg2rad(gamma))
            ly = np.sqrt(b**2 - xy**2)
            xz = c * np.cos(deg2rad(beta))
            yz = (b*c*np.cos(deg2rad(alpha)) - xy*xz)/ly
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
            print("Something is wrong with box_dim. Expected length 3, 6, or 9.")
            return atoms

    # Zero out near-zero tilt values
    # (avoiding floating precision issues by checking absolute < 1e-5)
    # Then check if it's basically orthogonal
    for i, val in enumerate([lx, ly, lz, xy, xz, yz]):
        if abs(val) < 1e-5:
            if i == 0: lx = 0.0
            elif i == 1: ly = 0.0
            elif i == 2: lz = 0.0
            elif i == 3: xy = 0.0
            elif i == 4: xz = 0.0
            elif i == 5: yz = 0.0

    # Compute cell volume using standard formula
    # volume = a*b*c * sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma)
    #                       + 2*cos(alpha)*cos(beta)*cos(gamma))
    alpha_rad = deg2rad(alpha)
    beta_rad  = deg2rad(beta)
    gamma_rad = deg2rad(gamma)
    cos_a, cos_b, cos_g = np.cos(alpha_rad), np.cos(beta_rad), np.cos(gamma_rad)
    volume_factor = np.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g)
    box_volume = a * b * c * volume_factor

    # Define transformation matrices
    fromFrac = np.array([
        [a, b*np.cos(gamma_rad),       c*np.cos(beta_rad)],
        [0, b*np.sin(gamma_rad),       c*(cos_a - cos_b*cos_g)/np.sin(gamma_rad) if np.sin(gamma_rad) != 0 else 0],
        [0, 0,                         c*volume_factor/np.sin(gamma_rad) if np.sin(gamma_rad) != 0 else c]
    ])

    toFrac = np.linalg.inv(fromFrac)  # Inverse of fromFrac

    # Transform each atom
    for atom in atoms:
        x, y, z = atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)
        # Convert to fractional
        frac = toFrac.dot([x, y, z])
        # Convert back to orthonormal coords using only lx, ly, lz
        x_ortho = frac[0] * lx
        y_ortho = frac[1] * ly
        z_ortho = frac[2] * lz
        atom['xfrac'] = round(float(frac[0]), 4)
        atom['yfrac'] = round(float(frac[1]), 4)
        atom['zfrac'] = round(float(frac[2]), 4)
        atom['x'] = round(x_ortho, 4)
        atom['y'] = round(y_ortho, 4)
        atom['z'] = round(z_ortho, 4)

    # Optionally store orthonormal box info if desired
    orto_box_dim = [lx, ly, lz]
    # box_volume is also available as a local variable

    return atoms
