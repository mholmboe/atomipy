import math

def triclinic_atom(atom, box_dim, angleparam, angletype):
    """
    Transforms an orthogonal atom list to a triclinic one using either angles
    (alpha, beta, gamma) or tilt factors (xy, xz, yz).

    Parameters
    ----------
    atom : list of dict
        Each dict should have 'x', 'y', 'z' coordinates (orthogonal).
    box_dim : list or tuple
        Orthogonal box dimensions: [lx, ly, lz].
    angleparam : list or tuple
        If angletype='angle', provide [alpha, beta, gamma] in degrees.
        If angletype='tilt',  provide [xy, xz, yz].
    angletype : str
        'angle' or 'tilt'.

    Returns
    -------
    atom : list of dict
        The same list of atoms, but with updated 'x', 'y', 'z' for triclinic coords.

    Notes
    -----
    - The final triclinic box parameters [lx, ly, lz, 0, 0, xy, 0, xz, yz]
      are not returned but are computed internally.
    - This closely mirrors the MATLAB triclinic_atom.m behavior.
    """

    def deg2rad(d):
        return math.radians(d)

    lx, ly, lz = box_dim[0], box_dim[1], box_dim[2]
    xy = xz = yz = 0.0

    if angletype.lower().startswith('angle'):
        # Using angles alpha, beta, gamma
        alpha, beta, gamma = angleparam
        a = lx
        # Estimate b, xy from gamma
        # (b is temporarily the distance in y-direction if xy used)
        b = ly / math.sqrt(1 - math.cos(deg2rad(gamma))**2)
        xy = b * math.cos(deg2rad(gamma))

        # Iteratively adjust c to match angles (mimicking the MATLAB approach)
        c = lz
        for _ in range(100):
            xz = c * math.cos(deg2rad(beta))
            yz = (b * c * math.cos(deg2rad(alpha)) - xy * xz) / ly
            c = math.sqrt(lz**2 + xz**2 + yz**2)

    else:
        # Using tilt factors xy, xz, yz
        xy, xz, yz = angleparam
        a = lx
        b = math.sqrt(ly**2 + xy**2)
        c = math.sqrt(lz**2 + xz**2 + yz**2)

        # Compute angles for reference (not strictly needed to do transforms):
        alpha = math.degrees(math.acos((ly * yz + xy * xz) / (b * c))) if b*c else 90.0
        beta  = math.degrees(math.acos(xz / c)) if c else 90.0
        gamma = math.degrees(math.acos(xy / b)) if b else 90.0

    # Volume factor for triclinic (1 - cos^2α - cos^2β - cos^2γ + 2cosαcosβcosγ)^0.5
    cos_a, cos_b, cos_g = math.cos(deg2rad(alpha)), math.cos(deg2rad(beta)), math.cos(deg2rad(gamma))
    v = math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g)

    # Matrix from fractional to triclinic coords
    from_frac = [
        [a,                 b*math.cos(deg2rad(gamma)),             c*math.cos(deg2rad(beta))],
        [0,                 b*math.sin(deg2rad(gamma)),             c*(cos_a - cos_b*cos_g) / (math.sin(deg2rad(gamma)) + 1e-15)],
        [0,                 0,                                      c*v / (math.sin(deg2rad(gamma)) + 1e-15)]
    ]

    # Update each atom coordinate
    if len(atom) > 0:
        for at in atom:
            # Convert orthonormal coords to fractional
            frac_x = at['x'] / lx
            frac_y = at['y'] / ly
            frac_z = at['z'] / lz
            # Apply from_frac matrix
            x_new = (from_frac[0][0]*frac_x +
                     from_frac[0][1]*frac_y +
                     from_frac[0][2]*frac_z)
            y_new = (from_frac[1][0]*frac_x +
                     from_frac[1][1]*frac_y +
                     from_frac[1][2]*frac_z)
            z_new = (from_frac[2][0]*frac_x +
                     from_frac[2][1]*frac_y +
                     from_frac[2][2]*frac_z)
            # Update atom
            at['x'] = x_new
            at['y'] = y_new
            at['z'] = z_new

    # Build final triclinic box dimensions
    final_box = [lx, ly, lz, 0, 0, xy, 0, xz, yz]
    # Zero out near-zero tilts
    final_box = [0 if abs(v) < 1e-5 else v for v in final_box]
    # If no tilt remains, reduce to orthogonal
    if sum(abs(v) for v in final_box[3:]) < 1e-7:
        final_box = final_box[:3]

    # In MATLAB, we do "assignin('caller','triclinic_Box_dim',Box_dim);" 
    # Here we simply do not return it by default.
    # If needed, you could return (atom, final_box).

    return atom
