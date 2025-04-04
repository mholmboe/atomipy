"""
This module provides functions for handling triclinic cells and transformations.

The main function converts atom coordinates between orthogonal and triclinic
coordinate systems using either angle parameters or tilt factors.
"""

import math
import numpy as np

def orthogonal_to_triclinic(atoms, box_dim, angleparam, angletype='angle', return_box=False):
    """
    Transforms an orthogonal atom list to a triclinic one using either angles
    (alpha, beta, gamma) or tilt factors (xy, xz, yz).

    Parameters
    ----------
    atoms : list of dict
        Each dict should have 'x', 'y', 'z' coordinates (orthogonal).
    box_dim : list or tuple
        Orthogonal box dimensions: [lx, ly, lz].
    angleparam : list or tuple
        If angletype='angle', provide [alpha, beta, gamma] in degrees.
        If angletype='tilt',  provide [xy, xz, yz].
    angletype : str, optional
        'angle' or 'tilt'. Default is 'angle'.
    return_box : bool, optional
        Whether to return the triclinic box dimensions. Default is False.

    Returns
    -------
    atoms : list of dict
        The same list of atoms, but with updated 'x', 'y', 'z' for triclinic coords.
    triclinic_box_dim : list, optional
        The triclinic box dimensions in the format [lx, ly, lz, 0, 0, xy, 0, xz, yz].
        Only returned if return_box is True.

    Examples
    --------
    # Convert using angles:
    atoms = ap.triclinic.orthogonal_to_triclinic(atoms, [10, 10, 10], [85, 95, 90])
    
    # Convert using tilt factors and get new box dimensions:
    atoms, new_box = ap.triclinic.orthogonal_to_triclinic(
        atoms, [10, 10, 10], [0.5, 0.2, 0.1], angletype='tilt', return_box=True
    )
    """
    # Helper function for angle conversion
    def deg2rad(d):
        return math.radians(d)

    # Extract box dimensions
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
    if atoms:
        # For large atom lists, vectorized approach can be faster
        if len(atoms) > 1000 and all('x' in atom and 'y' in atom and 'z' in atom for atom in atoms):
            # Extract coordinates as arrays
            x_vals = np.array([atom['x'] for atom in atoms])
            y_vals = np.array([atom['y'] for atom in atoms])
            z_vals = np.array([atom['z'] for atom in atoms])
            
            # Convert to fractional coordinates
            frac_x = x_vals / lx
            frac_y = y_vals / ly
            frac_z = z_vals / lz
            
            # Apply transformation matrix
            x_new = (from_frac[0][0]*frac_x + from_frac[0][1]*frac_y + from_frac[0][2]*frac_z)
            y_new = (from_frac[1][0]*frac_x + from_frac[1][1]*frac_y + from_frac[1][2]*frac_z)
            z_new = (from_frac[2][0]*frac_x + from_frac[2][1]*frac_y + from_frac[2][2]*frac_z)
            
            # Update atoms
            for i, atom in enumerate(atoms):
                atom['x'] = x_new[i]
                atom['y'] = y_new[i]
                atom['z'] = z_new[i]
        else:
            # Original approach for smaller atom lists
            for atom in atoms:
                # Convert orthonormal coords to fractional
                frac_x = atom['x'] / lx
                frac_y = atom['y'] / ly
                frac_z = atom['z'] / lz
                
                # Apply from_frac matrix
                x_new = (from_frac[0][0]*frac_x + from_frac[0][1]*frac_y + from_frac[0][2]*frac_z)
                y_new = (from_frac[1][0]*frac_x + from_frac[1][1]*frac_y + from_frac[1][2]*frac_z)
                z_new = (from_frac[2][0]*frac_x + from_frac[2][1]*frac_y + from_frac[2][2]*frac_z)
                
                # Update atom
                atom['x'] = x_new
                atom['y'] = y_new
                atom['z'] = z_new

    # Build final triclinic box dimensions
    triclinic_box_dim = [lx, ly, lz, 0, 0, xy, 0, xz, yz]
    
    # Zero out near-zero tilts for numerical stability
    triclinic_box_dim = [0 if abs(v) < 1e-5 else v for v in triclinic_box_dim]
    
    # If no tilt remains, reduce to orthogonal
    if sum(abs(v) for v in triclinic_box_dim[3:]) < 1e-7:
        triclinic_box_dim = triclinic_box_dim[:3]

    # Return atoms and optionally the box dimensions
    if return_box:
        return atoms, triclinic_box_dim
    return atoms


def triclinic_to_orthogonal(atoms, triclinic_box_dim, return_box=False):
    """
    Transforms a triclinic atom list to an orthogonal one.
    
    Parameters
    ----------
    atoms : list of dict
        Each dict should have 'x', 'y', 'z' coordinates (triclinic).
    triclinic_box_dim : list or tuple
        Triclinic box dimensions in format [lx, ly, lz, 0, 0, xy, 0, xz, yz].
    return_box : bool, optional
        Whether to return the orthogonal box dimensions. Default is False.
        
    Returns
    -------
    atoms : list of dict
        The same list of atoms, but with updated 'x', 'y', 'z' for orthogonal coords.
    orthogonal_box_dim : list, optional
        The orthogonal box dimensions in the format [lx, ly, lz].
        Only returned if return_box is True.
    
    Examples
    --------
    # Convert to orthogonal:
    atoms = ap.triclinic.triclinic_to_orthogonal(atoms, [10, 10, 10, 0, 0, 0.5, 0, 0.2, 0.1])
    """
    # Extract box dimensions
    if len(triclinic_box_dim) >= 9:
        lx, ly, lz, _, _, xy, _, xz, yz = triclinic_box_dim[:9]
    else:
        lx, ly, lz = triclinic_box_dim[:3]
        xy = xz = yz = 0.0
    
    # Calculate angles from tilt factors
    b = math.sqrt(ly**2 + xy**2)
    c = math.sqrt(lz**2 + xz**2 + yz**2)
    
    if b*c:
        alpha = math.degrees(math.acos((ly * yz + xy * xz) / (b * c)))
    else:
        alpha = 90.0
        
    if c:
        beta = math.degrees(math.acos(xz / c))
    else:
        beta = 90.0
        
    if b:
        gamma = math.degrees(math.acos(xy / b))
    else:
        gamma = 90.0
    
    # Helper function
    def deg2rad(d):
        return math.radians(d)
    
    # Calculate matrix for fractional coordinates
    cos_a, cos_b, cos_g = math.cos(deg2rad(alpha)), math.cos(deg2rad(beta)), math.cos(deg2rad(gamma))
    sin_g = math.sin(deg2rad(gamma))
    v = math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g)
    
    # Matrix from triclinic to fractional coords
    to_frac = [
        [1/lx, -cos_g/(lx*sin_g), (cos_a*cos_g - cos_b)/(lx*sin_g*v)],
        [0, 1/(ly*sin_g), (cos_b*cos_g - cos_a)/(ly*sin_g*v)],
        [0, 0, sin_g/(lz*v)]
    ]
    
    # Update each atom coordinate
    if atoms:
        for atom in atoms:
            # Get triclinic coordinates
            x, y, z = atom['x'], atom['y'], atom['z']
            
            # Apply to_frac matrix
            frac_x = to_frac[0][0]*x + to_frac[0][1]*y + to_frac[0][2]*z
            frac_y = to_frac[1][0]*x + to_frac[1][1]*y + to_frac[1][2]*z
            frac_z = to_frac[2][0]*x + to_frac[2][1]*y + to_frac[2][2]*z
            
            # Convert fractional to orthogonal
            atom['x'] = frac_x * lx
            atom['y'] = frac_y * ly
            atom['z'] = frac_z * lz
    
    # Calculate orthogonal box dimensions
    orthogonal_box_dim = [lx, ly, lz]
    
    # Return atoms and optionally the box dimensions
    if return_box:
        return atoms, orthogonal_box_dim
    return atoms

