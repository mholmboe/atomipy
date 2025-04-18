"""
This module provides functions for converting between triclinic and orthogonal coordinates.

The implementation follows the approach of the MATLAB function orto_MATLAB.m,
converting between triclinic and orthogonal coordinates using crystallographic transformations.
"""

import numpy as np
from .cell_utils import Box_dim2Cell, Cell2Box_dim

def orto_coordinates(atoms, box=None, add_to_atoms=True):
    """
    Transform triclinic atom coordinates to orthogonal coordinates.
    
    This function follows the approach of the MATLAB function orto_MATLAB.m.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with 'x', 'y', 'z' cartesian coordinates in the triclinic frame.
    box: a 1x3, 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
        a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
        a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells,
        or simple orthogonal box dimensions [lx, ly, lz].
    add_to_atoms : bool, optional
        If True, adds fractional and orthogonal coordinates to the atom dictionaries as 
        'xfrac', 'yfrac', 'zfrac' and 'x_ortho', 'y_ortho', 'z_ortho'. Default is True.
        
    Returns
    -------
    orto_atoms : list of dict
        The atoms list with orthogonalized coordinates.
    orto_box_dim : list
        The orthogonal box dimensions [lx, ly, lz].
    """
    # Initialize variables
    Cell = None
    Box_dim = None
    
    if box is None:
        raise ValueError("Box parameter must be provided")
    
    # Determine box format and convert as needed
    if len(box) == 9:
        Box_dim = box
        Cell = Box_dim2Cell(Box_dim)
    elif len(box) == 6:
        Cell = box
        Box_dim = Cell2Box_dim(Cell)
    elif len(box) == 3:  # Orthogonal box
        Box_dim = box
        Cell = list(box) + [90.0, 90.0, 90.0]
    else:
        raise ValueError("Box must be length 3, 6, or 9")
    
    # Extract box dimensions
    if len(Box_dim) == 3:
        # Orthogonal box
        lx, ly, lz = Box_dim
        xy, xz, yz = 0, 0, 0
    elif len(Box_dim) == 9:
        # Triclinic box in GROMACS format [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        lx, ly, lz = Box_dim[0], Box_dim[1], Box_dim[2]
        xy, xz, yz = Box_dim[5], Box_dim[7], Box_dim[8]
    
    # Extract cell parameters
    a, b, c, alpha, beta, gamma = Cell
    
    # Now create the box_dim for the orthogonal box
    # Calculate cell parameters from box dimensions if needed
    if len(Box_dim) == 3:  # Orthogonal box
        a = lx
        b = ly
        c = lz
        # For orthogonal box, angles are 90 degrees
        alpha, beta, gamma = 90.0, 90.0, 90.0
    elif len(Box_dim) == 9:  # GROMACS format
        a = lx
        b = np.sqrt(ly**2 + xy**2)
        c = np.sqrt(lz**2 + xz**2 + yz**2)
        alpha = np.degrees(np.arccos((ly * yz + xy * xz) / (b * c)))
        beta = np.degrees(np.arccos(xz / c))
        gamma = np.degrees(np.arccos(xy / b))
    
    # Set the box_dim for the orthogonal box
    orto_box_dim = [lx, ly, lz]
    
    # Clean up small values
    Box_dim = [x if abs(x) > 1e-5 else 0 for x in Box_dim]
    
    # Calculate volume term for transformation matrices
    v = np.sqrt(1 - np.cos(np.radians(alpha))**2 - np.cos(np.radians(beta))**2 - 
                 np.cos(np.radians(gamma))**2 + 2 * np.cos(np.radians(alpha)) * 
                 np.cos(np.radians(beta)) * np.cos(np.radians(gamma)))
    
    # From fractional to Cartesian coordinates (MATLAB's FromFrac matrix)
    from_frac = np.array([
        [a, b * np.cos(np.radians(gamma)), c * np.cos(np.radians(beta))],
        [0, b * np.sin(np.radians(gamma)), c * (np.cos(np.radians(alpha)) - 
                                                np.cos(np.radians(beta)) * 
                                                np.cos(np.radians(gamma))) / 
                                                np.sin(np.radians(gamma))],
        [0, 0, c * v / np.sin(np.radians(gamma))]
    ])
    
    # From Cartesian to fractional coordinates (MATLAB's ToFrac matrix)
    to_frac = np.array([
        [1/a, -np.cos(np.radians(gamma)) / (a * np.sin(np.radians(gamma))), 
         (np.cos(np.radians(alpha)) * np.cos(np.radians(gamma)) - 
          np.cos(np.radians(beta))) / (a * v * np.sin(np.radians(gamma)))],
        [0, 1 / (b * np.sin(np.radians(gamma))), 
         (np.cos(np.radians(beta)) * np.cos(np.radians(gamma)) - 
          np.cos(np.radians(alpha))) / (b * v * np.sin(np.radians(gamma)))],
        [0, 0, np.sin(np.radians(gamma)) / (c * v)]
    ])
    
    # Verify matrices are inverse of each other (to_frac * from_frac should be identity)
    # np.allclose(np.dot(to_frac, from_frac), np.identity(3))
    
    # Create output
    import copy
    orto_atoms = copy.deepcopy(atoms)
    
    # Process each atom
    for atom in orto_atoms:
        # Get cartesian coordinates
        xyz = np.array([atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)])
        
        # Convert to fractional coordinates
        frac_coords = np.dot(to_frac, xyz)
        
        # Convert to orthogonal coordinates
        ortho_coords = np.array([lx, ly, lz]) * frac_coords
        
        # Store fractional coordinates
        atom['xfrac'] = float(round(frac_coords[0], 4))
        atom['yfrac'] = float(round(frac_coords[1], 4))
        atom['zfrac'] = float(round(frac_coords[2], 4))
        
        if add_to_atoms:
            # Store orthogonal coordinates as additional fields
            atom['x_ortho'] = float(round(ortho_coords[0], 4))
            atom['y_ortho'] = float(round(ortho_coords[1], 4))
            atom['z_ortho'] = float(round(ortho_coords[2], 4))
        
        # Update original coordinates to orthogonal coordinates
        atom['x'] = float(round(ortho_coords[0], 4))
        atom['y'] = float(round(ortho_coords[1], 4))
        atom['z'] = float(round(ortho_coords[2], 4))
    
    # orto_box_dim was already defined above
    return orto_atoms, orto_box_dim

def cartesian_to_fractional(atoms, box, add_to_atoms=True):
    """
    Convert Cartesian coordinates to fractional coordinates.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with 'x', 'y', 'z' cartesian coordinates.
    box: a 1x3, 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
        a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
        a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells,
        or simple orthogonal box dimensions [lx, ly, lz].
        
    add_to_atoms : bool, optional
        If True, adds fractional coordinates to the atom dictionaries as 
        'xfrac', 'yfrac', 'zfrac'. Default is True.
        
    Returns
    -------
    frac_coords : numpy.ndarray
        Nx3 array of fractional coordinates, where N is the number of atoms.
    atoms : list of dict
        The original atoms list with added fractional coordinate fields
        if add_to_atoms is True.
    """
    # Initialize variables
    Cell = None
    Box_dim = None
    
    if box is None:
        raise ValueError("Box parameter must be provided")
    
    # Determine box format and convert as needed
    if len(box) == 9:
        Box_dim = box
        Cell = Box_dim2Cell(Box_dim)
    elif len(box) == 6:
        Cell = box
        Box_dim = Cell2Box_dim(Cell)
    elif len(box) == 3:  # Orthogonal box
        Box_dim = box
        Cell = list(box) + [90.0, 90.0, 90.0]
    else:
        raise ValueError("Box must be length 3, 6, or 9")
    
    # Extract box dimensions
    if len(Box_dim) == 3:
        # Orthogonal box
        lx, ly, lz = Box_dim
        xy, xz, yz = 0, 0, 0
    elif len(Box_dim) == 9:
        # Triclinic box in GROMACS format [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        lx, ly, lz = Box_dim[0], Box_dim[1], Box_dim[2]
        xy, xz, yz = Box_dim[5], Box_dim[7], Box_dim[8]
        
    # Extract cell parameters (we've already set Cell from the box parameter)
    a, b, c, alpha, beta, gamma = Cell
    
    # Calculate volume term
    v = np.sqrt(1 - np.cos(np.radians(alpha))**2 - np.cos(np.radians(beta))**2 - 
                np.cos(np.radians(gamma))**2 + 2 * np.cos(np.radians(alpha)) * 
                np.cos(np.radians(beta)) * np.cos(np.radians(gamma)))
    
    # Build transformation matrix from Cartesian to fractional
    to_frac = np.array([
        [1/a, -np.cos(np.radians(gamma)) / (a * np.sin(np.radians(gamma))), 
         (np.cos(np.radians(alpha)) * np.cos(np.radians(gamma)) - 
          np.cos(np.radians(beta))) / (a * v * np.sin(np.radians(gamma)))],
        [0, 1 / (b * np.sin(np.radians(gamma))), 
         (np.cos(np.radians(beta)) * np.cos(np.radians(gamma)) - 
          np.cos(np.radians(alpha))) / (b * v * np.sin(np.radians(gamma)))],
        [0, 0, np.sin(np.radians(gamma)) / (c * v)]
    ])
    
    # Extract cartesian coordinates from atoms
    cart_coords = np.array([[atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)] 
                           for atom in atoms])
    
    # Apply transformation to each point
    frac_coords = np.zeros_like(cart_coords)
    for i, xyz in enumerate(cart_coords):
        frac_coords[i] = np.dot(to_frac, xyz)
    
    # Add fractional coordinates to atoms if requested
    if add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['xfrac'] = float(round(frac_coords[i, 0], 4))
            atom['yfrac'] = float(round(frac_coords[i, 1], 4))
            atom['zfrac'] = float(round(frac_coords[i, 2], 4))
        return frac_coords, atoms
    
    return frac_coords

def fractional_to_cartesian(atoms=None, frac_coords=None, box=None, add_to_atoms=True):
    """
    Convert fractional coordinates to Cartesian coordinates.
    
    Parameters
    ----------
    atoms : list of dict, optional
        List of atom dictionaries with 'xfrac', 'yfrac', 'zfrac' fractional coordinates.
        If provided, these will be used for the conversion. Either atoms or frac_coords must be provided.
    frac_coords : numpy.ndarray, optional
        Nx3 array of fractional coordinates, where N is the number of atoms.
        Either atoms or frac_coords must be provided.
    box: a 1x3, 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
        a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
        a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells,
        or simple orthogonal box dimensions [lx, ly, lz].
    add_to_atoms : bool, optional
        If True and atoms is provided, adds cartesian coordinates to the atom dictionaries
        as 'x', 'y', 'z'. Default is True.
        
    Returns
    -------
    cart_coords : numpy.ndarray
        Nx3 array of cartesian coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with updated cartesian coordinate fields
        if add_to_atoms is True.
    """
    # Initialize variables
    Cell = None
    Box_dim = None
    
    if (atoms is None and frac_coords is None) or box is None:
        raise ValueError("Either (atoms or frac_coords) and box must be provided")
    
    # Determine box format and convert as needed
    if len(box) == 9:
        Box_dim = box
        Cell = Box_dim2Cell(Box_dim)
    elif len(box) == 6:
        Cell = box
        Box_dim = Cell2Box_dim(Cell)
    elif len(box) == 3:  # Orthogonal box
        Box_dim = box
        Cell = list(box) + [90.0, 90.0, 90.0]
    else:
        raise ValueError("Box must be length 3, 6, or 9")
    
    # Extract cell parameters
    a, b, c, alpha, beta, gamma = Cell
    
    # Calculate volume term
    v = np.sqrt(1 - np.cos(np.radians(alpha))**2 - np.cos(np.radians(beta))**2 - 
                np.cos(np.radians(gamma))**2 + 2 * np.cos(np.radians(alpha)) * 
                np.cos(np.radians(beta)) * np.cos(np.radians(gamma)))
    
    # Build transformation matrix from fractional to cartesian
    from_frac = np.array([
        [a, b * np.cos(np.radians(gamma)), c * np.cos(np.radians(beta))],
        [0, b * np.sin(np.radians(gamma)), c * (np.cos(np.radians(alpha)) - 
                                                np.cos(np.radians(beta)) * 
                                                np.cos(np.radians(gamma))) / 
                                                np.sin(np.radians(gamma))],
        [0, 0, c * v / np.sin(np.radians(gamma))]
    ])
    
    # Get fractional coordinates from atoms if needed
    if frac_coords is None:
        frac_coords = np.array([[atom.get('xfrac', 0.0), 
                                 atom.get('yfrac', 0.0), 
                                 atom.get('zfrac', 0.0)] 
                                for atom in atoms])
    
    # Apply transformation to each point
    cart_coords = np.zeros_like(frac_coords)
    for i, frac in enumerate(frac_coords):
        cart_coords[i] = np.dot(from_frac, frac)
    
    # Add cartesian coordinates to atoms if requested
    if atoms is not None and add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['x'] = float(round(cart_coords[i, 0], 4))
            atom['y'] = float(round(cart_coords[i, 1], 4))
            atom['z'] = float(round(cart_coords[i, 2], 4))
        return cart_coords, atoms
    
    return cart_coords

def wrap_coordinates(atoms, box=None, in_place=True):
    """
    Wrap atoms into the primary unit cell (0 ≤ frac < 1).
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with cartesian coordinates.
    box: a 1x3, 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
        a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
        a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells,
        or simple orthogonal box dimensions [lx, ly, lz].
    in_place : bool, optional
        If True, modifies the input atoms list in place. Default is True.
        
    Returns
    -------
    atoms : list of dict
        The atoms list with wrapped coordinates.
    """
    if box is None:
        raise ValueError("Box parameter must be provided")
    
    # Create a copy of atoms if not modifying in place
    if not in_place:
        import copy
        atoms = copy.deepcopy(atoms)
    
    # Convert to fractional coordinates
    frac_coords, atoms = cartesian_to_fractional(atoms, box=box, add_to_atoms=True)
    
    # Wrap fractional coordinates to [0, 1)
    for i, atom in enumerate(atoms):
        atom['xfrac'] = atom['xfrac'] % 1.0
        atom['yfrac'] = atom['yfrac'] % 1.0
        atom['zfrac'] = atom['zfrac'] % 1.0
    
    # Convert back to cartesian coordinates
    fractional_to_cartesian(atoms, box=box, add_to_atoms=True)
    
    return atoms

def triclinic_to_orthogonal(atoms, box=None, add_to_atoms=True):
    """
    Convert coordinates from a triclinic cell to an orthogonal cell.
    This is a wrapper around orto_coordinates for backward compatibility.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with 'x', 'y', 'z' cartesian coordinates in the triclinic frame.
    box : a 1x3, 1x6 or 1x9 list representing cell dimensions (in Angstroms):
        - For orthogonal boxes, a 1x3 list [lx, ly, lz] where box = Box_dim, and Cell would be [lx, ly, lz, 90, 90, 90]
        - For cell parameters, a 1x6 list [a, b, c, alpha, beta, gamma] (Cell format)
        - For triclinic boxes, a 1x9 list [lx, ly, lz, 0, 0, xy, 0, xz, yz] (GROMACS Box_dim format)
    add_to_atoms : bool, optional
        If True, adds orthogonal coordinates to the atom dictionaries as 
        'x_ortho', 'y_ortho', 'z_ortho'. Default is True.
        
    Returns
    -------
    ortho_box : list
        Orthogonal box dimensions [lx, ly, lz].
    ortho_coords : numpy.ndarray
        Nx3 array of orthogonal coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with added orthogonal coordinate fields
        if add_to_atoms is True.
    """
    orto_atoms, orto_box_dim = orto_coordinates(atoms, box=box, add_to_atoms=add_to_atoms)
    
    # Extract orthogonal coordinates
    ortho_coords = np.array([[atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)] 
                           for atom in orto_atoms])
    
    # Return orthogonal box, coordinates, and original atoms with added fields
    return orto_box_dim, ortho_coords, atoms

def get_cell_vectors(box_dim=None, cell=None):
    """
    Get the three cell vectors from box dimensions or cell parameters.
    
    Parameters
    ----------
    box: a 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
        a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
        a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
        Note that for orthogonal boxes Cell = Box_dim.
        
    Returns
    -------
    cell_vectors : numpy.ndarray
        3x3 array with the three cell vectors as rows.
    """
    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # If box_dim is provided but not cell, calculate cell parameters
    if cell is None:
        cell, _ = Box_dim2Cell(box_dim[:6] if len(box_dim) >= 6 else box_dim)
    
    # Extract cell parameters
    a, b, c, alpha, beta, gamma = cell
    
    # Convert angles to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Calculate trigonometric values
    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    sin_gamma = np.sin(gamma_rad)
    
    # Calculate volume term
    v = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 
                2 * cos_alpha * cos_beta * cos_gamma)
    
    # Build transformation matrix from fractional to cartesian
    from_frac = np.array([
        [a, b * cos_gamma, c * cos_beta],
        [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
        [0, 0, c * v / sin_gamma]
    ])
    
    # Cell vectors are columns of the from_frac matrix
    a_vec = from_frac[:, 0]
    b_vec = from_frac[:, 1]
    c_vec = from_frac[:, 2]
    
    # Return cell vectors as rows
    return np.array([a_vec, b_vec, c_vec])

# Add direct transformation functions to replace direct_fract module
def direct_cartesian_to_fractional(atoms, box_dim=None, cell=None, add_to_atoms=True):
    """
    Direct conversion from Cartesian coordinates to fractional coordinates.
    This function provides a direct implementation that follows the MATLAB approach
    without intermediate orthogonalization steps.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with 'x', 'y', 'z' cartesian coordinates.
    box: a 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
        a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
        a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
        Note that for orthogonal boxes Cell = Box_dim.
    add_to_atoms : bool, optional
        If True, adds fractional coordinates to the atom dictionaries as 
        'xfrac', 'yfrac', 'zfrac'. Default is True.
        
    Returns
    -------
    frac_coords : numpy.ndarray
        Nx3 array of fractional coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with added fractional coordinate fields
        if add_to_atoms is True.
    """
    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # If box_dim is provided but not cell, calculate cell parameters
    if cell is None:
        cell, _ = Box_dim2Cell(box_dim[:6] if len(box_dim) >= 6 else box_dim)
    
    # Extract cell parameters
    a, b, c, alpha, beta, gamma = cell
    
    # Convert angles to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Calculate trigonometric values
    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    sin_gamma = np.sin(gamma_rad)
    
    # Calculate volume term
    v = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 
                2 * cos_alpha * cos_beta * cos_gamma)
    
    # Build transformation matrix from Cartesian to fractional (ToFrac in MATLAB)
    to_frac = np.array([
        [1/a, -np.cos(gamma_rad) / (a * np.sin(gamma_rad)), 
         (np.cos(alpha_rad) * np.cos(gamma_rad) - 
          np.cos(beta_rad)) / (a * v * np.sin(gamma_rad))],
        [0, 1 / (b * np.sin(gamma_rad)), 
         (np.cos(beta_rad) * np.cos(gamma_rad) - 
          np.cos(alpha_rad)) / (b * v * np.sin(gamma_rad))],
        [0, 0, np.sin(gamma_rad) / (c * v)]
    ])
    
    # Extract cartesian coordinates from atoms
    cart_coords = np.array([[atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)] 
                             for atom in atoms])
    
    # Apply transformation to each point
    frac_coords = np.zeros_like(cart_coords)
    for i, xyz in enumerate(cart_coords):
        frac_coords[i] = np.dot(to_frac, xyz)
    
    # Add fractional coordinates to atoms if requested
    if add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['xfrac'] = float(round(frac_coords[i, 0], 4))
            atom['yfrac'] = float(round(frac_coords[i, 1], 4))
            atom['zfrac'] = float(round(frac_coords[i, 2], 4))
        return frac_coords, atoms
    
    return frac_coords

def direct_fractional_to_cartesian(atoms=None, frac_coords=None, box_dim=None, cell=None, add_to_atoms=True):
    """
    Direct conversion from fractional coordinates to Cartesian coordinates.
    This function provides a direct implementation that follows the MATLAB approach
    without intermediate orthogonalization steps.
    
    Parameters
    ----------
    atoms : list of dict, optional
        List of atom dictionaries with 'xfrac', 'yfrac', 'zfrac' fractional coordinates.
        If provided, these will be used for the conversion. Either atoms or frac_coords must be provided.
    frac_coords : numpy.ndarray, optional
        Nx3 array of fractional coordinates, where N is the number of atoms.
        Either atoms or frac_coords must be provided.
    box: a 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
        a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
        a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
        Note that for orthogonal boxes Cell = Box_dim.
    add_to_atoms : bool, optional
        If True and atoms is provided, adds cartesian coordinates to the atom dictionaries
        as 'x', 'y', 'z'. Default is True.
        
    Returns
    -------
    cart_coords : numpy.ndarray
        Nx3 array of cartesian coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with updated cartesian coordinate fields
        if add_to_atoms is True.
    """
    if (atoms is None and frac_coords is None) or (box_dim is None and cell is None):
        raise ValueError("Either (atoms or frac_coords) and (box_dim or cell) must be provided")
    
    # If box_dim is provided but not cell, calculate cell parameters
    if cell is None:
        cell, _ = Box_dim2Cell(box_dim[:6] if len(box_dim) >= 6 else box_dim)
    
    # Extract cell parameters
    a, b, c, alpha, beta, gamma = cell
    
    # Convert angles to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Calculate trigonometric values
    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)
    sin_gamma = np.sin(gamma_rad)
    
    # Calculate volume term
    v = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 
                2 * cos_alpha * cos_beta * cos_gamma)
    
    # Build transformation matrix from fractional to cartesian (FromFrac in MATLAB)
    from_frac = np.array([
        [a, b * cos_gamma, c * cos_beta],
        [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
        [0, 0, c * v / sin_gamma]
    ])
    
    # Get fractional coordinates from atoms if needed
    if frac_coords is None:
        frac_coords = np.array([[atom.get('xfrac', 0.0), 
                                 atom.get('yfrac', 0.0), 
                                 atom.get('zfrac', 0.0)] 
                                for atom in atoms])
    
    # Apply transformation to each point
    cart_coords = np.zeros_like(frac_coords)
    for i, frac in enumerate(frac_coords):
        cart_coords[i] = np.dot(from_frac, frac)
    
    # Add cartesian coordinates to atoms if requested
    if atoms is not None and add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['x'] = float(round(cart_coords[i, 0], 4))
            atom['y'] = float(round(cart_coords[i, 1], 4))
            atom['z'] = float(round(cart_coords[i, 2], 4))
        return cart_coords, atoms
    
    return cart_coords
