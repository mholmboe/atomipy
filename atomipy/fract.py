"""
This module provides functions for converting between cartesian and fractional coordinates.

Functions in this module allow conversion of atom coordinates between cartesian (real-space)
and fractional (unit cell) representations, using either Box_dim or Cell parameters.
"""

import numpy as np
from .cell_utils import Box_dim2Cell, Cell2Box_dim


def cartesian_to_fractional(atoms, box_dim=None, cell=None, add_to_atoms=True):
    """
    Convert cartesian coordinates to fractional coordinates.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with 'x', 'y', 'z' cartesian coordinates.
    box_dim : list or array, optional
        Box dimensions in one of the supported formats:
        - [lx, ly, lz] (orthogonal)
        - [lx, ly, lz, xy, xz, yz] (triclinic with tilt factors)
        - [lx, ly, lz, alpha, beta, gamma] (triclinic with angles)
        - [lx, ly, lz, 0, 0, xy, 0, xz, yz] (GROMACS format)
        Either box_dim or cell must be provided.
    cell : list or array, optional
        Cell parameters as [a, b, c, alpha, beta, gamma].
        Either box_dim or cell must be provided.
    add_to_atoms : bool, optional
        If True, adds fractional coordinates to the atom dictionaries as 
        'xfrac', 'yfrac', 'zfrac'. Default is True.
        
    Returns
    -------
    frac_coords : numpy.ndarray
        Nx3 array of fractional coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with added 'xfrac', 'yfrac', 'zfrac' fields
        if add_to_atoms is True.
        
    Examples
    --------
    # Using box dimensions:
    frac_coords = ap.fract.cartesian_to_fractional(atoms, box_dim=[10, 10, 10])
    
    # Using cell parameters:
    frac_coords = ap.fract.cartesian_to_fractional(
        atoms, cell=[10, 10, 10, 90, 90, 90]
    )
    
    # Get fractional coordinates without modifying atoms:
    frac_coords = ap.fract.cartesian_to_fractional(
        atoms, box_dim=[10, 10, 10], add_to_atoms=False
    )
    """
    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # If cell is provided but not box_dim, use Box_dim2Cell to get the box
    if box_dim is None:
        _, box = Box_dim2Cell(cell)
    else:
        # Get the box matrix from box_dim
        if len(box_dim) == 3:
            # Orthogonal box
            box = np.zeros((3, 3))
            box[0, 0] = box_dim[0]  # lx
            box[1, 1] = box_dim[1]  # ly
            box[2, 2] = box_dim[2]  # lz
        elif len(box_dim) == 6:
            # Could be [lx, ly, lz, xy, xz, yz] or [lx, ly, lz, alpha, beta, gamma]
            # Let's check if the last three values might be angles
            if all(angle > 0 and angle < 180 for angle in box_dim[3:6]):
                # Treat as [lx, ly, lz, alpha, beta, gamma]
                cell_params = box_dim
                _, box = Box_dim2Cell(cell_params)
            else:
                # Treat as [lx, ly, lz, xy, xz, yz]
                lx, ly, lz, xy, xz, yz = box_dim
                box = np.zeros((3, 3))
                box[0, 0] = lx
                box[1, 1] = ly
                box[2, 2] = lz
                box[0, 1] = xy
                box[0, 2] = xz
                box[1, 2] = yz
        elif len(box_dim) == 9:
            # GROMACS format: [lx, ly, lz, 0, 0, xy, 0, xz, yz]
            lx, ly, lz = box_dim[0], box_dim[1], box_dim[2]
            xy, xz, yz = box_dim[5], box_dim[7], box_dim[8]
            
            box = np.zeros((3, 3))
            box[0, 0] = lx
            box[1, 1] = ly
            box[2, 2] = lz
            box[0, 1] = xy
            box[0, 2] = xz
            box[1, 2] = yz
        else:
            raise ValueError(f"Invalid box_dim length: {len(box_dim)}. Expected 3, 6, or 9.")
    
    # Invert the box matrix to get the transformation matrix
    inv_box = np.linalg.inv(box)
    
    # Extract cartesian coordinates from atoms
    cart_coords = np.array([[atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)] 
                            for atom in atoms])
    
    # Convert to fractional coordinates using the transformation matrix
    frac_coords = np.dot(cart_coords, inv_box)
    
    # Add fractional coordinates to atoms if requested
    if add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['xfrac'] = float(frac_coords[i, 0])
            atom['yfrac'] = float(frac_coords[i, 1])
            atom['zfrac'] = float(frac_coords[i, 2])
        return frac_coords, atoms
    
    return frac_coords


def fractional_to_cartesian(atoms=None, frac_coords=None, box_dim=None, cell=None, add_to_atoms=True):
    """
    Convert fractional coordinates to cartesian coordinates.
    
    Parameters
    ----------
    atoms : list of dict, optional
        List of atom dictionaries with 'xfrac', 'yfrac', 'zfrac' fractional coordinates.
        If provided, these will be used for the conversion. Either atoms or frac_coords must be provided.
    frac_coords : numpy.ndarray, optional
        Nx3 array of fractional coordinates, where N is the number of atoms.
        Either atoms or frac_coords must be provided.
    box_dim : list or array, optional
        Box dimensions in one of the supported formats:
        - [lx, ly, lz] (orthogonal)
        - [lx, ly, lz, xy, xz, yz] (triclinic with tilt factors)
        - [lx, ly, lz, alpha, beta, gamma] (triclinic with angles)
        - [lx, ly, lz, 0, 0, xy, 0, xz, yz] (GROMACS format)
        Either box_dim or cell must be provided.
    cell : list or array, optional
        Cell parameters as [a, b, c, alpha, beta, gamma].
        Either box_dim or cell must be provided.
    add_to_atoms : bool, optional
        If True and atoms is provided, adds cartesian coordinates to the atom dictionaries
        as 'x', 'y', 'z'. Default is True.
        
    Returns
    -------
    cart_coords : numpy.ndarray
        Nx3 array of cartesian coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with updated 'x', 'y', 'z' fields
        if atoms was provided and add_to_atoms is True.
        
    Examples
    --------
    # Using box dimensions and atoms with fractional coordinates:
    cart_coords = ap.fract.fractional_to_cartesian(atoms, box_dim=[10, 10, 10])
    
    # Using cell parameters and explicit fractional coordinates:
    cart_coords = ap.fract.fractional_to_cartesian(
        frac_coords=frac_coords, cell=[10, 10, 10, 90, 90, 90]
    )
    
    # Get cartesian coordinates without modifying atoms:
    cart_coords = ap.fract.fractional_to_cartesian(
        atoms, box_dim=[10, 10, 10], add_to_atoms=False
    )
    """
    if (atoms is None and frac_coords is None) or (box_dim is None and cell is None):
        raise ValueError("Either (atoms or frac_coords) and (box_dim or cell) must be provided")
    
    # Get fractional coordinates from atoms if needed
    if frac_coords is None:
        frac_coords = np.array([[atom.get('xfrac', 0.0), atom.get('yfrac', 0.0), atom.get('zfrac', 0.0)] 
                               for atom in atoms])
    
    # If cell is provided but not box_dim, use Box_dim2Cell to get the box
    if box_dim is None:
        _, box = Box_dim2Cell(cell)
    else:
        # Get the box matrix from box_dim
        if len(box_dim) == 3:
            # Orthogonal box
            box = np.zeros((3, 3))
            box[0, 0] = box_dim[0]  # lx
            box[1, 1] = box_dim[1]  # ly
            box[2, 2] = box_dim[2]  # lz
        elif len(box_dim) == 6:
            # Could be [lx, ly, lz, xy, xz, yz] or [lx, ly, lz, alpha, beta, gamma]
            # Let's check if the last three values might be angles
            if all(angle > 0 and angle < 180 for angle in box_dim[3:6]):
                # Treat as [lx, ly, lz, alpha, beta, gamma]
                cell_params = box_dim
                _, box = Box_dim2Cell(cell_params)
            else:
                # Treat as [lx, ly, lz, xy, xz, yz]
                lx, ly, lz, xy, xz, yz = box_dim
                box = np.zeros((3, 3))
                box[0, 0] = lx
                box[1, 1] = ly
                box[2, 2] = lz
                box[0, 1] = xy
                box[0, 2] = xz
                box[1, 2] = yz
        elif len(box_dim) == 9:
            # GROMACS format: [lx, ly, lz, 0, 0, xy, 0, xz, yz]
            lx, ly, lz = box_dim[0], box_dim[1], box_dim[2]
            xy, xz, yz = box_dim[5], box_dim[7], box_dim[8]
            
            box = np.zeros((3, 3))
            box[0, 0] = lx
            box[1, 1] = ly
            box[2, 2] = lz
            box[0, 1] = xy
            box[0, 2] = xz
            box[1, 2] = yz
        else:
            raise ValueError(f"Invalid box_dim length: {len(box_dim)}. Expected 3, 6, or 9.")
    
    # Convert fractional to cartesian using the box matrix
    cart_coords = np.dot(frac_coords, box)
    
    # Add cartesian coordinates to atoms if requested
    if atoms is not None and add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['x'] = float(cart_coords[i, 0])
            atom['y'] = float(cart_coords[i, 1])
            atom['z'] = float(cart_coords[i, 2])
        return cart_coords, atoms
    
    return cart_coords


def wrap_coordinates(atoms, box_dim=None, cell=None, in_place=True):
    """
    Wrap atoms into the primary unit cell (0 ≤ frac < 1).
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with cartesian coordinates.
    box_dim : list or array, optional
        Box dimensions in one of the supported formats.
        Either box_dim or cell must be provided.
    cell : list or array, optional
        Cell parameters as [a, b, c, alpha, beta, gamma].
        Either box_dim or cell must be provided.
    in_place : bool, optional
        If True, modifies the input atoms list in place. Default is True.
        
    Returns
    -------
    atoms : list of dict
        The atoms list with wrapped coordinates.
        
    Examples
    --------
    # Wrap coordinates in place:
    ap.fract.wrap_coordinates(atoms, box_dim=[10, 10, 10])
    
    # Get a new list with wrapped coordinates:
    wrapped_atoms = ap.fract.wrap_coordinates(atoms, box_dim=[10, 10, 10], in_place=False)
    """
    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # Create a copy if not in_place
    if not in_place:
        atoms = [atom.copy() for atom in atoms]
    
    # Convert to fractional coordinates
    frac_coords, atoms = cartesian_to_fractional(atoms, box_dim, cell, add_to_atoms=True)
    
    # Wrap fractional coordinates to [0, 1)
    frac_coords = frac_coords % 1.0
    
    # Update atoms with wrapped fractional coordinates
    for i, atom in enumerate(atoms):
        atom['xfrac'] = float(frac_coords[i, 0])
        atom['yfrac'] = float(frac_coords[i, 1])
        atom['zfrac'] = float(frac_coords[i, 2])
    
    # Convert back to cartesian
    fractional_to_cartesian(atoms, box_dim=box_dim, cell=cell, add_to_atoms=True)
    
    return atoms


def shift_origin(atoms, shift_vec, box_dim=None, cell=None, in_place=True):
    """
    Shift all atoms by a given vector in fractional coordinates.
    
    This is useful for changing the origin of the coordinate system or
    recentering the atoms in the unit cell.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with cartesian coordinates.
    shift_vec : list or array of length 3
        The shift vector in fractional coordinates [dx, dy, dz].
    box_dim : list or array, optional
        Box dimensions in one of the supported formats.
        Either box_dim or cell must be provided.
    cell : list or array, optional
        Cell parameters as [a, b, c, alpha, beta, gamma].
        Either box_dim or cell must be provided.
    in_place : bool, optional
        If True, modifies the input atoms list in place. Default is True.
        
    Returns
    -------
    atoms : list of dict
        The atoms list with shifted coordinates.
        
    Examples
    --------
    # Shift atoms by half a unit cell in the x direction:
    ap.fract.shift_origin(atoms, [0.5, 0, 0], box_dim=[10, 10, 10])
    
    # Center atoms at the middle of the unit cell:
    ap.fract.shift_origin(atoms, [0.5, 0.5, 0.5], box_dim=[10, 10, 10])
    """
    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # Create a copy if not in_place
    if not in_place:
        atoms = [atom.copy() for atom in atoms]
    
    # Convert to fractional coordinates
    frac_coords, atoms = cartesian_to_fractional(atoms, box_dim, cell, add_to_atoms=True)
    
    # Apply shift
    shift_vec = np.array(shift_vec)
    frac_coords = frac_coords + shift_vec
    
    # Update atoms with shifted fractional coordinates
    for i, atom in enumerate(atoms):
        atom['xfrac'] = float(frac_coords[i, 0])
        atom['yfrac'] = float(frac_coords[i, 1])
        atom['zfrac'] = float(frac_coords[i, 2])
    
    # Convert back to cartesian
    fractional_to_cartesian(atoms, box_dim=box_dim, cell=cell, add_to_atoms=True)
    
    return atoms


def get_cell_vectors(box_dim=None, cell=None):
    """
    Get the three cell vectors from box dimensions or cell parameters.
    
    Parameters
    ----------
    box_dim : list or array, optional
        Box dimensions in one of the supported formats.
        Either box_dim or cell must be provided.
    cell : list or array, optional
        Cell parameters as [a, b, c, alpha, beta, gamma].
        Either box_dim or cell must be provided.
        
    Returns
    -------
    cell_vectors : numpy.ndarray
        3x3 array with the three cell vectors as rows.
        
    Examples
    --------
    # Get cell vectors from box dimensions:
    vectors = ap.fract.get_cell_vectors(box_dim=[10, 10, 10])
    
    # Get cell vectors from cell parameters:
    vectors = ap.fract.get_cell_vectors(cell=[10, 10, 10, 90, 90, 90])
    """
    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # If cell is provided but not box_dim, use Cell2Box_dim to get box_dim
    if box_dim is None:
        box_dim = Cell2Box_dim(cell)
    
    # Get the box matrix from box_dim
    if len(box_dim) == 3:
        # Orthogonal box
        box = np.zeros((3, 3))
        box[0, 0] = box_dim[0]  # lx
        box[1, 1] = box_dim[1]  # ly
        box[2, 2] = box_dim[2]  # lz
    elif len(box_dim) == 6:
        # Could be [lx, ly, lz, xy, xz, yz] or [lx, ly, lz, alpha, beta, gamma]
        # Let's check if the last three values might be angles
        if all(angle > 0 and angle < 180 for angle in box_dim[3:6]):
            # Treat as [lx, ly, lz, alpha, beta, gamma]
            cell_params = box_dim
            _, box = Box_dim2Cell(cell_params)
        else:
            # Treat as [lx, ly, lz, xy, xz, yz]
            lx, ly, lz, xy, xz, yz = box_dim
            box = np.zeros((3, 3))
            box[0, 0] = lx
            box[1, 1] = ly
            box[2, 2] = lz
            box[0, 1] = xy
            box[0, 2] = xz
            box[1, 2] = yz
    elif len(box_dim) == 9:
        # GROMACS format: [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        lx, ly, lz = box_dim[0], box_dim[1], box_dim[2]
        xy, xz, yz = box_dim[5], box_dim[7], box_dim[8]
        
        box = np.zeros((3, 3))
        box[0, 0] = lx
        box[1, 1] = ly
        box[2, 2] = lz
        box[0, 1] = xy
        box[0, 2] = xz
        box[1, 2] = yz
    else:
        raise ValueError(f"Invalid box_dim length: {len(box_dim)}. Expected 3, 6, or 9.")
    
    return box
