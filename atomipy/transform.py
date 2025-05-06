#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transform module for coordinate transformations in atomipy.

This module consolidates functionality from the previous fract.py, triclinic.py, and ortho.py
modules into a unified interface for handling all coordinate transformations in molecular systems.

Key functionality includes:

1. Cartesian-fractional coordinate conversions:
   - cartesian_to_fractional(): Convert Cartesian coordinates to fractional coordinates
   - fractional_to_cartesian(): Convert fractional coordinates to Cartesian coordinates
   - direct_cartesian_to_fractional(): Direct conversion using crystallographic matrices
   - direct_fractional_to_cartesian(): Direct conversion using crystallographic matrices

2. Triclinic-orthogonal transformations:
   - triclinic_to_orthogonal(): Convert triclinic coordinates to orthogonal system
   - orthogonal_to_triclinic(): Convert orthogonal coordinates to triclinic system

3. Utility functions:
   - wrap_coordinates(): Ensure coordinates are within the primary unit cell
   - get_orthogonal_box(): Get orthogonal box dimensions from triclinic parameters
   - get_cell_vectors(): Calculate cell vectors from box parameters

All functions support both direct atom dictionary input and numpy array input,
making them flexible for various use cases.
"""

import numpy as np
from .cell_utils import Box_dim2Cell, Cell2Box_dim


def cartesian_to_fractional(atoms=None, cart_coords=None, box_dim=None, box=None, add_to_atoms=True):
    """
    Convert cartesian coordinates to fractional coordinates.
    
    Parameters
    ----------
    atoms : list of dict, optional
        List of atom dictionaries with 'x', 'y', 'z' cartesian coordinates.
        If provided, these will be used for the conversion. Either atoms or cart_coords must be provided.
    cart_coords : numpy.ndarray, optional
        Nx3 array of cartesian coordinates, where N is the number of atoms.
        Either atoms or cart_coords must be provided.
    box_dim : list or array, optional
        Box dimensions with 3, 6, or 9 parameters. 3 parameters for orthogonal box [lx, ly, lz].
        6 parameters for triclinic box [lx, ly, lz, xy, xz, yz]. 9 parameters for full description
        [lx, ly, lz, 0, 0, xy, 0, xz, yz]. Either box_dim or box must be provided.
    box : list or array, optional
        Box parameters in the format [a, b, c, alpha, beta, gamma].
        Either box_dim or box must be provided.
    add_to_atoms : bool, optional
        If True and atoms is provided, adds fractional coordinates to the atom dictionaries
        as 'xfrac', 'yfrac', 'zfrac'. Default is True.
        
    Returns
    -------
    frac_coords : numpy.ndarray
        Nx3 array of fractional coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with updated fractional coordinate fields
        if add_to_atoms is True.
    """
    if (atoms is None and cart_coords is None) or (box_dim is None and box is None):
        raise ValueError("Either (atoms or cart_coords) and (box_dim or box) must be provided")
    
    # Get cartesian coordinates from atoms if needed
    if cart_coords is None:
        cart_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
    
    # Get box parameters if not provided
    if box is None:
        box, _ = Box_dim2Cell(box_dim[:6] if len(box_dim) >= 6 else box_dim)
    
    # Extract box parameters
    a, b, c, alpha, beta, gamma = box
    
    # Get cell vectors
    cell_vectors = get_cell_vectors(box)
    
    # Compute the transformation matrix
    # This is the inverse of the matrix used in fractional_to_cartesian
    det = np.dot(cell_vectors[0], np.cross(cell_vectors[1], cell_vectors[2]))
    inv_mat = np.zeros((3, 3))
    inv_mat[0] = np.cross(cell_vectors[1], cell_vectors[2]) / det
    inv_mat[1] = np.cross(cell_vectors[2], cell_vectors[0]) / det
    inv_mat[2] = np.cross(cell_vectors[0], cell_vectors[1]) / det
    
    # Apply transformation to each point
    frac_coords = np.zeros_like(cart_coords)
    for i, cart in enumerate(cart_coords):
        frac_coords[i] = np.dot(inv_mat, cart)
    
    # Add fractional coordinates to atoms if requested
    if atoms is not None and add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['xfrac'] = float(round(frac_coords[i, 0], 6))
            atom['yfrac'] = float(round(frac_coords[i, 1], 6))
            atom['zfrac'] = float(round(frac_coords[i, 2], 6))
        return frac_coords, atoms
    
    return frac_coords


def fractional_to_cartesian(atoms=None, frac_coords=None, box_dim=None, box=None, add_to_atoms=True):
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
        Box dimensions with 3, 6, or 9 parameters. 3 parameters for orthogonal box [lx, ly, lz].
        6 parameters for triclinic box [lx, ly, lz, xy, xz, yz]. 9 parameters for full description
        [lx, ly, lz, 0, 0, xy, 0, xz, yz]. Either box_dim or box must be provided.
    box : list or array, optional
        Box parameters in the format [a, b, c, alpha, beta, gamma].
        Either box_dim or box must be provided.
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
    if (atoms is None and frac_coords is None) or (box_dim is None and box is None):
        raise ValueError("Either (atoms or frac_coords) and (box_dim or box) must be provided")
    
    # Get box parameters if not provided
    if box is None:
        box, _ = Box_dim2Cell(box_dim[:6] if len(box_dim) >= 6 else box_dim)
    
    # Get fractional coordinates from atoms if needed
    if frac_coords is None:
        frac_coords = np.array([[atom.get('xfrac', 0.0), 
                                atom.get('yfrac', 0.0), 
                                atom.get('zfrac', 0.0)] 
                                for atom in atoms])
    
    # Get cell vectors
    cell_vectors = get_cell_vectors(box)
    
    # Apply transformation to each point
    cart_coords = np.zeros_like(frac_coords)
    for i, frac in enumerate(frac_coords):
        cart_coords[i] = (cell_vectors[0] * frac[0] +
                        cell_vectors[1] * frac[1] +
                        cell_vectors[2] * frac[2])
    
    # Add cartesian coordinates to atoms if requested
    if atoms is not None and add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['x'] = float(round(cart_coords[i, 0], 6))
            atom['y'] = float(round(cart_coords[i, 1], 6))
            atom['z'] = float(round(cart_coords[i, 2], 6))
        return cart_coords, atoms
    
    return cart_coords


def wrap_coordinates(atoms=None, coords=None, frac_coords=None, box_dim=None, box=None,
                    add_to_atoms=True, return_type='fractional'):
    """
    Wrap coordinates to ensure they are within the primary unit cell (0 to 1 in fractional coordinates).
    
    Parameters
    ----------
    atoms : list of dict, optional
        List of atom dictionaries with coordinate information.
        If provided along with add_to_atoms=True, the wrapped coordinates will be added to atoms.
    coords : numpy.ndarray, optional
        Nx3 array of cartesian coordinates to wrap.
    frac_coords : numpy.ndarray, optional
        Nx3 array of fractional coordinates to wrap. If provided, coords is ignored.
    box_dim : list or array, optional
        Box dimensions with 3, 6, or 9 parameters. Required if coords is provided.
    box : list or array, optional
        Box parameters in the format [a, b, c, alpha, beta, gamma]. Required if coords is provided.
    add_to_atoms : bool, optional
        If True and atoms is provided, updates the atom dictionaries with wrapped coordinates.
    return_type : str, optional
        Specifies the type of coordinates to return: 'fractional' (default) or 'cartesian'.
        
    Returns
    -------
    wrapped_coords : numpy.ndarray
        Nx3 array of wrapped coordinates in the specified return_type.
    atoms : list of dict, optional
        The original atoms list with updated coordinate fields if add_to_atoms is True.
    """
    if atoms is None and coords is None and frac_coords is None:
        raise ValueError("Either atoms, coords, or frac_coords must be provided")
    
    # Step 1: Get fractional coordinates
    if frac_coords is None:
        if coords is not None:
            # Convert cartesian to fractional
            frac_coords = cartesian_to_fractional(cart_coords=coords, box_dim=box_dim, box=box)
        elif atoms is not None:
            # Check if atoms already have fractional coordinates
            if all('xfrac' in atom for atom in atoms):
                frac_coords = np.array([[atom['xfrac'], atom['yfrac'], atom['zfrac']] for atom in atoms])
            else:
                # Extract cartesian coordinates and convert to fractional
                cart_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
                frac_coords = cartesian_to_fractional(cart_coords=cart_coords, box_dim=box_dim, box=box)
    
    # Step 2: Perform the wrapping operation on fractional coordinates
    wrapped_frac = frac_coords.copy()
    wrapped_frac = wrapped_frac % 1.0  # Simple modulo operation for primary cell wrapping
    
    # Step 3: Update atoms if requested
    if atoms is not None and add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['xfrac'] = float(round(wrapped_frac[i, 0], 6))
            atom['yfrac'] = float(round(wrapped_frac[i, 1], 6))
            atom['zfrac'] = float(round(wrapped_frac[i, 2], 6))
    
    # Step 4: Return coordinates in the requested format
    if return_type.lower() == 'cartesian':
        # Convert wrapped fractional coordinates back to cartesian
        wrapped_cart = fractional_to_cartesian(frac_coords=wrapped_frac, box_dim=box_dim, box=box)
        if atoms is not None and add_to_atoms:
            for i, atom in enumerate(atoms):
                atom['x'] = float(round(wrapped_cart[i, 0], 6))
                atom['y'] = float(round(wrapped_cart[i, 1], 6))
                atom['z'] = float(round(wrapped_cart[i, 2], 6))
            return wrapped_cart, atoms
        return wrapped_cart
    else:  # fractional
        if atoms is not None and add_to_atoms:
            return wrapped_frac, atoms
        return wrapped_frac


def triclinic_to_orthogonal(atoms=None, coords=None, box_dim=None, box=None, add_to_atoms=True):
    """
    Convert coordinates from triclinic to orthogonal representation.
    
    Parameters
    ----------
    atoms : list of dict, optional
        List of atom dictionaries with position information.
        If provided along with add_to_atoms=True, orthogonal coordinates
        will be added to atoms as 'x_ortho', 'y_ortho', 'z_ortho'.
    coords : numpy.ndarray, optional
        Nx3 array of triclinic coordinates.
    box_dim : list or array, optional
        Box dimensions as [Lx, Ly, Lz] or [Lx, Ly, Lz, xy, xz, yz].
    box : list or array, optional
        Box parameters as [a, b, c, alpha, beta, gamma].
    add_to_atoms : bool, optional
        If True and atoms is provided, adds orthogonal coordinates to the atom dictionaries.
        
    Returns
    -------
    ortho_coords : numpy.ndarray
        Nx3 array of orthogonal coordinates.
    atoms : list of dict, optional
        The original atoms list with added orthogonal coordinate fields
        if add_to_atoms is True.
    ortho_box : array
        The orthogonal box dimensions [a', b', c'].
    """
    if (atoms is None and coords is None) or (box_dim is None and box is None):
        raise ValueError("Either (atoms or coords) and (box_dim or box) must be provided")
    
    # Get the coordinate array from atoms if needed
    if coords is None and atoms is not None:
        coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
    
    # Convert box_dim to box parameters if needed
    if box is None:
        box, _ = Box_dim2Cell(box_dim)
    
    # Extract box parameters
    a, b, c, alpha, beta, gamma = box
    
    # Convert to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Calculate orthogonal box dimensions
    a_ortho = a
    b_ortho = b * np.sin(gamma_rad)
    c_ortho = c * np.sin(beta_rad)
    
    # Build transformation matrix
    transform_matrix = np.array([
        [1, np.cos(gamma_rad), np.cos(beta_rad)],
        [0, np.sin(gamma_rad), (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)],
        [0, 0, np.sqrt(1 - np.cos(beta_rad)**2 - ((np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad))**2)]
    ])
    
    # Scale transformation matrix by box dimensions
    scale_matrix = np.diag([a, b, c])
    transform_matrix = np.dot(transform_matrix, scale_matrix)
    
    # Apply transformation to each coordinate
    ortho_coords = np.zeros_like(coords)
    for i, coord in enumerate(coords):
        ortho_coords[i] = np.dot(transform_matrix, coord)
    
    # Update atoms if requested
    if atoms is not None and add_to_atoms:
        if ortho_coords.shape[0] == len(atoms):
            for i, atom in enumerate(atoms):
                atom['x_ortho'] = float(round(ortho_coords[i, 0], 6))
                atom['y_ortho'] = float(round(ortho_coords[i, 1], 6))
                atom['z_ortho'] = float(round(ortho_coords[i, 2], 6))
            return ortho_coords, atoms, np.array([a_ortho, b_ortho, c_ortho])
    
    return ortho_coords, np.array([a_ortho, b_ortho, c_ortho])


def orthogonal_to_triclinic(ortho_coords, box, atoms=None, add_to_atoms=True):
    """
    Convert coordinates from orthogonal to triclinic representation.
    
    Parameters
    ----------
    ortho_coords : numpy.ndarray
        Nx3 array of orthogonal coordinates.
    box : list or array
        Box parameters as [a, b, c, alpha, beta, gamma].
    atoms : list of dict, optional
        List of atom dictionaries, if provided and add_to_atoms is True, triclinic
        coordinates will be added to atoms as 'x', 'y', 'z'.
    add_to_atoms : bool, optional
        If True and atoms is provided, adds triclinic coordinates to the atom dictionaries.
        
    Returns
    -------
    tri_coords : numpy.ndarray
        Nx3 array of triclinic coordinates.
    atoms : list of dict, optional
        The original atoms list with updated triclinic coordinate fields
        if add_to_atoms is True.
    """
    if ortho_coords is None or box is None:
        raise ValueError("Both ortho_coords and box must be provided")
    
    # Extract box parameters
    a, b, c, alpha, beta, gamma = box
    
    # Convert to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Build inverse transformation matrix
    transform_matrix = np.array([
        [1, np.cos(gamma_rad), np.cos(beta_rad)],
        [0, np.sin(gamma_rad), (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)],
        [0, 0, np.sqrt(1 - np.cos(beta_rad)**2 - ((np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad))**2)]
    ])
    
    # Scale transformation matrix by box dimensions
    scale_matrix = np.diag([a, b, c])
    transform_matrix = np.dot(transform_matrix, scale_matrix)
    
    # Invert the transformation matrix
    inv_transform = np.linalg.inv(transform_matrix)
    
    # Apply inverse transformation to each coordinate
    cart_coords = np.zeros_like(ortho_coords)
    for i, coord in enumerate(ortho_coords):
        cart_coords[i] = np.dot(inv_transform, coord)
    
    # Update atoms if requested
    if atoms is not None and add_to_atoms:
      for i, atom in enumerate(atoms):
        atom['x'] = float(round(cart_coords[i, 0], 6))
        atom['y'] = float(round(cart_coords[i, 1], 6))
        atom['z'] = float(round(cart_coords[i, 2], 6))
      return cart_coords, atoms
    
    return cart_coords


def get_orthogonal_box(box_dim=None, box=None):
    """
    Get the dimensions of the orthogonal box representing 
    the triclinic cell.
    
    Parameters
    ----------
    box_dim : list or array, optional
        Box dimensions with 3, 6, or 9 parameters.
    box : list or array, optional
        Box parameters in the format [a, b, c, alpha, beta, gamma].
        
    Returns
    -------
    ortho_box : numpy.ndarray
        Orthogonal box dimensions [a', b', c'].
    """
    if box_dim is None and box is None:
        raise ValueError("Either box_dim or box must be provided")
    
    # Convert box_dim to box parameters if needed
    if box is None:
        box, _ = Box_dim2Cell(box_dim)
    
    # Extract box parameters
    a, b, c, alpha, beta, gamma = box
    
    # Convert to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Calculate orthogonal box dimensions
    a_ortho = a
    b_ortho = b * np.sin(gamma_rad)
    c_ortho = c * np.sin(beta_rad)
    
    return np.array([a_ortho, b_ortho, c_ortho])


def get_cell_vectors(box):
    """
    Calculate cell vectors from box parameters.
    
    Parameters
    ----------
    box : list or array
        Box parameters in the format [a, b, c, alpha, beta, gamma].
        
    Returns
    -------
    cell_vectors : numpy.ndarray
        3x3 array of cell vectors, where each row is a unit cell vector.
    """
    a, b, c, alpha, beta, gamma = box
    
    # Convert angles to radians
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Calculate the cell vectors
    v1 = np.array([a, 0, 0])
    v2 = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
    
    # Calculate the third vector components
    cx = c * np.cos(beta_rad)
    cy = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
    cz = c * np.sqrt(1 - np.cos(beta_rad)**2 - ((np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad))**2)
    
    v3 = np.array([cx, cy, cz])
    
    return np.array([v1, v2, v3])


def direct_cartesian_to_fractional(atoms=None, cart_coords=None, box_dim=None, cell=None, add_to_atoms=True):
    """
    Direct conversion from cartesian coordinates to fractional coordinates.
    This function provides a direct implementation that follows the MATLAB approach
    without intermediate orthogonalization steps.
    
    Parameters
    ----------
    atoms : list of dict, optional
        List of atom dictionaries with 'x', 'y', 'z' cartesian coordinates.
        If provided, these will be used for the conversion. Either atoms or cart_coords must be provided.
    cart_coords : numpy.ndarray, optional
        Nx3 array of cartesian coordinates, where N is the number of atoms.
        Either atoms or cart_coords must be provided.
    box_dim : list or array, optional
        Box dimensions in the format [lx, ly, lz] or [lx, ly, lz, 0, 0, xy, 0, xz, yz].
        Either box_dim or cell must be provided.
    cell : list or array, optional
        Cell parameters in the format [a, b, c, alpha, beta, gamma].
        Either box_dim or cell must be provided.
    add_to_atoms : bool, optional
        If True and atoms is provided, adds fractional coordinates to the atom dictionaries
        as 'xfrac', 'yfrac', 'zfrac'. Default is True.
        
    Returns
    -------
    frac_coords : numpy.ndarray
        Nx3 array of fractional coordinates, where N is the number of atoms.
    atoms : list of dict, optional
        The original atoms list with updated fractional coordinate fields
        if add_to_atoms is True.
    """
    if (atoms is None and cart_coords is None) or (box_dim is None and cell is None):
        raise ValueError("Either (atoms or cart_coords) and (box_dim or cell) must be provided")
    
    # If cell is provided but not box_dim, calculate box_dim
    if box_dim is None:
        box_dim, _ = Cell2Box_dim(cell)
    
    # If cell is not provided, calculate from box_dim
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
    
    # Build transformation matrix from cartesian to fractional (ToFrac in MATLAB)
    to_frac = np.array([
        [1/a, -cos_gamma / (a * sin_gamma), 
         (cos_alpha * cos_gamma - cos_beta) / (a * v * sin_gamma)],
        [0, 1 / (b * sin_gamma), 
         (cos_beta * cos_gamma - cos_alpha) / (b * v * sin_gamma)],
        [0, 0, sin_gamma / (c * v)]
    ])
    
    # Get cartesian coordinates from atoms if needed
    if cart_coords is None:
        cart_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
    
    # Apply transformation to each point
    frac_coords = np.zeros_like(cart_coords)
    for i, cart in enumerate(cart_coords):
        frac_coords[i] = np.dot(to_frac, cart)
    
    # Add fractional coordinates to atoms if requested
    if atoms is not None and add_to_atoms:
        for i, atom in enumerate(atoms):
            atom['xfrac'] = float(round(frac_coords[i, 0], 6))
            atom['yfrac'] = float(round(frac_coords[i, 1], 6))
            atom['zfrac'] = float(round(frac_coords[i, 2], 6))
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
    box_dim : list or array, optional
        Box dimensions in the format [lx, ly, lz] or [lx, ly, lz, 0, 0, xy, 0, xz, yz].
        Either box_dim or cell must be provided.
    cell : list or array, optional
        Cell parameters in the format [a, b, c, alpha, beta, gamma].
        Either box_dim or cell must be provided.
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
