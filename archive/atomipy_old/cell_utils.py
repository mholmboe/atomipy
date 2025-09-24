"""
Utility functions for converting between different cell representations.

This module provides functions to convert between Box_dim (a 1D array of
the simulation box dimensions) and Cell (a 3×3 matrix representation of
the unit cell).
"""

import numpy as np


def Box_dim2Cell(box_dim):
    """
    Convert Box_dim to box matrix and cell parameters.
    
    Box_dim is a 1D array or list of box dimensions, typically in Angstroms.
    For an orthogonal box, Box_dim is [Lx, Ly, Lz].
    For a triclinic box, Box_dim is [Lx, Ly, Lz, xy, xz, yz] or
    [Lx, Ly, Lz, alpha, beta, gamma] (angles in degrees).
    
    Args:
        box_dim: A list or numpy array with box dimensions.
            Length 3: [Lx, Ly, Lz] - orthogonal box
            Length 6: [Lx, Ly, Lz, xy, xz, yz] - triclinic box with box vectors
            Length 6: [Lx, Ly, Lz, alpha, beta, gamma] - triclinic box with angles (degrees)
            Length 9: [xx, xy, xz, yx, yy, yz, zx, zy, zz] - full 3×3 matrix in row-major order
    
    Returns:
        box: A 3×3 numpy array representing the box matrix
        cell: A 1×6 numpy array with [a, b, c, alfa, beta, gamma]
    """
    box_dim = np.array(box_dim, dtype=float)
    
    if len(box_dim) == 3:
        # Orthogonal box: [Lx, Ly, Lz]
        lx, ly, lz = box_dim
        box = np.array([
            [lx, 0.0, 0.0],
            [0.0, ly, 0.0],
            [0.0, 0.0, lz]
        ])
        # Create cell parameters for orthogonal box
        cell = np.array([lx, ly, lz, 90.0, 90.0, 90.0])
    elif len(box_dim) == 6:
        # Check if the values are angles or box vectors
        if all(angle > 0 and angle < 180 for angle in box_dim[3:6]):
            # Triclinic box with angles: [Lx, Ly, Lz, alpha, beta, gamma]
            lx, ly, lz, alpha, beta, gamma = box_dim
            
            # Store cell parameters directly
            cell = np.array([lx, ly, lz, alpha, beta, gamma])
            
            # Convert angles from degrees to radians
            alpha_rad = np.radians(alpha)
            beta_rad = np.radians(beta)
            gamma_rad = np.radians(gamma)
            
            # Calculate box vectors
            cos_alpha = np.cos(alpha_rad)
            cos_beta = np.cos(beta_rad)
            cos_gamma = np.cos(gamma_rad)
            sin_gamma = np.sin(gamma_rad)
            
            box = np.zeros((3, 3))
            box[0, 0] = lx
            box[1, 0] = 0.0
            box[2, 0] = 0.0
            
            box[0, 1] = ly * cos_gamma
            box[1, 1] = ly * sin_gamma
            box[2, 1] = 0.0
            
            box[0, 2] = lz * cos_beta
            box[1, 2] = lz * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
            box[2, 2] = lz * np.sqrt(1.0 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 
                                     2.0 * cos_alpha * cos_beta * cos_gamma) / sin_gamma
        else:
            # Triclinic box with box vectors: [Lx, Ly, Lz, xy, xz, yz]
            lx, ly, lz, xy, xz, yz = box_dim
            
            cell = np.array([
                [lx, xy, xz],
                [0.0, ly, yz],
                [0.0, 0.0, lz]
            ])
    elif len(box_dim) == 9:
        # GRO 9-component format:
        # Box_dim = [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        lx = box_dim[0]
        ly = box_dim[1]
        lz = box_dim[2]
        xy = box_dim[5]  # different index than documented - based on your MATLAB code
        xz = box_dim[7]
        yz = box_dim[8]
        
        # Construct the unit cell vectors to match your MATLAB implementation
        box = np.zeros((3, 3))
        box[0, 0] = lx            # xx
        box[0, 1] = 0.0           # xy
        box[0, 2] = 0.0           # xz
        box[1, 0] = xy            # yx
        box[1, 1] = ly            # yy
        box[1, 2] = 0.0           # yz
        box[2, 0] = xz            # zx
        box[2, 1] = yz            # zy
        box[2, 2] = lz            # zz
        
        # Calculate cell parameters from box dimensions using MATLAB formulas
        a = lx
        b = np.sqrt(ly**2 + xy**2)
        c = np.sqrt(lz**2 + xz**2 + yz**2)
        
        # Calculate angles using formulas from MATLAB implementation
        cos_alfa = (ly*yz + xy*xz)/(b*c)
        cos_beta = xz/c
        cos_gamma = xy/b
        
        # Convert to degrees
        alfa = np.degrees(np.arccos(np.clip(cos_alfa, -1.0, 1.0)))
        beta = np.degrees(np.arccos(np.clip(cos_beta, -1.0, 1.0)))
        gamma = np.degrees(np.arccos(np.clip(cos_gamma, -1.0, 1.0)))
        
        cell = np.array([a, b, c, alfa, beta, gamma])
    else:
        raise ValueError(f"Invalid Box_dim length: {len(box_dim)}. Expected 3, 6, or a 9.")
        
    return cell


def Cell2Box_dim(cell, original_box_dim=None):
    """
    Convert cell parameters [a, b, c, alfa, beta, gamma] to Box_dim.
    
    Args:
        cell: A 1×6 numpy array with cell parameters [a, b, c, alfa, beta, gamma]
                    where a, b, c are lengths and alfa, beta, gamma are angles in degrees.
        original_box_dim: Optional, the original Box_dim from which the cell parameters were derived.
                        If provided, it will be used to ensure consistency in triclinic parameters.
    
    Returns:
        box_dim: A numpy array with box dimensions
    """
    cell = np.array(cell, dtype=float)
    
    if len(cell) == 3:
        cell = list(cell) + [90.0, 90.0, 90.0]
    
    if len(cell) != 6:
        raise ValueError(f"Expected 6 cell parameters, got {len(cell)}")
    
    a, b, c, alfa, beta, gamma = cell
    
    # Check if this is an orthogonal box (all angles ~90 degrees)
    if (np.isclose(alfa, 90.0) and np.isclose(beta, 90.0) and np.isclose(gamma, 90.0)):
        # Simple orthogonal box - return only [lx, ly, lz] as a 1x3 array
        return np.array([a, b, c], dtype=float)
    
    # Convert angles from degrees to radians for non-orthogonal calculations
    alfa_rad = np.radians(alfa)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    
    # Calculate box vectors according to MATLAB implementation
    lx = a
    xy = b * np.cos(gamma_rad)  # Use direct cosine, not offset by π/2
    ly = np.sqrt(b**2 - xy**2)
    xz = c * np.cos(beta_rad)   # Use direct cosine, not offset by π/2
    
    # Calculate yz using the formula from MATLAB implementation
    yz = (b * c * np.cos(alfa_rad) - xy * xz) / ly
    
    # Calculate lz
    lz = np.sqrt(c**2 - xz**2 - yz**2)
    
    # For non-orthogonal box with original dimensions provided
    if original_box_dim is not None and len(original_box_dim) == 9:
        # If original box dimensions are provided, maintain the triclinic parameters
        # but update the box lengths to match cell parameters a, b, c
        return np.array([lx, ly, lz, 0.0, 0.0, original_box_dim[5], 0.0, original_box_dim[7], original_box_dim[8]])

    # Return in GRO format: [lx, ly, lz, 0, 0, xy, 0, xz, yz]
    return np.array([lx, ly, lz, 0.0, 0.0, xy, 0.0, xz, yz])
