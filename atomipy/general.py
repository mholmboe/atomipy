"""
General module for atomipy - provides utility functions for structure manipulation.
"""

import numpy as np
from .transform import (
    cartesian_to_fractional,
    fractional_to_cartesian,
    triclinic_to_orthogonal,
    orthogonal_to_triclinic
)
from .cell_utils import normalize_box, Cell2Box_dim


def scale(atoms, Box, scale_factors, resname=None):
    """
    Scale atom coordinates and box dimensions.
    
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with coordinates
    Box : list of float
        Box dimensions (1x3, 1x6, or 1x9).
    scale_factors : list of float or float
        Scaling factors for x, y, z dimensions. If a single value is provided,
        it will be applied to all dimensions.
    resname : str, optional
        Only scale atoms with this residue name, if provided. Default is None (all atoms).
        
    Returns
    -------
    tuple
        (scaled_atoms, scaled_box) - Atoms with scaled coordinates and new box dimensions
    
    Notes
    -----
    - For triclinic cells, the function first converts to orthogonal coordinates, 
      applies scaling, and then converts back to triclinic.
    - The function preserves the original cell angles for triclinic boxes.

    Examples
    --------
    atoms, Box_dim = ap.import_gro("structure.gro")
    atoms_scaled, Box_scaled = scale(atoms, Box_dim, 2.0)
    atoms_scaled, Box_scaled = scale(atoms, Box_dim, [1.0, 1.0, 2.0], resname="WAT")
    """
    # Make a deep copy to avoid modifying the original atoms
    atoms = [atom.copy() for atom in atoms]
    
    # Handle single scale factor
    if isinstance(scale_factors, (int, float)):
        scale_factors = [scale_factors, scale_factors, scale_factors]
    
    # Ensure scale_factors has 3 elements
    if len(scale_factors) != 3:
        raise ValueError("scale_factors must be a single number or list of 3 numbers")
    
    # Determine which atoms to scale
    if resname is not None:
        indices_to_scale = [i for i, atom in enumerate(atoms) 
                           if atom.get('resname', '').upper() == resname.upper()]
    else:
        indices_to_scale = list(range(len(atoms)))
    
    # Normalize box
    Box_dim, Cell = normalize_box(Box)
    
    # Handle triclinic box if needed
    is_triclinic = len(Box_dim) > 3 and any(Box_dim[i] != 0 for i in [5, 7, 8])
    
    if is_triclinic:
        # Convert to orthogonal coordinates
        cell_params = Cell
        ortho_atoms = triclinic_to_orthogonal(atoms=atoms, Box=Box_dim)
        
        # Scale orthogonal coordinates
        for i in indices_to_scale:
            ortho_atoms[i]['x'] *= scale_factors[0]
            ortho_atoms[i]['y'] *= scale_factors[1]
            ortho_atoms[i]['z'] *= scale_factors[2]
        
        # Scale box dimensions while preserving angles
        new_cell = cell_params.copy()
        new_cell[0] *= scale_factors[0]  # a
        new_cell[1] *= scale_factors[1]  # b
        new_cell[2] *= scale_factors[2]  # c
        new_box = Cell2Box_dim(new_cell)
        
        # Convert back to triclinic coordinates
        scaled_atoms = orthogonal_to_triclinic(atoms=ortho_atoms, Box=new_box)
    else:
        # Direct scaling for orthogonal box
        for i in indices_to_scale:
            atoms[i]['x'] *= scale_factors[0]
            atoms[i]['y'] *= scale_factors[1]
            atoms[i]['z'] *= scale_factors[2]
        
        # Scale box dimensions
        new_box = Box_dim.copy() if hasattr(Box_dim, 'copy') else list(Box_dim)
        for i in range(min(3, len(Box_dim))):
            new_box[i] *= scale_factors[i]
            
        scaled_atoms = atoms
    
    return scaled_atoms, new_box
