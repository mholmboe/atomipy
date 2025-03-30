"""
This module provides functions for replicating atomic structures along unit cell dimensions.

The main function converts atoms to fractional coordinates, replicates the structure
along the unit cell vectors, and converts back to cartesian coordinates with updated
box dimensions.
"""

import copy
import numpy as np
from .fract import cartesian_to_fractional, fractional_to_cartesian, get_cell_vectors
from .cell_utils import Box_dim2Cell, Cell2Box_dim


def replicate_cell(atoms, box_dim=None, cell=None, replicate=[1, 1, 1], keep_molid=True, 
                  keep_resname=True, renumber_index=True):
    """
    Replicates a unit cell along specified directions.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with cartesian coordinates.
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
    replicate : list or array of length 3, optional
        Number of replications in the a, b, c directions. Default is [1, 1, 1].
    keep_molid : bool, optional
        If True, keeps the original molecule IDs for all replicas. Default is True.
    keep_resname : bool, optional
        If True, keeps the original residue names for all replicas. Default is True.
    renumber_index : bool, optional
        If True, renumbers atom indices sequentially. Default is True.
        
    Returns
    -------
    new_atoms : list of dict
        The replicated atoms list with updated coordinates.
    new_box_dim : list
        The updated box dimensions for the replicated cell.
        
    Examples
    --------
    # Replicate 2x2x1:
    new_atoms, new_box = ap.replicate.replicate_cell(atoms, box_dim=[10, 10, 10], replicate=[2, 2, 1])
    
    # Replicate 2x2x2 and assign new molecule IDs:
    new_atoms, new_box = ap.replicate.replicate_cell(
        atoms, box_dim=[10, 10, 10], replicate=[2, 2, 2], keep_molid=False
    )
    """
    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # If only cell is provided, derive box_dim from it
    if box_dim is None:
        box_dim = Cell2Box_dim(cell)
    
    # Handle integer input for replicate
    if isinstance(replicate, int):
        replicate = [replicate, replicate, replicate]
    
    # Ensure replicate is a list/array of length 3
    if len(replicate) != 3:
        raise ValueError("replicate must be a list/array of length 3")
    
    # Convert to fractional coordinates
    frac_coords, atoms_with_frac = cartesian_to_fractional(atoms, box_dim=box_dim, add_to_atoms=True)
    
    # Create list to store all replicated atoms
    replicated_atoms = []
    
    # Generate all replicas
    max_index = 0
    max_molid = 0
    max_resid = 0
    
    # Find maximum values of index, molid, resid for renumbering
    for atom in atoms:
        if 'index' in atom and atom['index'] > max_index:
            max_index = atom['index']
        if 'molid' in atom and atom['molid'] > max_molid:
            max_molid = atom['molid']
        if 'resid' in atom and atom['resid'] > max_resid:
            max_resid = atom['resid']
    
    # Generate replicas
    for i in range(replicate[0]):
        for j in range(replicate[1]):
            for k in range(replicate[2]):
                # Skip the original cell if all offsets are 0
                if i == 0 and j == 0 and k == 0:
                    replica = copy.deepcopy(atoms_with_frac)
                else:
                    replica = copy.deepcopy(atoms_with_frac)
                    
                    # Update indices and IDs if needed
                    for atom in replica:
                        # Update fractional coordinates
                        atom['xfrac'] += i
                        atom['yfrac'] += j
                        atom['zfrac'] += k
                        
                        # Update atom index
                        if renumber_index and 'index' in atom:
                            atom['index'] += (i * replicate[1] * replicate[2] + 
                                             j * replicate[2] + k) * (max_index + 1)
                        
                        # Update molecule ID if not keeping original
                        if not keep_molid and 'molid' in atom:
                            atom['molid'] += (i * replicate[1] * replicate[2] + 
                                             j * replicate[2] + k) * (max_molid + 1)
                        
                        # Update residue ID (always update for proper PDB/GRO format)
                        if 'resid' in atom:
                            atom['resid'] += (i * replicate[1] * replicate[2] + 
                                             j * replicate[2] + k) * (max_resid + 1)
                        
                        # Residue name remains unchanged if keep_resname is True (default)
                
                # Add replica to the list
                replicated_atoms.extend(replica)
    
    # Scale the box dimensions
    if len(box_dim) == 3:
        # Orthogonal box
        new_box_dim = [
            box_dim[0] * replicate[0],
            box_dim[1] * replicate[1],
            box_dim[2] * replicate[2]
        ]
    elif len(box_dim) == 6:
        # Could be [lx, ly, lz, xy, xz, yz] or [lx, ly, lz, alpha, beta, gamma]
        if all(angle > 0 and angle < 180 for angle in box_dim[3:6]):
            # [lx, ly, lz, alpha, beta, gamma]
            new_box_dim = [
                box_dim[0] * replicate[0],
                box_dim[1] * replicate[1],
                box_dim[2] * replicate[2],
                box_dim[3], box_dim[4], box_dim[5]  # Angles stay the same
            ]
        else:
            # [lx, ly, lz, xy, xz, yz]
            # For triclinic boxes with tilt factors, we need to scale both lengths and tilt factors
            new_box_dim = [
                box_dim[0] * replicate[0],
                box_dim[1] * replicate[1],
                box_dim[2] * replicate[2],
                box_dim[3] * replicate[0],  # xy scales with x
                box_dim[4] * replicate[0],  # xz scales with x
                box_dim[5] * replicate[1]   # yz scales with y
            ]
    elif len(box_dim) == 9:
        # GROMACS format: [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        new_box_dim = [
            box_dim[0] * replicate[0],  # lx
            box_dim[1] * replicate[1],  # ly
            box_dim[2] * replicate[2],  # lz
            box_dim[3],                 # 0
            box_dim[4],                 # 0
            box_dim[5] * replicate[0],  # xy scales with x
            box_dim[6],                 # 0
            box_dim[7] * replicate[0],  # xz scales with x
            box_dim[8] * replicate[1]   # yz scales with y
        ]
    else:
        raise ValueError(f"Invalid box_dim length: {len(box_dim)}. Expected 3, 6, or 9.")
    
    # Convert back to cartesian coordinates
    cart_coords = fractional_to_cartesian(replicated_atoms, box_dim=new_box_dim, add_to_atoms=True)
    
    # Return the replicated atoms and new box dimensions
    return replicated_atoms, new_box_dim


def replicate_atom(atoms, box_dim=None, cell=None, replicate=[1, 1, 1], dim_order='xyz', 
                  add_molid=False, renumber_index=True):
    """
    Legacy function for compatibility with old replicate_atom functionality.
    
    This function calls replicate_cell internally but maintains the same interface
    as the original replicate_atom function.
    
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
    replicate : list or array of length 3, optional
        Number of replications in the a, b, c directions. Default is [1, 1, 1].
    dim_order : str, optional
        Order of replication (e.g., 'xyz', 'yxz'). Default is 'xyz'.
        Note: This parameter is accepted for backward compatibility but is ignored.
        The replication is always performed in all three dimensions simultaneously.
    add_molid : bool, optional
        If True, adds new molecule IDs for each replica. Default is False.
    renumber_index : bool, optional
        If True, renumbers atom indices sequentially. Default is True.
        
    Returns
    -------
    new_atoms : list of dict
        The replicated atoms list with updated coordinates.
        
    Notes
    -----
    The dim_order parameter is accepted for backward compatibility but is ignored.
    The replication is always performed in all three dimensions simultaneously.
    """
    # Call replicate_cell with the appropriate parameters
    replicated_atoms, new_box_dim = replicate_cell(
        atoms=atoms,
        box_dim=box_dim,
        cell=cell,
        replicate=replicate,
        keep_molid=not add_molid,
        keep_resname=True,
        renumber_index=renumber_index
    )
    
    # For backward compatibility, only return the atoms
    return replicated_atoms


def update_atom_indices(atoms):
    """
    Update atom indices to be sequential.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
        
    Returns
    -------
    atoms : list of dict
        The atoms list with updated indices.
    """
    for i, atom in enumerate(atoms):
        atom['index'] = i + 1
    return atoms
