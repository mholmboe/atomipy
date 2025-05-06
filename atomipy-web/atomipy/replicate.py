"""
This module provides functions for replicating atomic structures along unit cell dimensions.

The main function converts atoms to fractional coordinates, replicates the structure
along the unit cell vectors, and converts back to cartesian coordinates with updated
box dimensions.
"""

import copy
import numpy as np
from .transform import (
    cartesian_to_fractional, fractional_to_cartesian, get_cell_vectors,
    direct_fractional_to_cartesian
)
from .cell_utils import Box_dim2Cell, Cell2Box_dim
from . import write_conf # Import for debug writing

def replicate_system(atoms, box, replicate=[1, 1, 1], keep_molid=True, 
                   keep_resname=True, renumber_index=True):
    """
    Replicates a unit cell along specified directions using crystallographic approach.
    
    This function implements the standard crystallographic replication process:
    1. Convert input coordinates to fractional coordinates (unit cube)
    2. Stack the unit cell n1×n2×n3 times in a,b,c directions
    3. Scale the cell parameters accordingly (preserving angles)
    4. Convert fractional coordinates back to cartesian
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with cartesian coordinates.
    box : a 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
            a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
            a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
            Note that for orthogonal boxes Cell = Box_dim.
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
    new_cell : list
        The replicated cell parameters as [a, b, c, alpha, beta, gamma].
        
    Examples
    --------
    # Replicate 2x2x1 using cell parameters:
    new_atoms, new_box, new_cell = ap.replicate.replicate_system(
        atoms, box=[10, 10, 10, 90, 90, 90], replicate=[2, 2, 1]
    )
    
    # Replicate 2x2x2 with box dimensions and assign new molecule IDs:
    new_atoms, new_box, new_cell = ap.replicate.replicate_system(
        atoms, box=[10, 10, 10], replicate=[2, 2, 2], keep_molid=False
    )
    """

        # Possibly convert Box_dim into [a,b,c,alpha,beta,gamma] form
    if box is not None:
        if len(box) == 9:
            box_dim = box 
            # Convert from Box_dim format to Cell format
            cell = Box_dim2Cell(box_dim)
        elif len(box) == 6:
            cell = box
        elif len(box) == 3:
            # Orthogonal box
            cell = list(box) + [90.0, 90.0, 90.0]
            box_dim = box

    if box_dim is None and cell is None:
        raise ValueError("Either box_dim or cell must be provided")
    
    # Handle integer input for replicate
    if isinstance(replicate, int):
        replicate = [replicate, replicate, replicate]
    
    # Ensure replicate is a list/array of length 3
    if len(replicate) != 3:
        raise ValueError("replicate must be a list/array of length 3")
    
    # Calculate box dimensions from cell parameters (for coordinate conversion)
    box_dim = Cell2Box_dim(cell)
    
    # Step 2: Convert to fractional coordinates (unit cube)
    # Using the new fract.py module which converts through orthogonal coordinates
    frac_coords, atoms_with_frac = cartesian_to_fractional(atoms, box=cell, add_to_atoms=True)
    
    # Step 3: Create storage for replicated atoms
    replicated_atoms = []
    
    # Find maximum values for renumbering
    max_index = max([atom.get('index', 0) for atom in atoms], default=0)
    max_molid = max([atom.get('molid', 0) for atom in atoms], default=0)
    max_resid = max([atom.get('resid', 0) for atom in atoms], default=0)
    
    # Store original cartesian coordinates for each atom (to handle special cases)
    original_coords = {}
    for i, atom in enumerate(atoms):
        original_coords[i] = {
            'x': atom.get('x', 0.0),
            'y': atom.get('y', 0.0),
            'z': atom.get('z', 0.0)
        }
    
    # --- Step 4: Sequential Replication --- 
    current_replication = copy.deepcopy(atoms_with_frac)
    num_atoms_base = len(atoms)
    max_idx_base = max_index
    max_molid_base = max_molid
    max_resid_base = max_resid
    
    # Stage 1: Replicate along X (a-vector)
    if replicate[0] > 1:
        replicated_x = []
        atoms_to_replicate_x = copy.deepcopy(current_replication)
        num_atoms_stage = len(atoms_to_replicate_x)
        max_idx_stage = max((atom.get('index', 0) for atom in atoms_to_replicate_x), default=-1)
        max_molid_stage = max((atom.get('molid', 0) for atom in atoms_to_replicate_x), default=-1)
        max_resid_stage = max((atom.get('resid', 0) for atom in atoms_to_replicate_x), default=-1)
        
        for i in range(replicate[0]):
            replica = copy.deepcopy(atoms_to_replicate_x)
            offset_idx = i * (max_idx_stage + 1)
            offset_molid = i * (max_molid_stage + 1)
            offset_resid = i * (max_resid_stage + 1)
            
            for atom in replica:
                atom['xfrac'] += i
                if renumber_index and 'index' in atom:
                    atom['index'] += offset_idx
                if not keep_molid and 'molid' in atom:
                    atom['molid'] += offset_molid
                if 'resid' in atom:
                    atom['resid'] += offset_resid
            replicated_x.extend(replica)
        current_replication = replicated_x

    # Stage 2: Replicate along Y (b-vector)
    if replicate[1] > 1:
        replicated_xy = []
        atoms_to_replicate_y = copy.deepcopy(current_replication) # Result from Stage 1
        num_atoms_stage = len(atoms_to_replicate_y)
        max_idx_stage = max((atom.get('index', 0) for atom in atoms_to_replicate_y), default=-1)
        max_molid_stage = max((atom.get('molid', 0) for atom in atoms_to_replicate_y), default=-1)
        max_resid_stage = max((atom.get('resid', 0) for atom in atoms_to_replicate_y), default=-1)

        for j in range(replicate[1]):
            replica = copy.deepcopy(atoms_to_replicate_y)
            offset_idx = j * (max_idx_stage + 1)
            offset_molid = j * (max_molid_stage + 1)
            offset_resid = j * (max_resid_stage + 1)
            
            for atom in replica:
                atom['yfrac'] += j
                if renumber_index and 'index' in atom:
                    atom['index'] += offset_idx
                if not keep_molid and 'molid' in atom:
                    atom['molid'] += offset_molid
                if 'resid' in atom:
                    atom['resid'] += offset_resid
            replicated_xy.extend(replica)
        current_replication = replicated_xy

    # Stage 3: Replicate along Z (c-vector)
    if replicate[2] > 1:
        replicated_xyz = []
        atoms_to_replicate_z = copy.deepcopy(current_replication) # Result from Stage 2
        num_atoms_stage = len(atoms_to_replicate_z)
        max_idx_stage = max((atom.get('index', 0) for atom in atoms_to_replicate_z), default=-1)
        max_molid_stage = max((atom.get('molid', 0) for atom in atoms_to_replicate_z), default=-1)
        max_resid_stage = max((atom.get('resid', 0) for atom in atoms_to_replicate_z), default=-1)

        for k in range(replicate[2]):
            replica = copy.deepcopy(atoms_to_replicate_z)
            offset_idx = k * (max_idx_stage + 1)
            offset_molid = k * (max_molid_stage + 1)
            offset_resid = k * (max_resid_stage + 1)
            
            for atom in replica:
                atom['zfrac'] += k
                if renumber_index and 'index' in atom:
                    atom['index'] += offset_idx
                if not keep_molid and 'molid' in atom:
                    atom['molid'] += offset_molid
                if 'resid' in atom:
                    atom['resid'] += offset_resid
            replicated_xyz.extend(replica)
        current_replication = replicated_xyz

    replicated_atoms = current_replication # Final result after all stages
    # --------------------------------------
    
    # Step 5: Create replicated cell by scaling the original cell parameters
    new_cell = [
        cell[0] * replicate[0],  # a - scale by replication in x
        cell[1] * replicate[1],  # b - scale by replication in y
        cell[2] * replicate[2],  # c - scale by replication in z
        cell[3],                 # alpha - angles remain unchanged
        cell[4],                 # beta - angles remain unchanged
        cell[5]                  # gamma - angles remain unchanged
    ]
    
    # Step 6: Generate new box dimensions from the replicated cell
    new_box_dim = Cell2Box_dim(new_cell)
    
    # Step 7: For triclinic cells, we need to be careful with the coordinate transformation
    # to preserve atomic planes
    
    # First, convert the fractional coordinates to integers (which unit cell they belong to)
    # and offsets within the unit cell (0-1 range)
    for atom in replicated_atoms:
        # Calculate which unit cell this atom belongs to in each direction
        ix = int(atom['xfrac'])
        iy = int(atom['yfrac'])
        iz = int(atom['zfrac'])
        
        # Calculate the fractional offset within that unit cell (0-1 range)
        x_offset = atom['xfrac'] - ix
        y_offset = atom['yfrac'] - iy
        z_offset = atom['zfrac'] - iz
        
        # For proper replication that preserves planes, we calculate the new position
        # using the unit cell indices and the original fractional coordinates
        atom['xfrac'] = (ix + x_offset) / replicate[0]
        atom['yfrac'] = (iy + y_offset) / replicate[1] 
        atom['zfrac'] = (iz + z_offset) / replicate[2]
    
    # Step 8: Convert replicated atoms back to cartesian coordinates
    # using the new (scaled) cell parameters directly for crystallographic accuracy
    cart_coords = direct_fractional_to_cartesian(replicated_atoms, cell=new_cell, add_to_atoms=True)
    
    # Clean up temporary attributes used during replication
    for atom in replicated_atoms:
        if '_original_idx' in atom:
            del atom['_original_idx']
    
    # As a final check, recalculate cell parameters from box dimensions to ensure consistency
    # This is redundant but serves as a sanity check
    new_cell = Box_dim2Cell(new_box_dim)
    
    # Round box dimensions to 5 decimal places
    new_box_dim = [round(float(val), 5) for val in new_box_dim]
    
    # Calculate new cell parameters from new box dimensions
    new_cell = Box_dim2Cell(new_box_dim)
    
    # Return the replicated atoms, new box dimensions, and new cell parameters
    return replicated_atoms, new_box_dim, new_cell


def replicate_atom(atoms, box_dim=None, cell=None, replicate=[1, 1, 1], dim_order='xyz', 
                  add_molid=False, renumber_index=True):
    """
    Legacy function for compatibility with old replicate_atom functionality.
    
    This function calls replicate_system internally but maintains the same interface
    as the original replicate_atom function.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with cartesian coordinates.
    box : a 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
            a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
            a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
            Note that for orthogonal boxes Cell = Box_dim.
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
    # Call replicate_system with the appropriate parameters
    replicated_atoms, new_box_dim = replicate_system(
        atoms=atoms,
        box_dim=box_dim,
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
