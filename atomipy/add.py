"""
Add module for atomipy - provides functions for combining and updating atom structures.


This module contains functions to update atom indices and combine multiple atom structures.
"""

import copy
import numpy as np


def update(*atoms_list, molid=None, use_resname=True):
    """
    Update atom indices and optionally combine multiple atom structures.
    
    This function serves several purposes:
    1. When called with a single atoms structure, it updates all indices to be consecutive
       and assigns molecule IDs based on both molid and resname boundaries
    2. When called with multiple atoms structures, it combines them into one structure 
       with consecutive indices and molecule IDs
    3. It maintains field/attribute consistency across all atom dictionaries
    
    Parameters
    ----------
    *atoms_list : variable length argument list of atom structures
        One or more lists of atom dictionaries
    molid : int, optional
        If provided, sets all molecule IDs to this value
    use_resname : bool, optional
        If True, molecule boundaries are also determined by changes in residue name.
        Default is True.
        
    Returns
    -------
    atoms : list of dict
        Combined atoms list with updated indices and molecule IDs, and consistently 
        ordered attributes
        
    Examples
    --------
    # Update indices of a single structure:
    new_atoms = ap.update(atoms)
    
    # Combine multiple structures:
    new_atoms = ap.update(atoms1, atoms2, atoms3)
    
    # Combine structures and set specific molecule ID:
    new_atoms = ap.update(atoms1, atoms2, molid=5)
    
    # Update structure without using residue names for molecule boundaries:
    new_atoms = ap.update(atoms, use_resname=False)
    """
    # Make deep copies to avoid modifying originals
    atoms_copies = [copy.deepcopy(atoms) for atoms in atoms_list if atoms]
    
    # Handle case with no input or all empty inputs
    if not atoms_copies:
        return []
    
    # Ensure field consistency across all atom structures
    # Find common fields across all structures
    all_fields = set(atoms_copies[0][0].keys()) if atoms_copies[0] else set()
    for atoms in atoms_copies[1:]:
        if atoms:  # Skip empty structures
            all_fields = all_fields.intersection(atoms[0].keys())
    
    # Remove fields not common to all structures
    for i, atoms in enumerate(atoms_copies):
        if not atoms:  # Skip empty structures
            continue
        for j, atom in enumerate(atoms):
            # Keep only common fields
            atoms_copies[i][j] = {k: v for k, v in atom.items() if k in all_fields}
    
    # Handle single atom structure case
    if len(atoms_copies) == 1:
        return _update_single_structure(atoms_copies[0], molid, use_resname)
    
    # Combine multiple atom structures
    result_atoms = []
    current_molid = 1
    
    # Process and combine all structures
    for i, atoms in enumerate(atoms_copies):
        if not atoms:  # Skip empty structures
            continue
            
        # For the first structure, just update it
        if not result_atoms:
            result_atoms = _update_single_structure(atoms, current_molid, use_resname)
            # Get the highest molid for the next structure
            current_molid = max(atom['molid'] for atom in result_atoms) + 1
            continue
        
        # Update the current structure
        updated_atoms = _update_single_structure(atoms, None, use_resname)
        
        # Adjust molids to avoid conflicts
        # Check if the structure has a single molid or multiple molids
        unique_molids = set(atom['molid'] for atom in updated_atoms)
        
        if len(unique_molids) == 1:
            # If just one molid, simply set all to the current molid
            for atom in updated_atoms:
                atom['molid'] = current_molid
        else:
            # If multiple molids, preserve their relative relationships
            min_molid = min(unique_molids)
            offset = current_molid - min_molid
            
            for atom in updated_atoms:
                atom['molid'] += offset
        
        # Append the updated structure
        result_atoms.extend(updated_atoms)
        
        # Update the current molid for the next structure
        current_molid = max(atom['molid'] for atom in result_atoms) + 1
    
    # Update indices to be consecutive
    for i, atom in enumerate(result_atoms):
        atom['index'] = i + 1
    
    # If a specific molid was provided, set all to that value
    if molid is not None:
        for atom in result_atoms:
            atom['molid'] = molid
    
    # Order attributes consistently
    result_atoms = order_attributes(result_atoms)
    
    return result_atoms


def _update_single_structure(atoms, molid=None, use_resname=True):
    """
    Helper function to update a single atom structure.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries to update
    molid : int, optional
        If provided, sets all molecule IDs to this value
    use_resname : bool, optional
        If True, molecule boundaries are also determined by changes in residue name
        
    Returns
    -------
    atoms : list of dict
        Updated atoms list
    """
    if not atoms:
        return []
    
    # Make a deep copy to avoid modifying the original
    atoms = copy.deepcopy(atoms)
    
    # Update indices
    for i, atom in enumerate(atoms):
        atom['index'] = i + 1
    
    # If a specific molid was provided, just set all to that value
    if molid is not None:
        for atom in atoms:
            atom['molid'] = molid
        return atoms
    
    # Make sure all atoms have a molid
    for i, atom in enumerate(atoms):
        if 'molid' not in atom:
            # If no molids exist, assign sequential molids
            atom['molid'] = i + 1
    
    # Update molids based on boundaries
    current_molid = 1
    atoms[0]['molid'] = current_molid
    
    for i in range(1, len(atoms)):
        # Check for molecule boundary
        new_molecule = False
        
        # Different molid indicates a boundary
        if atoms[i]['molid'] != atoms[i-1]['molid']:
            new_molecule = True
        
        # Different residue name can also indicate a boundary if use_resname is True
        if use_resname and 'resname' in atoms[i] and 'resname' in atoms[i-1]:
            if atoms[i]['resname'] != atoms[i-1]['resname']:
                new_molecule = True
        
        # Set the molid
        if new_molecule:
            current_molid += 1
        
        atoms[i]['molid'] = current_molid
    
    return atoms


def order_attributes(atoms):
    """
    Order all attributes alphabetically in each atom dictionary.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
        
    Returns
    -------
    atoms : list of dict
        The atoms list with attributes ordered.
    """
    ordered_atoms = []
    
    for atom in atoms:
        # Create a new ordered dictionary by sorting keys
        ordered_dict = {key: atom[key] for key in sorted(atom.keys())}
        ordered_atoms.append(ordered_dict)
        
    return ordered_atoms
