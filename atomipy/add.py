"""
Add module for atomipy - provides functions for combining and updating atom structures.

This module contains functions to update atom indices and combine multiple atom structures.
"""

import copy
import numpy as np


def update(*atoms_list, molid=None):
    """
    Update atom indices and optionally combine multiple atom structures.
    
    This function serves two main purposes:
    1. When called with a single atoms structure, it updates all indices to be consecutive
    2. When called with multiple atoms structures, it combines them into one structure 
       with consecutive indices and molecule IDs
    
    Parameters
    ----------
    *atoms_list : variable length argument list of atom structures
        One or more lists of atom dictionaries
    molid : int, optional
        If provided, sets all molecule IDs to this value
        
    Returns
    -------
    atoms : list of dict
        Combined atoms list with updated indices and molecule IDs
        
    Examples
    --------
    # Update indices of a single structure:
    new_atoms = ap.update(atoms)
    
    # Combine multiple structures:
    new_atoms = ap.update(atoms1, atoms2, atoms3)
    
    # Combine structures and set specific molecule ID:
    new_atoms = ap.update(atoms1, atoms2, molid=5)
    """
    # Make deep copies to avoid modifying originals
    atoms_copies = [copy.deepcopy(atoms) for atoms in atoms_list]
    
    # Handle case with no input
    if not atoms_copies:
        return []
    
    # Handle single atom structure case
    if len(atoms_copies) == 1:
        atoms = atoms_copies[0]
        # Update indices
        for i, atom in enumerate(atoms):
            atom['index'] = i + 1
        
        # Update molids if specified
        if molid is not None:
            for atom in atoms:
                atom['molid'] = molid
                
        return atoms
    
    # Combine multiple atom structures
    result_atoms = []
    current_index = 1
    current_molid = 1
    
    # Process each atom structure
    for atoms in atoms_copies:
        if not atoms:  # Skip empty structures
            continue
            
        # Find the maximum molid in the current result if any
        if result_atoms:
            max_molid = max(atom['molid'] for atom in result_atoms if 'molid' in atom)
            current_molid = max_molid + 1
        
        # Process each atom in the current structure
        for atom in atoms:
            atom_copy = copy.deepcopy(atom)
            atom_copy['index'] = current_index
            
            # Update molid if not specified by parameter
            if molid is None:
                # Keep consistent molid within each structure
                if 'molid' in atom:
                    # Preserve relative molid differences within each structure
                    # but offset by the current_molid to ensure uniqueness
                    if len(result_atoms) > 0 and 'molid' in result_atoms[-1]:
                        last_molid = atoms[0]['molid'] if atoms else 0
                        offset = current_molid - last_molid
                        atom_copy['molid'] = atom['molid'] + offset
                    else:
                        atom_copy['molid'] = current_molid
                else:
                    atom_copy['molid'] = current_molid
            else:
                atom_copy['molid'] = molid
                
            result_atoms.append(atom_copy)
            current_index += 1
            
        # Increment molid for the next structure if molid not specified
        if molid is None:
            current_molid += 1
    
    return result_atoms


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
