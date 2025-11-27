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
    
    # Normalize field consistency across all structures:
    # Use the union of all fields and fill missing ones with None to avoid dropping keys
    all_fields = set()
    for atoms in atoms_copies:
        if atoms:
            all_fields.update(atoms[0].keys())
    for i, atoms in enumerate(atoms_copies):
        if not atoms:
            continue
        for j, atom in enumerate(atoms):
            atoms_copies[i][j] = {k: atom.get(k) for k in all_fields}
    
    # If only one structure, just update indices/molids as requested
    if len(atoms_copies) == 1:
        return _update_single_structure(atoms_copies[0], molid, use_resname)

    result_atoms = []
    current_molid = 1

    for atoms in atoms_copies:
        if not atoms:
            continue

        # If all atoms already have molids, preserve them and just adjust to avoid clashes
        if all('molid' in a for a in atoms):
            # Deep copy already done; adjust offset if needed
            unique_molids = set(a['molid'] for a in atoms)
            min_molid = min(unique_molids)
            # If this is the first appended chunk, keep molids; otherwise, offset to continue sequence
            if result_atoms:
                offset = current_molid - min_molid
                for a in atoms:
                    a['molid'] = a['molid'] + offset
            result_atoms.extend(atoms)
            current_molid = max(a['molid'] for a in result_atoms) + 1
            continue

        # Otherwise, compute molids based on boundaries
        updated_atoms = _update_single_structure(atoms, None, use_resname)

        unique_molids = set(atom['molid'] for atom in updated_atoms)
        if len(unique_molids) == 1:
            for atom in updated_atoms:
                atom['molid'] = current_molid
        else:
            min_molid = min(unique_molids)
            offset = current_molid - min_molid
            for atom in updated_atoms:
                atom['molid'] += offset

        result_atoms.extend(updated_atoms)
        current_molid = max(atom['molid'] for atom in result_atoms) + 1

    # Update indices to be consecutive
    for i, atom in enumerate(result_atoms):
        atom['index'] = i + 1

    # If a specific molid was provided, set all to that value
    if molid is not None:
        for atom in result_atoms:
            atom['molid'] = molid

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

    # If all atoms already have a molid, preserve them as-is (no regrouping)
    if all('molid' in atom for atom in atoms):
        return atoms

    # Make sure all atoms have a molid
    for i, atom in enumerate(atoms):
        if 'molid' not in atom:
            atom['molid'] = i + 1
    
    # Update molids based on boundaries, starting from the first atom's molid
    current_molid = atoms[0].get('molid', 1)
    atoms[0]['molid'] = current_molid
    
    for i in range(1, len(atoms)):
        new_molecule = False
        if atoms[i]['molid'] != atoms[i-1]['molid']:
            new_molecule = True
        if use_resname and 'resname' in atoms[i] and 'resname' in atoms[i-1]:
            if atoms[i]['resname'] != atoms[i-1]['resname']:
                new_molecule = True
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
