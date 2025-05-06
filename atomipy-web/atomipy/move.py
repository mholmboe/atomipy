"""
Move module for atomipy - provides functions for translating atoms.

This module contains functions to translate atomic coordinates in atomipy data structures.
"""

import copy
import numpy as np


def translate(atoms, trans_vec, resname="all"):
    """
    Translate atom coordinates by a specified vector.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with cartesian coordinates.
    trans_vec : list or array-like
        Translation vector [x, y, z] in Angstroms.
    resname : str, optional
        Residue name to translate. Default is "all" which translates all atoms.
        
    Returns
    -------
    atoms : list of dict
        The atoms list with updated coordinates.
        
    Examples
    --------
    # Translate all atoms by [1, 2, 3]:
    new_atoms = ap.translate(atoms, [1, 2, 3])
    
    # Only translate water molecules:
    new_atoms = ap.translate(atoms, [1, 2, 3], "SOL")
    """
    # Make a deep copy to avoid modifying the original
    atoms_copy = copy.deepcopy(atoms)
    
    # Convert translation vector to numpy array for consistency
    trans_vec = np.array(trans_vec, dtype=float)
    
    # Ensure we're working with a 3D vector
    if len(trans_vec) != 3:
        raise ValueError("Translation vector must have exactly 3 components (x, y, z).")
    
    # Determine which atoms to translate based on resname
    if resname.lower() == "all":
        # Translate all atoms
        for atom in atoms_copy:
            atom['x'] += trans_vec[0]
            atom['y'] += trans_vec[1]
            atom['z'] += trans_vec[2]
    else:
        # Translate only atoms with matching resname
        for atom in atoms_copy:
            if atom['resname'] == resname:
                atom['x'] += trans_vec[0]
                atom['y'] += trans_vec[1]
                atom['z'] += trans_vec[2]
    
    return atoms_copy
