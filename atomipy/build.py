"""
Build module for atomipy - provides functions for isomorphous substitution, solvation,
and other structure building operations.
"""


import copy
import os
import random
import numpy as np
from .dist_matrix import dist_matrix
from .cell_list_dist_matrix import cell_list_dist_matrix
from .move import translate
from .transform import cartesian_to_fractional, fractional_to_cartesian
from .cell_utils import Cell2Box_dim


def is_centrosymmetric_along_z(atoms, tolerance=0.1):
    """
    Check if a structure is approximately centrosymmetric along the z-axis.
    
    
    This function analyzes the distribution of z-coordinates to determine if the
    structure is centrosymmetric, without checking atom by atom.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with coordinates.
    tolerance : float, optional
        Tolerance for considering z-coordinates as symmetric. Default is 0.1 Å.
        
    Returns
    -------
    bool
        True if the structure appears to be centrosymmetric along z, False otherwise.

    Examples
    --------
    from atomipy import is_centrosymmetric_along_z
    atoms, Box = ap.import_gro("structure.gro")
    is_centrosymmetric_along_z(atoms)
    """
    # Extract all z-coordinates
    z_coords = [atom['z'] for atom in atoms]
    
    # Find center along z
    z_center = np.mean(z_coords)
    
    # Calculate distances from center
    z_dists = [abs(z - z_center) for z in z_coords]
    
    # Sort distances for comparison
    z_dists.sort()
    
    # For a centrosymmetric structure, we expect pairs of atoms at similar distances
    # from the center (on opposite sides)
    # Group distances that are similar
    grouped_dists = {}
    for dist in z_dists:
        found_match = False
        for key in grouped_dists:
            if abs(dist - key) < tolerance:
                grouped_dists[key] += 1
                found_match = True
                break
        
        if not found_match:
            grouped_dists[dist] = 1
    
    # Check if most distances appear in pairs (even counts)
    even_count = sum(1 for count in grouped_dists.values() if count % 2 == 0)
    odd_count = len(grouped_dists) - even_count
    
    # If most groups have even counts, structure is likely centrosymmetric
    return even_count > odd_count


def substitute(atoms, Box, num_oct_subst, o1_type, o2_type, min_o2o2_dist,
               num_tet_subst=0, t1_type=None, t2_type=None, min_t2t2_dist=5.5,
               lo_limit=None, hi_limit=None, dimension=3):
    """
    Perform isomorphous substitution by replacing atom types.
    
    This function performs isomorphous substitution, replacing O1->O2 atom types
    (octahedral sites) and optionally T1->T2 atom types (tetrahedral sites).
    Substitutions are distributed equally between the top and bottom halves of
    the structure along the specified dimension.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries with coordinates and types.
    Box : list
        Box dimensions. Can be:
        - 1x3 list [lx, ly, lz] for orthogonal boxes
        - 1x6 list [a, b, c, alpha, beta, gamma] for cell parameters
        - 1x9 list for triclinic boxes (GROMACS format)
    num_oct_subst : int
        Number of octahedral substitutions to perform (O1->O2).
    o1_type : str
        Atom type to be replaced in octahedral sites.
    o2_type : str
        Replacement atom type for octahedral sites.
    min_o2o2_dist : float
        Minimum distance between O2 atoms in Angstroms.
    num_tet_subst : int, optional
        Number of tetrahedral substitutions to perform (T1->T2). Default is 0.
    t1_type : str, optional
        Atom type to be replaced in tetrahedral sites.
    t2_type : str, optional
        Replacement atom type for tetrahedral sites.
    min_t2t2_dist : float, optional
        Minimum distance between T2 atoms in Angstroms. Default is 5.5.
    lo_limit : float, optional
        Lower limit for substitution sites along the specified dimension.
        Default is -1e9 (no limit).
    hi_limit : float, optional
        Upper limit for substitution sites along the specified dimension.
        Default is 1e9 (no limit).
    dimension : int, optional
        Dimension for spatial limits: 0=x, 1=y, 2=z. Default is 2 (z).
        
    Returns
    -------
    atoms : list of dict
        Updated atoms list with substituted atom types.
        
    Examples
    --------
    # Basic octahedral substitution (5 Al->Mgo replacements):
    atoms = ap.substitute(atoms, Box, 5, 'Al', 'Mgo', 5.5)
    
    # Both octahedral and tetrahedral substitutions:
    atoms = ap.substitute(atoms, Box, 5, 'Al', 'Mgo', 5.5, 
                         num_tet_subst=2, t1_type='Si', t2_type='Alt', min_t2t2_dist=5.5)
    
    # With spatial limits (only substitute between z=-2.5 and z=12.5):
    atoms = ap.substitute(atoms, Box, 5, 'Al', 'Mgo', 5.5, 
                         num_tet_subst=2, t1_type='Si', t2_type='Alt', min_t2t2_dist=5.5,
                         lo_limit=-2.5, hi_limit=12.5, dimension=2)
    
    Notes
    -----
    The algorithm ensures that substitutions are distributed equally between the top and
    bottom halves of the structure and maintains minimum separation distances between
    substituted atoms.
    """
    # Make a deep copy to avoid modifying the original
    atoms = copy.deepcopy(atoms)
    
    # Set default limits if not provided
    if lo_limit is None:
        lo_limit = -1e9
    if hi_limit is None:
        hi_limit = 1e9
    
    # Map dimension name to index (support both MATLAB-style 1-3 and Python-style 0-2)
    if dimension == 1:
        dim_index = 0  # x
        dim_name = 'x'
    elif dimension == 2:
        dim_index = 1  # y
        dim_name = 'y'
    elif dimension == 3:
        dim_index = 2  # z
        dim_name = 'z'
    else:
        dim_index = dimension  # Assume Python-style 0-based
        dim_name = ['x', 'y', 'z'][dim_index]
    
    # Calculate average position along the specified dimension for centering
    avg_pos = np.mean([atom[dim_name] for atom in atoms])
    shift_z = 0
    
    # Check if structure is centrosymmetric along the specified dimension
    is_symmetric = False
    if dimension == 2 or dimension == 3:  # If working with z-dimension
        is_symmetric = is_centrosymmetric_along_z(atoms)
        if is_symmetric:
            print("Structure appears to be centrosymmetric along z. Centering will be applied.")
        else:
            print("Structure does not appear to be centrosymmetric along z.")
    
    # Center the structure around 0 along the specified dimension if needed or if centrosymmetric
    if avg_pos > 1 or is_symmetric:
        shift_z = avg_pos
        atoms = translate(atoms, [0, 0, -shift_z] if dim_index == 2 else 
                         ([-shift_z, 0, 0] if dim_index == 0 else [0, -shift_z, 0]))
        print(f"Translated structure by {-shift_z:.3f} Å along {dim_name} to center around 0")
    
    # Total number of substitutions
    num_total_subst = num_oct_subst + num_tet_subst
    
    # Initialize list to store O2 atoms
    o2_atoms = []
    o2_indices = []
    
    # ========== OCTAHEDRAL SUBSTITUTION ==========
    if num_oct_subst > 0:
        print(f"\n=== Performing octahedral substitution: {o1_type} -> {o2_type} ===")
        
        # Determine element for o2_type to ensure consistent updates
        from .element import element
        dummy_atom = {'type': o2_type}
        element([dummy_atom])
        o2_element = dummy_atom.get('element', '')

        # Find all O1 and O2 atoms
        ind_o1 = [i for i, atom in enumerate(atoms) if atom['type'] == o1_type or atom['type'] == o2_type]
        
        # If no atoms found by type, try to look by element
        if not ind_o1:
            print(f"No atoms with type '{o1_type}' found. Trying to match by element instead...")
            # Ensure element field is populated for all atoms
            atoms = element(atoms)
            
            # Try to match by element field
            ind_o1 = [i for i, atom in enumerate(atoms) if atom.get('element', '').upper() == o1_type.upper() 
                      or atom.get('element', '').upper() == o2_type.upper()]
            
            if ind_o1:
                print(f"Found {len(ind_o1)} atoms by element match.")
            else:
                raise ValueError(f"No atoms of type or element '{o1_type}' or '{o2_type}' found for octahedral substitution")
        
        # Extract O1 atoms
        o1_atoms = [atoms[i] for i in ind_o1]
        o1_data = np.array([[atom['x'], atom['y'], atom['z']] for atom in o1_atoms])
        ave_oct_z = np.mean(o1_data[:, dim_index])
        
        # Create random permutation for O1 atoms
        rand_o1_index = np.random.permutation(len(o1_atoms))
        
        # Select indices for substitution
        oct_subst_index = []
        
        # Calculate distance matrix for O1 atoms
        print(f"Calculating distance matrix for {len(o1_atoms)} octahedral sites...")
        o1_dist_matrix, _, _, _ = dist_matrix(o1_atoms, Box)
        
        # Perform substitutions
        i = 0
        n_oct_lo = 0
        n_oct_hi = 0
        n_oct_mid = 0
        
        while (n_oct_lo + n_oct_hi + n_oct_mid) < num_oct_subst and i < len(o1_atoms):
            # Find existing O2 atoms
            ind_o2_local = [j for j, atom in enumerate(o1_atoms) if atom['type'] == o2_type]
            
            # Check if current candidate is too close to existing O2 atoms
            o2_distances = o1_dist_matrix[rand_o1_index[i], ind_o2_local] if ind_o2_local else []
            too_close = np.any(np.array(o2_distances) < min_o2o2_dist) if len(o2_distances) > 0 else False
            
            # Get position along specified dimension
            current_pos = o1_data[rand_o1_index[i], dim_index]
            
            # Check if within limits and not too close to existing O2
            if not too_close and lo_limit < current_pos < hi_limit:
                # Distribute substitutions between low and high halves
                if n_oct_lo < num_oct_subst / 2 and current_pos < ave_oct_z:
                    oct_subst_index.append(rand_o1_index[i])
                    n_oct_lo += 1
                    o1_atoms[rand_o1_index[i]]['type'] = o2_type
                    o1_atoms[rand_o1_index[i]]['element'] = o2_element
                elif n_oct_hi < num_oct_subst / 2 and current_pos > ave_oct_z:
                    oct_subst_index.append(rand_o1_index[i])
                    n_oct_hi += 1
                    o1_atoms[rand_o1_index[i]]['type'] = o2_type
                    o1_atoms[rand_o1_index[i]]['element'] = o2_element
                elif (n_oct_lo + n_oct_hi + n_oct_mid) < num_oct_subst and current_pos == ave_oct_z:
                    oct_subst_index.append(rand_o1_index[i])
                    n_oct_mid += 1
                    o1_atoms[rand_o1_index[i]]['type'] = o2_type
                    o1_atoms[rand_o1_index[i]]['element'] = o2_element
            
            i += 1
        
        # Update main atoms list with octahedral substitutions
        for idx, global_idx in enumerate(ind_o1):
            atoms[global_idx]['type'] = o1_atoms[idx]['type']
            if 'element' in o1_atoms[idx]:
                atoms[global_idx]['element'] = o1_atoms[idx]['element']
        
        # Store O2 atoms for later use
        o2_indices = [ind_o1[idx] for idx in oct_subst_index]
        o2_atoms = [atoms[idx] for idx in o2_indices]
        
        # Report octahedral substitution results
        print(f"Octahedral substitutions: {n_oct_lo} in lower half, {n_oct_hi} in upper half, {n_oct_mid} at center")
        
        if i >= len(o1_atoms):
            print("Warning: Stopped the loop - ran out of candidate atoms")
        
        if (n_oct_lo == n_oct_hi and (n_oct_lo + n_oct_hi) == num_oct_subst) or n_oct_mid == num_oct_subst:
            print("✓ Octahedral substitution successful!")
        else:
            print("⚠ Octahedral substitution not optimal!")
    
    # ========== TETRAHEDRAL SUBSTITUTION ==========
    if num_tet_subst > 0:
        if t1_type is None or t2_type is None:
            raise ValueError("t1_type and t2_type must be specified for tetrahedral substitution")
        
        print(f"\n=== Performing tetrahedral substitution: {t1_type} -> {t2_type} ===")
        
        # Find all T1 atoms
        ind_t1 = [i for i, atom in enumerate(atoms) if atom['type'] == t1_type]
        
        # If no atoms found by type, try to look by element
        if not ind_t1:
            print(f"No atoms with type '{t1_type}' found. Trying to match by element instead...")
            # Element function should already be called from octahedral section if needed
            # Make sure element field is populated for all atoms
            if not any('element' in atom for atom in atoms[:10]):
                from .element import element
                atoms = element(atoms)
            
            # Try to match by element field
            ind_t1 = [i for i, atom in enumerate(atoms) if atom.get('element', '').upper() == t1_type.upper()]
            
            if ind_t1:
                print(f"Found {len(ind_t1)} atoms by element match.")
            else:
                raise ValueError(f"No atoms of type or element '{t1_type}' found for tetrahedral substitution")
        
        # Extract T1 atoms
        t1_atoms = [atoms[i] for i in ind_t1]
        t1_data = np.array([[atom['x'], atom['y'], atom['z']] for atom in t1_atoms])
        ave_tet_z = np.mean(t1_data[:, dim_index])
        
        # Create random permutation for T1 atoms
        rand_t1_index = np.random.permutation(len(t1_atoms))
        
        # Calculate distance matrix for T1 atoms
        print(f"Calculating distance matrix for {len(t1_atoms)} tetrahedral sites...")
        t1_dist_matrix, _, _, _ = dist_matrix(t1_atoms, Box)
        
        # Calculate distance matrix between T1 and O2 atoms if octahedral subst was done
        t1o2_dist_matrix = None
        if num_oct_subst > 0 and o2_atoms:
            print(f"Calculating T1-O2 distance matrix...")
            # Create combined list with T1 atoms first, then O2 atoms
            combined_atoms = t1_atoms + o2_atoms
            combined_dist_matrix, _, _, _ = dist_matrix(combined_atoms, Box)
            # Extract the T1-O2 block
            t1o2_dist_matrix = combined_dist_matrix[:len(t1_atoms), len(t1_atoms):]
        
        # Perform substitutions
        i = 0
        n_tet_lo = 0
        n_tet_hi = 0
        tet_subst_index = []
        
        while (n_tet_lo + n_tet_hi) < num_tet_subst and i < len(t1_atoms):
            # Find existing T2 atoms
            ind_t2_local = [j for j, atom in enumerate(t1_atoms) if atom['type'] == t2_type]
            
            # Check if current candidate is too close to existing T2 atoms
            t2_distances = t1_dist_matrix[rand_t1_index[i], ind_t2_local] if ind_t2_local else []
            too_close_t2 = np.any(np.array(t2_distances) < min_t2t2_dist) if len(t2_distances) > 0 else False
            
            # Check if too close to O2 atoms
            too_close_o2 = False
            if t1o2_dist_matrix is not None:
                to_distances = t1o2_dist_matrix[rand_t1_index[i], :]
                too_close_o2 = np.any(to_distances < min_t2t2_dist)
            
            # Get position along specified dimension
            current_pos = t1_data[rand_t1_index[i], dim_index]
            
            # Check if within limits and not too close to existing T2 or O2
            if not too_close_t2 and not too_close_o2 and lo_limit < current_pos < hi_limit:
                # Distribute substitutions between low and high halves
                if n_tet_lo < num_tet_subst / 2 and current_pos <= ave_tet_z:
                    tet_subst_index.append(rand_t1_index[i])
                    n_tet_lo += 1
                    t1_atoms[rand_t1_index[i]]['type'] = t2_type
                elif n_tet_hi < num_tet_subst / 2 and current_pos >= ave_tet_z:
                    tet_subst_index.append(rand_t1_index[i])
                    n_tet_hi += 1
                    t1_atoms[rand_t1_index[i]]['type'] = t2_type
            
            i += 1
        
        # Update main atoms list with tetrahedral substitutions
        for idx, global_idx in enumerate(ind_t1):
            atoms[global_idx]['type'] = t1_atoms[idx]['type']
        
        # Report tetrahedral substitution results
        print(f"Tetrahedral substitutions: {n_tet_lo} in lower half, {n_tet_hi} in upper half")
        
        if i >= len(t1_atoms):
            print("Warning: Stopped the loop - ran out of candidate atoms")
        
        if n_tet_lo == n_tet_hi and (n_tet_lo + n_tet_hi) == num_tet_subst:
            print("✓ Tetrahedral substitution successful!")
        else:
            print("⚠ Tetrahedral substitution not optimal!")
    
    # Shift structure back to original position if it was translated
    if abs(shift_z) > 0:
        atoms = translate(atoms, [0, 0, shift_z] if dim_index == 2 else 
                         ([shift_z, 0, 0] if dim_index == 0 else [0, shift_z, 0]))
    
    # ========== FINAL DISTANCE CHECKS ==========
    print("\n=== Distance verification ===")
    
    if num_oct_subst > 0:
        # Check minimum O2-O2 distance
        o2_atoms_final = [atom for atom in atoms if atom['type'] == o2_type]
        if len(o2_atoms_final) > 1:
            o2_dist_matrix, _, _, _ = dist_matrix(o2_atoms_final, Box)
            # Get off-diagonal minimum (exclude diagonal zeros)
            np.fill_diagonal(o2_dist_matrix, np.inf)
            min_o2_dist = np.min(o2_dist_matrix)
            print(f"Minimum {o2_type}-{o2_type} distance: {min_o2_dist:.3f} Å")
    
    if num_tet_subst > 0:
        # Check minimum T2-T2 distance
        t2_atoms_final = [atom for atom in atoms if atom['type'] == t2_type]
        if len(t2_atoms_final) > 1:
            t2_dist_matrix, _, _, _ = dist_matrix(t2_atoms_final, Box)
            np.fill_diagonal(t2_dist_matrix, np.inf)
            min_t2_dist = np.min(t2_dist_matrix)
            print(f"Minimum {t2_type}-{t2_type} distance: {min_t2_dist:.3f} Å")
        
        # Check minimum T2-O2 distance if both substitutions were done
        if num_oct_subst > 0 and o2_atoms_final and t2_atoms_final:
            combined_atoms = t2_atoms_final + o2_atoms_final
            combined_dist_matrix, _, _, _ = dist_matrix(combined_atoms, Box)
            t2o2_block = combined_dist_matrix[:len(t2_atoms_final), len(t2_atoms_final):]
            if t2o2_block.size > 0:
                min_t2o2_dist = np.min(t2o2_block)
                print(f"Minimum {t2_type}-{o2_type} distance: {min_t2o2_dist:.3f} Å")
    # Print composition
    print("\n=== Composition ===")
    _print_composition(atoms)
    
    return atoms, Box, None


def merge(atoms1, atoms2, Box, type_mode='molid', atom_label=None, min_distance=None):
    """
    Merge two atom lists by removing atoms from atoms2 that are too close to atoms1.
    
    Parameters
    ----------
    atoms1 : list of dict
        First list of atoms (usually solute)
    atoms2 : list of dict
        Second list of atoms (usually solvent)
    Box : list of float
        Box dimensions for the system
    type_mode : str, optional
        Mode for identifying atoms to check distances: 'molid' or 'index'
        If 'molid', removes entire molecules that are too close. Default is 'molid'.
        If 'index', removes individual atoms that are too close.
    atom_label : str or list of str, optional
        Atom type(s) in atoms2 for applying the smaller min_distance cutoff.
        Default is None (same cutoff for all atoms).
    min_distance : float or list of float, optional
        Minimum distance threshold. Can be a single value or [standard, small]
        where 'small' applies to atom_label types and 'standard' to other types.
        Default is None (error if not provided).
        
    Returns
    -------
    list of dict
        Copy of atoms2 with atoms/molecules too close to atoms1 removed
    """
    # Validate inputs
    if min_distance is None:
        raise ValueError("min_distance must be provided")
    
    if isinstance(min_distance, (int, float)):
        standard_dist = min_distance
        small_dist = min_distance
    elif len(min_distance) == 2:
        standard_dist, small_dist = min_distance
    else:
        raise ValueError("min_distance must be a single value or a list of two values")
    
    # Make copies to avoid modifying originals
    atoms1 = copy.deepcopy(atoms1)
    atoms2 = copy.deepcopy(atoms2)
    
    # Calculate distances between atoms1 and atoms2
    if len(atoms1) + len(atoms2) > 20000:
        # For large systems, use cell list method
        combined = atoms1 + atoms2
        dist_matrix_result = cell_list_dist_matrix(combined, Box)
        # Extract only the atoms1-atoms2 part
        dist_matrix_result = dist_matrix_result[:len(atoms1), len(atoms1):]
        if isinstance(dist_matrix_result, tuple):
            # If dist_matrix returns a tuple, use the first element (the actual distance matrix)
            dist_matrix_result = dist_matrix_result[0][:len(atoms1), len(atoms1):]
    else:
        # For smaller systems, use standard distance matrix
        xs1 = np.array([atom['x'] for atom in atoms1])
        ys1 = np.array([atom['y'] for atom in atoms1])
        zs1 = np.array([atom['z'] for atom in atoms1])
        coords1 = np.column_stack((xs1, ys1, zs1))
        
        xs2 = np.array([atom['x'] for atom in atoms2])
        ys2 = np.array([atom['y'] for atom in atoms2])
        zs2 = np.array([atom['z'] for atom in atoms2])
        coords2 = np.column_stack((xs2, ys2, zs2))
        
        # Calculate distances
        dist_matrix_result = dist_matrix(atoms1 + atoms2, Box)
        if isinstance(dist_matrix_result, tuple):
            # If dist_matrix returns a tuple, use the first element (the actual distance matrix)
            dist_matrix_result = dist_matrix_result[0][:len(atoms1), len(atoms1):]
        else:
            dist_matrix_result = dist_matrix_result[:len(atoms1), len(atoms1):]
    
    # Create a mask for atoms to remove
    to_remove = []
    for i in range(len(atoms2)):
        atom = atoms2[i]
        # Check if this atom is of the type that should use small_dist
        use_small_dist = False
        if atom_label is not None:
            if isinstance(atom_label, str) and atom.get('type', '') == atom_label:
                use_small_dist = True
            elif isinstance(atom_label, (list, tuple)) and atom.get('type', '') in atom_label:
                use_small_dist = True
        
        dist_threshold = small_dist if use_small_dist else standard_dist
        
        # Check if any atom in atoms1 is too close
        for j in range(len(atoms1)):
            if dist_matrix_result[j, i] < dist_threshold:
                to_remove.append(i)
                break
    
    # Handle removal based on type_mode
    if type_mode.lower() == 'molid':
        # Get unique molids to remove
        molids_to_remove = set(atoms2[i]['molid'] for i in to_remove)
        result = [atom for atom in atoms2 if atom['molid'] not in molids_to_remove]
    else:  # index mode
        # Remove individual atoms
        indices_to_keep = [i for i in range(len(atoms2)) if i not in to_remove]
        result = [atoms2[i] for i in indices_to_keep]
    
    return result


def _print_composition(atoms):
    """Print the composition of the structure."""
    # Count atom types
    type_counts = {}
    for atom in atoms:
        atom_type = atom.get('type', 'Unknown')
        type_counts[atom_type] = type_counts.get(atom_type, 0) + 1
    
    # Sort by type name
    sorted_types = sorted(type_counts.items())
    
    # Print composition
    total_atoms = len(atoms)
    print(f"Total atoms: {total_atoms}")
    for atom_type, count in sorted_types:
        percentage = (count / total_atoms) * 100
        print(f"  {atom_type}: {count} ({percentage:.1f}%)")


def molecule(atoms, molid=1, resname=None):
    """
    Assign molecule ID (and optionally residue name) to atoms.

    This function assigns the same molid to all atoms in the list, treating
    them as a single molecule. Similar to MATLAB's molecule_atom function.

    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
    molid : int, optional
        Molecule ID to assign to all atoms (default: 1).
    resname : str, optional
        Residue name to assign to all atoms. If None, existing resnames 
        are preserved (default: None).

    Returns
    -------
    list of dict
        Updated atoms with molid (and optionally resname) assigned.

    Notes
    -----
    - This function modifies atoms in-place and also returns the modified list.
    - Use this to group atoms as a single molecular unit before topology generation.

    Examples
    --------
    # Assign all mineral atoms to molecule 1 with resname 'MIN':
    MIN = ap.molecule(MIN, molid=1, resname='MIN')
    
    # Assign ions to molecule 2:
    IONS = ap.molecule(IONS, molid=2)
    """
    atoms = copy.deepcopy(atoms)
    
    for atom in atoms:
        atom['molid'] = molid
        if resname is not None:
            atom['resname'] = resname
    
    n_atoms = len(atoms)
    if resname:
        print(f"Assigned molid={molid} and resname='{resname}' to {n_atoms} atoms")
    else:
        print(f"Assigned molid={molid} to {n_atoms} atoms")
    
    return atoms


def slice(atoms, limits, remove_partial_molecules=True):
    """
    Extract atoms within a region defined by limits.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries
    limits : list of float
        Region limits [xlo, ylo, zlo, xhi, yhi, zhi] or [xhi, yhi, zhi] 
        (the latter assumes xlo=ylo=zlo=0)
    remove_partial_molecules : bool, optional
        If True, remove molecules that are partially outside the region.
        Default is True.
        
    Returns
    -------
    list of dict
        Atoms within the specified region
    """
    # Standardize limits to [xlo, ylo, zlo, xhi, yhi, zhi] format
    if len(limits) == 3:
        xlo, ylo, zlo = 0, 0, 0
        xhi, yhi, zhi = limits
    elif len(limits) == 6:
        xlo, ylo, zlo, xhi, yhi, zhi = limits
    else:
        raise ValueError("Limits must be a list of length 3 [xhi, yhi, zhi] "
                         "or 6 [xlo, ylo, zlo, xhi, yhi, zhi]")
    
    # Check which atoms are within the region
    in_region = []
    for i, atom in enumerate(atoms):
        if (xlo <= atom['x'] <= xhi and 
            ylo <= atom['y'] <= yhi and 
            zlo <= atom['z'] <= zhi):
            in_region.append(i)
    
    if remove_partial_molecules:
        # Get molids of all atoms in region
        molids_in_region = set(atoms[i]['molid'] for i in in_region)
        
        # Check for each molid if all its atoms are in the region
        complete_molids = set()
        for molid in molids_in_region:
            # Find all atoms with this molid
            molid_atoms = [i for i, atom in enumerate(atoms) if atom['molid'] == molid]
            # Check if all atoms of this molecule are in the region
            if all(i in in_region for i in molid_atoms):
                complete_molids.add(molid)
        
        # Select only atoms from complete molecules
        selected_atoms = [copy.deepcopy(atom) for atom in atoms 
                         if atom['molid'] in complete_molids]
    else:
        # Select all atoms in region
        selected_atoms = [copy.deepcopy(atoms[i]) for i in in_region]
    
    return selected_atoms


def _get_surface_atoms(atoms, distance_threshold=2.5):
    """
    Identify atoms at the surface of a structure.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries
    distance_threshold : float, optional
        Distance threshold for considering an atom to be at the surface.
        Default is 2.5 Å.
        
    Returns
    -------
    list of dict
        List of surface atoms
    """
    # First compute the convex hull using atom coordinates
    coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
    
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(coords)
        
        # Get indices of atoms on the hull
        surface_indices = np.unique(hull.simplices.flatten())
        surface_atoms = [atoms[i] for i in surface_indices]
        
        return surface_atoms
    except ImportError:
        # If scipy is not available, use a distance-based approach
        surface_atoms = []
        
        # Use cell list distance matrix for efficiency
        dist_matrix = cell_list_dist_matrix(atoms, [1000, 1000, 1000])  # Large Box to ignore PBC
        
        for i, atom in enumerate(atoms):
            # An atom is on the surface if it has empty space around it
            # We check if any direction has no atoms within distance_threshold
            neighbors = 0
            for j in range(len(atoms)):
                if i != j and dist_matrix[i, j] < distance_threshold:
                    neighbors += 1
            
            # If atom has few neighbors, it's likely on the surface
            if neighbors < 12:  # Typical number for dense packing is 12-14
                surface_atoms.append(atom)
        
        return surface_atoms


def ionize(ion_type, resname, limits, num_ions, min_distance=None, solute_atoms=None, 
           placement='random', direction=None, direction_value=None):
    """
    Add ions to a system within specified region limits.
    
    Parameters
    ----------
    ion_type : str
        Type of ion to add (e.g., 'Na', 'Cl')
    resname : str
        Residue name for the ions
    limits : list of float
        Region limits [xlo, ylo, zlo, xhi, yhi, zhi] or [xhi, yhi, zhi]
    num_ions : int
        Number of ions to add
    min_distance : float, optional
        Minimum distance between ions and existing atoms.
        If None, uses twice the sum of ionic radii
    solute_atoms : list of dict, optional
        Existing solute atoms to avoid overlaps with
    placement : str, optional
        Where to place ions: 'random', 'surface', or 'bulk'
    direction : str or float, optional
        Direction constraint for surface/bulk placement ('x', 'y', 'z' or value)
    direction_value : float, optional
        Value constraint for the specified direction
        
    Returns
    -------
    list of dict
        List of added ion atoms
        
    Examples
    --------
    # Add 10 Na+ ions randomly in a Box:
    ions = ap.ionize('Na', 'NA', [0, 0, 0, 50, 50, 50], 10)
    
    # Add 8 Cl- ions near a protein surface with minimum distance 3 Å:
    ions = ap.ionize('Cl', 'CL', [0, 0, 0, 60, 60, 60], 8, 3.0, protein_atoms, 'surface')
    """
    from .radius import ionic_radius
    
    # Standardize limits to [xlo, ylo, zlo, xhi, yhi, zhi] format
    if len(limits) == 3:
        xlo, ylo, zlo = 0, 0, 0
        xhi, yhi, zhi = limits
    elif len(limits) == 6:
        xlo, ylo, zlo, xhi, yhi, zhi = limits
    else:
        raise ValueError("Limits must be a list of length 3 [xhi, yhi, zhi] "
                         "or 6 [xlo, ylo, zlo, xhi, yhi, zhi]")
    
    # Calculate Box dimensions for the region
    box_dim = [xhi - xlo, yhi - ylo, zhi - zlo]
    
    # If minimum distance not specified, use ionic radii
    if min_distance is None:
        ion_radius = ionic_radius().get(ion_type, 1.0)  # Default to 1.0 if not found
        if solute_atoms:
            # Use twice the sum of radii as a safe default
            min_distance = 2 * ion_radius
        else:
            # If no solute, just use twice the ionic radius
            min_distance = 2 * ion_radius
    
    # Initialize result
    ions = []
    
    # If using surface or bulk placement, need to process solute
    surface_atoms = None
    if solute_atoms and placement.lower() in ['surface', 'bulk']:
        surface_atoms = _get_surface_atoms(solute_atoms)
    
    # Generate candidate positions
    max_attempts = num_ions * 100  # Limit total attempts
    attempts = 0
    
    while len(ions) < num_ions and attempts < max_attempts:
        # Generate random position within limits
        x = xlo + np.random.random() * (xhi - xlo)
        y = ylo + np.random.random() * (yhi - ylo)
        z = zlo + np.random.random() * (zhi - zlo)
        
        # Apply direction constraint if specified
        if direction is not None and direction_value is not None:
            if direction.lower() == 'x':
                x = direction_value
            elif direction.lower() == 'y':
                y = direction_value
            elif direction.lower() == 'z':
                z = direction_value
        
        # Create ion atom
        ion = {
            'index': len(ions) + 1,
            'molid': len(ions) + 1,  # Each ion is its own molecule
            'resname': resname,
            'type': ion_type,
            'element': ion_type,  # Assume type is element
            'x': x,
            'y': y,
            'z': z
        }
        
        # Check for overlaps with solute atoms
        overlap = False
        dist_to_surface = float('inf')
        
        if solute_atoms:
            for solute_atom in solute_atoms:
                dx = ion['x'] - solute_atom['x']
                dy = ion['y'] - solute_atom['y']
                dz = ion['z'] - solute_atom['z']
                dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                
                if dist < min_distance:
                    overlap = True
                    break
                
                # For surface placement, track distance to surface atoms
                if placement.lower() == 'surface' and surface_atoms:
                    if solute_atom in surface_atoms:
                        dist_to_surface = min(dist_to_surface, dist)
        
        # Check for overlaps with already placed ions
        if not overlap and ions:
            for placed_ion in ions:
                dx = ion['x'] - placed_ion['x']
                dy = ion['y'] - placed_ion['y']
                dz = ion['z'] - placed_ion['z']
                dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                
                if dist < min_distance:
                    overlap = True
                    break
        
        # Add ion if no overlaps and it meets placement criteria
        if not overlap:
            # For surface placement, only accept ions near the surface
            if placement.lower() == 'surface' and surface_atoms:
                # Only accept if close to the surface but not overlapping
                if dist_to_surface < 5.0:
                    ions.append(ion)
            # For bulk placement, only accept ions away from the surface
            elif placement.lower() == 'bulk' and surface_atoms:
                # Only accept if far from the surface
                if dist_to_surface > 8.0:
                    ions.append(ion)
            # For random placement, accept any valid position
            else:
                ions.append(ion)
        
        attempts += 1
    
    if len(ions) < num_ions:
        print(f"Warning: Could only place {len(ions)} out of {num_ions} requested ions")
    
    return ions


def insert(molecule_atoms, limits, rotate='random', min_distance=2.0, 
           num_molecules=1, solute_atoms=None, type_constraints=None, z_diff=None):
    """
    Insert molecules into a system within specified region limits.
    
    Parameters
    ----------
    molecule_atoms : list of dict
        Template molecule atoms to insert
    limits : list of float
        Region limits [xlo, ylo, zlo, xhi, yhi, zhi] or [xhi, yhi, zhi]
    rotate : str or list of float, optional
        Rotation to apply: 'random' or specific [alpha, beta, gamma] angles
    min_distance : float, optional
        Minimum distance between inserted molecules and existing atoms
    num_molecules : int, optional
        Number of molecules to insert
    solute_atoms : list of dict, optional
        Existing solute atoms to avoid overlaps with
    type_constraints : list or tuple, optional
        Atom types (type1, type2) with z-positioning constraints
    z_diff : float, optional
        Minimum z-difference between the constrained atom types
        
    Returns
    -------
    list of dict
        List of inserted molecule atoms
        
    Examples
    --------
    # Insert 5 copies of a molecule randomly in a Box:
    new_atoms = ap.insert(molecule_atoms, [0, 0, 0, 50, 50, 50], num_molecules=5)
    
    # Insert a molecule with specific orientation and position constraints:
    new_atoms = ap.insert(molecule_atoms, [0, 0, 0, 50, 50, 50], 
                         rotate=[0, 0, 90], type_constraints=['C', 'N'], z_diff=3.0)
    """
    from .move import rotate as rotate_atoms
    from .move import place
    
    # Standardize limits to [xlo, ylo, zlo, xhi, yhi, zhi] format
    if len(limits) == 3:
        xlo, ylo, zlo = 0, 0, 0
        xhi, yhi, zhi = limits
    elif len(limits) == 6:
        xlo, ylo, zlo, xhi, yhi, zhi = limits
    else:
        raise ValueError("Limits must be a list of length 3 [xhi, yhi, zhi] "
                         "or 6 [xlo, ylo, zlo, xhi, yhi, zhi]")
    
    # Calculate Box dimensions for the region
    box_dim = [xhi - xlo, yhi - ylo, zhi - zlo]
    
    # Initialize result
    inserted_atoms = []
    
    # Try to insert each molecule
    molecule_counter = 0
    max_attempts_per_molecule = 100
    max_total_attempts = num_molecules * max_attempts_per_molecule
    total_attempts = 0
    
    while molecule_counter < num_molecules and total_attempts < max_total_attempts:
        # Make a copy of the template molecule
        temp_molecule = copy.deepcopy(molecule_atoms)
        
        # Rotate the molecule
        temp_molecule = rotate_atoms(temp_molecule, box_dim, rotate)
        
        # Generate random position within limits
        pos_x = xlo + np.random.random() * (xhi - xlo)
        pos_y = ylo + np.random.random() * (yhi - ylo)
        pos_z = zlo + np.random.random() * (zhi - zlo)
        
        # Place the molecule at the random position
        temp_molecule = place(temp_molecule, [pos_x, pos_y, pos_z])
        
        # Make sure the whole molecule is within the limits
        in_bounds = True
        for atom in temp_molecule:
            if (atom['x'] < xlo or atom['x'] > xhi or
                atom['y'] < ylo or atom['y'] > yhi or
                atom['z'] < zlo or atom['z'] > zhi):
                in_bounds = False
                break
        
        if not in_bounds:
            total_attempts += 1
            continue
        
        # Check for overlaps with solute atoms
        overlap = False
        if solute_atoms:
            for solute_atom in solute_atoms:
                for mol_atom in temp_molecule:
                    dx = mol_atom['x'] - solute_atom['x']
                    dy = mol_atom['y'] - solute_atom['y']
                    dz = mol_atom['z'] - solute_atom['z']
                    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    if dist < min_distance:
                        overlap = True
                        break
                if overlap:
                    break
        
        # Check for overlaps with already inserted molecules
        if not overlap and inserted_atoms:
            for existing_atom in inserted_atoms:
                for mol_atom in temp_molecule:
                    dx = mol_atom['x'] - existing_atom['x']
                    dy = mol_atom['y'] - existing_atom['y']
                    dz = mol_atom['z'] - existing_atom['z']
                    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    if dist < min_distance:
                        overlap = True
                        break
                if overlap:
                    break
        
        # Check type constraints if specified
        if not overlap and type_constraints and z_diff is not None:
            type1, type2 = type_constraints
            
            # Find atoms of the constrained types
            type1_atoms = [atom for atom in temp_molecule if atom['type'] == type1]
            type2_atoms = [atom for atom in temp_molecule if atom['type'] == type2]
            
            if type1_atoms and type2_atoms:
                # Check z-positioning constraint
                z1 = np.median([atom['z'] for atom in type1_atoms])
                z2 = np.median([atom['z'] for atom in type2_atoms])
                
                if z1 < z2 + z_diff:  # type1 should be above type2 by at least z_diff
                    overlap = True  # Not actually overlap, but constraint not met
        
        # Add molecule if no overlaps and constraints satisfied
        if not overlap:
            # Update molid for the molecule
            current_molid = max([atom['molid'] for atom in inserted_atoms]) + 1 if inserted_atoms else 1
            for atom in temp_molecule:
                atom['molid'] = current_molid
                atom['index'] = len(inserted_atoms) + 1
                inserted_atoms.append(atom)
            
            molecule_counter += 1
        
        total_attempts += 1
    
    if molecule_counter < num_molecules:
        print(f"Warning: Could only insert {molecule_counter} out of {num_molecules} requested molecules")
    
    return inserted_atoms


# =====================================================
# Hydrogen Manipulation Functions
# =====================================================

def add_H_atom(atoms, Box, target_type, h_type='H', bond_length=0.96, coordination=1, max_h_per_atom=1):
    """
    Add hydrogen atoms to under-coordinated atoms of a specific type.
    
    This is commonly used to protonate edge oxygen atoms or other sites.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
    Box : list
        Box dimensions [lx, ly, lz, ...].
    target_type : str
        Atom type to check for protonation (e.g., 'Oh', 'Ob').
    h_type : str, optional
        Atom type for the new hydrogen atoms. Default is 'H'.
    bond_length : float, optional
        Distance for the new O-H bond. Default is 0.96 A.
    coordination : int, optional
        Target coordination number. If current CN < coordination, H is added. Default is 1.
        NOTE: This checks ALL neighbors. For specificity (e.g. only metal neighbors), simpler logic 
        is usually sufficient: if an atom has fewer than X total neighbors, add H.
    max_h_per_atom : int, optional
        Maximum number of H atoms to add per target atom. Default is 1.
        
    Returns
    -------
    list of dict
        Updated atoms list with new hydrogens.
    """
    print(f"Adding H atoms to '{target_type}' (target CN={coordination})...")
    
    # Needs neighbor list
    if not atoms or 'neigh' not in atoms[0]:
        print("  Calculating neighbor list...")
        from .bond_angle import bond_angle
        atoms, _, _ = bond_angle(atoms, Box, rmaxM=2.45, rmaxH=1.2, same_molecule_only=False) # Use typical metal-oxygen cutoff
        
    new_atoms = []
    count_added = 0
    
    # Iterate through existing atoms containing neighbors
    n_original = len(atoms)
    
    for i in range(n_original):
        atom = atoms[i]
        
        # Check if this is a target atom
        if atom.get('type') != target_type:
            continue
            
        neighbors = atom.get('neigh', [])
        current_cn = len(neighbors)
        
        # Check if under-coordinated
        if current_cn < coordination:
            n_needed = min(coordination - current_cn, max_h_per_atom)
            
            if n_needed <= 0:
                continue
                
            origin = np.array([atom['x'], atom['y'], atom['z']])
            
            # Get vectors to existing neighbors
            neighbor_vectors = []
            for n_idx in neighbors:
                n_atom = atoms[n_idx]
                n_pos = np.array([n_atom['x'], n_atom['y'], n_atom['z']])
                vec = n_pos - origin
                # Minimal simple PBC check assuming orthogonal Box for vector calculation
                if Box is not None and len(Box) >= 3:
                     for d in range(3):
                         L = Box[d]
                         if vec[d] > L/2: vec[d] -= L
                         elif vec[d] < -L/2: vec[d] += L
                neighbor_vectors.append(vec / np.linalg.norm(vec))
            
            # Add H atoms
            for _ in range(n_needed):
                # Generate random unit vector
                for attempt in range(20):
                    # Random point on sphere
                    costheta = random.uniform(-1, 1)
                    phi = random.uniform(0, 2*np.pi)
                    theta = np.arccos(costheta)
                    
                    hx = np.sin(theta) * np.cos(phi)
                    hy = np.sin(theta) * np.sin(phi)
                    hz = np.cos(theta)
                    h_vec = np.array([hx, hy, hz])
                    
                    # Check angle with existing neighbors (avoid overlap)
                    too_close = False
                    for n_vec in neighbor_vectors:
                        dot = np.dot(h_vec, n_vec)
                        if dot > 0.9: # Very close overlap (~25 deg)
                            too_close = True
                            break
                    
                    if not too_close or attempt == 19:
                        # Accept position
                        h_pos = origin + h_vec * bond_length
                        
                        # Create new atom
                        new_h = {
                            'type': h_type,
                            'element': 'H',
                            'resname': atom.get('resname', 'MIN'),
                            'molid': atom.get('molid', 1),
                            'x': h_pos[0],
                            'y': h_pos[1],
                            'z': h_pos[2],
                            'charge': 0.4, # Default placeholder
                            'mass': 1.008
                        }
                        new_atoms.append(new_h)
                        
                        # Treat this new H as a neighbor for subsequent H additions to same atom
                        neighbor_vectors.append(h_vec)
                        break
            
            count_added += len(new_atoms) - count_added # Track added count

    print(f"  Added {len(new_atoms)} new '{h_type}' atoms.")
    
    # Combine lists
    full_list = atoms + new_atoms
    
    # Re-index
    for i, atom in enumerate(full_list):
        atom['index'] = i + 1
        # Clear neighbor lists as they are now stale/incomplete
        if 'neigh' in atom:
            del atom['neigh']
        
    return full_list


def adjust_H_atom(atoms, Box, h_type='H', neighbor_type='O', distance=0.96):
    """
    Adjust bond lengths of hydrogen atoms involved in specified bonds.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
    Box : list
        Box dimensions.
    h_type : str
        Type of hydrogen atom to adjust.
    neighbor_type : str
        Type of atom the hydrogen is bonded to.
    distance : float
        Target bond length.
        
    Returns
    -------
    list of dict
        Updated atoms list.
    """
    print(f"Adjusting '{h_type}'-'{neighbor_type}' bonds to {distance} A...")
    
    # Needs neighbor list
    if not atoms or 'neigh' not in atoms[0]:
        from .bond_angle import bond_angle
        atoms, _, _ = bond_angle(atoms, Box, rmaxH=2.0, rmaxM=2.5, same_molecule_only=False) # Larger cutoff for potentially distorted bonds
    
    adj_count = 0
    xyz_box = Box[:3] if Box else [1000, 1000, 1000]
    
    for i, atom in enumerate(atoms):
        if atom.get('type') != h_type:
            continue
            
        neighbors = atom.get('neigh', [])
        
        # Find the bonded neighbor of correct type
        bonded_idx = -1
        for n_idx in neighbors:
            if atoms[n_idx].get('type') == neighbor_type:
                bonded_idx = n_idx
                break
        
        if bonded_idx != -1:
            neighbor = atoms[bonded_idx]
            
            # Vector neighbor -> H
            dx = atom['x'] - neighbor['x']
            dy = atom['y'] - neighbor['y']
            dz = atom['z'] - neighbor['z']
            
            # MIC
            if dx > xyz_box[0]/2: dx -= xyz_box[0]
            elif dx < -xyz_box[0]/2: dx += xyz_box[0]
            if dy > xyz_box[1]/2: dy -= xyz_box[1]
            elif dy < -xyz_box[1]/2: dy += xyz_box[1]
            if dz > xyz_box[2]/2: dz -= xyz_box[2]
            elif dz < -xyz_box[2]/2: dz += xyz_box[2]
            
            current_dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if current_dist > 0.01: # Avoid division by zero
                scale = distance / current_dist
                
                # Update H position relative to neighbor
                atom['x'] = neighbor['x'] + dx * scale
                atom['y'] = neighbor['y'] + dy * scale
                atom['z'] = neighbor['z'] + dz * scale
                
                adj_count += 1
                
    print(f"  Adjusted {adj_count} bonds.")
    return atoms


def adjust_Hw_atom(atoms, Box, water_resname='SOL', water_model='OPC3'):
    """
    Repair water molecules: add missing hydrogens and fix geometry.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries.
    Box : list
        Box dimensions.
    water_resname : str
        Residue name for water (e.g. 'SOL', 'WAT').
    water_model : str
        Water model to determine ideal geometry ('SPC', 'TIP3P', 'OPC3', etc.).
        Default 'OPC3': r(OH)=0.9572, angle=104.52
        
    Returns
    -------
    list of dict
        Updated atoms list.
    """
    print(f"Adjusting water molecules ('{water_resname}', model={water_model})...")
    
    # Model parameters
    if 'SPC' in water_model.upper():
        r_oh = 1.0
        angle_hoh_deg = 109.47
        qh = 0.41 if 'E' not in water_model else 0.4238
    else: # Default/OPC3/TIP3P
        r_oh = 0.9572
        angle_hoh_deg = 104.52
        qh = 0.447585 if 'OPC' in water_model else 0.417
        
    angle_rad = np.deg2rad(angle_hoh_deg)
    
    # Find water atoms
    water_indices = [i for i, a in enumerate(atoms) if a.get('resname', '') == water_resname]
    if not water_indices:
        print("  No water atoms found.")
        return atoms
        
    # Organize by molid (or residue number if molid missing)
    molecules = {}
    for idx in water_indices:
        atom = atoms[idx]
        mol_id = atom.get('molid', atom.get('resnr', -1))
        if mol_id not in molecules:
            molecules[mol_id] = []
        molecules[mol_id].append(idx)
        
    new_atoms = []
    
    for mol_id, indices in molecules.items():
        oxy_indices = [i for i in indices if atoms[i]['type'].upper().startswith('O')]
        hyd_indices = [i for i in indices if atoms[i]['type'].upper().startswith('H')]
        
        if not oxy_indices:
            continue
            
        oxy_idx = oxy_indices[0] # Assume 1 oxygen per water
        oxy = atoms[oxy_idx]
        ox_pos = np.array([oxy['x'], oxy['y'], oxy['z']])
        
        current_h_count = len(hyd_indices)
        
        if current_h_count < 2:
            # Generate vectors
            costheta = random.uniform(-1, 1)
            phi = random.uniform(0, 2*np.pi)
            theta = np.arccos(costheta)
            h1_vec = np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta)
            ])
            
            # Orthogonal vector for rotation axis
            tmp_vec = np.array([1, 0, 0])
            if abs(np.dot(h1_vec, tmp_vec)) > 0.9: tmp_vec = np.array([0, 1, 0])
            axis = np.cross(h1_vec, tmp_vec)
            axis /= np.linalg.norm(axis)
            
            # H2 vector rotated by angle
            h2_vec = h1_vec * np.cos(angle_rad) + np.cross(axis, h1_vec) * np.sin(angle_rad)
            
            # ... Logic for partial existence (1 H) skipped for brevity of initial impl, 
            # assume simplest case of adding missing ones from scratch relative to O or random
            # Just create new Hs
            
            h_positions_to_add = []
            if current_h_count == 0:
                h_positions_to_add = [ox_pos + h1_vec * r_oh, ox_pos + h2_vec * r_oh]
            elif current_h_count == 1:
                # Keep H1, add H2 relative to O-H1
                h_idx = hyd_indices[0]
                h1_vec_exist = np.array([atoms[h_idx]['x'], atoms[h_idx]['y'], atoms[h_idx]['z']]) - ox_pos
                h1_vec_exist /= np.linalg.norm(h1_vec_exist)
                
                # New axis perpendicular to existing bond
                tmp_vec = np.array([1, 0, 0])
                if abs(np.dot(h1_vec_exist, tmp_vec)) > 0.9: tmp_vec = np.array([0, 1, 0])
                axis = np.cross(h1_vec_exist, tmp_vec)
                axis /= np.linalg.norm(axis)
                
                h2_vec_new = h1_vec_exist * np.cos(angle_rad) + np.cross(axis, h1_vec_exist) * np.sin(angle_rad)
                h_positions_to_add = [ox_pos + h2_vec_new * r_oh]
                
            for h_pos in h_positions_to_add:
                new_h = {
                    'type': 'HW',
                    'element': 'H',
                    'resname': oxy['resname'],
                    'molid': oxy['molid'],
                    'x': h_pos[0],
                    'y': h_pos[1],
                    'z': h_pos[2],
                    'charge': qh,
                    'mass': 1.008
                }
                new_atoms.append(new_h)
                    
        elif current_h_count == 2:
            # Fix geometry of existing 2 H
            h1_idx, h2_idx = hyd_indices[0], hyd_indices[1]
            h1, h2 = atoms[h1_idx], atoms[h2_idx]
            
            v1 = np.array([h1['x'], h1['y'], h1['z']]) - ox_pos
            v2 = np.array([h2['x'], h2['y'], h2['z']]) - ox_pos
            
            # Bisector
            bisector = v1 + v2
            if np.linalg.norm(bisector) < 0.1: bisector = np.array([1, 0, 0])
            bisector /= np.linalg.norm(bisector)
            
            # Plane normal
            normal = np.cross(v1, v2)
            if np.linalg.norm(normal) < 0.1: normal = np.array([0, 0, 1])
            normal /= np.linalg.norm(normal)
            
            # New vectors
            half_angle = angle_rad / 2
            v1_new = bisector * np.cos(half_angle) + np.cross(normal, bisector) * np.sin(half_angle)
            v2_new = bisector * np.cos(-half_angle) + np.cross(normal, bisector) * np.sin(-half_angle)
            
            p1 = ox_pos + v1_new * r_oh
            p2 = ox_pos + v2_new * r_oh
            
            atoms[h1_idx].update({'x': p1[0], 'y': p1[1], 'z': p1[2]})
            atoms[h2_idx].update({'x': p2[0], 'y': p2[1], 'z': p2[2]})

    print(f"  Added {len(new_atoms)} missing water hydrogens.")
    
    full_list = atoms + new_atoms
    for i, atom in enumerate(full_list):
        atom['index'] = i + 1
        # Clear neighbor lists as they are now stale/incomplete
        if 'neigh' in atom:
            del atom['neigh']
            
    return full_list
