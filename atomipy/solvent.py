"""
Solvent module for atomipy - provides functions for handling solvent molecules.

This module provides functions for:
- Identifying solvent molecules (e.g., water) in molecular systems
- Solvating structures with water or other solvents
"""

import copy
import os
import random
import numpy as np


def find_H2O(atoms, Box_dim=None, rmin=1.25):
    """
    Find water molecules in a molecular system using distance-based detection.

    Identifies water molecules by computing distances and finding oxygen atoms
    that have exactly 2 hydrogen atoms within the O-H bond distance cutoff.
    This works regardless of atom order or existing molid assignments.

    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries. Should contain 'element' or 'type' fields
        for atom identification, and 'x', 'y', 'z' coordinates.
    Box_dim : list of float, optional
        Box dimensions for periodic boundary conditions. If None, non-periodic
        distances are used (default: None).
    rmin : float, optional
        Maximum O-H bond distance in Angstroms for water detection
        (default: 1.25, matches typical water O-H bond length).

    Returns
    -------
    water_atoms : list of dict
        List of atom dictionaries belonging to water molecules. Atoms are
        updated with new molids and types ('Ow', 'Hw').
    non_water_atoms : list of dict
        Remaining atoms that are not part of water molecules.

    Notes
    -----
    - Uses distance matrix to find O atoms with exactly 2 H neighbors
    - Assigns new sequential molids to each water molecule
    - Changes atom types to 'Ow' and 'Hw' for water oxygens and hydrogens
    - Sets resname to 'SOL' for all water atoms
    - Works regardless of input atom order or existing molid values

    Examples
    --------
    atoms, Box_dim = ap.import_gro("system.gro")
    atoms = ap.element(atoms)
    SOL, noSOL = ap.find_H2O(atoms, Box_dim)
    print(f"Found {len(SOL)//3} water molecules")

    See Also
    --------
    assign_resname : Assign residue names including 'SOL' for water
    """
    from .dist_matrix import dist_matrix
    
    if not atoms:
        return [], []
    
    atoms = copy.deepcopy(atoms)
    
    # Find all O and H atom indices based on element/type
    O_indices = []
    H_indices = []
    
    for i, atom in enumerate(atoms):
        element = atom.get('element', '')
        atom_type = atom.get('type', atom.get('atom_name', ''))
        
        # Check element first (can be 'O', 'Ow', 'Oh', etc), then type
        is_oxygen = (element and element[0].upper() == 'O') or \
                    (not element and atom_type and atom_type[0].upper() == 'O')
        is_hydrogen = (element and element[0].upper() == 'H') or \
                      (not element and atom_type and atom_type[0].upper() == 'H')
        
        if is_oxygen:
            O_indices.append(i)
        elif is_hydrogen:
            H_indices.append(i)
    
    if not O_indices or len(H_indices) < 2:
        print("Found 0 water molecules (0 atoms)")
        return [], atoms
    
    # Compute distance matrix
    if Box_dim is not None:
        dist_mat, _, _, _ = dist_matrix(atoms, Box_dim)
    else:
        # Non-periodic distance calculation
        coords = np.array([[a['x'], a['y'], a['z']] for a in atoms])
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dist_mat = np.sqrt(np.sum(diff**2, axis=2))
    
    # Find water molecules: O atoms with exactly 2 H neighbors within rmin
    water_atoms = []
    water_indices_set = set()
    molid = 1
    
    for o_idx in O_indices:
        # Find H atoms within rmin of this O
        h_neighbors = []
        for h_idx in H_indices:
            if dist_mat[o_idx, h_idx] < rmin:
                h_neighbors.append(h_idx)
        
        # Water molecule: O with exactly 2 H neighbors
        if len(h_neighbors) == 2:
            # Check these atoms haven't already been assigned to water
            if o_idx in water_indices_set:
                continue
            if h_neighbors[0] in water_indices_set or h_neighbors[1] in water_indices_set:
                continue
            
            # Create water molecule atoms with new molid and types
            o_atom = copy.deepcopy(atoms[o_idx])
            o_atom['molid'] = molid
            o_atom['type'] = 'Ow'
            o_atom['resname'] = 'SOL'
            o_atom['_orig_index'] = o_idx
            
            h1_atom = copy.deepcopy(atoms[h_neighbors[0]])
            h1_atom['molid'] = molid
            h1_atom['type'] = 'Hw'
            h1_atom['resname'] = 'SOL'
            h1_atom['_orig_index'] = h_neighbors[0]
            
            h2_atom = copy.deepcopy(atoms[h_neighbors[1]])
            h2_atom['molid'] = molid
            h2_atom['type'] = 'Hw'
            h2_atom['resname'] = 'SOL'
            h2_atom['_orig_index'] = h_neighbors[1]
            
            water_atoms.extend([o_atom, h1_atom, h2_atom])
            water_indices_set.add(o_idx)
            water_indices_set.add(h_neighbors[0])
            water_indices_set.add(h_neighbors[1])
            molid += 1
    
    # Create non-water atoms list
    non_water_atoms = [copy.deepcopy(atoms[i]) for i in range(len(atoms)) 
                       if i not in water_indices_set]
    
    # Update indices for non-water atoms
    for i, atom in enumerate(non_water_atoms):
        atom['index'] = i + 1
    
    # Update indices for water atoms
    for i, atom in enumerate(water_atoms):
        atom['index'] = i + 1
        # Remove temporary original index marker
        if '_orig_index' in atom:
            del atom['_orig_index']
    
    n_water = len(water_atoms) // 3
    print(f"Found {n_water} water molecules ({len(water_atoms)} atoms)")
    
    return water_atoms, non_water_atoms


def solvate(limits, density=1000.0, min_distance=2.0, max_solvent='max', 
           solute_atoms=None, solvent_type='spce', custom_solvent=None, custom_box=None,
           include_solute=False):
    """
    Solvate a structure or region with water or other solvent molecules.
    
    Parameters
    ----------
    limits : list of float
        Region limits [xlo, ylo, zlo, xhi, yhi, zhi] or [xhi, yhi, zhi] 
        (the latter assumes xlo=ylo=zlo=0)
    density : float, optional
        Solvent density in kg/m³, default is 1000.0 (water)
    min_distance : float or list of float, optional
        Minimum distance between solute and solvent atoms. 
        Can be a single value or [larger, smaller]
        where smaller applies to hydrogens and larger to other atoms.
        Default is 2.0 Å.
    max_solvent : int or str, optional
        Maximum number of solvent molecules, or 'max' for all possible,
        or 'shell' for a shell around solute. If 'shell', can be one of
        'shell10', 'shell15', 'shell20', 'shell25', or 'shell30' for different
        shell thicknesses in Å. Default is 'max'.
    solute_atoms : list of dict, optional
        Solute atoms to be solvated. Default is None.
    solvent_type : str, optional
        Type of solvent ('spce', 'tip3p', etc.). Default is 'spce'.
    custom_solvent : list of dict, optional
        Custom solvent atoms structure. Default is None.
    custom_box : list of float, optional
        Box dimensions for the custom solvent. Default is None.
    include_solute : bool, optional
        If True, return solute + solvent; if False (default), return solvent only.
        
    Returns
    -------
    list of dict
        Solvated structure (solvent only or solute + solvent)
        
    Notes
    -----
    - If solute_atoms is provided, atoms from the solvent that overlap with 
      the solute will be removed.
    - The 'shell' option creates a solvent shell of specified thickness around the solute.
    - Atom types in water molecules are typically 'OW', 'HW1', and 'HW2'.
    
    Examples
    --------
    # Solvate a box region:
    water = ap.solvate([0, 0, 0, 30, 30, 30])
    
    # Solvate around a protein with shell of 10 Å:
    water = ap.solvate([0, 0, 0, 60, 60, 60], solute_atoms=protein, max_solvent='shell10')
    """
    # Import functions here to avoid circular imports
    from .build import merge, slice as build_slice
    from .dist_matrix import dist_matrix
    from .cell_utils import Cell2Box_dim
    from . import import_conf
    
    # Standardize limits to [xlo, ylo, zlo, xhi, yhi, zhi] format
    if len(limits) == 3:
        xlo, ylo, zlo = 0, 0, 0
        xhi, yhi, zhi = limits
    elif len(limits) == 6:
        xlo, ylo, zlo, xhi, yhi, zhi = limits
    else:
        raise ValueError("Limits must be a list of length 3 [xhi, yhi, zhi] "
                         "or 6 [xlo, ylo, zlo, xhi, yhi, zhi]")
    
    # Calculate Box dimensions
    box_dim = [xhi - xlo, yhi - ylo, zhi - zlo]
    
    # Parse shell thickness if provided
    shell_thickness = None
    if isinstance(max_solvent, str) and max_solvent.startswith('shell'):
        if max_solvent == 'shell10':
            shell_thickness = 10.0
        elif max_solvent == 'shell15':
            shell_thickness = 15.0
        elif max_solvent == 'shell20':
            shell_thickness = 20.0
        elif max_solvent == 'shell25':
            shell_thickness = 25.0
        elif max_solvent == 'shell30':
            shell_thickness = 30.0
        elif max_solvent == 'shell':
            shell_thickness = 10.0  # Default
        else:
            try:
                # Try to extract thickness from string like 'shell12.5'
                shell_thickness = float(max_solvent[5:])
            except ValueError:
                shell_thickness = 10.0  # Default if parsing fails
    
    # Load solvent structure
    if custom_solvent is not None and custom_box is not None:
        solvent_atoms = copy.deepcopy(custom_solvent)
        solvent_box = custom_box
    else:
        # Load standard solvent
        solvent_atoms, solvent_box = _load_solvent(solvent_type)
    
    # Find how many atoms per solvent molecule (e.g., 3 for SPC water)
    unique_molids = set(atom['molid'] for atom in solvent_atoms)
    atoms_per_molecule = len(solvent_atoms) / len(unique_molids)
    
    # Calculate how many solvent molecules are needed to fill the Box
    # For water, 1000 kg/m³ is about 33.3 molecules per nm³
    volume_nm3 = (box_dim[0] / 10) * (box_dim[1] / 10) * (box_dim[2] / 10)
    molecules_per_nm3 = density / 30  # Approximate for water
    n_molecules_needed = int(volume_nm3 * molecules_per_nm3)
    
    if isinstance(max_solvent, int):
        n_molecules_needed = min(n_molecules_needed, max_solvent)
    
    # Calculate how many times to replicate the solvent Box
    nx = int(np.ceil((xhi - xlo) / solvent_box[0]))
    ny = int(np.ceil((yhi - ylo) / solvent_box[1]))
    nz = int(np.ceil((zhi - zlo) / solvent_box[2]))
    
    # Create the full solvent Box by replication
    full_solvent = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for atom in solvent_atoms:
                    new_atom = atom.copy()
                    new_atom['x'] = atom['x'] + ix * solvent_box[0] + xlo
                    new_atom['y'] = atom['y'] + iy * solvent_box[1] + ylo
                    new_atom['z'] = atom['z'] + iz * solvent_box[2] + zlo
                    # Adjust molid to keep track of different molecules
                    new_atom['molid'] = atom['molid'] + (ix * ny * nz + iy * nz + iz) * len(unique_molids)
                    full_solvent.append(new_atom)
    
    # Slice to get only atoms within the target region
    sliced_solvent = build_slice(full_solvent, [xlo, ylo, zlo, xhi, yhi, zhi])
    
    # Randomize the order of molecules for unbiased selection
    unique_molids = set(atom['molid'] for atom in sliced_solvent)
    molid_list = list(unique_molids)
    random.shuffle(molid_list)
    
    if shell_thickness is None and solute_atoms is None:
        # Just fill the Box
        if n_molecules_needed < len(molid_list):
            # Take only the required number of molecules
            selected_molids = set(molid_list[:n_molecules_needed])
            solvent_result = [atom for atom in sliced_solvent 
                             if atom['molid'] in selected_molids]
        else:
            solvent_result = sliced_solvent
    elif shell_thickness is not None and solute_atoms is not None:
        # Create a shell of solvent around the solute
        # First, merge solvent with solute, removing overlaps
        non_overlapping = merge(solute_atoms, sliced_solvent, box_dim, 
                               atom_label=['HW1', 'HW2'], 
                               min_distance=[min_distance, min_distance/2])
        
        # Now create a shell by keeping only molecules within the shell distance
        # but not closer than min_distance
        combined = solute_atoms + non_overlapping
        distances = dist_matrix(combined, box_dim)[0]
        
        # Extract solute-solvent distances
        solute_solvent_dist = distances[:len(solute_atoms), len(solute_atoms):]
        
        # Identify molecules in the shell
        shell_molids = set()
        for i, atom in enumerate(non_overlapping):
            dist_idx = i + len(solute_atoms)  # Index in the combined distance matrix
            min_dist = min(distances[:len(solute_atoms), dist_idx])
            if min_dist <= shell_thickness:
                shell_molids.add(atom['molid'])
        
        # Keep only molecules in the shell
        solvent_result = [atom for atom in non_overlapping 
                         if atom['molid'] in shell_molids]
        
        # If max_solvent is a number, limit the molecules
        if isinstance(max_solvent, int):
            unique_shell_molids = list(shell_molids)
            random.shuffle(unique_shell_molids)
            if max_solvent < len(unique_shell_molids):
                selected_molids = set(unique_shell_molids[:max_solvent])
                solvent_result = [atom for atom in solvent_result 
                                 if atom['molid'] in selected_molids]
    else:
        # Just remove overlapping solvent molecules
        solvent_result = merge(solute_atoms, sliced_solvent, box_dim,
                              atom_label=['HW1', 'HW2'],
                              min_distance=[min_distance, min_distance/2])
        
        # If max_solvent is a number, limit the molecules
        if isinstance(max_solvent, int):
            unique_molids = set(atom['molid'] for atom in solvent_result)
            molid_list = list(unique_molids)
            random.shuffle(molid_list)
            if max_solvent < len(molid_list):
                selected_molids = set(molid_list[:max_solvent])
                solvent_result = [atom for atom in solvent_result 
                                 if atom['molid'] in selected_molids]
    
    # Calculate statistics for output
    n_solvent_molecules = len(set(atom['molid'] for atom in solvent_result))
    n_solvent_atoms = len(solvent_result)
    
    # Calculate subvolume
    volume_angstrom3 = box_dim[0] * box_dim[1] * box_dim[2]
    
    # Print the solvation statistics
    print(f"  Subvolume: {box_dim[0]:.2f} x {box_dim[1]:.2f} x {box_dim[2]:.2f} Å³ = {volume_angstrom3:.2f} Å³")
    print(f"  Added {n_solvent_molecules} water molecules ({n_solvent_atoms} atoms)")
    
    # Combine solute and solvent if solute is provided
    if solute_atoms is not None:
        # Update molids to avoid conflicts
        max_solute_molid = max(atom['molid'] for atom in solute_atoms) if solute_atoms else 0
        for atom in solvent_result:
            atom['molid'] += max_solute_molid + 1
        
        # Combine
        result = solute_atoms + solvent_result if include_solute else solvent_result
    else:
        result = solvent_result
    
    return result


def _load_solvent(solvent_type='spce'):
    """
    Load a pre-equilibrated solvent Box.
    
    Parameters
    ----------
    solvent_type : str, optional
        Type of solvent ('spce', 'tip3p', etc.)
        
    Returns
    -------
    tuple
        (solvent_atoms, solvent_box)
    """
    from . import import_conf
    from .cell_utils import Cell2Box_dim
    
    # Define path to solvent structures
    base_path = os.path.dirname(os.path.abspath(__file__))
    structures_path = os.path.join(base_path, 'structures', 'water')
    
    # Map solvent type to filename
    solvent_files = {
        'spce': '864_spce.pdb',
        'spc': '864_spce.pdb',  # Alias
        'tip3p': '864_tip3p.pdb',
        'tip4p': '864_tip4p.pdb',
    }
    
    # Get the appropriate solvent file
    if solvent_type.lower() in solvent_files:
        solvent_file = os.path.join(structures_path, solvent_files[solvent_type.lower()])
    else:
        raise ValueError(f"Unsupported solvent type: {solvent_type}")
    
    # Import the solvent structure
    atoms, Cell = import_conf.pdb(solvent_file)
    box_dim = Cell2Box_dim(Cell)
    
    return atoms, box_dim

