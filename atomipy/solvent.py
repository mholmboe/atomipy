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
from . import config


def find_H2O(atoms, Box_dim=None, rmin=1.30):
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
    from .distances import neighbor_list_fast
    
    if not atoms:
        return [], []
    
    atoms = copy.deepcopy(atoms)
    
    # 1. Identify types/elements for all atoms once
    N = len(atoms)
    is_oxygen = np.zeros(N, dtype=bool)
    is_hydrogen = np.zeros(N, dtype=bool)
    
    for i, atom in enumerate(atoms):
        element = atom.get('element', '')
        atom_type = atom.get('type', atom.get('atom_name', ''))
        
        # Check element first (can be 'O', 'Ow', 'Oh', etc), then type
        is_o = (element and element[0].upper() == 'O') or \
               (not element and atom_type and atom_type[0].upper() == 'O')
        is_h = (element and element[0].upper() == 'H') or \
               (not element and atom_type and atom_type[0].upper() == 'H')
        
        if is_o:
            is_oxygen[i] = True
        elif is_h:
            is_hydrogen[i] = True
    
    if not np.any(is_oxygen) or np.sum(is_hydrogen) < 2:
        print("Found 0 water molecules (0 atoms)")
        return [], atoms
    
    # 2. Use the central dispatcher to find all pairs within rmin
    # It automatically selects between Direct and Sparse based on config.SPARSE_THRESHOLD
    from .distances import get_neighbor_list
    i_idx, j_idx, dists, _, _, _ = get_neighbor_list(atoms, Box_dim, cutoff=rmin, rmaxH=rmin)
    
    # 3. Filter for O-H pairs
    # Since neighbor_list_fast returns i < j, we check both (i=O, j=H) and (i=H, j=O)
    is_oh_pair = (is_oxygen[i_idx] & is_hydrogen[j_idx]) | (is_hydrogen[i_idx] & is_oxygen[j_idx])
    
    oh_i = i_idx[is_oh_pair]
    oh_j = j_idx[is_oh_pair]
    oh_dists = dists[is_oh_pair]
    
    # Map H indices to their closest O neighbor
    h_to_o = {}
    for k in range(len(oh_i)):
        idx1, idx2 = oh_i[k], oh_j[k]
        o_idx = idx1 if is_oxygen[idx1] else idx2
        h_idx = idx2 if is_oxygen[idx1] else idx1
        
        dist = oh_dists[k]
        if h_idx not in h_to_o or dist < h_to_o[h_idx][1]:
            h_to_o[h_idx] = (o_idx, dist)
    
    # Map O indices back to their assigned H neighbors
    o_to_h = {}
    for h_idx, (o_idx, dist) in h_to_o.items():
        if o_idx not in o_to_h:
            o_to_h[o_idx] = []
        o_to_h[o_idx].append(h_idx)
    
    # 4. Find water molecules: O with exactly 2 H neighbors
    water_atoms = []
    water_indices_set = set()
    molid = 1
    
    # Sort O indices for deterministic results
    for o_idx in sorted(o_to_h.keys()):
        h_neighbors = o_to_h[o_idx]
        
        if len(h_neighbors) == 2:
            h1_idx, h2_idx = h_neighbors
            
            # Check these atoms haven't already been assigned to water
            if o_idx in water_indices_set or h1_idx in water_indices_set or h2_idx in water_indices_set:
                continue
            
            # Create water molecule atoms
            for idx, role in [(o_idx, 'Ow'), (h1_idx, 'Hw'), (h2_idx, 'Hw')]:
                atom = copy.deepcopy(atoms[idx])
                atom['molid'] = molid
                atom['type'] = role
                atom['resname'] = 'SOL'
                water_atoms.append(atom)
                water_indices_set.add(idx)
            
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
           solute_atoms=None, Box=None, solvent_type='spce', custom_solvent=None, custom_box=None,
           include_solute=False):
    """
    Solvate a structure or region with water or other solvent molecules.
    
    Performance Note
    ----------------
    This function uses the central dispatcher `get_neighbor_list` for shell distance 
    calculations, automatically switching to a memory-efficient sparse neighbor list 
    for systems larger than `config.SPARSE_THRESHOLD` to prevent O(N^2) memory bottlenecks.
    
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
    from .distances import dist_matrix, get_neighbor_list
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
        non_overlapping = merge(solute_atoms, sliced_solvent, Box if Box is not None else box_dim, 
                               atom_label=['HW1', 'HW2'], 
                               min_distance=[min_distance, min_distance/2])
        
        # Now create a shell by keeping only molecules within the shell distance
        # but not closer than min_distance
        combined = solute_atoms + non_overlapping
        
        # Use the central dispatcher to get neighbors and distances
        # It automatically handles Direct vs Sparse based on config.SPARSE_THRESHOLD
        n_solute = len(solute_atoms)
        i_idx, j_idx, dists, _, _, _ = get_neighbor_list(
            combined, 
            box_dim, 
            cutoff=shell_thickness, 
            dm_method=dm_method if 'dm_method' in locals() else None
        )
        
        # Filter for solute-solvent pairs
        # i_idx < j_idx holds in the dispatcher. We need i < n_solute and j >= n_solute
        mask = (i_idx < n_solute) & (j_idx >= n_solute)
        
        # Map which solvent indices are within reach
        shell_indices = set(j_idx[mask] - n_solute)
        shell_molids = set(non_overlapping[idx]['molid'] for idx in shell_indices)
        
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
        solvent_result = merge(solute_atoms, sliced_solvent, Box if Box is not None else box_dim,
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

def spc2tip4p(atoms, Box=None, om_dist=0.15):
    """
    Convert SPC water molecules to TIP4P model by adding an M dummy site.
    
    Parameters
    ----------
    atoms : list of dict
        List of atom dictionaries (should already be identified as water).
    Box : list of float, optional
        Box dimensions for PBC handling.
    om_dist : float, optional
        Distance from oxygen to the M site along the bisector (default 0.15 A).
        
    Returns
    -------
    list of dict
        New list of atoms with M sites added.
    """
    # Group by molid
    molids = np.array([a['molid'] for a in atoms])
    unique_molids = np.unique(molids)
    
    new_atoms = []
    
    for mid in unique_molids:
        mol_indices = np.where(molids == mid)[0]
        if len(mol_indices) != 3:
            # Not a 3-site water molecule, copy as is
            for idx in mol_indices:
                new_atoms.append(copy.deepcopy(atoms[idx]))
            continue
            
        # Find O and Hs
        o_idx = -1
        h_indices = []
        for idx in mol_indices:
            atom = atoms[idx]
            atype = atom.get('type', '').upper()
            element = atom.get('element', '').upper()
            if 'O' in atype or element == 'O':
                o_idx = idx
            elif 'H' in atype or element == 'H':
                h_indices.append(idx)
        
        if o_idx == -1 or len(h_indices) != 2:
            # Not a standard water structure
            for idx in mol_indices:
                new_atoms.append(copy.deepcopy(atoms[idx]))
            continue
            
        o_atom = copy.deepcopy(atoms[o_idx])
        h1_atom = copy.deepcopy(atoms[h_indices[0]])
        h2_atom = copy.deepcopy(atoms[h_indices[1]])
        
        # Add O and Hs to new list
        new_atoms.extend([o_atom, h1_atom, h2_atom])
        
        # Calculation of M site
        o_pos = np.array([o_atom['x'], o_atom['y'], o_atom['z']])
        h1_pos = np.array([h1_atom['x'], h1_atom['y'], h1_atom['z']])
        h2_pos = np.array([h2_atom['x'], h2_atom['y'], h2_atom['z']])
        
        v1 = h1_pos - o_pos
        v2 = h2_pos - o_pos
        
        if Box is not None:
            # PBC correction for vectors
            if len(Box) == 3:
                v1 -= np.round(v1 / np.array(Box)) * np.array(Box)
                v2 -= np.round(v2 / np.array(Box)) * np.array(Box)
            elif len(Box) == 9:
                lx, ly, lz = Box[0], Box[1], Box[2]
                v1[0] -= lx * np.round(v1[0] / lx)
                v1[1] -= ly * np.round(v1[1] / ly)
                v1[2] -= lz * np.round(v1[2] / lz)
                v2[0] -= lx * np.round(v2[0] / lx)
                v2[1] -= ly * np.round(v2[1] / ly)
                v2[2] -= lz * np.round(v2[2] / lz)

        # Normalized vectors
        v1_u = v1 / np.linalg.norm(v1)
        v2_u = v2 / np.linalg.norm(v2)
        
        # Bisector
        bisector = v1_u + v2_u
        bisector_u = bisector / np.linalg.norm(bisector)
        
        m_pos = o_pos + om_dist * bisector_u
        
        # Create M atom
        m_atom = copy.deepcopy(o_atom)
        m_atom['x'], m_atom['y'], m_atom['z'] = m_pos
        m_atom['type'] = 'MW'
        m_atom['element'] = 'M'
        m_atom['charge'] = -1.0484  # Example charge for TIP4P
        m_atom['index'] = 0 # Will be updated by update()
        
        # Adjust O and H charges for TIP4P if needed, 
        # but usually the user handles forcefield assignment separately.
        # However, for convenience we can set some defaults.
        o_atom['charge'] = 0.0
        h1_atom['charge'] = 0.5242
        h2_atom['charge'] = 0.5242
        
        new_atoms.append(m_atom)
        
    from .build import update
    return update(new_atoms)

def tip3p2tip4p(atoms, Box=None):
    """Alias for spc2tip4p using default TIP4P distance."""
    return spc2tip4p(atoms, Box, om_dist=0.15)

