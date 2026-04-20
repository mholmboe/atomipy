#!/usr/bin/env python3
"""Final comparative log generation script."""
import sys
import os
import copy
import time
import numpy as np

# Setup paths
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

import atomipy as ap
from atomipy.dist_matrix import dist_matrix as dm_direct
from atomipy.cell_list_dist_matrix import cell_list_dist_matrix as dm_old
from atomipy.cell_list_dist_matrix_fast import cell_list_dist_matrix_fast as dm_fast

# Load the system
atoms_clay, Box_clay = ap.import_auto("Clay_system.gro")
atoms_clay = ap.element(atoms_clay)

# To produce the stats, we need to run the structural typing (minff) once
# Since the typing depends on the topology, we just need to ensure the matrix used is correct.
def generate_log_for_method(method_label, dm_func, log_name, needs_6=False):
    print(f"Generating stats for {method_label} -> {log_name}...")
    atoms = copy.deepcopy(atoms_clay)
    
    # Standard cutoffs
    rmaxH, rmaxM = 1.2, 2.45
    
    # 1. Calculate matrices
    if needs_6:
        dmat, dx, dy, dz, _, _ = dm_func(atoms, Box_clay, cutoff=max(rmaxH, rmaxM), rmaxH=rmaxH)
    else:
        dmat, dx, dy, dz = dm_func(atoms, Box_clay)
        
    # 2. Build topology using these EXACT matrices
    # (Extracting the core logic from bond_angle.py)
    types = np.array([atom.get('type', atom.get('name', '')) for atom in atoms])
    is_h = np.array([bool(t and t[0].upper() == 'H') for t in types])
    cutoff_matrix = np.where(is_h[:, np.newaxis] | is_h[np.newaxis, :], rmaxH, rmaxM)
    mask = (dmat > 0) & (dmat <= cutoff_matrix)
    
    # Populate the atom dictionaries
    for i in range(len(atoms)):
        atoms[i]['neigh'] = []
        atoms[i]['bonds'] = []
        atoms[i]['angles'] = []
        
        # Find neighbors
        neigh_indices = np.where(mask[i])[0]
        atoms[i]['neigh'] = list(neigh_indices)
        
        # Populate bonds
        for n_idx in neigh_indices:
            atoms[i]['bonds'].append((n_idx, dmat[i, n_idx]))
            
        # Populate angles
        for idx_a in range(len(neigh_indices)):
            for idx_b in range(idx_a + 1, len(neigh_indices)):
                n1, n2 = neigh_indices[idx_a], neigh_indices[idx_b]
                v1 = [dx[i, n1], dy[i, n1], dz[i, n1]]
                v2 = [dx[i, n2], dy[i, n2], dz[i, n2]]
                cos_val = np.dot(v1, v2) / (dmat[i, n1] * dmat[i, n2])
                angle = np.degrees(np.arccos(np.clip(cos_val, -1.0, 1.0)))
                atoms[i]['angles'].append(((n1, n2), angle))

    # 3. Generate stats log
    ap.get_structure_stats(atoms, Box_clay, log_file=log_name)
    print(f"  Done: {log_name}")

# Generate all 3
generate_log_for_method("Direct", dm_direct, "stats_direct.log")
generate_log_for_method("Old Cell-List", dm_old, "stats_old_cl.log", needs_6=True)
generate_log_for_method("New Fast Cell-List", dm_fast, "stats_fast_cl.log", needs_6=True)

print("\nVerification complete. All 3 logs are ready for comparison.")
