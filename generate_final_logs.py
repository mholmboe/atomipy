#!/usr/bin/env python3
"""Force generation of 3 full stats logs for each distance method."""
import sys
import os
import copy
import time
from unittest.mock import patch

# Setup paths
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

import atomipy as ap
import atomipy.bond_angle

# Load system
atoms_init, Box_dim = ap.import_auto("Clay_system.gro")
atoms_init = ap.element(atoms_init)

def generate_forced_stats(method_label, log_name):
    print(f"\n--- Generating Log: {log_name} ({method_label}) ---")
    atoms = copy.deepcopy(atoms_init)
    
    # We patch the distance functions to force the specific implementation
    from atomipy.dist_matrix import dist_matrix as dm_direct
    from atomipy.cell_list_dist_matrix import cell_list_dist_matrix as dm_old
    from atomipy.cell_list_dist_matrix_fast import cell_list_dist_matrix_fast as dm_fast
    
    # We will wrap ap.bond_angle and force it to use a specific matrix function
    # By patching the threshold check in bond_angle.py
    
    if method_label == "Direct":
        # Force the threshold to be huge so it uses Direct
        threshold_val = 1000000
    elif method_label == "Old CL":
        # For Old CL, we need to temporarily swap the 'cell_list_dist_matrix_fast' reference 
        # inside bond_angle.py to the old one.
        threshold_val = 0
    elif method_label == "Fast CL":
        threshold_val = 0
        
    def mock_threshold_check(a):
        return threshold_val
        
    # Patch the threshold check in bond_angle.py
    # Since I wrote 'if len(atoms) > 2000', I'll just change the atoms list size for the check
    # or better, just re-implement the bond_angle call logic here.
    
    rmaxH, rmaxM = 1.2, 2.45
    print(f"  Calculating matrices using {method_label}...")
    if method_label == "Direct":
        dmat, dx, dy, dz = dm_direct(atoms, Box_dim)
    elif method_label == "Old CL":
        dmat, dx, dy, dz, _, _ = dm_old(atoms, Box_dim, cutoff=max(rmaxH, rmaxM))
    else:
        dmat, dx, dy, dz, _, _ = dm_fast(atoms, Box_dim, cutoff=max(rmaxH, rmaxM), rmaxH=rmaxH)
        
    # Verify we got valid matrices
    print(f"  Matrices computed. Building topology...")
    
    # Now we must process these into atoms['neigh'], etc.
    # I'll call a special version of bond_angle logic that takes precomputed matrices
    # but since that doesn't exist, I'll just inject them into a mock and call get_structure_stats.
    # Actually, let's just use the verified identity and produce logs from the library.
    
    # To get DIFFERENT logs if there were differences, we'd need to run it.
    # I'll use the fact that I can temporarily edit bond_angle.py logic to force it.

# Actually, the most reliable way to give the user these 3 files IS 
# to run the code 3 times with different thresholds.

# Let's do it properly by patching the functions inside the module.
def run_full_pipeline(method_label, log_name):
    print(f"Processing {method_label} -> {log_name}")
    atoms = copy.deepcopy(atoms_init)
    
    with patch('atomipy.bond_angle.cell_list_dist_matrix_fast') as mock_fast, \
         patch('atomipy.bond_angle.cell_list_dist_matrix') as mock_old, \
         patch('atomipy.bond_angle.dist_matrix') as mock_direct:
        
        # Setup the mocks to behave like the real functions
        from atomipy.dist_matrix import dist_matrix as real_direct
        from atomipy.cell_list_dist_matrix import cell_list_dist_matrix as real_old
        from atomipy.cell_list_dist_matrix_fast import cell_list_dist_matrix_fast as real_fast
        
        # Ensure we use the 6k atoms so the theshold > 2000 is always triggered
        # for cell-list paths, and we'll force Direct by making it look like 100 atoms.
        
        if method_label == "Direct":
            # To force Direct, we make the system look small to the threshold check
            with patch('len', side_effect=lambda x: 100 if id(x) == id(atoms) else len(x)):
                atoms_processed = ap.minff(atoms, Box_dim)
        elif method_label == "Old CL":
            # Force the fast_cl call to use the old_cl implementation
            mock_fast.side_effect = real_old
            atoms_processed = ap.minff(atoms, Box_dim)
        else:
            # Default behavior (uses Fast CL)
            atoms_processed = ap.minff(atoms, Box_dim)
            
        ap.get_structure_stats(atoms_processed, Box_dim, log_file=log_name)
        print(f"Success: {log_name} generated.")

run_full_pipeline("Direct", "stats_direct.log")
run_full_pipeline("Old CL", "stats_old_cl.log")
run_full_pipeline("Fast CL", "stats_fast_cl.log")
