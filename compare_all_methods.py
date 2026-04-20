#!/usr/bin/env python3
"""Comprehensive comparison of atomipy distance matrix methods and their impact on bonds/angles."""
import sys
import os
import numpy as np
import time

# Setup paths
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

import atomipy as ap
from atomipy.dist_matrix import dist_matrix as direct_dm
from atomipy.cell_list_dist_matrix import cell_list_dist_matrix as old_cl
from atomipy.cell_list_dist_matrix_fast import cell_list_dist_matrix_fast as new_cl

def compare_topology(atoms, Box, label):
    print(f"\n{'='*60}")
    print(f" BENCHMARK: {label} ({len(atoms)} atoms)")
    print(f"{'='*60}\n")

    results = {}
    
    # helper to run full bond_angle
    def run_full_pipeline(method_name, threshold_val):
        print(f"Running {method_name} (forcing threshold to {threshold_val})...")
        t0 = time.perf_counter()
        
        # We patch the threshold in bond_angle briefly just for this call
        # Or better yet, just mock the logic
        # For this test, I will manually call bond_angle and ensure it uses the specific method
        
        # Re-importing inside to ensure we can manipulate
        from atomipy.bond_angle import bond_angle as ba
        from atomipy.dist_matrix import dist_matrix as dm
        from atomipy.cell_list_dist_matrix import cell_list_dist_matrix as cl_old
        from atomipy.cell_list_dist_matrix_fast import cell_list_dist_matrix_fast as cl_new
        
        # Manual bypass logic
        rmaxH, rmaxM = 1.2, 2.45
        if method_name == "Direct":
            dmat, dx, dy, dz = dm(atoms, Box)
        elif method_name == "Old Cell-List":
            dmat, dx, dy, dz, _, _ = cl_old(atoms, Box, cutoff=max(rmaxH, rmaxM))
        else:
            dmat, dx, dy, dz, _, _ = cl_new(atoms, Box, cutoff=max(rmaxH, rmaxM), rmaxH=rmaxH)
            
        # Re-using the logic from bond_angle.py but forcing our matrices
        # (This avoids re-patching the library file)
        # However, it's easier to just call bond_angle and let it use the thresholds
        # with different atom counts if we wanted. 
        # But here we want to compare THE SAME 6k system.
        
        # We'll run the full topology extraction
        # We need the bond_angle logic to find angles accurately
        # I'll use a modified version of the bond_angle logic here
        
        # ... (simplified topology extraction same for all) ...
        types = np.array([atom.get('type', atom.get('name', '')) for atom in atoms])
        is_h = np.array([bool(t and t[0].upper() == 'H') for t in types])
        cutoff_matrix = np.where(is_h[:, np.newaxis] | is_h[np.newaxis, :], rmaxH, rmaxM)
        mask = (dmat > 0) & (dmat <= cutoff_matrix)
        ii, jj = np.where(np.triu(mask, k=1))
        bonds = np.column_stack((ii, jj))
        bdists = dmat[ii, jj]
        
        # To calculate angles:
        # We only really need to check the math of the angle calculation
        # Pick 500 random angles to compare
        angles = []
        for bidx in range(min(500, len(bonds))):
            i, j = bonds[bidx]
            # find neighbors of i
            neighs = np.where(mask[i])[0]
            for n in neighs:
                if n != j:
                    # Angle between j-i-n
                    v1 = [dx[i, j], dy[i, j], dz[i, j]]
                    v2 = [dx[i, n], dy[i, n], dz[i, n]]
                    cos_val = np.dot(v1, v2) / (dmat[i, j] * dmat[i, n])
                    ang = np.degrees(np.arccos(np.clip(cos_val, -1.0, 1.0)))
                    angles.append(ang)

        t1 = time.perf_counter()
        return {
            "time": t1-t0,
            "num_bonds": len(bonds),
            "bonds": bonds,
            "bdists": bdists,
            "num_angles": len(angles),
            "angles": np.array(angles)
        }

    res_direct = run_full_pipeline("Direct", 0)
    res_old = run_full_pipeline("Old Cell-List", 100000)
    res_new = run_full_pipeline("New Cell-List", 100000)

    def diff(name1, r1, name2, r2):
        print(f"\n--- Comparison: {name1} vs {name2} ---")
        print(f"  Bond Count: {r1['num_bonds']} vs {r2['num_bonds']}")
        
        if r1['num_bonds'] == r2['num_bonds']:
            dist_err = np.max(np.abs(r1['bdists'] - r2['bdists']))
            print(f"  Max Bond Dist Error: {dist_err:.2e}")
        else:
            print("  BOND COUNT MISMATCH - Skipping detail")
            
        if r1['num_angles'] == r2['num_angles'] and r1['num_angles'] > 0:
            ang_err = np.max(np.abs(r1['angles'] - r2['angles']))
            print(f"  Max Angle Value Error: {ang_err:.2e}")
        
        match = (r1['num_bonds'] == r2['num_bonds']) and (dist_err < 1e-5)
        print(f"  RESULT: {'MATCH' if match else 'MISMATCH'}")

    diff("Direct", res_direct, "Old Cell-List", res_old)
    diff("Direct", res_direct, "New Cell-List", res_new)
    diff("Old Cell-List", res_old, "New Cell-List", res_new)

    diff("Direct", "Old Cell-List")
    diff("Direct", "New Cell-List")
    diff("Old Cell-List", "New Cell-List")

# 1. Orthogonal
atoms_clay, Box_clay = ap.import_auto("Clay_system.gro")
atoms_clay = ap.element(atoms_clay)
# Sub-sample to 1000 atoms for faster matrix comparison if needed, but 6k is okay
compare_topology(atoms_clay, Box_clay, "Clay System (Orthogonal)")

# 2. Triclinic
kaol_file = "Kaolinite_GII_0.0487.gro"
if os.path.exists(kaol_file):
    atoms_k, Box_k = ap.import_auto(kaol_file)
    atoms_k = ap.element(atoms_k)
    compare_topology(atoms_k, Box_k, "Kaolinite (Triclinic)")
