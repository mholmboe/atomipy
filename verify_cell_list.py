#!/usr/bin/env python3
"""Compare old and new cell-list implementations for accuracy and speed."""
import time
import sys
import os
import numpy as np

# Setup paths
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

import atomipy as ap
from atomipy.cell_list_dist_matrix import cell_list_dist_matrix as old_cl
from atomipy.cell_list_dist_matrix_fast import cell_list_dist_matrix_fast as new_cl

def run_comparison(label, atoms, Box):
    print(f"\n--- {label} ({len(atoms)} atoms) ---")
    
    # Old Cell-List
    t0 = time.perf_counter()
    d_old, x_old, y_old, z_old, b_old, dist_old = old_cl(atoms, Box, cutoff=2.45)
    t1 = time.perf_counter()
    print(f"Old Cell-List: {t1-t0:.3f}s")
    
    # New Cell-List
    t2 = time.perf_counter()
    d_new, x_new, y_new, z_new, bonds, dists = new_cl(atoms, Box, cutoff=2.45)
    t3 = time.perf_counter()
    print(f"New Cell-List: {t3-t2:.3f}s")
    
    # Accuracy Check
    max_err = np.max(np.abs(d_old - d_new))
    print(f"Max distance error: {max_err:.2e}")
    
    # Check bond counts
    old_bonds = np.sum(d_old > 0) // 2
    new_bonds = len(bonds)
    print(f"Bond counts: Old={old_bonds}, New={new_bonds}")
    
    if max_err > 1e-5:
        print("FAILED: Accuracy mismatch!")
    else:
        print("PASSED: Accuracy check.")

# 1. Orthogonal System
atoms_clay, Box_clay = ap.import_auto("Clay_system.gro")
atoms_clay = ap.element(atoms_clay)
run_comparison("Orthogonal (Clay)", atoms_clay, Box_clay)

# 2. Triclinic System (Fake it if needed, or use a known one)
# Kaolinite might be triclinic
kaol_file = "/Users/miho0052/Dropbox/Coding/Windsurf/atominpython/Kaolinite_GII_0.0487.gro"
if os.path.exists(kaol_file):
    atoms_k, Box_k = ap.import_auto(kaol_file)
    atoms_k = ap.element(atoms_k)
    run_comparison("Triclinic (Kaolinite)", atoms_k, Box_k)
else:
    # Create synthetic triclinic box
    print("\nCreating synthetic triclinic system...")
    pos = np.random.rand(1000, 3) * 20.0
    syn_atoms = [{'x': p[0], 'y': p[1], 'z': p[2], 'type': 'O'} for p in pos]
    # Box: a, b, c, alpha, beta, gamma
    Box_syn = [30, 30, 30, 75, 85, 95]
    run_comparison("Synthetic Triclinic", syn_atoms, Box_syn)

