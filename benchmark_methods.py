#!/usr/bin/env python3
"""Benchmark Direct Matrix vs Cell-List methods for a 6000 atom system."""
import time
import tracemalloc
import sys
import os
import numpy as np

# Setup paths
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

import atomipy as ap
from atomipy.dist_matrix import dist_matrix
from atomipy.cell_list_dist_matrix import cell_list_dist_matrix

input_file = "Clay_system.gro"
atoms, Box_dim = ap.import_auto(input_file)
atoms = ap.element(atoms)

print(f"Benchmarking {len(atoms)} atoms...")

def profile_call(label, func, *args, **kwargs):
    tracemalloc.start()
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    t1 = time.perf_counter()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    elapsed = t1 - t0
    peak_mb = peak / (1024 * 1024)
    print(f"[{label}] Time: {elapsed:.3f}s, Peak RAM: {peak_mb:.1f} MB")
    return elapsed, peak_mb

# Method 1: Direct Matrix (N^2)
profile_call("Direct Matrix (N^2)", dist_matrix, atoms, Box_dim)

# Method 2: Cell-List (N)
# Note: cell_list_dist_matrix also returns bond_list and dist_list
profile_call("Cell-List (O(N))", cell_list_dist_matrix, atoms, Box_dim, cutoff=2.45)
