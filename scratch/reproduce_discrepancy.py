import os
import sys
import numpy as np

# Add repo root to sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import atomipy as ap

file_path = 'DifferentSystemSizes/System10000.gro'
atoms, Box_dim = ap.import_auto(file_path)
print(f"Loaded {len(atoms)} atoms.")

# Test find_H2O with Direct method
print("\nTesting Direct method:")
ap.config.SPARSE_THRESHOLD = 99999  # Force Direct
SOL_d, noSOL_d = ap.find_H2O(atoms, Box_dim)

# Test find_H2O with Sparse method
print("\nTesting Sparse method:")
ap.config.SPARSE_THRESHOLD = 1000   # Force Sparse
SOL_s, noSOL_s = ap.find_H2O(atoms, Box_dim)

print(f"\nResults Comparison:")
print(f"Direct: {len(SOL_d)} SOL atoms")
print(f"Sparse: {len(SOL_s)} SOL atoms")

if len(SOL_d) != len(SOL_s):
    print("!!! DISCREPANCY DETECTED between Direct and Sparse methods !!!")
    
    # Identify which atoms are in SOL_d but not in SOL_s
    indices_d = set(a['index'] for a in SOL_d)
    indices_s = set(a['index'] for a in SOL_s)
    
    diff = indices_d - indices_s
    if diff:
        print(f"Indices in Direct but NOT in Sparse: {sorted(list(diff))[:30]}...")
        # Check their coordinates
        for idx in sorted(list(diff))[:6]:
            atom = next(a for a in SOL_d if a['index'] == idx)
            print(f"  Atom {idx}: {atom['type']} at ({atom['x']:.3f}, {atom['y']:.3f}, {atom['z']:.3f})")
    
    diff_inv = indices_s - indices_d
    if diff_inv:
        print(f"Indices in Sparse but NOT in Direct: {sorted(list(diff_inv))[:30]}...")
else:
    print("Matches perfectly.")
