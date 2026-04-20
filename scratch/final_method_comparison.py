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

# Targeted Ow 
targeted = [6949, 6952, 6955, 6958, 6961, 6964]
indices = [t-1 for t in targeted]

print(f"Checking get_neighbor_list results for targeted Ow (cutoff=1.30)...")
i, j, d, dx, dy, dz = ap.dist_matrix.get_neighbor_list(atoms, Box_dim, cutoff=1.30)

for idx in indices:
    found = []
    for k in range(len(i)):
        if i[k] == idx:
            found.append((j[k]+1, d[k]))
        elif j[k] == idx:
            found.append((i[k]+1, d[k]))
    print(f"Ow {idx+1} neighbors: {found}")

# Total count check
print(f"\nTotal pairs found by Sparse: {len(i)}")

# Check with Direct for comparison
print("Checking with Direct...")
ap.config.SPARSE_THRESHOLD = 99999
i_d, j_d, d_d, dx_d, dy_d, dz_d = ap.dist_matrix.get_neighbor_list(atoms, Box_dim, cutoff=1.30)
print(f"Total pairs found by Direct: {len(i_d)}")

for idx in indices:
    found_d = []
    for k in range(len(i_d)):
        if i_d[k] == idx:
            found_d.append((j_d[k]+1, d_d[k]))
        elif j_d[k] == idx:
            found_d.append((i_d[k]+1, d_d[k]))
    print(f"Ow {idx+1} neighbors (Direct): {found_d}")
