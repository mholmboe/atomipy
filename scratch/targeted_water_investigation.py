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

targeted_o = [6949, 6952, 6955, 6958, 6961, 6964] # Indices (1-based) from previous log
indices = [i-1 for i in targeted_o]

print(f"Targeted Oxygen indices (0-based): {indices}")

# Test with various cutoffs
for cutoff in [1.25, 1.30, 2.0, 5.0]:
    print(f"\nTesting cutoff = {cutoff}:")
    from atomipy.dist_matrix import get_neighbor_list
    i_idx, j_idx, dists, _, _, _ = get_neighbor_list(atoms, Box_dim, cutoff=cutoff)
    
    for o_idx in indices:
        found = False
        for k in range(len(i_idx)):
            if i_idx[k] == o_idx or j_idx[k] == o_idx:
                neigh_idx = j_idx[k] if i_idx[k] == o_idx else i_idx[k]
                neigh_type = atoms[neigh_idx]['type']
                if neigh_type.startswith('H'):
                    print(f"  Ow {o_idx+1} found H {neigh_idx+1} at distance {dists[k]:.4f}")
                    found = True
        if not found:
            print(f"  Ow {o_idx+1} NOT FOUND any H neighbors!")
