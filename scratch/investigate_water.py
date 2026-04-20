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

# Identify atoms named Ow and Hw (or similar)
is_o = [i for i, a in enumerate(atoms) if a['type'].startswith('Ow')]
is_h = [i for i, a in enumerate(atoms) if a['type'].startswith('Hw') or a['type'] == 'H']

print(f"Total potential water atoms: {len(is_o)} O, {len(is_h)} H")

# Run find_H2O to see what's excluded
SOL, noSOL = ap.find_H2O(atoms, Box_dim, rmin=1.30)
sol_indices = set(a['index'] for a in SOL)

# Find Ow atoms that are NOT in SOL
missed_o = [i for i in is_o if (i+1) not in sol_indices]
print(f"\nMissed {len(missed_o)} Ow atoms from find_H2O:")

# For each missed Oxygen, find its closest Hydrogens
from atomipy.dist_matrix import get_neighbor_list
i_idx, j_idx, dists, dx_s, dy_s, dz_s = get_neighbor_list(atoms, Box_dim, cutoff=5.0)

for o_idx in missed_o:
    atom_o = atoms[o_idx]
    print(f"\nOw {o_idx+1} at ({atom_o['x']:.3f}, {atom_o['y']:.3f}, {atom_o['z']:.3f})")
    
    # Find all hydrogens close to this oxygen
    neighbors = []
    for k in range(len(i_idx)):
        if i_idx[k] == o_idx:
            neighbor_idx = j_idx[k]
            if neighbor_idx in is_h:
                neighbors.append((neighbor_idx, dists[k]))
        elif j_idx[k] == o_idx:
            neighbor_idx = i_idx[k]
            if neighbor_idx in is_h:
                neighbors.append((neighbor_idx, dists[k]))
    
    neighbors.sort(key=lambda x: x[1])
    print(f"  Closest hydrogens:")
    for h_idx, d in neighbors[:5]:
        atom_h = atoms[h_idx]
        print(f"    H {h_idx+1} ({atom_h['type']}) at distance {d:.4f} A")
