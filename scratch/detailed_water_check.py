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

# Identification logic from solvent.py
N = len(atoms)
is_oxygen = np.zeros(N, dtype=bool)
is_hydrogen = np.zeros(N, dtype=bool)
for i, atom in enumerate(atoms):
    element = atom.get('element', '')
    atom_type = atom.get('type', atom.get('atom_name', ''))
    is_o = (element and element[0].upper() == 'O') or (not element and atom_type and atom_type[0].upper() == 'O')
    is_h = (element and element[0].upper() == 'H') or (not element and atom_type and atom_type[0].upper() == 'H')
    if is_o: is_oxygen[i] = True
    elif is_h: is_hydrogen[i] = True

rmin = 1.30
from atomipy.dist_matrix import get_neighbor_list
i_idx, j_idx, dists, _, _, _ = get_neighbor_list(atoms, Box_dim, cutoff=rmin, rmaxH=rmin)

is_oh_pair = (is_oxygen[i_idx] & is_hydrogen[j_idx]) | (is_hydrogen[i_idx] & is_oxygen[j_idx])
oh_i = i_idx[is_oh_pair]
oh_j = j_idx[is_oh_pair]

o_to_h = {}
for k in range(len(oh_i)):
    idx1, idx2 = oh_i[k], oh_j[k]
    o_idx = idx1 if is_oxygen[idx1] else idx2
    h_idx = idx2 if is_oxygen[idx1] else idx1
    if o_idx not in o_to_h: o_to_h[o_idx] = []
    o_to_h[o_idx].append(h_idx)

print("\nOxygen atoms with != 2 Hydrogen neighbors within rmin=1.30:")
for o_idx in sorted(o_to_h.keys()):
    if not is_oxygen[o_idx]: continue # Should be True
    h_neighbors = o_to_h[o_idx]
    if len(h_neighbors) != 2:
        print(f"Ow {o_idx+1} has {len(h_neighbors)} H neighbors: {h_neighbors}")
        for h_idx in h_neighbors:
             # Find distance
             mask = ((i_idx == o_idx) & (j_idx == h_idx)) | ((i_idx == h_idx) & (j_idx == o_idx))
             d = dists[mask][0]
             print(f"  - H {h_idx+1} at distance {d:.4f} A")

# Check if any Ow are missing from o_to_h entirely
all_o_indices = np.where(is_oxygen)[0]
for o_idx in all_o_indices:
    if o_idx not in o_to_h:
        # Check if it has ANY H neighbors (maybe name matching failed?)
        print(f"Ow {o_idx+1} has NO H neighbors in o_to_h.")
