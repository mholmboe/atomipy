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

# Find H atoms that are neighbors to multiple O atoms
h_to_o = {}
for o_idx, h_list in o_to_h.items():
    for h_idx in h_list:
        if h_idx not in h_to_o: h_to_o[h_idx] = []
        h_to_o[h_idx].append(o_idx)

shared_h = {h: o_list for h, o_list in h_to_o.items() if len(o_list) > 1}
if shared_h:
    print(f"Found {len(shared_h)} shared hydrogens!")
    for h, o_list in shared_h.items():
        print(f"  H {h+1} is neighbor to O: {[o+1 for o in o_list]}")
else:
    print("No shared hydrogens found.")

# Re-run find_H2O exact logic
water_found = []
water_indices_set = set()
for o_idx in sorted(o_to_h.keys()):
    h_neighbors = o_to_h[o_idx]
    if len(h_neighbors) == 2:
        h1_idx, h2_idx = h_neighbors
        if o_idx in water_indices_set or h1_idx in water_indices_set or h2_idx in water_indices_set:
            print(f"Skipping Ow {o_idx+1} because some atoms were already taken.")
            continue
        water_indices_set.add(o_idx); water_indices_set.add(h1_idx); water_indices_set.add(h2_idx)
        water_found.append(o_idx)

print(f"Logic result: found {len(water_found)} waters.")
