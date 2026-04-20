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
targeted = [6949, 6950, 6951]
for idx_p1 in targeted:
    i = idx_p1 - 1
    atom = atoms[i]
    element = atom.get('element', '')
    atom_type = atom.get('type', atom.get('atom_name', ''))
    is_o = (element and element[0].upper() == 'O') or (not element and atom_type and atom_type[0].upper() == 'O')
    is_h = (element and element[0].upper() == 'H') or (not element and atom_type and atom_type[0].upper() == 'H')
    print(f"Atom {idx_p1}: type={atom_type}, element={element}, is_o={is_o}, is_h={is_h}")

# Check find_H2O internal state
rmin = 1.30
from atomipy.dist_matrix import get_neighbor_list
i_idx, j_idx, dists, _, _, _ = get_neighbor_list(atoms, Box_dim, cutoff=rmin, rmaxH=rmin)

# Count H neighbors for 6949
o_idx = 6949 - 1
h_count = 0
for k in range(len(i_idx)):
    u, v = i_idx[k], j_idx[k]
    if u == o_idx or v == o_idx:
        neigh = v if u == o_idx else u
        # Check if neigh is H
        n_atom = atoms[neigh]
        n_element = n_atom.get('element', '')
        n_type = n_atom.get('type', n_atom.get('atom_name', ''))
        n_is_h = (n_element and n_element[0].upper() == 'H') or (not n_element and n_type and n_type[0].upper() == 'H')
        if n_is_h:
            h_count += 1
            print(f"  Found H neighbor {neigh+1} at distance {dists[k]:.4f}")

print(f"Total H neighbors for Ow 6949: {h_count}")
