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

# Targeted Ow 6949 and its Hydrogen neighbors 6951, 6950
o_idx = 6949 - 1
h1_idx = 6951 - 1
h2_idx = 6950 - 1

print(f"Investigating Ow {o_idx+1} and H {h1_idx+1}, {h2_idx+1}")

# Manual Cell calculation mirroring cell_list_dist_matrix.py
cutoff = 1.30
rmaxH = 1.30
positions = np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=np.float32)

from atomipy.cell_utils import normalize_box
Box_dim_v, Cell = normalize_box(Box_dim)
a, b, c = Cell[0], Cell[1], Cell[2]
alpha, beta, gamma = Cell[3], Cell[4], Cell[5]
ar, br, gr = np.radians([alpha, beta, gamma])
ax = a
bx = b * np.cos(gr)
by = b * np.sin(gr)
cx = c * np.cos(br)
cy = c * (np.cos(ar) - np.cos(br) * np.cos(gr)) / np.sin(gr)
cz = np.sqrt(max(0, c**2 - cx**2 - cy**2))
H = np.array([[ax, bx, cx], [0, by, cy], [0, 0, cz]], dtype=np.float32)
Hinv = np.linalg.inv(H)

frac_coords = (Hinv @ positions.T).T
frac_coords = frac_coords % 1.0

n_cells = np.maximum(np.floor(np.array([H[0,0], H[1,1], H[2,2]]) / 1.30), 1).astype(int)
cell_idx = np.floor(frac_coords * n_cells).astype(int)
cell_idx = np.clip(cell_idx, 0, n_cells - 1)

print(f"Matrix H diagonal: {np.diag(H)}")
print(f"n_cells: {n_cells}")

o_cell = cell_idx[o_idx]
h1_cell = cell_idx[h1_idx]
h2_cell = cell_idx[h2_idx]

print(f"Ow {o_idx+1} cell: {o_cell} (frac: {frac_coords[o_idx]})")
print(f"H {h1_idx+1} cell: {h1_cell} (frac: {frac_coords[h1_idx]})")
print(f"H {h2_idx+1} cell: {h2_cell} (frac: {frac_coords[h2_idx]})")

# Check if H cells are within 27 neighbors of O cell
def is_neighbor(c1, c2, n_cells):
    diff = (c2 - c1 + n_cells // 2) % n_cells - n_cells // 2 # WRONG wrapping
    # Better:
    d = []
    for i in range(3):
        delta = c2[i] - c1[i]
        # Handle PBC
        if delta > n_cells[i] // 2: delta -= n_cells[i]
        elif delta < -n_cells[i] // 2: delta += n_cells[i]
        d.append(delta)
    return all(abs(x) <= 1 for x in d), d

is_h1_neigh, d1 = is_neighbor(o_cell, h1_cell, n_cells)
is_h2_neigh, d2 = is_neighbor(o_cell, h2_cell, n_cells)

print(f"H 6951 is neighbor of Ow 6949? {is_h1_neigh} (delta: {d1})")
print(f"H 6950 is neighbor of Ow 6949? {is_h2_neigh} (delta: {d2})")
