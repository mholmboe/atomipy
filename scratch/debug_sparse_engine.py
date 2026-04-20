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

# Manual Sparse Engine logic for debugging
N = len(atoms)
positions = np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=np.float32)
types = np.array([atom.get('type', '') for atom in atoms])
is_h = np.array([bool(t and t[0].upper() == 'H') for t in types])

Box_dim_v, Cell = ap.cell_utils.normalize_box(Box_dim)
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
frac_coords = (Hinv @ positions.T).T % 1.0

max_cutoff = 1.30
n_cells = np.maximum(np.floor(np.diag(H) / max_cutoff), 1).astype(int)
cell_idx = np.floor(frac_coords * n_cells).astype(int)
cell_idx = np.clip(cell_idx, 0, n_cells - 1)
flat_idx = cell_idx[:, 0] * (n_cells[1] * n_cells[2]) + cell_idx[:, 1] * n_cells[2] + cell_idx[:, 2]

num_cells_total = np.prod(n_cells)
head = np.full(num_cells_total, -1, dtype=np.int32)
next_atom = np.full(N, -1, dtype=np.int32)
for i in range(N):
    c = flat_idx[i]
    next_atom[i] = head[c]
    head[c] = i

offsets = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)

target_o = 6949-1
target_h = 6951-1

print(f"Target: O {target_o+1}, H {target_h+1}")

for cx in range(n_cells[0]):
    for cy in range(n_cells[1]):
        for cz_idx in range(n_cells[2]):
            c1 = cx * (n_cells[1] * n_cells[2]) + cy * n_cells[2] + cz_idx
            i_head = head[c1]
            if i_head == -1: continue
            
            atoms1 = []
            curr = i_head
            while curr != -1:
                atoms1.append(curr)
                curr = next_atom[curr]
            atoms1 = np.array(atoms1, dtype=np.int32)
            
            if target_o not in atoms1 and target_h not in atoms1: continue
            
            for off in offsets:
                nx = (cx + off[0]) % n_cells[0]
                ny = (cy + off[1]) % n_cells[1]
                nz = (cz_idx + off[2]) % n_cells[2]
                
                c2 = nx * (n_cells[1] * n_cells[2]) + ny * n_cells[2] + nz
                j_head = head[c2]
                if j_head == -1: continue
                
                atoms2 = []
                curr = j_head
                while curr != -1:
                    atoms2.append(curr)
                    curr = next_atom[curr]
                atoms2 = np.array(atoms2, dtype=np.int32)
                
                if (target_o in atoms1 and target_h in atoms2) or (target_h in atoms1 and target_o in atoms2):
                    print(f"\nProcessing cell pair {c1} vs {c2} (offset {off})")
                    p1 = positions[atoms1]
                    p2 = positions[atoms2]
                    diff = p2[np.newaxis, :, :] - p1[:, np.newaxis, :]
                    diff_frac = (Hinv @ diff.reshape(-1, 3).T).T
                    diff_frac = diff_frac - np.round(diff_frac)
                    diff_cart = (H @ diff_frac.T).T.reshape(len(atoms1), len(atoms2), 3)
                    d = np.sqrt(np.sum(diff_cart**2, axis=2))
                    
                    ii, jj = np.where(d <= max_cutoff)
                    for k in range(len(ii)):
                        idx1 = atoms1[ii[k]]
                        idx2 = atoms2[jj[k]]
                        if (idx1 == target_o and idx2 == target_h) or (idx1 == target_h and idx2 == target_o):
                            print(f"  FOUND pair ({idx1+1}, {idx2+1}) with d={d[ii[k], jj[k]]:.4f}")
                            if not (idx1 < idx2):
                                print(f"  REJECTED because idx1 ({idx1+1}) >= idx2 ({idx2+1})")
                            else:
                                print(f"  ACCEPTED!")
