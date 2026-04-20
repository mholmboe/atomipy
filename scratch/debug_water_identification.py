import atomipy as ap
import numpy as np

# Load the generated system
atoms, Box_dim = ap.import_gro("DifferentSystemSizes_v2/System2000.gro")
# Box_dim is already 1x9 for triclinic

print(f"Box_dim: {Box_dim}")
print(f"Total atoms: {len(atoms)}")

# Filter for waters (assuming resname SOL or types OW, HW1, HW2)
ow_atoms = [i for i, a in enumerate(atoms) if a.get('type', '').upper() == 'OW']
hw_atoms = [i for i, a in enumerate(atoms) if a.get('type', '').upper() in ['HW1', 'HW2', 'HW']]

print(f"Found {len(ow_atoms)} OW and {len(hw_atoms)} HW")

# Check distance of the first OW to all atoms with PBC
from atomipy.dist_matrix import get_neighbor_list
i_idx, j_idx, dists, _, _, _ = get_neighbor_list(atoms, Box_dim, cutoff=1.30)

# Count how many H neighbors each O has
o_neighbors = {i: [] for i in ow_atoms}
for k in range(len(i_idx)):
    i, j = i_idx[k], j_idx[k]
    d = dists[k]
    
    if i in ow_atoms and j in hw_atoms:
        o_neighbors[i].append(j)
    elif j in ow_atoms and i in hw_atoms:
        o_neighbors[j].append(i)

counts = [len(v) for v in o_neighbors.values()]
from collections import Counter
print(f"H-neighbor counts for OW atoms: {Counter(counts)}")

if 0 in o_neighbors:
    idx = 0
    ow = atoms[ow_atoms[0]]
    print(f"\nExample OW[0]: {ow}")
    print(f"Neighbors found: {o_neighbors[ow_atoms[0]]}")
    
    # Manually check distance to the first few HW
    for h_idx in hw_atoms[:10]:
        hw = atoms[h_idx]
        dx = ow['x'] - hw['x']
        dy = ow['y'] - hw['y']
        dz = ow['z'] - hw['z']
        
        # Simple MIC for orthogonal just to see
        if len(Box_dim) == 3:
            L = Box_dim
            dx -= L[0] * np.round(dx / L[0])
            dy -= L[1] * np.round(dy / L[1])
            dz -= L[2] * np.round(dz / L[2])
            dist = np.sqrt(dx**2 + dy**2 + dz**2)
            print(f"  Dist to HW[{h_idx}]: {dist:.3f}")
