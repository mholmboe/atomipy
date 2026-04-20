import os
import sys

# Add repo root to sys.path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import atomipy as ap

# Set threshold to ensure we use Sparse if needed
ap.config.SPARSE_THRESHOLD = 5000

file_path = 'DifferentSystemSizes/System10000.gro'
atoms, Box_dim = ap.import_auto(file_path)
print(f"Loaded {len(atoms)} atoms from file.")

SOL, noSOL = ap.find_H2O(atoms, Box_dim)
print(f"find_H2O found {len(SOL)} SOL atoms and {len(noSOL)} noSOL atoms.")
print(f"Total: {len(SOL) + len(noSOL)}")

noSOL = ap.assign_resname(noSOL)
IONS = [a for a in noSOL if a.get('resname') == 'ION']
MIN = [a for a in noSOL if a.get('resname') == 'MIN']
print(f"IONS: {len(IONS)}, MIN: {len(MIN)}")

System = ap.update(MIN, IONS, SOL)
print(f"Final System: {len(System)} atoms.")

# Check for lost atoms
lost = [a['index'] for a in noSOL if a.get('resname') not in ['ION', 'MIN']]
if lost:
    print(f"Lost {len(lost)} atoms from noSOL!")
    print(f"Lost atom resnames: {set(a.get('resname') for a in noSOL if a.get('resname') not in ['ION', 'MIN'])}")
    print(f"Lost atom types: {set(a.get('type') for a in noSOL if a.get('resname') not in ['ION', 'MIN'])}")
