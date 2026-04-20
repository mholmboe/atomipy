import atomipy as ap
import numpy as np

def create_dummy_system(n):
    atoms = []
    for i in range(n):
        atoms.append({
            'index': i + 1,
            'type': 'O',
            'element': 'O',
            'x': np.random.random() * 20,
            'y': np.random.random() * 20,
            'z': np.random.random() * 20,
            'molid': 1
        })
    return atoms, [20.0, 20.0, 20.0]

print(f"Current SPARSE_THRESHOLD: {ap.config.SPARSE_THRESHOLD}")

# Test 1: System just below threshold
n_small = ap.config.SPARSE_THRESHOLD - 1
print(f"\nTesting system with {n_small} atoms (expected: Direct)...")
atoms_small, Box = create_dummy_system(n_small)
ap.bond_angle(atoms_small, Box)

# Test 2: System at threshold
n_large = ap.config.SPARSE_THRESHOLD
print(f"\nTesting system with {n_large} atoms (expected: Sparse)...")
atoms_large, Box = create_dummy_system(n_large)
ap.bond_angle(atoms_large, Box)

# Test 3: Changing threshold at runtime
print(f"\nChanging SPARSE_THRESHOLD to 1000...")
ap.config.SPARSE_THRESHOLD = 1000
print(f"New SPARSE_THRESHOLD: {ap.config.SPARSE_THRESHOLD}")

print(f"\nTesting system with 1500 atoms (expected: Sparse now)...")
atoms_new, Box = create_dummy_system(1500)
ap.bond_angle(atoms_new, Box)
