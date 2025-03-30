import sys
import os

# Add the project directory to Python path for importing atomipy
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import atomipy as ap
import numpy as np

# Import the GRO file
print("Importing preem.gro...")
atoms, box_dim = ap.import_gro('preem.gro')

# Print some information about the atoms
print(f"Successfully imported {len(atoms)} atoms")
print(f"Box dimensions: {box_dim}")

# Print the first few atoms
print("\nFirst 5 atoms:")
for i, atom in enumerate(atoms[:5]):
    print(f"Atom {i+1}: {atom['resname']} {atom.get('element', atom.get('type', '?'))} at ({atom['x']:.3f}, {atom['y']:.3f}, {atom['z']:.3f})")

# Try to guess elements for all atoms
print("\nGuessing elements...")
for i in range(len(atoms)):
    atoms[i] = ap.element(atoms[i])

# Print the first few atoms again with elements
print("\nFirst 5 atoms with elements:")
for i, atom in enumerate(atoms[:5]):
    print(f"Atom {i+1}: {atom['resname']} {atom.get('element', '?')} at ({atom['x']:.3f}, {atom['y']:.3f}, {atom['z']:.3f})")
