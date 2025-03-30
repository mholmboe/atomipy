# Atomipy: The atom Toolbox in Python

A modular Python toolbox for handling and analyzing molecular structures, particularly for mineral slabs with periodic boundary conditions. This toolbox is inspired by the MATLAB [atom Toolbox](https://github.com/mholmboe/atom) and provides similar functionality for working with PDB and Gromacs GRO files but with a streamlined Python interface.

## Overview

This toolbox is designed to import, export, and analyze molecular structures with a focus on mineral slabs containing different atom types (Si, Al, Fe, Mg, Ti, Li, F, O, H). It handles periodic and triclinic simulation cells, and provides functions for calculating bonds, angles, and distances while respecting periodic boundary conditions.

The molecular structure information is stored in dictionaries where each atom has fields for coordinates, neighbors, bonds, angles, element type, and more.

## Requirements

- NumPy
- SciPy (optional for some functionality)

## Key Features

- Import/export PDB and Gromacs GRO files
- Handle triclinic simulation cells with periodic boundary conditions
- Calculate bond distances and angles
- Element type assignment
- Coordination number analysis
- MINFF forcefield atom typing
- Distance matrices with PBC corrections
- Formal charge assignment for ions and water
- Coordinate transformations (orthogonal/triclinic, fractional/cartesian)
- Supercell replication with proper handling of triclinic cells

## Function Documentation

### File I/O

- `import_pdb(file_path)`: Import a PDB file, returning a list of atom dictionaries and cell parameters
- `import_gro(file_path)`: Import a Gromacs GRO file, including velocities if present
- `write_pdb(atoms, cell, file_path)`: Write atoms to a PDB file
- `write_gro(atoms, Box_dim, file_path)`: Write atoms to a Gromacs GRO file, including velocities if present

### Atom Properties

- `element(atom)`: Guess the chemical element for an atom based on its residue name
- `mass()`: Returns a dictionary of atomic masses for different elements
- `radius()`: Returns a dictionary of van der Waals radii for different elements

### Structure Analysis

- `dist_matrix(atoms, Box_dim)`: Calculate a full distance matrix between all atoms with PBC
- `bond_angle(atoms, Box_dim, rmaxH=1.2, rmaxM=2.45)`: Compute bonds and angles for a given atomic structure


### Coordinate Transformations

- `triclinic.orthogonal_to_triclinic(atoms, box_dim, angleparam)`: Convert atom coordinates from orthogonal to triclinic space
- `ortho.triclinic_to_orthogonal(atoms, box_dim)`: Convert atom coordinates from triclinic to orthogonal space
- `fract.cartesian_to_fractional(atoms, box_dim)`: Convert atom coordinates to fractional coordinates
- `fract.fractional_to_cartesian(atoms, box_dim)`: Convert fractional coordinates to cartesian
- `replicate.replicate_cell(atoms, box_dim, replicate)`: Create supercells by replicating in a, b, c directions

### Cell Utilities

- `cell_utils.Box_dim2Cell(box_dim)`: Convert box dimensions to cell parameters
- `cell_utils.Cell2Box_dim(cell)`: Convert cell parameters to box dimensions

### Charges

- `charge_formal.assign_formal_charges(atoms)`: Assign formal charges to ions and water molecules
- `charge_minff.charge_minff(atoms)`: Assign MINFF charges to atoms

### Force Field

- `minff(atom)`: Assign MINFF forcefield specific atom types to each atom

## Data Structure

All atomic information is stored in a list of dictionaries called `atoms`. Each atom dictionary contains the following fields:

- `molid`: Molecule ID
- `index`: Atom index
- `resname`: Residue name
- `x`, `y`, `z`: Coordinates
- `vx`, `vy`, `vz`: Velocities (if present)
- `neigh`: List of neighbor indices
- `bonds`: List of pairs `(j, distance)` where j is the index of a bonded atom
- `angles`: List of pairs `((j, k), angle)` where j,k are indices of atoms forming an angle with the central atom
- `element`: Chemical element
- `type`: Atom type
- `fftype`: Force field specific atom type
- `cn`: Coordination number

The simulation cell is represented in two ways:
- `Box_dim`: A 1x9 array used in Gromacs GRO files for triclinic cells
- `cell`: A 1x6 array [a, b, c, alpha, beta, gamma] used in PDB files

## Usage Examples

```python
# Import the entire package
import atomipy as ap
import numpy as np

# Import a PDB file
atoms, cell = ap.import_pdb("structure.pdb")

# Guess elements for all atoms
for i in range(len(atoms)):
    atoms[i] = ap.element(atoms[i])

# Convert cell parameters to box dimensions using the utility function
Box_dim = ap.cell_utils.Cell2Box_dim(cell)

# Calculate bonds and angles
atoms = ap.bond_angle(atoms, Box_dim)

# Assign formal charges for ions and water
atoms = ap.charge_formal.assign_formal_charges(atoms)

# Convert to fractional coordinates
frac_coords, atoms = ap.fract.cartesian_to_fractional(atoms, Box_dim)

# Create a 2x2x1 supercell
replicated_atoms, new_Box_dim = ap.replicate.replicate_cell(atoms, Box_dim, replicate=[2, 2, 1])

# Convert from triclinic to orthogonal coordinates if needed
ortho_atoms = ap.ortho.triclinic_to_orthogonal(replicated_atoms, new_Box_dim)


# Export to GRO format
ap.write_gro(replicated_atoms, new_Box_dim, "structure.gro")
```



## Differences from MATLAB atom Toolbox

This Python implementation is designed to provide similar functionality to the MATLAB atom Toolbox while following Python's conventions and making use of NumPy for efficient numerical operations. The data structure is dictionary-based rather than struct-based, and the function interfaces are designed for Python's style.

## License

This project is released under the MIT License.
