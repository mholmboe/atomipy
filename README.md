# Atomipy: The atom Toolbox in Python

A modular Python toolbox for handling and analyzing molecular structures, particularly for mineral slabs with periodic boundary conditions. This toolbox is a light version of the MATLAB [atom Toolbox](https://github.com/mholmboe/atom) and can mainly be used to generate molecular topology files for the [**MINFF**](https://github.com/mholmboe/minff) forcefield with a streamlined Python interface.

## Overview

This toolbox is designed to import, export, and analyze molecular structures with a focus on mineral slabs containing different elements (Si, Al, Fe, Mg, Ca, Ti, Li, F, O, H). It handles periodic and triclinic simulation cells, and provides functions for calculating bonds, angles, and distances while respecting periodic boundary conditions, and hence is ideal for generating molecular topology files for mineral bulk/slab systems that can be modelled using the [**MINFF**](https://github.com/mholmboe/minff) forcefield.

The molecular structure information is stored in dictionaries where each atom has fields for coordinates, neighbors, bonds, angles, element type, and more.

## Requirements

- NumPy (>=1.18.0)
- tqdm (>=4.45.0) - for progress bars
- Numba (>=0.50.0, optional) - for performance optimization via JIT compilation

## Key Features

- MINFF forcefield atom typing
- Import/export PDB and Gromacs GRO files
- Generating topology files for MINFF forcefield, for Gromacs (.itp), NAMD (.psf) and LAMMPS (.data)
- Handle both orthogonal and triclinic simulation cells with periodic boundary conditions
- Calculate bond distances and angles
- Element type assignment
- Coordination number analysis
- Distance matrices with PBC corrections (using both full matrix and efficient cell-list algorithms)
- Progress tracking for computationally intensive calculations
- Formal charge assignment for ions and water
- Coordinate transformations (orthogonal/triclinic, fractional/cartesian)
- Supercell replication with proper handling of triclinic cells
- Support for case-insensitive element matching in charge assignments

## Function Documentation

### File I/O

- `import_pdb(file_path)`: Import a PDB file, returning a list of atom dictionaries and cell parameters
- `import_gro(file_path)`: Import a Gromacs GRO file, including velocities if present
- `write_pdb(atoms, cell, file_path)`: Write atoms to a PDB file
- `write_gro(atoms, Box_dim, file_path)`: Write atoms to a Gromacs GRO file, including velocities if present

### Force Field

- `minff(atom)`: Assign MINFF forcefield specific atom types to each atom

### Molecular Topology

- `write_itp(atoms, Box_dim, file_path)`: Write a Gromacs topology file
- `write_psf(atoms, Box_dim, file_path)`: Write a NAMD topology file
- `write_data(atoms, Box_dim, file_path)`: Write a LAMMPS topology file

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

## Example Scripts

### run_atomi.py

The `run_atomi.py` script demonstrates a comprehensive workflow for processing mineral structures using atomipy. This script serves as an excellent starting point for users new to the package.

#### What the script does:

1. **Imports a GRO structure file** (Kaolinite_GII_0.0487.gro) containing mineral coordinates and box dimensions
2. **Assigns chemical elements** to each atom using chemical knowledge-based rules
3. **Creates a supercell** by replicating the unit cell to reach a target size (40 Å in each dimension)
4. **Saves the replicated structure** in both GRO and PDB formats
5. **Calculates bonds and angles** based on distance criteria with periodic boundary awareness
6. **Assigns specialized MINFF atom types** based on their chemical environment and coordination
7. **Generates a molecular topology file** (ITP) for use in molecular dynamics simulations
8. **Writes the final structure** with all assigned properties to output files

#### Running the script:

Simply execute `python run_atomi.py` from the command line. The script requires a GRO file named "Kaolinite_GII_0.0487.gro" in the same directory.

#### Output files:

- `replicated_structure.gro` - The enlarged supercell in GROMACS format
- `replicated_structure.pdb` - The enlarged supercell in PDB format
- `molecular_topology.itp` - GROMACS topology file with bond and angle definitions
- `preem.gro` - Final structure with MINFF typing and calculated properties

### Other Examples

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
