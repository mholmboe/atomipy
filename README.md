# atomipy: The atom library in Python

A modular Python toolbox for handling and analyzing molecular structures, particularly for mineral slabs with periodic boundary conditions. This toolbox is a light version of the MATLAB [**atom**](https://github.com/mholmboe/atom) library and can mainly be used to generate molecular topology files for the [**MINFF**](https://github.com/mholmboe/minff) forcefield with a streamlined Python interface. For small bulk systems, a simple online server running [**www.atomipy.io**](https://www.atomipy.io) is now available. Test cases for hydrated montmorillonite using the general and tailored MINFF parameters (angle force constant 500 kJ/mol/rad²) can be found in the [**example cases of the atom Toolbox**](https://github.com/mholmboe/atom/tree/master/ATOM_scripts_lecture/MINFF).

The package now supports generating GROMACS n2t (atom name to type) files for both MINFF and CLAYFF forcefields, enabling seamless integration with GROMACS utilities like gmx x2top for enhanced topology handling.

## Overview

This toolbox is designed to import, export, and analyze molecular structures with a focus on mineral slabs containing the elements Si, Al, Fe, Mg, Ca, Ti, Li, F, O, H. It handles periodic and triclinic simulation cells, and provides functions for calculating bonds, angles, and distances while taking periodic boundary conditions into account, and hence is ideal for generating molecular topology files for mineral bulk/slab systems that can be modelled using the [**MINFF**](https://github.com/mholmboe/minff) forcefield. However it also has the capability to handle clay minerals, hydroxides, and oxyhydroxides using CLAYFF (Cygan, R.T.; Liang, J.J.; Kalinichev, A.G. Molecular Models of Hydroxide, Oxyhydroxide, and Clay Phases and the Development of a General Force Field. *J. Phys. Chem. B* **2004**, *108*, 1255-1266).

The molecular structure information is stored in dictionaries where each atom has fields for coordinates, neighbors, bonds, angles, element type, and more.

## Key Features

- Support for multiple forcefields: MINFF and CLAYFF atom typing and parameter assignment
- Import/export PDB, Gromacs GRO, and XYZ files
- Generation of GROMACS n2t (atom name to type) files for use with gmx x2top
  - Neighbour distances honour periodic boundary conditions when the box is provided
  - Nearly identical environments are merged to avoid duplicate entries caused by floating-point noise
  - All file formats need to contain system size information
  - For .pdb files, system dimensions should be in the CRYST1 record
  - For .gro files, box dimensions should be in the last line
  - For .xyz files, system dimensions should be on the second line after a # character
- Generating topology files for MINFF and CLAYFF forcefields, for Gromacs (.itp), NAMD (.psf) and LAMMPS (.data)
- Handle both orthogonal and triclinic simulation cells with periodic boundary conditions
- Calculate bond distances and angles
- Element type assignment
- Coordination number analysis
- Distance matrices with PBC corrections (using both full matrix and efficient cell-list algorithms)
- Progress tracking for computationally intensive calculations
- Consolidated charge module with support for formal, MINFF, and CLAYFF charge assignments
- Unified coordinate transformation system:
  - Cartesian-to-fractional and fractional-to-cartesian conversions
  - Orthogonal-to-triclinic and triclinic-to-orthogonal transformations
  - Direct transformation methods for crystallographic calculations
  - Coordinate wrapping for periodic boundary conditions
- Supercell replication with proper handling of triclinic cells
- Support for case-insensitive element matching in charge assignments

## Requirements

- NumPy (>=1.18.0)
- tqdm (>=4.45.0) - for progress bars
- Numba (>=0.50.0, optional) - for performance optimization via JIT compilation

## Installation

If you're new to Python, follow these simple steps to get started with atomipy:

### Step 1: Install Python

1. Download and install Python from [python.org](https://www.python.org/downloads/) (version 3.7 or newer recommended)
2. During installation on Windows, make sure to check "Add Python to PATH"

### Step 2: Install atomipy

#### Method 1: Install from GitHub (recommended)

1. Open a terminal or command prompt
2. Install Git if you don't have it already: [git-scm.com](https://git-scm.com/downloads)
3. Clone the repository and install:
   ```bash
   git clone https://github.com/mholmboe/atomipy.git
   cd atomipy
   pip install -e .
   ```

#### Method 2: Manual Installation

1. Download this repository as a ZIP file (click the green "Code" button on GitHub and select "Download ZIP")
2. Extract the ZIP file to a folder on your computer
3. Open a terminal or command prompt
4. Navigate to the extracted folder:
   ```bash
   cd path/to/extracted/atomipy
   ```
5. Install the package and its dependencies:
   ```bash
   pip install -e .
   ```

### Step 3: Verify Installation

To verify that atomipy is installed correctly, run this simple test:

```python
# Create a file named test_atomipy.py with these contents:
import atomipy as ap
print("atomipy installed successfully!")
```

Then run it with:
```bash
python test_atomipy.py
```

## Getting Started for Python Beginners

Here's a simple example to help you get started with atomipy:

```python
# Create a file named my_first_atomipy.py with these contents:
import atomipy as ap

# Step 1: Load a structure file
print("Loading a GRO file...")
atoms, box_dim = ap.import_gro("example.gro")  # Replace with your GRO file
print(f"Loaded {len(atoms)} atoms")

# Step 2: Assign elements based on atom names
print("Assigning elements to atoms...")
for i in range(len(atoms)):
    atoms[i] = ap.element(atoms[i])

# Step 3: Calculate bonds and angles
print("Calculating bonds and angles...")
atoms = ap.bond_angle(atoms, box=box_dim)

# Step 4: Save as a new file
print("Saving processed structure...")
ap.write_gro(atoms, box=box_dim, "processed.gro")
print("Done!")
```

Run this script with:
```bash
python my_first_atomipy.py
```

### Gromacs .n2t File Generation Example

Here's an example of generating a GROMACS n2t file for use with gmx x2top:

```python
import atomipy as ap

# Load structure
atoms, cell, box = ap.import_auto("structure.gro")

# Process with forcefield (optional)
atoms, _ = ap.minff(atoms, box)

# Generate n2t file (box argument now precedes the optional output path)
n2t_path = ap.write_n2t(atoms, box=box, n2t_file="minff_atomtypes.n2t")
print(f"N2T file saved to: {n2t_path}")
```

The `box` argument accepts orthogonal `[lx, ly, lz]` vectors, 1×6 cell parameter lists, or the 1×9 `Box_dim` layout and will be normalised internally.

Alternatively, you can use the included MINFF helper script:

```bash
python minff2n2t.py structure.gro --output minff_atomtypes.n2t
```

For the simplest possible conversion you can call the structure-only helper:

```bash
python struct2n2t.py structure.gro
```

This script auto-detects the input format, forwards the box dimensions, and writes `<structure>.n2t`.

### Common Issues for Beginners

- **ModuleNotFoundError**: Make sure you've installed all required packages.
- **File not found errors**: Check that your file paths are correct and that the files exist.
- **No module named 'atomipy'**: Make sure you've installed the package correctly.

## Function Documentation

### File I/O

- `import_pdb(file_path)`: Import a PDB file, returning a list of atom dictionaries and box parameters
- `import_gro(file_path)`: Import a Gromacs GRO file, including velocities if present
- `write_pdb(atoms, box, file_path)`: Write atoms to a PDB file
- `write_gro(atoms, box, file_path)`: Write atoms to a Gromacs GRO file, including velocities if present

### Force Field

- `minff(atoms, box, ffname='minff', rmaxlong=2.45, rmaxH=1.2)`: Assign MINFF forcefield specific atom types to each atom
- `clayff(atoms, box, ffname='clayff', rmaxlong=2.45, rmaxH=1.2)`: Assign CLAYFF forcefield specific atom types to each atom
- `write_n2t(atoms, box=None, n2t_file=None, verbose=True)`: Generate a GROMACS n2t (atom name to type) file based on structural analysis, honouring periodic boundary conditions when a 1×3 box, 1×6 cell, or 1×9 ``Box_dim`` array is supplied and merging nearly identical environments

### Molecular Topology

- `write_itp(atoms, box, file_path)`: Write a Gromacs topology file
- `write_psf(atoms, box, file_path)`: Write a NAMD topology file
- `write_data(atoms, box, file_path)`: Write a LAMMPS topology file

### Atom Properties

- `element(atom)`: Guess the chemical element for an atom based on its residue name
- `mass()`: Returns a dictionary of atomic masses for different elements
- `radius()`: Returns a dictionary of van der Waals radii for different elements

### Structure Analysis

- `bond_angle(atoms, box, rmaxH=1.2, rmaxM=2.45)`: Compute bonds and angles for a given atomic structure
- `dist_matrix(atoms, box)`: Calculate a full distance matrix between all atoms with PBC
- `cell_list_dist_matrix(atoms, box)`: Calculate a sparse distance matrix between all atoms with PBC


### Coordinate Transformations

- `transform.orthogonal_to_triclinic(ortho_coords, box, atoms=None)`: Convert coordinates from orthogonal to triclinic space
- `transform.triclinic_to_orthogonal(atoms=None, coords=None, box_dim=None, box=None)`: Convert coordinates from triclinic to orthogonal space
- `transform.cartesian_to_fractional(atoms=None, cart_coords=None, box_dim=None, box=None)`: Convert coordinates to fractional coordinates
- `transform.fractional_to_cartesian(atoms=None, frac_coords=None, box_dim=None, box=None)`: Convert fractional coordinates to cartesian
- `transform.direct_cartesian_to_fractional(atoms=None, cart_coords=None, box_dim=None, cell=None)`: Direct conversion using crystallographic matrices
- `transform.direct_fractional_to_cartesian(atoms=None, frac_coords=None, box_dim=None, cell=None)`: Direct conversion using crystallographic matrices
- `transform.wrap_coordinates(atoms=None, coords=None, frac_coords=None, box_dim=None, box=None)`: Wrap coordinates to ensure they're within the primary unit cell
- `transform.get_cell_vectors(box)`: Calculate cell vectors from box parameters
- `transform.get_orthogonal_box(box_dim=None, box=None)`: Get orthogonal box dimensions from triclinic parameters
- `replicate.replicate_system(atoms, box, replicate=[1, 1, 1])`: Create supercells by replicating in a, b, c directions

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

The simulation cell can be represented in three ways with the unified `box` parameter:
- **Orthogonal box**: A 1x3 array `[lx, ly, lz]` for simple rectangular boxes
- **Cell parameters**: A 1x6 array `[a, b, c, alpha, beta, gamma]` (used in PDB files)
- **Triclinic box**: A 1x9 array `[lx, ly, lz, 0, 0, xy, 0, xz, yz]` (used in Gromacs GRO files for triclinic cells)

Conversion utilities `Box_dim2Cell()` and `Cell2Box_dim()` can be used to convert between these formats.

## Common Workflow Examples

Here are some common workflows that demonstrate how to use Atomipy for specific tasks:

## Box Parameter Handling

Atomipy now uses a standardized approach for handling simulation box parameters across all functions:

- All functions accept a generalized `box` parameter (replacing the previous `Box_dim` parameter in most functions)
- The `box` parameter supports three formats:
  1. **Orthogonal box**: A 1x3 array `[lx, ly, lz]` for simple rectangular boxes
  2. **Cell parameters**: A 1x6 array `[a, b, c, alpha, beta, gamma]` for crystallographic notation
  3. **Triclinic box**: A 1x9 array `[lx, ly, lz, 0, 0, xy, 0, xz, yz]` using GROMACS triclinic box format

- Conversion utilities are automatically applied internally based on the format provided
- Backward compatibility is maintained for functions still using `Box_dim` and `Cell` parameters

Example usage with the different formats:

```python
# Using orthogonal box format
box_ortho = [30.0, 30.0, 30.0]  # lx, ly, lz in Angstroms
atoms = ap.bond_angle(atoms, box=box_ortho)

# Using cell parameters format
box_cell = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]  # a, b, c, alpha, beta, gamma
atoms = ap.bond_angle(atoms, box=box_cell)

# Using triclinic box format
box_triclinic = [30.0, 30.0, 30.0, 0.0, 0.0, 5.0, 0.0, 5.0, 2.0]  # GROMACS format
atoms = ap.bond_angle(atoms, box=box_triclinic)
```

### Basic Structure Processing

```python
import atomipy as ap

# Load a structure
atoms, box_dim = ap.import_gro("my_structure.gro")

# Assign elements
for atom in atoms:
    atom = ap.element(atom)

# Calculate bonds and angles
atoms, bonds, angles = ap.bond_angle(atoms, box_dim)

# Save processed structure
ap.write_gro(atoms, box_dim, "processed.gro")
```

### Creating a Topology File for Molecular Dynamics

```python
import atomipy as ap

# Load structure
atoms, box_dim = ap.import_gro("my_mineral.gro")

# Process the structure (elements, bonds, etc.)
for atom in atoms:
    atom = ap.element(atom)

# Assign forcefield atom types (choose either MINFF or CLAYFF)
# For MINFF:
ap.minff(atoms, box_dim)
# OR for CLAYFF:
# ap.clayff(atoms, box_dim)

# Write a topology file for different simulation programs
# For GROMACS:
ap.write_itp(atoms, box_dim, "topology.itp")

# For LAMMPS:
ap.write_lmp(atoms, box_dim, "topology.lmp") 

# For NAMD:
ap.write_psf(atoms, box_dim, "topology.psf")
```

### Creating a Supercell

```python
import atomipy as ap

# Load structure
atoms, box_dim = ap.import_gro("unit_cell.gro")

# Process atoms
for atom in atoms:
    atom = ap.element(atom)

# Create a 2x2x2 supercell
replicated_atoms, new_box_dim = ap.replicate.replicate_cell(atoms, box_dim, replicate=[2, 2, 2])

# Save the supercell
ap.write_gro(replicated_atoms, new_box_dim, "supercell.gro")
```

## Example Scripts

### Example Scripts

The package includes example scripts that demonstrate comprehensive workflows for processing mineral structures using atomipy:

#### run_minff_atomi.py and run_clayff_atomi.py

These scripts demonstrate workflows for using MINFF and CLAYFF forcefields respectively. Both serve as excellent starting points for users new to the package.

#### What the script does:

1. **Imports a GRO structure file** (Kaolinite_GII_0.0487.gro) containing mineral coordinates and box dimensions
2. **Assigns chemical elements** to each atom using chemical knowledge-based rules
3. **Creates a supercell** by replicating the unit cell to reach a target size (40 Å in each dimension)
4. **Saves the replicated structure** in both GRO and PDB formats
5. **Calculates bonds and angles** based on distance criteria with periodic boundary awareness
6. **Assigns specialized MINFF atom types** based on their chemical environment and coordination
7. **Generates a molecular topology file** (ITP) for use in molecular dynamics simulations
8. **Writes the final structure** with all assigned properties to output files

#### Running the scripts:

Simply execute `python run_minff_atomi.py` or `python run_clayff_atomi.py` from the command line. The scripts require a GRO file named "Kaolinite_GII_0.0487.gro" in the same directory.

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
atoms, bonds, angles = ap.bond_angle(atoms, Box_dim)

# Assign charges using one of the available methods
# For formal charges (ions and water):
atoms = ap.assign_formal_charges(atoms)
# For MINFF charges:
# atoms = ap.charge_minff(atoms, Box_dim)
# For CLAYFF charges:
# atoms = ap.charge_clayff(atoms, Box_dim)

# Convert to fractional coordinates
frac_coords, atoms = ap.fract.cartesian_to_fractional(atoms, Box_dim)

# Create a 2x2x1 supercell
replicated_atoms, new_Box_dim = ap.replicate.replicate_cell(atoms, Box_dim, replicate=[2, 2, 1])

# Convert from triclinic to orthogonal coordinates if needed
ortho_atoms = ap.ortho.triclinic_to_orthogonal(replicated_atoms, new_Box_dim)


# Export to GRO format
ap.write_gro(replicated_atoms, new_Box_dim, "structure.gro")
```


## Differences from atom MATLAB library

This Python implementation is designed to provide similar functionality to the MATLAB atom library while following Python's conventions and making use of NumPy for efficient numerical operations. The data structure is dictionary-based rather than struct-based, and the function interfaces are designed for Python's style.

## License

This project is released under the MIT License.
