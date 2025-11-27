# atomipy: The atom library in Python

A modular Python toolbox for handling and analyzing molecular structures, particularly for mineral slabs with periodic boundary conditions, PBC. This toolbox is a light version of the MATLAB [**atom**](https://github.com/mholmboe/atom) library and can in particular be used to generate molecular topology files for the [**MINFF**](https://github.com/mholmboe/minff) forcefield with a streamlined Python interface. For small bulk systems, a simple online server running [**www.atomipy.io**](https://www.atomipy.io) is now available. Test cases for hydrated montmorillonite using the general and tailored MINFF parameters (angle force constant 500 kJ/mol/rad²) can be found in the [**example cases of the atom Toolbox**](https://github.com/mholmboe/atom/tree/master/ATOM_scripts_lecture/MINFF).

The package now supports generating GROMACS n2t (atom name to type) files for both MINFF and CLAYFF forcefields, enabling seamless integration with GROMACS utilities like gmx x2top for enhanced topology handling.

## Contents
- [Overview](#overview)
- [Common Variables](#common-variables)
  - [Structure Containers](#structure-containers)
  - [Atom dictionary fields](#atom-dictionary-fields)
- [Key Features](#key-features)
- [Requirements](#requirements)
- [Installation](#installation)
- [MINFF atom typing & topology generation](#minff-atom-typing--topology-generation)
- [XRD pattern simulation](#xrd-pattern-simulation)
- [Bond valence analysis](#bond-valence-analysis)
- [Ionic/crystal radii and bond estimates (Shannon radii)](#ioniccrystal-radii-and-bond-estimates-shannon-radii)

## Overview

This toolbox is designed to import, export, and analyze molecular structures with a focus on mineral slabs containing the elements Si, Al, Fe, Mg, Ca, Ti, Li, F, O, H. It handles periodic and triclinic simulation cells, and provides functions for calculating bonds, angles, and distances while taking periodic boundary conditions into account, and hence is ideal for generating molecular topology files for mineral bulk/slab systems that can be modelled using the [**MINFF**](https://github.com/mholmboe/minff) forcefield. However it also has the capability to handle clay minerals, hydroxides, and oxyhydroxides using CLAYFF (Cygan, R.T.; Liang, J.J.; Kalinichev, A.G. Molecular Models of Hydroxide, Oxyhydroxide, and Clay Phases and the Development of a General Force Field. *J. Phys. Chem. B* **2004**, *108*, 1255-1266).

Built-in atom typing for MINFF and CLAYFF is a core feature: you can assign atom types, charges, and generate ITP/PSF/LAMMPS topologies directly from structures. For quick runs without local setup, use the web server at [www.atomipy.io](https://www.atomipy.io).

The molecular structure information is stored in dictionaries where each atom has fields for coordinates, neighbors, bonds, angles, element type, and more.

## Common Variables

Understanding the core containers and fields used in atomipy makes it easier to follow the examples and extend the code. The tables below summarise the names, expected types, and roles of the variables you will encounter most often.

### Structure Containers

| Variable | Type | Description |
| --- | --- | --- |
| `atoms` | `list[dict]` | Primary structure container; each dict holds coordinates, typing, charge, and connectivity for one atom. |
| `Box` | `list[float]` (len 3/6/9) | Canonical Cell dimensions accepted by most APIs; supports `[lx, ly, lz]`, `[a, b, c, alpha, beta, gamma]`, or the GROMACS-style `[lx, ly, lz, 0, 0, xy, 0, xz, yz]` in Angstrom and degrees. |
| `Box_dim` | `list[float]` (len 3/9) | Box dimensions returned by GRO import (1x3 orthogonal or 1x9 triclinic) and used for GRO/ITP/LAMMPS writers. |
| `Cell` | `list[float]` (len 6) | Unit-cell parameters `[a, b, c, alpha, beta, gamma]` in Angstrom and degrees. |
| `Bond_index` | `numpy.ndarray` or `list` | Bond records `[i, j, distance]` returned by `bond_angle`; indices are 0- or 1-based depending on export step. |
| `Angle_index` | `numpy.ndarray` or `list` | Angle records `[i, j, k, angle, dx12, dy12, dz12, dx23, dy23, dz23]` returned by `bond_angle`. |
| `distances` | `numpy.ndarray` (N, N) | Full pairwise distance matrix produced by `dist_matrix`. |
| `dx`, `dy`, `dz` | `numpy.ndarray` (N, N) | Cartesian component offsets that accompany `distances`. |
| `cart_coords` | `numpy.ndarray` (N, 3) | Cartesian coordinates created by transformation helpers such as `fractional_to_cartesian`. |
| `frac_coords` | `numpy.ndarray` (N, 3) | Fractional coordinates created by `cartesian_to_fractional` and wrapping utilities. |
| `total_charge` | `float` | Net system charge computed during forcefield or charge analysis. |
| `atom_types` | `list[str]` | Unique atom or forcefield type labels gathered from the atom dictionaries. |

### Atom dictionary fields

| Key | Type | Description |
| --- | --- | --- |
| `index` | `int` | 1-based atom index assigned during import or when structures are combined. |
| `molid` | `int` | Molecule or residue identifier used to group atoms. |
| `resname` | `str` | Residue or component name (e.g. `MMT`, `ION`, `SOL`). |
| `resid` | `int` | Optional residue sequence number preserved from input files. |
| `type` | `str` | Working atom name used for matching templates and neighbour searches. |
| `fftype` | `str` | Forcefield-specific atom type assigned by `minff`, `clayff`, or related routines. |
| `element` | `str` | Chemical element symbol inferred by `element`. |
| `x`, `y`, `z` | `float` | Cartesian coordinates in Angstrom. |
| `vx`, `vy`, `vz` | `float` or `None` | Optional velocities in Angstrom per picosecond carried over from GRO files. |
| `xfrac`, `yfrac`, `zfrac` | `float` | Fractional coordinates stored when using the transformation helpers. |
| `neigh` | `list[int]` | Indices of bonded neighbours (0-based) populated by `bond_angle`. |
| `bonds` | `list[tuple[int, float]]` | `(neighbour_index, distance_angstrom)` pairs generated by `bond_angle`. |
| `angles` | `list[tuple[tuple[int, int], float]]` | `((neighbour1, neighbour2), angle_deg)` entries describing local angles. |
| `cn` | `int` | Coordination number saved when `calculate_coordination=True`. |
| `charge` | `float` | Partial charge in elementary charge units; initial imports may store strings before reassignment. |
| `mass` | `float` | Atomic mass in atomic mass units assigned via `set_atomic_masses()` or forcefield setup. |
| `com` | `dict` | Center of mass coordinates `{'x': x, 'y': y, 'z': z}` added by the `com()` function. |
| `occupancy` | `float` | PDB occupancy value (defaults to 1.0 when unavailable). |
| `temp_factor` | `float` | PDB temperature factor or B-factor (defaults to 0.0). |
| `chain_id` | `str` | PDB chain identifier used during export. |
| `icode` | `str` | PDB insertion code, if present. |
| `is_nm` | `bool` | Flag noting whether coordinates originated in nanometres before conversion to Angstrom. |

## Key Features

- Support for multiple forcefields: MINFF and CLAYFF atom typing and parameter assignment
- Import/export PDB, Gromacs GRO, and XYZ files
- Generation of GROMACS n2t (atom name to type) files for use with gmx x2top
  - Neighbour distances honour periodic boundary conditions when the Box is provided
  - Nearly identical environments are merged to avoid duplicate entries caused by floating-point noise
  - All file formats need to contain system size information
  - For .pdb files, system dimensions should be in the CRYST1 record
  - For .gro files, Box dimensions should be in the last line
  - For .xyz files, system dimensions should be on the second line after a # character
- Generating topology files for MINFF and CLAYFF forcefields, for Gromacs (.itp), NAMD (.psf) and LAMMPS (.data)
- High-performance X-ray diffraction pattern simulation and plotting via `diffraction.xrd`, with optional CLI helper `run_xrd_example.py`
- Handle both orthogonal and triclinic simulation cells with periodic boundary conditions
- Calculate bond distances, angles, dihedrals, and 1–4 pairs (`bond_angle_dihedral`)
- Bond valence analysis utilities (bond valence sums and Global Instability Index)
- Ionic/crystal radii utilities from Revised Shannon radii, including bond-distance estimates
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
- Isomorphous substitution for creating defects in mineral structures:
  - Replace atom types in octahedral and/or tetrahedral sites (e.g., Al→Mg, Si→Al)
  - Ensure minimum separation distances between substituted atoms
  - Distribute substitutions evenly across the structure
  - Optional spatial limits for targeted substitution regions

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
atoms, Box_dim = ap.import_gro("example.gro")  # Replace with your GRO file
print(f"Loaded {len(atoms)} atoms")

# Step 2: Assign elements based on atom names
print("Assigning elements to atoms...")
atoms = ap.element(atoms)

# Step 3: Calculate bonds and angles
print("Calculating bonds and angles...")
atoms = ap.bond_angle(atoms, Box=Box_dim)

# Step 4: Save as a new file
print("Saving processed structure...")
ap.write_gro(atoms, Box=Box_dim, file_path="processed.gro")
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

# Load structure (auto-detect format)
auto_result = ap.import_auto("structure.gro")
if len(auto_result) == 2:
    atoms, box_like = auto_result
else:
    atoms, box_like = auto_result[0], auto_result[-1]

# Normalise Box (Cell -> Box_dim if needed)
if len(box_like) == 6:
    Box = ap.Cell2Box_dim(box_like)
else:
    Box = box_like

# Process with forcefield (optional)
atoms, _ = ap.minff(atoms, Box)

# Generate n2t file (Box argument now precedes the optional output path)
n2t_path = ap.write_n2t(atoms, Box=Box, n2t_file="minff_atomtypes.n2t")
print(f"N2T file saved to: {n2t_path}")
```

The `Box` argument accepts orthogonal `[lx, ly, lz]` vectors, 1×6 Cell parameter lists, or the 1×9 `Box_dim` layout and will be normalised internally.

### Dihedrals and 1–4 pair detection

Use `bond_angle_dihedral` to build bonds/angles and derive dihedrals plus 1–4 pairs in one call:

```python
atoms, Box_dim = ap.import_gro("NMA_element.gro")
ap.element(atoms)
atoms, Bond_index, Angle_index, Dihedral_index, Pairlist = ap.bond_angle_dihedral(
    atoms, Box_dim, same_molecule_only=False, same_element_bonds=False
)
print(len(Bond_index), len(Angle_index), len(Dihedral_index), len(Pairlist))
```

### Import existing .itp files

To read a GROMACS topology and inspect interaction sections:

```python
from atomipy import import_itp_topology
itp = import_itp_topology("NMA.itp")
print("Pairs:", len(itp.get("pairs", {}).get("ai", [])))
print("Dihedrals:", len(itp.get("dihedrals", {}).get("ai", [])))
```

Alternatively, you can use the included MINFF helper script:

```bash
python minff2n2t.py structure.gro --output minff_atomtypes.n2t
```

For the simplest possible conversion you can call the structure-only helper:

```bash
python struct2n2t.py structure.gro
```

This script auto-detects the input format, forwards the Box dimensions, and writes `<structure>.n2t`.

### MINFF atom typing & topology generation

Assign MINFF atom types/charges and write topology files directly from a structure:

```python
import atomipy as ap

atoms, cell = ap.import_pdb("Kaolinite_GII_0.0487.pdb")
typed = ap.minff(atoms, Box=cell, log=True)

# Write out topologies
ap.write_itp(typed, Box=cell, file_path="minff_Kao.itp")
ap.write_pdb(typed, Box=cell, file_path="minff_Kao.pdb")
ap.write_lmp(typed, Box=cell, file_path="minff_Kao.data")
```

If your input is GRO, pass the returned `Box_dim` to `minff`/writers:

```python
atoms, Box_dim = ap.import_gro("structure.gro")
typed = ap.minff(atoms, Box=Box_dim)
ap.write_itp(typed, Box=Box_dim, file_path="structure.itp")
```

Convenient helpers:
- `run_create_itp_example.py` shows a minimal end-to-end MINFF + ITP + typed PDB workflow.
- `minff2n2t.py` generates a MINFF-typed `.n2t` mapping for gmx x2top.

### XRD pattern simulation

atomipy includes a fast X-ray diffraction module (`atomipy.diffraction.xrd`) that can turn a PDB/GRO/XYZ structure plus its box/cell into a calculated powder pattern with optional plotting and data export.

Command-line helper:
```bash
python run_xrd_example.py Kaolinite_GII_0.0487.pdb --two-theta 5 70 --save-output
```

Minimal API example:
```python
from atomipy import import_auto, xrd

atoms, box = import_auto("Kaolinite_GII_0.0487.pdb")
two_theta, intensity, fig = xrd(
    atoms=atoms,
    Box=box,
    wavelength=1.54187,      # Cu K-alpha
    two_theta_range=(5, 70), # degrees
    angle_step=0.02,
    save_output=False,
)
fig.show()
```

### Common Issues for Beginners

- **ModuleNotFoundError**: Make sure you've installed all required packages.
- **File not found errors**: Check that your file paths are correct and that the files exist.
- **No module named 'atomipy'**: Make sure you've installed the package correctly.

## Function Documentation

### File I/O

- `import_pdb(file_path)`: Import a PDB file, returning `(atoms, Cell)` - a list of atom dictionaries and Cell parameters
- `import_gro(file_path)`: Import a Gromacs GRO file, returning `(atoms, Box_dim)` (1x3 orthogonal or 1x9 triclinic), including velocities if present
- `import_xyz(file_path)`: Import an XYZ file, returning `(atoms, Cell)`
- `import_auto(file_path)`: Auto-detect format and import, returning `(atoms, Cell)` for PDB/XYZ or `(atoms, Box_dim)` for GRO
- `write_pdb(atoms, Box, file_path)`: Write atoms to a PDB file
- `write_gro(atoms, Box, file_path)`: Write atoms to a Gromacs GRO file, including velocities if present
- `write_xyz(atoms, Box, file_path)`: Write atoms to an XYZ file

### Force Field

- `minff(atoms, Box, ffname='minff', rmaxlong=2.45, rmaxH=1.2)`: Assign MINFF forcefield specific atom types to each atom
- `clayff(atoms, Box, ffname='clayff', rmaxlong=2.45, rmaxH=1.2)`: Assign CLAYFF forcefield specific atom types to each atom
- `write_n2t(atoms, Box=None, n2t_file=None, verbose=True)`: Generate a GROMACS n2t (atom name to type) file based on structural analysis, honouring periodic boundary conditions when a 1×3 Box, 1×6 Cell, or 1×9 ``Box_dim`` array is supplied and merging nearly identical environments

### Molecular Topology

- `write_itp(atoms, Box, file_path)`: Write a Gromacs topology file
- `write_psf(atoms, Box, file_path)`: Write a NAMD topology file
- `write_data(atoms, Box, file_path)`: Write a LAMMPS topology file

### Atom Properties

- `element(atoms)`: Assign chemical elements to atoms based on their atom types and residue names
- `mass()`: Returns a comprehensive dictionary of atomic masses for different elements (includes all elements from the periodic table)
- `radius()`: Returns a dictionary of van der Waals radii for different elements
- `set_atomic_masses(atoms)`: Set mass attributes for each atom in the list based on element type
- `com(atoms, add_to_atoms=True)`: Calculate the center of mass of a molecule or slab without PBC wrapping, optionally adding COM coordinates to each atom dictionary

### Structure Analysis

- `bond_angle(atoms, Box, rmaxH=1.2, rmaxM=2.45)`: Compute bonds and angles for a given atomic structure
- `dist_matrix(atoms, Box)`: Calculate a full distance matrix between all atoms with PBC
- `cell_list_dist_matrix(atoms, Box)`: Calculate a sparse distance matrix between all atoms with PBC


### Coordinate Transformations

- `orthogonal_to_triclinic(ortho_coords, Box, atoms=None)`: Convert coordinates from orthogonal to triclinic space
- `triclinic_to_orthogonal(atoms=None, coords=None, Box)`: Convert coordinates from triclinic to orthogonal space
- `cartesian_to_fractional(atoms=None, cart_coords=None, Box)`: Convert Cartesian coordinates to fractional coordinates
- `fractional_to_cartesian(atoms=None, frac_coords=None, Box)`: Convert fractional coordinates to Cartesian
- `direct_cartesian_to_fractional(atoms=None, cart_coords=None, Box)`: Direct conversion using crystallographic matrices
- `direct_fractional_to_cartesian(atoms=None, frac_coords=None, Box)`: Direct conversion using crystallographic matrices
- `wrap(atoms, Box, return_type='cartesian')`: Wrap atom coordinates into the primary simulation cell. Simple interface supporting orthogonal (1x3), triclinic (1x9), and Cell parameter (1x6) formats. Returns atoms with updated coordinates and fractional coordinates
- `wrap_coordinates(atoms=None, coords=None, frac_coords=None, Box=None, add_to_atoms=True, return_type='fractional')`: Advanced wrapping function with more control over input/output formats
- `get_cell_vectors(Box)`: Calculate Cell vectors from Box parameters

### Diffraction

- `xrd(atoms, Box, wavelength=1.54187, angle_step=0.02, two_theta_range=(2, 90), ...)`: Calculate and optionally plot/save an XRD powder pattern using the high-performance diffraction module (see `run_xrd_example.py` for CLI usage)
- `get_orthogonal_box(Box)`: Get orthogonal Box dimensions from triclinic parameters
- `replicate_system(atoms, Box, replicate=[1, 1, 1])`: Create supercells by replicating in a, b, c directions

### Structure Building

- `substitute(atoms, box, num_oct_subst, o1_type, o2_type, min_o2o2_dist, ...)`: Perform isomorphous substitution by replacing atom types (e.g., Al→Mg in octahedral sites). Advanced features include:
  - Automatic detection of centrosymmetric structures for proper centering during substitution
  - Fallback to element-based matching when atom type names aren't found (using the element module)
  - Even distribution of substitutions across the structure
  - Support for both octahedral and tetrahedral substitutions with minimum distance constraints
  - Optional spatial limits for targeting substitutions to specific regions
- `is_centrosymmetric_along_z(atoms, tolerance=0.1)`: Check if a structure is approximately centrosymmetric along the z-axis by analyzing the distribution of z-coordinates

### Cell Utilities

- `cell_utils.Box_dim2Cell(Box_dim)`: Convert Box dimensions to Cell parameters
- `cell_utils.Cell2Box_dim(Cell)`: Convert Cell parameters to Box dimensions

### Charges

- `assign_formal_charges(atoms)`: Assign formal charges to ions and water molecules based on residue names (SOL, ION) and element types
- `charge_minff(atoms, Box, atom_labels=None, charges=None, resname=None)`: Assign MINFF charges to atoms based on coordination environment
- `charge_clayff(atoms, Box, atom_labels=None, charges=None, resname=None)`: Assign CLAYFF charges to atoms based on coordination environment
- `balance_charges(atoms, resname=None)`: Balance the total charge of the system by adjusting oxygen charges

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

The simulation Cell can be represented in three ways with the unified `Box` parameter:
- **Orthogonal box**: A 1x3 array `[lx, ly, lz]` for simple rectangular boxes
- **Cell parameters**: A 1x6 array `[a, b, c, alpha, beta, gamma]` (used in PDB files)
- **Triclinic box**: A 1x9 array `[lx, ly, lz, 0, 0, xy, 0, xz, yz]` (used in Gromacs GRO files for triclinic cells)

Conversion utilities `Box_dim2Cell()` and `Cell2Box_dim()` can be used to convert between these formats.

## Common Workflow Examples

Here are some common workflows that demonstrate how to use Atomipy for specific tasks:

## Box Parameter Handling

Atomipy now uses a standardized approach for handling simulation Box parameters across all functions:

- All functions accept a generalized `Box` parameter and normalise internally
- The `Box` parameter supports three formats:
  1. **Orthogonal box**: A 1x3 array `[lx, ly, lz]` for simple rectangular boxes
  2. **Cell parameters**: A 1x6 array `[a, b, c, alpha, beta, gamma]` for crystallographic notation
  3. **Triclinic box**: A 1x9 array `[lx, ly, lz, 0, 0, xy, 0, xz, yz]` using GROMACS triclinic Box format

- Conversion utilities are automatically applied internally based on the format provided

Example usage with the different formats:

```python
# Using orthogonal Box format
Box_ortho = [30.0, 30.0, 30.0]  # lx, ly, lz in Angstroms
atoms = ap.bond_angle(atoms, Box=Box_ortho)

# Using Cell parameters format
Box_cell = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]  # a, b, c, alpha, beta, gamma
atoms = ap.bond_angle(atoms, Box=Box_cell)

# Using triclinic Box format
Box_triclinic = [30.0, 30.0, 30.0, 0.0, 0.0, 5.0, 0.0, 5.0, 2.0]  # GROMACS format
atoms = ap.bond_angle(atoms, Box=Box_triclinic)
```

### Basic Structure Processing

```python
import atomipy as ap

# Load structure
atoms, Box_dim = ap.import_gro("my_structure.gro")

# Assign elements
atoms = ap.element(atoms)

# Calculate bonds and angles
atoms, bonds, angles = ap.bond_angle(atoms, Box=Box_dim)

# Save processed structure
ap.write_gro(atoms, Box=Box_dim, file_path="processed.gro")
```

### Calculating Center of Mass

```python
import atomipy as ap

# Load structure
atoms, Box_dim = ap.import_gro("my_molecule.gro")

# Ensure atoms have element information (needed for mass assignment)
atoms = ap.element(atoms)

# Calculate center of mass and add it to each atom
com_coords = ap.com(atoms, add_to_atoms=True)
print(f"Center of mass coordinates: x={com_coords[0]:.3f}, y={com_coords[1]:.3f}, z={com_coords[2]:.3f} Å")

# Access the COM stored in each atom
print(f"COM from first atom: {atoms[0]['com']}")

# Calculate COM without adding to atoms
com_coords_only = ap.com(atoms, add_to_atoms=False)
```

### Creating a Topology File for Molecular Dynamics

```python
import atomipy as ap

# Load structure
atoms, Box = ap.import_gro("my_mineral.gro")

# Process the structure (elements, bonds, etc.)
atoms = ap.element(atoms)

# Assign forcefield atom types (choose either MINFF or CLAYFF)
# For MINFF:
ap.minff(atoms, Box=Box)
# OR for CLAYFF:
# ap.clayff(atoms, Box=Box)

# Write a topology file for different simulation programs
# For GROMACS:
ap.write_itp(atoms, Box=Box, file_path="topology.itp")

# For LAMMPS:
ap.write_lmp(atoms, Box=Box, file_path="topology.lmp") 

# For NAMD:
ap.write_psf(atoms, Box=Box, file_path="topology.psf")
```

### Creating a Supercell

```python
import atomipy as ap

# Load structure
atoms, Box = ap.import_gro("unit_cell.gro")

# Process atoms
atoms = ap.element(atoms)

# Create a 2x2x2 supercell
replicated_atoms, new_Box, new_Cell = ap.replicate_system(atoms, Box, replicate=[2, 2, 2])

# Save the supercell
ap.write_gro(replicated_atoms, Box=new_Box, file_path="supercell.gro")
```

### Isomorphous Substitution in Clay Minerals

```python
import atomipy as ap

# Load a clay mineral structure (e.g., pyrophyllite or montmorillonite)
atoms, Cell = ap.import_pdb("Pyrophyllite_GII_0.071.pdb")

# The substitute function now automatically detects if a structure is centrosymmetric
# and centers it appropriately during substitution (important for pyrophyllite)

# Perform octahedral substitution: Replace 16 Al atoms with Mg atoms
# Using the mineral-specific atom types (will automatically look up by element if type not found)
atoms = ap.substitute(
    atoms, 
    Cell, 
    num_oct_subst=16, 
    o1_type='Alo',   # Uses 'Alo' atom type for octahedral Al
    o2_type='Mgo',   # Will become Mg in octahedral position
    min_o2o2_dist=5.2  # Minimum distance between Mg atoms in Angstroms
)

# Perform both octahedral and tetrahedral substitutions
# Replace 4 Al with Mg (octahedral) and 8 Si with Al (tetrahedral)
# If atom types aren't found, it will try to match by element instead
atoms = ap.substitute(
    atoms, 
    Cell, 
    num_oct_subst=4, 
    o1_type='Al',    # Will find atoms by element if 'Al' type not found
    o2_type='Mgo', 
    min_o2o2_dist=5.5,
    num_tet_subst=8,
    t1_type='Si',
    t2_type='Alt',   # Will become Al in tetrahedral position
    min_t2t2_dist=5.5
)

# Limit substitutions to specific region (e.g., only z > 10 Å)
# The function ensures even distribution between upper and lower parts
atoms = ap.substitute(
    atoms, 
    Box, 
    num_oct_subst=4, 
    o1_type='Al', 
    o2_type='Mgo', 
    min_o2o2_dist=5.5,
    lo_limit=10.0,
    hi_limit=50.0,
    dimension=2  # 0=x, 1=y, 2=z (Python-style indexing)
)

# Check if a structure is centrosymmetric along z (can be used independently)
is_symmetric = ap.is_centrosymmetric_along_z(atoms)
print(f"Is the structure centrosymmetric along z? {is_symmetric}")

# Save the substituted structure
ap.write_pdb(atoms, Cell=Box_dim, file_path="substituted_clay.pdb")
```

## Example Scripts

### Example Scripts

The package includes example scripts that demonstrate comprehensive workflows for processing mineral structures using atomipy:

#### run_minff_atomi.py and run_clayff_atomi.py

These scripts demonstrate workflows for using MINFF and CLAYFF forcefields respectively. Both serve as excellent starting points for users new to the package.

#### What the script does:

1. **Imports a GRO structure file** (Kaolinite_GII_0.0487.gro) containing mineral coordinates and Box dimensions
2. **Assigns chemical elements** to each atom using chemical knowledge-based rules
3. **Creates a supercell** by replicating the unit Cell to reach a target size (40 Å in each dimension)
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

### Bond valence analysis

Compute bond valence sums (BVS) and the Global Instability Index (GII) using the IUCr parameter table (`bvparm2020.cif`):

```python
import atomipy as ap

atoms, Box_dim = ap.import_gro("structure.gro")
params = ap.load_bv_params()  # optional, uses bvparm2020.cif by default
results, gii = ap.compute_bvs(atoms, Box_dim, params=params)
print("GII =", gii)
for atom in results[:5]:
    print(atom["index"], atom["element"], atom["bvs"], atom["expected_ox"], atom["delta"])
```

`results` contains per-atom bond valence sums, the inferred oxidation state, and the deviation (`delta`). Use `ap.summarize_bvs(results)` to get quick per-element averages and the worst-offending sites.

### Ionic/crystal radii and bond estimates (Shannon radii)

Use the Revised Shannon radii table to fetch ionic/crystal radii and estimate bond distances:

```python
import atomipy as ap

# Radii lookups
r_al = ap.get_radius("Al", 3, 6)                     # ionic radius
r_si_cryst = ap.get_radius("Si", 4, 4, prefer="crystal")

# Bond distance estimates (radius sums)
d_al_o = ap.bond_distance("Al", 3, 6, "O", -2, 4)
d_si_o = ap.bond_distance("Si", 4, 4, "O", -2, 4, use_crystal=True)
```

### Other Examples

```python
# Import the entire package
import atomipy as ap
import numpy as np

# Import a PDB file
atoms, Cell = ap.import_pdb("structure.pdb")

# Guess elements for all atoms
atoms = ap.element(atoms)

# Calculate bonds and angles
atoms, bonds, angles = ap.bond_angle(atoms, Box=Cell)

# Assign charges using one of the available methods
# For formal charges (ions and water):
atoms = ap.assign_formal_charges(atoms)
# For MINFF charges:
# atoms = ap.charge_minff(atoms, Box=Cell)
# For CLAYFF charges:
# atoms = ap.charge_clayff(atoms, Box=Cell)

# Convert to fractional coordinates
frac_coords, atoms = ap.cartesian_to_fractional(atoms, Box=Cell)

# Create a 2x2x1 supercell
replicated_atoms, new_Box, new_Cell = ap.replicate_system(atoms, Cell, replicate=[2, 2, 1])

# Convert from triclinic to orthogonal coordinates if needed
ortho_atoms, _, _ = ap.triclinic_to_orthogonal(replicated_atoms, Box=new_Box)


# Export to GRO format
ap.write_gro(replicated_atoms, Box=new_Box, file_path="structure.gro")
```


## Differences from atom MATLAB library

This Python implementation is designed to provide similar functionality to the MATLAB atom library while following Python's conventions and making use of NumPy for efficient numerical operations. The data structure is dictionary-based rather than struct-based, and the function interfaces are designed for Python's style.

## License

This project is released under the MIT License.
