# Release Notes: atomipy v0.98 (from v0.97)

We are thrilled to release **atomipy v0.98**! This release introduces robust topology merging, advanced composition analysis, and a suite of interactive visual MD improvements in the builder interface.

---

## 🧬 Core Package & Engine Improvements (`atomipy`)

### 1. 🔀 Multi-Component Topology Merging & System Composition
* **Advanced Topology Merging (`atomipy.merge_top`)**: Added `merge_top.py` for merging multiple complex topologies and itp file linkages cleanly.
* **System Composition Analysis (`atomipy.composition`)**: Added `composition.py` to analyze molecular composition, count residues/components, and calculate system masses.
* **OpenMM Force Unit Compatibility**: Refactored OpenMM force calculation to utilize native `unit.md_unit_system` compatibility, avoiding platform-specific unit dictionary lookup crashes during energy calculations.

### 2. 🌀 Molecular Dynamics & Energy Minimization Refinements
* **Step-by-Step Energy Minimization Trajectories**: Implemented iterative L-BFGS energy minimization chunks allowing real-time trajectory recording and energy relaxation tracking.
* **Max Force Norm Auditing**: Enhanced energy minimization stdout logging to dynamically compute and print maximum Euclidean force norm vectors in standard GROMACS/MD units (`kJ/mol/nm`).

---

## 🌐 Web Visual Builder Improvements (`atomipy-web-module`)

### 1. 🎛️ Unified Configurable Reporting & Visualizing
* **Decoupled Log and PDB Frequencies**: Re-engineered the UI to allow configuring **Log frequency** and **PDB trajectory frequency** independently for both Energy Minimization and MD simulations.
* **Fluctuating Bounding Box Rendering**: Enabled dynamic unit cell/periodic box sizing in the 3Dmol trajectory viewer (supporting NPT simulation box expansion and contraction).
* **Clean Log Header Interceptor**: Injected an output stream interceptor wrapper (`CleanHeaderStream`) that cleanly filters quotation marks (`"`) and comment hashes (`#`) from CSV output logs.
* **Top-Left Warnings Position**: Repositioned the Workflow Warnings panel to the top-left corner of the canvas to optimize editor space and avoid overlapping bottom docks.

---

# Release Notes: atomipy v0.97 (from v0.96)

We are excited to announce **atomipy v0.97**! This major release brings **mixed organic-mineral simulations** to the platform. By introducing advanced bridges to ParmEd and the OpenFF ecosystem, users can now parametrize organic molecules via SMILES and seamlessly merge them into complex mineral models.

---

## 🧬 Mixed Organic/Mineral Forcefields & Topologies

### 1. ParmEd Bridge (`atomipy.parmed_bridge`)
* **Dynamic Organic Parametrization**: Parametrize any organic molecule from a SMILES string directly using GAFF-2.11 or CGenFF forcefields without leaving Python.
* **Mixed System Topology Engine**: Safely merge complex CLAYFF/MINFF mineral structures with parametrized organic molecules. The engine handles combination rule alignment and intelligently merges structures into a unified configuration.
* **Universal Export formats**: Export the resulting mixed topology to AMBER (`.prmtop`), NAMD (`.psf`), GROMACS (`.itp`), and LAMMPS (`.data`).

### 2. Standalone OpenFF Worker Microservice
* **Dockerized SMIRNOFF Engine**: Introduced an isolated FastAPI microservice for OpenFF Sage, Parsley, and Rosemary parametrizations to avoid dependency conflicts with ParmEd.

---

## 🌐 Web Visual Builder Improvements (`atomipy-web-module`)

### 1. Organic Molecular Generation Node
* **SMILES Node Integration**: A brand new interactive `OrganicNode` added to the Visual Builder. Users can inject organic components into their mineral workflows, preview conformations on-the-fly, and select target forcefields visually.
* **Advanced Code Generation (`graphExecution.ts`)**: The visual node compiler was re-engineered to recognize organic and mixed topologies and automatically configure downstream export tools for unified formats like AMBER prmtop.

---

# Release Notes: atomipy v0.96 (from v0.95)

## 🛠️ Core Package & Engine Improvements (`atomipy`)

### 1. 🧪 Integrated Molecular Dynamics Simulations with OpenMM
* **GROMACS-to-OpenMM Bridge**: Added the `atomipy.openmm_interface` module containing the `load_minff_into_openmm` parser. It dynamically links GROMACS `.top`/`.itp` files directly into OpenMM topology and system objects, allowing seamless MD execution from Python.
* **Auto-Platform Selector**: Implemented an automated hardware platform detection selector. Simulations will run on the fastest available platform (prioritizing **CUDA** or **OpenCL** GPU platforms for maximum performance, then seamlessly falling back to CPU or Reference platforms).
* **Dynamic Parameter Injection**: Automatically handles forcefield constraints (e.g. flexible water models OPC3/SPC/TIP3P, custom mineral bond/angle parameters, and ion corrections) dynamically at load time.
* **Unicode/Encoding Safety**: Updated all trajectory and PDB reading/writing interfaces to use strict UTF-8 mode, resolving encoding crashes on systems default-locale-configured as ASCII when writing unicode ion labels (`Cl−`).

### 2. 🛡️ Robust Edge Case Handling & Empty Box Safety
* **Dummy Atom Instantiation**: Introduced a safe-fallback protocol for empty box configurations. If a box has zero atoms, `atomipy` instantiates a single `Dummy` atom at `(0,0,0)` with zero charge and mass to guarantee compatibility with visualization tools and molecular physics engines.

---

## 🌐 Web Visual Builder Improvements (`atomipy-web-module`)

### 1. 🎛️ Node-based Interactive MD Simulations
* **Simulation Node**: Added a draggable visual `SimulationNode` allowing researchers to set time steps, temperatures, and run NVT molecular dynamics simulations directly from the workflow graph.
* **Live Progress Plots**: Integrated real-time potential energy and temperature trajectory charts that update dynamically as the simulation runs in the background.
* **Auto-Unzipping Download Bundle**: Refactored the bundle export feature to package completed trajectories, GROMACS topologies, and system coordinates into a single download that auto-unzips locally.

---

*Thank you to all our users and contributors! For the full source code and visual builder guides, visit the [atomipy GitHub Repository](https://github.com/mholmboe/atomipy).*
