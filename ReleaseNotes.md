# Release Notes: atomipy v0.96 (from v0.95)

We are thrilled to announce the release of **atomipy v0.96**! This version introduces complete **Molecular Dynamics (MD) Simulation capabilities** integrated directly with **OpenMM**, an intelligent GPU/CPU acceleration auto-selector, and major interactive visual workflow refinements.

---

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
