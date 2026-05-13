# Release Notes: atomipy v0.95 (from v0.94)

We are thrilled to announce the release of **atomipy v0.95**! This version introduces substantial performance robustness, automated structural verification pipelines, and significant interactive user experience (UX) refinements in the companion Web Visual Builder dashboard.

---

## 🛠️ Core Package & Engine Improvements (`atomipy`)

### 1. Robust Periodic Boundary Conditions (PBC) & Triclinic Normalization
* **Triclinic Box Robustness**: Refactored `get_neighbor_list` and coordinate normalization routines in `distances.py` to parse length-9 arrays (GROMACS/matrix formats) and length-6 dimensions gracefully without throwing ValueErrors.
* **Coordinate Wrapping Fix**: Resolved an edge-case coordinate wrapping issue in PBC-aware sparse neighbor lists, ensuring high-accuracy atom connectivity maps even in narrow-cell (e.g. Portlandite) or highly oblique crystal systems.

### 2. Bond Valence Sum (BVS) & Chemistry Refinements
* **Dynamic Refinement**: Implemented a self-consistent oxidation state assignment loop in `bond_valence.py` to dynamically adjust transition metal valences (e.g., Fe²⁺ vs Fe³⁺ in Wüstite), driving the global instability index (GII) to its mathematical optimum.
* **Symmetry & CIF Parsing**: Fixed a deepcopy bug in `fuse_atoms` and refined the space-group symmetry expander to handle complex unit cell boundaries correctly.

### 3. Powder X-Ray Diffraction (XRD) Simulation
* **Element-Type Prioritization**: Standardized element mapping in scattering factor routines to prioritize actual element characters over forcefield type overrides (e.g., `Mgo` -> `Mg`), preventing coordinate lookup failures during simulations.

### 4. 🧪 Automated Testing Suite
* Fully migrated `atomipy` from script-based verification to a professional, high-performance **`pytest`** test suite. Running `pytest tests/` now provides 100% green passing coverage over:
  * Asymmetric CIF parsing and expansion.
  * PBC-aware coordination and skewed distance metrics.
  * BVS dynamic convergence.
  * Solvation, clay layers replication, and forcefield parameterization pipelines.

---

## 🌐 Web Visual Builder Improvements (`atomipy-web-module`)

### 1. Seamless 3Dmol vs JSmol Toggling
* Resolved an interactive layout lock where switching between renderers resulted in empty views. Toggling between WebGL (3Dmol) and HTML5 (JSmol) views now displays structures **instantly** without requiring system rebuilds, while keeping memory extremely lean by cleanly destroying inactive instances.

### 2. Spawn-Layering and Focus Protection
* Newly added nodes in the visual workflow dashboard are now automatically selected, bringing them to the **front of the visual stack** (z-index elevation). This prevents newly created nodes from spawning hidden underneath large elements like the 500x500px `ViewerNode`.

---

## 📦 Package Distribution Security

* **`MANIFEST.in` Integration**: Added a package manifest file ensuring all non-Python data presets, forcefield `.itp` definition files, and CIF template configurations are correctly bundled during PyPI builds and source installs (`pip install`).
* **Sleek Root Cleanup**: Cleared non-essential scratch files, local workflows, and guide templates (`verify_bvs_h.py`, `.github/workflows/`, `CONTRIBUTING.md`) to deliver an ultra-clean public distribution.

---

*Thank you to all our users and contributors! For the full source code and visual builder guides, visit the [atomipy GitHub Repository](https://github.com/mholmboe/atomipy).*
