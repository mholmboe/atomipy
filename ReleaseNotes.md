# Release Notes: atomipy v0.97 (from v0.96)

*Released 2026-06-24*

atomipy **v0.97** adds a local **GROMACS simulation engine** and a suite of
**trajectory-analysis** tools, and hardens solvation.

---

## ⚙️ Local GROMACS engine (`atomipy.gromacs`)
* Run **grompp + mdrun** directly from an atomipy system, reusing the MINFF/CLAYFF
  topology writers: `detect_gmx`, `mdp`/`build_defines`, `stage_minff`, `run_stage`,
  `run_pipeline`, `run_local_gmx`, `trjconv_to_pdb`, `energy_timeseries`.
* Single-stage runs (EM **or** NVT **or** NPT) — chain them in any order; accepts a
  verbatim `mdp_text` to override the generated `.mdp`. EM writes a `.trr` (steepest
  descent ignores `.xtc`) via `nstxout`.

## 📈 Trajectory analysis (`atomipy.analysis`)
* `density_profile` / `density_frames` (number / mass / charge, along x/y/z).
* `rdf_frames` and `calculate_rdf(return_cn=True)` — ensemble g(r) **plus** running
  coordination *n(r)*.
* `msd` (PBC-unwrapped, multi-origin; 3D / 2D-xy / 1D-z) + `displacement_distribution`
  (van Hove self-part).
* `vacf` — finite-difference velocity autocorrelation, Green–Kubo *D*, and the
  vibrational power spectrum / DOS.
* `find_hbonds` / `hbonds_frames` — geometric donor–H···acceptor analysis
  (gmx-hbond convention), selectable by atom type and/or residue name.

## 💧 Solvation & transforms
* `solvate`: isotropically scales the water template to the target **density**,
  removes solvent–solvent overlaps, renumbers molids contiguously, and treats an
  integer `max_solvent` as an exact count (raises if the box can't hold it).
* `scale(atoms, Box, scale_factors)` reimplemented via fractional coordinates —
  scalar = isotropic, length-3 = per-axis; correct for triclinic cells.

## 🚀 Use it without local setup
* This release powers the web-module's **GROMACS** Simulate node, including
  **GPU GROMACS on Google Colab** via the app launcher's optional GROMACS cell.
  See `atomipy-web-module`.

---

# Release Notes: atomipy v0.96 (from v0.95)

*Released 2026-06-06*

atomipy **v0.96** is a large consolidation release gathering everything developed
since v0.95. Highlights: an OpenMM simulation bridge, a hub-and-spoke topology
interchange and multi-component merging engine, mixed organic/mineral systems via
GAFF (ACPYPE) and OpenFF Sage, a new **Dummy FF** for non-MINFF inorganics, a
rules-based oxidation-state guesser, Chemical-JSON I/O, bundled organic (428) and
inorganic-crystal (517) libraries, and provenance headers on every exported file.

> **Note:** ParmEd is no longer a dependency. Mixed organic/mineral parametrization
> now runs entirely through ACPYPE (GAFF) and the pure-Python OpenFF Interchange stack.

---

## 🧪 Molecular Dynamics with OpenMM
* **GROMACS → OpenMM bridge** (`atomipy.openmm_interface`, `load_minff_into_openmm`): links GROMACS `.top`/`.itp` files directly into OpenMM `System`/`Topology` objects for MD from Python. Falls back to a PDB coordinate file when no GRO is present.
* **Auto platform selection**: prefers CUDA/OpenCL GPUs, falling back to CPU/Reference.
* **Dynamic parameter injection**: flexible water models (OPC3/SPC(/E)/TIP3P), custom mineral bond/angle parameters, and ion corrections handled at load time.
* **Iterative L-BFGS energy minimization** with real-time trajectory recording and maximum-force-norm logging (kJ/mol/nm).
* **UTF-8-safe** trajectory/PDB I/O (handles unicode ion labels such as `Cl⁻`).

## 🔀 Topology Interchange & Multi-Component Merging
* **Hub-and-spoke topology interchange** package: a `Topology` object exports to GROMACS `.itp`/`.top`, LAMMPS `.data`, CHARMM `.psf`/`.prm`, and JSON/XML.
* **`atomipy.merge_top`**: merge multiple CLAYFF/MINFF mineral + organic/ion/water topologies into one consistent system, with combination-rule alignment, water-model `#define` emission, and ion-model `#define` resolution.
* **`atomipy.composition`** engine: classify atoms (water/ion/organic/mineral), count residues/components and system mass; recognize amino-acid, DNA/RNA-nucleotide, and zeolite residue names.
* **Nameable minerals**: mineral `moleculetype` follows the residue name, with uniquified names (`MIN`, `MIN_1`, …).
* **Explicit O–M–O / M–O–H angles** generated in merged mineral topologies (with a flag to omit `[angles]`).
* API: the GROMACS `.top` writer is exposed as **`write_gmx_top`**.
* Ion-spelling tolerance (`Na` vs `Na+`, `Ca` vs `Ca2+`).

## 🧬 Mixed Organic / Mineral Systems
* Parametrize organics from **SMILES or file** via **GAFF-2.11 / GAFF-1** (ACPYPE, bundles antechamber — no separate AmberTools install) or **OpenFF Sage / Parsley** (pure-Python OpenFF Interchange).
* Standalone, dockerized **OpenFF worker** microservice that isolates the OpenFF stack.
* Merge parametrized organics into mineral models seamlessly via `merge_top`.

## 🧱 Dummy FF — qualitative model for non-MINFF inorganics
* **`assign_dummy_mineral_params`**: a frozen-framework model that lets MINFF-unsupported inorganics interact with water/solutes electrostatically (EM/NVT only).
* **Pauling effective charges** `q_eff = q_formal·[1 − exp(−¼(χ_O − χ_M)²)]`; H fixed at **+0.4**; **F** supported.
* **MINFF coordination-resolved oxygen charges** `q_O = −2.0 + Σⱼ (Formalⱼ − Partialⱼ)/CNⱼ`.
* **Self-calculated per-element LJ** from vdW radii (UFF) by default; element-appropriate LJ for pure metals (Heinz).
* Exposes MINFF global cutoffs (`rmaxlong`, `rmaxH`); coordination typing is purely geometric (ignores molids).
* Self-contained dummy+water topology writer; supports organics/ions; **water/ions/organics are left untouched** (identical to a mixed MINFF run).

## ⚛️ Oxidation States & Charges
* **`guess_oxidation_states`**: a rules-based oxidation/formal-charge guesser (H +1, O −2, Si +4, Al +3, Mg +2, Ti +4, alkali +1, alkaline-earth +2, halogens −1, …) with ionic and electronegativity engines — a fast alternative to BVS.

## 📁 File Formats & I/O
* **Chemical JSON (`.cjson`)** import *and* export.
* **Bundled organic molecule library** (428 molecules: amino acids, nucleotides, sugars, alcohols, …) with **L-before-D** enantiomer ordering.
* **Bundled inorganic crystal library** (517 CIFs) with a gemmi-based CIF/mmCIF reader (symmetry expansion). **gemmi is now a core dependency.**
* `write_sdf` emits a proper **V2000 bonds block** with bond orders.
* **`fit_box(atoms, padding, cubic, center)`**: size an orthogonal box to the molecule + margin (equivalent to `gmx editconf -d`).
* Clarified/expanded import support via `import_auto`: `.pdb`, `.gro`, `.xyz`, `.cif`/`.mmcif`, `.poscar`/`.contcar`, `.pqr`, `.cjson`.
* **Provenance headers**: every written structure (`.pdb`/`.gro`/`.xyz`/`.cif`/`.poscar`/`.sdf`/`.pqr`) and topology (`.itp`/`.top`/`.psf`/`.data`/`.prm`) file now carries a `Generated by atomipy 0.96 on <date>` comment (a `_generator` key for JSON/cjson, which cannot hold comments).

## 🛡️ Robustness
* **Empty-box safety**: a single zero-charge/zero-mass `Dummy` atom is instantiated for empty configurations to keep visualization/physics engines happy.
* Robust CIF parsing; PBC/MIC neighbour-list fixes for skewed cells; replication off-by-one fix; `fuse_atoms` deepcopy fix.

## 🌐 Ecosystem
* See the dedicated release notes for the web **Visual Builder** (`atomipy-web-module` v0.4.0) and the **Topology Generator** (`atomipy-topology-generator` v0.1.0), both of which embed atomipy 0.96.

---

*Thank you to all our users and contributors! For the full source code and guides, visit the [atomipy GitHub repository](https://github.com/mholmboe/atomipy).*
