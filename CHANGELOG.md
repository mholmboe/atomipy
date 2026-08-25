# Changelog

All notable changes to **atomipy** are documented here.

## 0.98

### Force field (MINFF / CLAYFF)
- **Corrected MINFF parameters for Ca-minerals.** Fixed the Lennard-Jones and
  charge parameters for the calcium atom types `Cao`/`Cah` and the fluorine
  type `Fs` — i.e. minerals such as lime (CaO), portlandite (Ca(OH)₂) and
  fluorite (CaF₂). The Ca-mineral atom-type ordering was corrected and the
  changes were synced across the `.itp`, GMINFF/TMINFF `.json`, and
  `par_minff.prm` files, with refreshed `UC_conf` reference structures.
- Synced `min.ff` to the canonical [`mholmboe/minff`](https://github.com/mholmboe/minff):
  re-optimised Hectorite-F, added `Fee3` in `GMINFF_k1500`.
- **Signed ion types throughout.** All ions carry an ASCII sign uniformly across
  `min.ff` and the JSON parameter files (`Na+`, `Cl-`, `Ca2+`, `Al3+`, …), with
  unsigned moleculetype/residue names. `minff()`/`clayff()` now assign the full
  formal charge to every monatomic ion (previously several were left at 0).
- CHARMM/NAMD PSF+PRM writers use sanitized, sign-free ion type names.

### Performance
- **XRD ~14× faster, ~44× less memory.** `diffraction.xrd` now factors the
  phase term into per-axis exponential tables reduced to a BLAS matrix product,
  applies Friedel's law, and accumulates peak profiles slice-wise. Results are
  numerically identical to 0.97 (float roundoff only).

### Features
- `reduce_supercell` — the inverse of `replicate_system` (fold a replicated
  system back to a single unit cell).
- Trajectory readers: pure-Python DCD reader, optional `libxdrfile` wrapper for
  GROMACS `.xtc`/`.trr`, optional mdtraj backend, and a general `trjconv()`.
- Ditrigonal distortion (tetrahedral rotation angle) analysis.
- GROMACS Dummy-FF: framework freezing (OpenMM parity) and atomtype hoisting so
  organics work alongside the dummy mineral; Shannon r_min LJ mode.

### Fixes
- `move.center()` guards against an empty atom list instead of dividing by zero.
- Assorted fixes: element mis-assignment, RDF range, coordination filter,
  `Box_dim2Cell`, centre-of-mass, PDB residue handling, solvate water–water
  declash for pure-water boxes.

## 0.97
- Baseline release (2026-08-02): local GROMACS engine + trajectory analysis.
