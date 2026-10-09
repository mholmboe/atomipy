# Changelog

All notable changes to **atomipy** are documented here.

## Unreleased

### Force field (MINFF v1.0 define names)
MINFF v1.0 selects parameter sets with different GROMACS defines: the angle force
constant is a define of its own (`-DMINFF_k500`), and a tailored (TMINFF) set needs
the mineral as well (`-DMontmorillonite -DMINFF_k500`). The old `-DGMINFF_k500` and
`-DMontmorillonite_k500` no longer select anything (there is no shim in `min.ff`).
- **Bundled `min.ff` and parameter files synced wholesale** to the canonical
  [`mholmboe/minff`](https://github.com/mholmboe/minff) v1.0: the per-k and
  all-k `ffnonbonded_tminff*.itp` files are replaced by `ffnonbonded_tminff.itp`,
  `ffbonded_tminff.itp`, `forcefield_tminff.itp`, `ffbonded_gminff.itp` and
  `forcefield_gminff.itp`; the GMINFF/TMINFF `.json` files carry the corrected
  Ca-mineral parameters, `Fee3` in k1500, the `MW` sites of 4-site waters and ASCII ion signs.
  TMINFF users should include `forcefield_tminff.itp`.
- **JSON block keys:** the general keys are renamed `GMINFF_k500` -> `MINFF_k500` (all four
  force constants); the tailored keys are unchanged (`Montmorillonite_k500`, ...).
  `load_forcefield(..., blocks=['GMINFF_k500'])` must become `blocks=['MINFF_k500']`.
- **Defines:** new `atomipy.minff_defines.minff_defines()` builds the defines for a
  general or tailored set, and `gromacs.build_defines()` gains `mineral=`
  (`build_defines(mineral="Montmorillonite")` -> `-DMontmorillonite -DMINFF_k500 -DFLEXIBLE`).
  The default `minff_variant` of `build_defines`, `write_merged_top` and `merge_and_write`
  is now `'MINFF_k500'`. The old spellings `'GMINFF_k500'` and `'<Mineral>_k500'` are still
  accepted for one release and are converted, with a `FutureWarning`; a tailored set that
  lacks the force constant (or the mineral) raises `ValueError` instead of selecting nothing.
  `write_merged_top` writes the general sets only and rejects a tailored `minff_variant`.
- Docs, examples and scripts use the new names.

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
