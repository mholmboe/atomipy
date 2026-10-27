# Changelog

All notable changes to **atomipy** are documented here.

## 1.1

### Diffraction / scattering factors (`xrd`)
- **Virtual / dummy sites now contribute zero scattering.** TIP4P `MW`, TIP5P
  `LP*`, and `Du`/`DUM`/`MY1` sites carry charge but no electrons; `xrd()`
  previously let them fall through to the oxygen fallback and silently add
  spurious scattering (and they could mis-size the structure-factor sum). They
  are now given f = 0, so e.g. a TIP4P/TIP5P water box gives the correct pattern.
- **Explicit tabulated ion labels are honored.** `xrd()` ran `element()` first,
  which collapsed an atom typed `Fe2+` to `Fe` (→ Fe1.5+), `Hf4+` to `H`, etc.
  It now uses the atom's original label verbatim when it is a tabulated
  Waasmaier–Kirfel row (e.g. `Fe2+`, `Fe3+`, `Hf4+`, `K1+`) before falling back to
  the element collapse — so explicit oxidation states scatter correctly, matching
  the MATLAB toolbox. Clay forcefield types (`Feo`, `Si`, …), which are not table
  rows, are unchanged.
- **Lowercase aqueous-ion labels resolve correctly.** The bare-ion remap in
  `atomic_scattering_factors` is now case-insensitive, so `cu`/`na`/`mg` give the
  aqueous ion instead of the oxygen fallback.

## 1.0

### Diffraction / scattering factors
- **Comprehensive clay-forcefield atom-type remap in `atomic_scattering_factors`.**
  All CLAYFF/MINFF structural site types now resolve to their intended
  partial-charge scattering factor (case-insensitively) instead of silently
  falling through to the neutral element or oxygen:
  - `Al / Alt / Alo / Ale / Alr / Am / ao / ae / at → Al1.5+`
  - `Si / Sit / Sio / st → Si2+`
  - `Mgo / Mgh / Mge → Mg1+`
  - `Fe / Feo / Fet / Feo2 / Feo3 / Fet2 / Fet3 / Fee / Fee2 / Fee3 / Fe2 / DFeo → Fe1.5+`
  - `Lio → Li1+`, `Tio → Ti2+`, `Cao / Cah → Ca2+`
  Bare `Mg / Ca / Sr / Ba / Cu` stay aqueous 2+; `Mn / Mn2 / Mn3 / Mn4` are left
  neutral on purpose (their CLAYFF digits are site labels, not oxidation states,
  and not monotonic in charge). Previously only bare `Al` / `Si` / `Fe` (and
  `Feo` / `Fet`) were handled — the rest gave the wrong element (e.g. `Sit` read
  as neutral Si, +17%; many types fell back to oxygen). Only direct
  `atomic_scattering_factors` calls are affected; `xrd()` already collapsed these
  via `element()`. The intended partial charge is force-field dependent (MINFF and
  CLAYFF disagree per label); a single fixed mapping works only because the
  tabulated Waasmaier–Kirfel rows are coarse enough that both charges round to the
  same row. Settled jointly with the canonical `mholmboe/atom` MATLAB toolbox.

## 0.99

### Diffraction / scattering factors
- **Corrected the Fe1.5+ scattering-factor row.** The CLAYFF partial-charge
  `Fe1.5+` row (octahedral Fe carries Al's 1.575 charge) had f(0)=23.976 e,
  ~2% below the correct 24.5, and was the sole outlier among the curated
  partial-charge rows (Si2+, Al1.5+, Mg1+, O1-). Re-fit to the same recipe as
  its siblings — a Waasmaier–Kirfel fit to mean(Fe⁰, Fe³⁺) with f(0) pinned to
  24.5 — so f(0)=24.500 and the XRD intensities of Fe-bearing clays are no
  longer biased low at low angle. Coordinated with the canonical `mholmboe/atom`
  MATLAB toolbox (same corrected row).
- **`atomic_scattering_factors` now remaps bare `Fe`/`Feo`/`Fet` → `Fe1.5+`**
  (was neutral Fe), consistent with the existing `Al→Al1.5+` / `Si→Si2+`
  clay-cation remaps and the MATLAB toolbox. `xrd()` already applied this remap,
  so only direct `atomic_scattering_factors` calls change.

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
