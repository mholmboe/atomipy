# MINFF for OpenMM via GROMACS Topology: Implementation Plan

**Audience:** AI IDE / coding assistant tasked with adding OpenMM support to the [atomipy](https://github.com/mholmboe/atomipy) package and the [MINFF](https://github.com/mholmboe/minff) force-field distribution.

**Status:** Architecture decided. Implementation is small (~documentation + thin helper + tests). No new force-field format, no parallel parameter database, no per-instance Python topology construction.

**Goal:** OpenMM users should be able to run MINFF simulations by feeding atomipy-generated GROMACS topology files (`.top` + `.gro`) directly to OpenMM's `GromacsTopFile` / `GromacsGroFile` loaders, with one well-documented helper function on the atomipy side to manage the `defines` translation. This is a thin wrapper, not a port.

---

## 1. Why This Approach

MINFF is already correctly implemented in GROMACS. The MINFF topology files atomipy produces contain:

- An `[ atomtypes ]` block listing all MINFF types with σ, ε, mass, default charge
- A `[ bonds ]` section with **per-instance** equilibrium bond lengths (Case B convention — explicit numerical `r₀` and `k` per atom-pair line)
- An `[ angles ]` section with **per-instance** equilibrium angles (explicit numerical `θ₀` and `k` per atom-triplet line)
- `#ifdef` blocks selecting between MINFF variants (e.g., `GMINFF_k500`, `OPC3`, `OPC3_IOD_LM`), activated at runtime via the GROMACS `.mdp` `define` directive

OpenMM's `GromacsTopFile` reader natively supports all of these:

- Per-instance bonded parameters are read line-by-line from the `[ bonds ]` and `[ angles ]` sections — exactly as GROMACS itself reads them.
- The OpenMM parser includes a built-in C-preprocessor that handles `#ifdef`, `#ifndef`, `#else`, `#endif`, `#define`, and `#include`. The `defines` argument to `GromacsTopFile.__init__` accepts a Python dict that plays the role of GROMACS's `-DVAR` flags.
- Combination rule 2 (Lorentz-Berthelot, MINFF's convention) is OpenMM's default.

**Therefore: no new XML format is needed.** The work reduces to (a) confirming atomipy writes a clean, self-contained `.top` file, (b) shipping a small helper that handles the `defines` translation cleanly, and (c) documenting the workflow.

---

## 2. atomipy Data and Files Already in Place

| Item | Location | Status |
|------|----------|--------|
| MINFF JSON parameter files | `atomipy/ffparams/GMINFF/*.json`, `atomipy/ffparams/TMINFF/*.json` | Existing |
| GROMACS `.itp` writers | `atomipy/write_*.py` | Existing |
| Full `.top` file writer | New in unreleased atomipy version (per user) | Exists, not yet on GitHub |
| `bond_angle` topology builder | `atomipy/bond_angle.py` | Existing |
| `.gro` file writer | `atomipy/write_gro.py` | Existing |
| `#ifdef`-based variant selection | Used by atomipy `.top` writer | Existing |

**Pre-implementation check:** The new `.top` writer should produce a file containing all of these sections in this order:

```
[ defaults ]              ; nbfunc, comb-rule, gen-pairs, fudgeLJ, fudgeQQ
[ atomtypes ]             ; all MINFF types
#ifdef GMINFF_k500
  [ angletypes ]          ; default per-triplet k (overridden per-instance below)
#endif
[ moleculetype ]          ; e.g. "Mineral"
[ atoms ]                 ; per-atom: nr, type, resnr, residu, atom, cgnr, charge, mass
[ bonds ]                 ; per-pair: ai aj funct r0 k        (numerical values)
[ angles ]                ; per-triplet: ai aj ak funct theta0 k   (numerical values)
[ system ]
[ molecules ]             ; molecule_name count
```

If any of these are missing, OpenMM's `GromacsTopFile` will raise a parse error pointing to the issue. The file must be parseable by `gmx grompp` first — if GROMACS likes it, OpenMM almost certainly will too.

---

## 3. What to Implement

The deliverable is small. Five items, in order of priority:

### 3.1 New module: `atomipy/openmm_interface.py`

A single-file module containing one user-facing function and one helper. No new dependencies beyond OpenMM itself (which is `pip install openmm`).

**Function signature:**

```python
def load_minff_into_openmm(
    top_path,
    gro_path,
    defines,
    include_dir=None,
    nonbonded_method=None,
    nonbonded_cutoff_nm=1.0,
    constraints=None,
    rigid_water=False,
    ewald_error_tolerance=5e-4,
    use_dispersion_correction=True,
):
    """
    Build an OpenMM Topology, System, and positions from an atomipy-generated
    GROMACS topology (.top) and coordinate (.gro) pair.

    This is the recommended OpenMM entry point for MINFF simulations. The
    function is a thin wrapper around openmm.app.GromacsTopFile and
    openmm.app.GromacsGroFile that handles the standard MINFF defaults:
    flexible bonds and angles on the mineral body, OPC3-style water at the
    user's discretion, Particle Mesh Ewald for electrostatics, and the
    Lorentz-Berthelot combination rule (MINFF's convention, also OpenMM's
    default).

    Parameters
    ----------
    top_path : str
        Path to a GROMACS .top file written by atomipy (e.g. via
        ap.write_top(atoms, Box=Box, file_path='kao.top', forcefield='minff')).
        Must be self-contained: include #include directives and #ifdef blocks
        for variant selection, but resolve to a complete topology once
        preprocessed.
    gro_path : str
        Path to the matching .gro coordinate file.
    defines : dict[str, str] or list[str]
        Preprocessor variables to activate, equivalent to the GROMACS .mdp
        directive `define = -DGMINFF_k500 -DOPC3_IOD_LM -DOPC3`. Either:
          - dict form: {'GMINFF_k500': '', 'OPC3_IOD_LM': '', 'OPC3': ''}
          - list form: ['GMINFF_k500', 'OPC3_IOD_LM', 'OPC3']  (auto-converted)
        The empty-string values are sufficient for plain #ifdef checks.
    include_dir : str or None
        Directory containing the MINFF .itp files referenced by `#include`
        directives in top_path. If None, OpenMM falls back to its default
        (typically /usr/local/gromacs/share/gromacs/top), which is almost
        certainly wrong for MINFF. Provide this explicitly in production.
    nonbonded_method : openmm.app constant, optional
        Defaults to app.PME. Other valid values: app.CutoffPeriodic,
        app.CutoffNonPeriodic, app.NoCutoff, app.Ewald, app.LJPME.
    nonbonded_cutoff_nm : float
        Real-space cutoff in nanometers (default 1.0).
    constraints : openmm.app constant or None
        Default None — MINFF requires *flexible* bonds and angles on the
        mineral body for the explicit angle terms to function. Do not set
        to AllBonds or HBonds for mineral simulations.
    rigid_water : bool
        Default False (flexible water). Set True if running with a rigid
        water model (TIP3P-rigid, SPC-rigid, OPC3-rigid) — OpenMM will then
        apply SETTLE constraints to water molecules automatically.
    ewald_error_tolerance : float
        PME tolerance. Default 5e-4 matches OpenMM's standard.
    use_dispersion_correction : bool
        Long-range LJ correction. Default True matches GROMACS default.

    Returns
    -------
    topology : openmm.app.Topology
    system : openmm.System
    positions : openmm.unit.Quantity (N, 3) array in nanometers
    """
```

**Implementation (full body, ~30 lines):**

```python
import openmm as mm
import openmm.app as app
from openmm import unit


def _normalize_defines(defines):
    """Accept either a dict or a list of names; return a dict."""
    if isinstance(defines, dict):
        return dict(defines)
    return {name: '' for name in defines}


def load_minff_into_openmm(
    top_path,
    gro_path,
    defines,
    include_dir=None,
    nonbonded_method=None,
    nonbonded_cutoff_nm=1.0,
    constraints=None,
    rigid_water=False,
    ewald_error_tolerance=5e-4,
    use_dispersion_correction=True,
):
    """[docstring as above]"""
    if nonbonded_method is None:
        nonbonded_method = app.PME
    defines_dict = _normalize_defines(defines)

    gro = app.GromacsGroFile(gro_path)
    top = app.GromacsTopFile(
        top_path,
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=include_dir,
        defines=defines_dict,
    )
    system = top.createSystem(
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=nonbonded_cutoff_nm * unit.nanometer,
        constraints=constraints,
        rigidWater=rigid_water,
        ewaldErrorTolerance=ewald_error_tolerance,
        useDispersionCorrection=use_dispersion_correction,
    )
    return top.topology, system, gro.positions
```

That is the entire core implementation. Anything more elaborate is over-engineering.

### 3.2 Expose the function in atomipy's public API

In `atomipy/__init__.py`, add:

```python
try:
    from .openmm_interface import load_minff_into_openmm
except ImportError:
    # OpenMM not installed; the function is unavailable but the rest of
    # atomipy still works. This is the same pattern atomipy uses for
    # other optional integrations.
    pass
```

The `try/except` is so that users without OpenMM don't get an `ImportError` when they `import atomipy`. The function only becomes available if OpenMM is importable.

### 3.3 Demo script: `scripts/run_openmm_minff.py`

A single end-to-end example script, parallel in style to your existing GROMACS demo scripts (e.g., `run_minff_atomi.py`). It should do the full pipeline on a known test mineral (kaolinite):

```python
"""
run_openmm_minff.py
===================

End-to-end demo: take a CIF or GRO of a clay mineral, type it with MINFF
in atomipy, write GROMACS topology and coordinate files, then load them
into OpenMM and run a short simulation.

This is the recommended OpenMM workflow for MINFF.
"""
import atomipy as ap
import openmm as mm
import openmm.app as app
from openmm import unit

# 1. Standard atomipy preprocessing (unchanged from GROMACS workflow)
atoms, Box = ap.import_gro('Kaolinite_GII_0.0487.gro')
atoms = ap.element(atoms)
atoms = ap.minff(atoms, Box=Box)
atoms, _, _ = ap.bond_angle(atoms, Box=Box)

# 2. Write GROMACS topology + coords (existing atomipy functionality)
ap.write_gro(atoms, Box=Box, file_path='kao.gro')
ap.write_top(atoms, Box=Box, file_path='kao.top', forcefield='minff')

# 3. Load into OpenMM via the helper. The `defines` dict is the exact
#    equivalent of the GROMACS .mdp directive:
#       define = -DGMINFF_k500 -DOPC3_IOD_LM -DOPC3
topology, system, positions = ap.load_minff_into_openmm(
    top_path='kao.top',
    gro_path='kao.gro',
    defines=['GMINFF_k500', 'OPC3_IOD_LM', 'OPC3'],
    include_dir='atomipy/ffparams',
    nonbonded_cutoff_nm=1.0,
)

# 4. Single-point energy as a sanity check
sp_integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
sp_context = mm.Context(system, sp_integrator)
sp_context.setPositions(positions)
E0 = sp_context.getState(getEnergy=True).getPotentialEnergy()
print(f'Single-point potential energy: {E0}')

# 5. Energy minimization
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.001 * unit.picosecond,
)
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy(maxIterations=200)
E1 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f'After minimization:           {E1}')

# 6. Short production run
simulation.reporters.append(
    app.StateDataReporter('out.csv', 100, step=True,
                          potentialEnergy=True, temperature=True))
simulation.reporters.append(app.DCDReporter('out.dcd', 100))
simulation.step(1000)
print('Done.')
```

This file serves as both a runnable demo and the canonical reference for users figuring out the API.

### 3.4 Unit tests: `tests/test_openmm_interface.py`

Three small tests, all runnable in CI provided OpenMM is installed:

**Test 1: parsing works.**

```python
def test_top_file_parses():
    """The atomipy-generated .top must load via GromacsTopFile without errors."""
    # Use a small test mineral checked into the repo
    topology, system, positions = ap.load_minff_into_openmm(
        top_path='tests/data/kaolinite.top',
        gro_path='tests/data/kaolinite.gro',
        defines=['GMINFF_k500', 'OPC3'],
        include_dir='atomipy/ffparams',
    )
    n_atoms = system.getNumParticles()
    assert n_atoms > 0
    # The topology and system must agree on atom count
    assert sum(1 for _ in topology.atoms()) == n_atoms
```

**Test 2: forces are correctly registered.**

```python
def test_forces_present():
    """The system must contain exactly one each of the three expected forces."""
    _, system, _ = ap.load_minff_into_openmm(
        top_path='tests/data/kaolinite.top',
        gro_path='tests/data/kaolinite.gro',
        defines=['GMINFF_k500', 'OPC3'],
        include_dir='atomipy/ffparams',
    )
    force_types = [type(f).__name__ for f in system.getForces()]
    assert 'HarmonicBondForce' in force_types
    assert 'HarmonicAngleForce' in force_types
    assert 'NonbondedForce' in force_types
```

**Test 3: angle parameters are per-instance, not collapsed by type.**

This is the critical correctness check. If OpenMM's parser failed to read per-instance values, all Si–O–Si angles would have the same θ₀ — and the test would catch it.

```python
def test_per_instance_angles_preserved():
    """
    Verify that per-instance angle equilibrium values from the [ angles ]
    section are read into the HarmonicAngleForce individually, not collapsed
    to a single per-triplet-type default.

    For a real clay mineral with multiple distinct Si-O-Si sites, the θ₀
    values for those sites should differ measurably (>0.5° spread).
    """
    import math
    _, system, _ = ap.load_minff_into_openmm(
        top_path='tests/data/kaolinite.top',
        gro_path='tests/data/kaolinite.gro',
        defines=['GMINFF_k500', 'OPC3'],
        include_dir='atomipy/ffparams',
    )
    angle_force = [f for f in system.getForces()
                   if type(f).__name__ == 'HarmonicAngleForce'][0]
    # Collect all theta0 values for angles among heavy atoms
    theta0_values = []
    for i in range(angle_force.getNumAngles()):
        _, _, _, theta0, _ = angle_force.getAngleParameters(i)
        theta0_values.append(theta0.value_in_unit_system(
            __import__('openmm').unit.md_unit_system))
    # In a real kaolinite there are many unique geometries; spread must be > 1°
    spread_deg = (max(theta0_values) - min(theta0_values)) * 180.0 / math.pi
    assert spread_deg > 1.0, (
        f'Angle θ₀ spread is only {spread_deg:.3f}°, suggesting '
        'per-instance values were collapsed. Check the .top file format.'
    )
```

Test 3 is the one that catches the most likely failure mode (silently degraded parsing). It's worth keeping even if it looks paranoid.

### 3.5 Documentation update

Add a section to the atomipy `README.md`, slotted in alongside the existing GROMACS/LAMMPS/NAMD sections:

```markdown
### Running MINFF simulations in OpenMM

atomipy generates GROMACS topology files that OpenMM can load directly via
its built-in `GromacsTopFile` parser. The recommended workflow is:

```python
import atomipy as ap

# Type structure and write GROMACS .top + .gro (same as for native GROMACS)
atoms, Box = ap.import_gro('kao.gro')
atoms = ap.element(atoms)
atoms = ap.minff(atoms, Box=Box)
ap.write_gro(atoms, Box=Box, file_path='kao_out.gro')
ap.write_top(atoms, Box=Box, file_path='kao.top', forcefield='minff')

# Load into OpenMM. The `defines` argument is the OpenMM equivalent of the
# GROMACS .mdp directive `define = -DGMINFF_k500 -DOPC3_IOD_LM -DOPC3`.
topology, system, positions = ap.load_minff_into_openmm(
    top_path='kao.top',
    gro_path='kao_out.gro',
    defines=['GMINFF_k500', 'OPC3_IOD_LM', 'OPC3'],
    include_dir='atomipy/ffparams',
)

# Simulate as usual
import openmm as mm
import openmm.app as app
from openmm import unit

integrator = mm.LangevinMiddleIntegrator(
    300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()
simulation.step(10000)
```

**Important notes:**
- Set `constraints=None` and `rigid_water=False` (the defaults) for mineral
  body flexibility. MINFF's explicit angle terms require flexible bonds and
  angles to function correctly.
- Use `rigid_water=True` only if you are running a constrained water model
  (e.g., SETTLE-constrained OPC3) — OpenMM will then apply the constraints
  automatically.
- The `include_dir` parameter must point to wherever the MINFF `.itp` files
  referenced by `#include` directives in your `.top` live. For the default
  atomipy distribution that is `atomipy/ffparams/`.
```

That's the entire documentation addition.

---

## 4. Sanity Validation (Manual, One-Time)

Before declaring v1 done, run this comparison once:

1. Take a small MINFF structure (kaolinite is canonical).
2. In GROMACS:
   ```bash
   gmx grompp -f minim.mdp -c kao.gro -p kao.top -o sp.tpr
   gmx mdrun -deffnm sp -nsteps 0
   ```
   Record the total potential energy from `md.log`, broken down by component:
   bonds, angles, LJ-SR, Coulomb-SR, Coulomb-recip (PME), LJ-recip if
   `LJ-PME` is active.

3. In OpenMM, using the same input geometry:
   ```python
   topology, system, positions = ap.load_minff_into_openmm(
       'kao.top', 'kao.gro',
       defines=['GMINFF_k500', 'OPC3_IOD_LM', 'OPC3'],
       include_dir='atomipy/ffparams',
   )
   # Force-group decomposition for component-wise comparison
   for i, force in enumerate(system.getForces()):
       force.setForceGroup(i)
   integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
   ctx = mm.Context(system, integrator)
   ctx.setPositions(positions)
   for i, force in enumerate(system.getForces()):
       E = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy()
       print(f'{type(force).__name__}: {E}')
   ```

4. **Expected agreement: within ~0.1% on the total potential energy** for
   identical input geometry and matched PME settings. Component-by-component
   agreement should be even tighter for the bonded forces and a few percent
   for the electrostatic components (PME settings differ slightly between
   codes).

5. If a single component disagrees by more than a few percent, investigate:
   - Bond/angle force constant convention (should be identical, factor-of-½)
   - PME parameters (`ewaldErrorTolerance`, real-space cutoff vs GROMACS
     `fourierspacing`)
   - Combination rule (must be 2 for MINFF)
   - 1–4 exclusions (`fudgeLJ`, `fudgeQQ` in `[ defaults ]`)
   - Long-range LJ correction (`useDispersionCorrection=True` vs GROMACS
     `DispCorr = EnerPres`)

This manual validation is not automatable without a GROMACS installation in
CI, but it only needs to be done once per major MINFF release.

---

## 5. Pre-Implementation Checks

Confirm these before starting the implementation work. Each is a 5-minute task.

### 5.1 `ap.write_top` produces a self-contained `.top` file

Open a generated `.top` and verify it has, in order:
- `[ defaults ]` (with `nbfunc=1`, `comb-rule=2`, `gen-pairs=no`, `fudgeLJ`, `fudgeQQ` — MINFF's standard combination)
- `[ atomtypes ]` covering every type used in the structure
- One or more `[ moleculetype ]` ... `[ atoms ]` ... `[ bonds ]` ... `[ angles ]` blocks
- `[ system ]`
- `[ molecules ]`

If `[ defaults ]` is missing or the `#include` directives don't resolve, `GromacsTopFile` will error out with a clear message. Fix in the `.top` writer, not in the OpenMM loader.

### 5.2 The `[ angles ]` section contains per-instance numerical values

A line should look like:

```
;  ai   aj   ak  funct    theta0       k
   1    7   13     1      150.234    500.0
```

with explicit numerical θ₀ on every line (Case B). If lines look like:

```
   1    7   13     1
```

with no parameters, the topology relies on `[ angletypes ]` defaults — which would be Case A, and a different (simpler) story. The user has confirmed Case B, but verifying once is cheap.

### 5.3 `#include` paths in the `.top` are relative or known

The `#include` directives in atomipy's `.top` should reference paths that the user can resolve via the `include_dir` parameter. Convention: `#include "ffnonbonded.itp"` (relative, found via `include_dir`), not `#include "/absolute/path/ffnonbonded.itp"`.

If any include paths are absolute, change the `.top` writer to use relative paths.

### 5.4 OpenMM is at least version 8.0 (recommended)

OpenMM ≥ 7.5 supports the `defines` argument, but the `useDispersionCorrection` argument to `createSystem` was added later. The implementation should target OpenMM ≥ 8.0 to keep things simple. State the minimum version in the module docstring and in the atomipy README.

---

## 6. What Is Deliberately NOT in V1

These were considered and rejected for the v1 OpenMM deliverable. The IDE should not introduce them by accident.

- **A standalone OpenMM XML force-field file.** Not needed. The GROMACS `.top` carries all the parameters and OpenMM reads it natively.
- **A Python-API path that returns an `openmm.System` built atom-by-atom.** Not needed. `GromacsTopFile.createSystem(...)` does this.
- **A new JSON-to-XML converter.** Not needed.
- **A custom preprocessor for `#ifdef`.** Not needed; OpenMM has one.
- **A standalone `minff_openmm.py` for users without atomipy.** Not needed at v1 — users without atomipy can still write a 4-line OpenMM script around `GromacsTopFile` directly, since the `.top` files are standard. Documentation in the MINFF repo README is sufficient.

If demand emerges later for any of these (e.g., users who want pure-OpenMM XML for SystemGenerator integration, or in-memory System construction for ML/MD coupling), they become v2 features. Do not preemptively build them.

---

## 7. Acceptance Criteria

V1 is complete when:

- [ ] `atomipy/openmm_interface.py` exists, implements `load_minff_into_openmm()`, is ~50 lines of code.
- [ ] `ap.load_minff_into_openmm` is importable from the atomipy top-level namespace when OpenMM is installed.
- [ ] The function gracefully degrades (no ImportError on `import atomipy`) when OpenMM is not installed.
- [ ] `scripts/run_openmm_minff.py` runs end-to-end on `Kaolinite_GII_0.0487.gro` (or the canonical test mineral) without errors.
- [ ] The three unit tests in `tests/test_openmm_interface.py` pass. Test 3 (per-instance angle preservation) is mandatory.
- [ ] The atomipy README has a working "Running MINFF simulations in OpenMM" section.
- [ ] The manual GROMACS vs OpenMM single-point energy comparison on kaolinite agrees within ~0.1% (recorded in the commit message or PR description).
- [ ] OpenMM is listed as an optional dependency in `setup.py` / `pyproject.toml` with extras_require: `pip install atomipy[openmm]`.

---

## 8. File Layout to Produce

```
atomipy/
├── __init__.py                          (UPDATE: conditional import)
├── openmm_interface.py                  (NEW: ~50 lines)
└── ...

scripts/
└── run_openmm_minff.py                  (NEW: ~40 lines end-to-end demo)

tests/
├── data/
│   ├── kaolinite.top                    (NEW: checked-in test fixture)
│   └── kaolinite.gro                    (NEW: checked-in test fixture)
└── test_openmm_interface.py             (NEW: 3 small tests)

README.md                                (UPDATE: add "Running MINFF in OpenMM")
pyproject.toml                           (UPDATE: openmm as optional dep)
```

Total estimated work: a few hundred lines including tests and documentation, half a day at most.

---

## 9. Future Extensions (Out of Scope for V1)

For reference only — do not implement now.

- **Native OpenMM XML output**: a system-specific XML enumerating all bonds and angles with their explicit r₀ / θ₀. Useful for SystemGenerator-style workflows. ~200 lines of code if ever needed.
- **In-memory System builder**: `to_openmm_system(atoms, Box, ...)` that returns a `System` object directly without round-tripping through files. Useful for ML/MD coupling (e.g., Δ-learning with `openmm-torch`). ~150 lines.
- **OpenMM-Setup integration**: a guide for using MINFF in OpenMM-Setup's web GUI. Documentation only; no code.
- **OpenMM XTC reporter aliasing**: convenience helpers for matching GROMACS-compatible trajectory output. ~20 lines.

These are all worth considering once v1 is in users' hands and you have feedback. None is a prerequisite for OpenMM support.

---

## 10. Quick Reference: The Entire User-Facing API

After v1, this is what an OpenMM-curious atomipy user needs to know:

```python
import atomipy as ap

# Standard atomipy preprocessing
atoms, Box = ap.import_gro('mineral.gro')
atoms = ap.element(atoms)
atoms = ap.minff(atoms, Box=Box)

# Write GROMACS topology (existing functionality)
ap.write_gro(atoms, Box=Box, file_path='out.gro')
ap.write_top(atoms, Box=Box, file_path='out.top', forcefield='minff')

# Load into OpenMM (new functionality)
topology, system, positions = ap.load_minff_into_openmm(
    top_path='out.top',
    gro_path='out.gro',
    defines=['GMINFF_k500', 'OPC3'],
    include_dir='atomipy/ffparams',
)

# Standard OpenMM from here on
```

That's the whole story. Everything else is standard OpenMM that the user already knows or can learn from the OpenMM docs.

---

## 11. References

- atomipy: <https://github.com/mholmboe/atomipy>
- MINFF: <https://github.com/mholmboe/minff>
- OpenMM `GromacsTopFile` API: <https://docs.openmm.org/latest/api-python/generated/openmm.app.gromacstopfile.GromacsTopFile.html>
- OpenMM user guide, "Using GROMACS Files": <https://docs.openmm.org/latest/userguide/application/02_running_sims.html#using-gromacs-files>
- GROMACS preprocessor reference (`#ifdef`, `define`): <https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html>

End of implementation plan.
