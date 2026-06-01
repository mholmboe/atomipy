# atomipy.topology — interchange layer

Hub-and-spoke topology interchange built around one per-instance-first
`Topology` data model. Every format is a reader/writer that converts to/from it,
so cross-format conversion is two hops through the hub (no N×N converters).

Implements the spec `atomipy_topology_interchange_plan.md` (resolved decisions §14).

## Quick start
```python
from atomipy.import_topology import read_itp, read_json, read_data, read_psf, from_atoms_box
from atomipy.write_topology  import write_itp, write_json, write_data, write_psf, to_atoms_box

top = read_itp("system.itp", defines=["GMINFF_k500"])   # -> Topology (params inline)
write_data(top, "system.data", units="real", emit_input_snippet="in.settings")  # unique terms -> unique types
write_json(top, "system.json")                            # canonical artifact

top2 = from_atoms_box(atoms, box)                         # bridge atomipy state
atoms2, box2 = to_atoms_box(top2)
```

## Key decisions (see plan §14)
- **Per-instance-first.** Every bond/angle/dihedral stores its own `params`
  inline; shared type tables are *derived on demand* (`typing.extract_types`).
  This is what gives unique terms per atom *site* (MINFF), which SMIRNOFF cannot.
- **Canonical units: nm / kJ·mol⁻¹ / degrees / e / amu** (never Å). The
  factor-of-2 prefactor (GROMACS ½·k·x² vs CHARMM/LAMMPS K·x²) is handled in
  `functional_forms`, separately from units.
- **Zero-dependency core**: stdlib dataclasses + hand-written validator +
  `schema/topology.schema.json`. Optional `jsonschema`/`pint` are import-guarded.

## Modules
| Module | Role |
|---|---|
| `topology/model.py` | `Topology, Atom, AtomType, Bond, Angle, Dihedral, Improper, Pair, Box, …` + JSON (de)serialization, triclinic box math |
| `topology/units.py` | canonical unit system + conversions |
| `topology/functional_forms.py` | form registry: param kinds, GROMACS funct ↔ form, LAMMPS styles, **units + ½ prefactor** |
| `topology/typing.py` | `extract_types` / `expand_types` (per-instance ↔ typed) |
| `topology/reduce.py` | `harmonize` (cluster-mean per type-tuple, bimodal) + `average_by_type` — legacy reduction parity |
| `topology/forcefield.py` | `ForceField` + `from_itp`/`from_json` + `apply` — parametrize a bare structure from GROMACS `ffbonded/ffnonbonded.itp` **or** `ffparams/*.json` (interchangeable) |
| `topology/validate.py` | semantic validation + optional JSON-Schema |
| `topology/gmx_preprocessor.py` | minimal GROMACS `#include`/`#ifdef`/`#define` |
| `topology/xml_io.py` | optional type-tagged XML |
| `topology/elements.py` | element → mass / Z |
| `import_topology.py` | `read_itp/json/data/psf/xml`, `from_atoms_box` |
| `write_topology.py` | `write_itp/json/data/psf/xml`, `to_atoms_box` |

## Status
**Done + tested** (`tests/test_topology_interchange.py`, 9 tests): model, units +
prefactor, functional forms, typing, validation, JSON & XML round-trip, GROMACS
`.itp` read/write (inline **and** `[bondtypes]`/`[angletypes]` resolution),
LAMMPS `.data` read/write with per-term unique typing (the §16 inversion of the
legacy average-by-type writer), CHARMM `.psf`+`.prm` write / `.psf` read,
`from_atoms_box`/`to_atoms_box`, the preprocessor on real `#ifdef` libraries, and
**the critical per-site-uniqueness test** (two same-type bonds with different
params survive to `.itp` inline and `.data` as two types).

**Legacy parity** (`reduce.py`, benchmarked against `write_top.py` on Kaolinite
×[2,2,1], 504 angles, in `tests/benchmark_legacy_vs_topology.py`):
- default = fully **explicit** per-instance (103 distinct θ0) — matches legacy
  `explicit_angles=1` exactly, and extends explicit terms to LAMMPS `.data`.
- `reduce.harmonize(detect_bimodal=True)` → **103 → 12** distinct θ0 — matches
  legacy `harmonize_angles` exactly.
- `reduce.average_by_type()` → **103 → 10** angle types — matches the legacy
  `lmp()` type-averaged default exactly.

**Force-field sources** (`forcefield.py`): params come from GROMACS
`ffnonbonded.itp`+`ffbonded.itp` (`[atomtypes]`/`[bondtypes]`/`[angletypes]`, via
the preprocessor) **or** the bundled `ffparams/*.json` — the two are
interchangeable (a test parametrizes a bare structure from each and asserts the
bond/angle/LJ values match exactly). Usage:
```python
ff = read_forcefield_itp(["ffnonbonded_gminff.itp", "ffbonded.itp"], defines=["GMINFF_k500"])
# or: ff = read_forcefield_json("GMINFF/gminff_all.json", variant="GMINFF_k500")
apply_forcefield(structure_topology, ff)   # fills AtomType LJ/mass/charge + per-instance bond/angle params
```

**Known gaps / next steps**
- `[dihedraltypes]` resolution from external/unshipped `#include`d FF files
  (inline dihedral params already work; in-file `[bondtypes]`/`[angletypes]` do).
- RB ↔ Fourier ↔ periodic dihedral conversions (energy units convert; the
  algebraic re-expression is a TODO with a lossiness warning).
- CHARMM `.rtf` generation (deferred per §14.6 — `.psf`+`.prm` only).
- Energy-validation gold-standard test (§10.5) needs GROMACS + LAMMPS binaries;
  structural round-trips are in place as the proxy.
- Legacy `write_top`/`import_top` are left intact (not converted to deprecation
  shims) to avoid breaking current callers; migrate incrementally via
  `from_atoms_box`.
