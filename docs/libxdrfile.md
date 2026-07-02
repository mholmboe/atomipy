# Reading GROMACS `.xtc` / `.trr` trajectories

`ap.import_traj()` can read GROMACS binary trajectories (`.xtc`, `.trr`) in addition to
`.pdb`/`.gro`. Because those formats store **only coordinates + box** (no atom names or
types), pass a companion structure with `top=` so the frames carry per-atom metadata:

```python
frames = ap.import_traj("traj.xtc", top="conf.gro")   # or conf.pdb / .cif
D = ap.msd(frames, atom_types=["OW"], dt=2.0)          # now typed analysis works
```

`stride`, `start`, `stop` subsample frames at read time:

```python
frames = ap.import_traj("traj.xtc", top="conf.gro", start=100, stride=10)
```

## How atomipy reads them

Two backends, tried in order:

1. **`libxdrfile`** (optional, fast, no subprocess) via a small `ctypes` wrapper
   (`atomipy/xdrfile.py`). atomipy compiles no C of its own — it just loads the prebuilt
   `libxdrfile` shared library at runtime.
2. **`gmx trjconv`** fallback — if `libxdrfile` isn't found but a GROMACS install is
   available and you passed `top=`, atomipy converts the trajectory to a temporary
   multi-frame PDB and parses that.

If neither is available, `import_traj` raises a clear error.

`.xtc`/`.trr` store nm; atomipy converts to **Ångström** automatically.

## Getting `libxdrfile`

atomipy locates the library from, in order:

1. `$ATOMIPY_XDRFILE` — an explicit path to the shared library,
2. `ctypes.util.find_library("xdrfile")`,
3. common SONAMEs (`libxdrfile.so`, `libxdrfile.dylib`, `xdrfile.dll`).

Note: modern GROMACS (2020+) folds these routines into `libgromacs` and does **not**
install a standalone `libxdrfile`, so you typically build the small v1.1.4 library once:

```bash
# download xdrfile-1.1.4.tar.gz (GROMACS contrib), then:
tar xzf xdrfile-1.1.4.tar.gz && cd xdrfile-1.1.4

# option A — one-line shared build (Linux .so / macOS .dylib):
cc -shared -fPIC -O2 -Iinclude -o libxdrfile.so \
   src/xdrfile.c src/xdrfile_xtc.c src/xdrfile_trr.c      # use .dylib on macOS

# option B — autotools:
./configure --enable-shared && make

# then point atomipy at it:
export ATOMIPY_XDRFILE=/abs/path/to/libxdrfile.so
```

`libxdrfile` is LGPL-2.1; it is an **optional** runtime dependency, not bundled with
atomipy. If you already run GROMACS, the `gmx trjconv` fallback needs no extra build.
