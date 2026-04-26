# Contributing to atomipy

Thank you for contributing to **atomipy**! This guide describes how to collaborate safely and consistently on this scientific Python package.

---

## Table of Contents

- [Getting Started](#getting-started)
- [Branch Strategy](#branch-strategy)
- [Making Changes](#making-changes)
- [Commit Message Style](#commit-message-style)
- [Pull Requests](#pull-requests)
- [Running Tests](#running-tests)
- [Code Style](#code-style)
- [Versioning](#versioning)

---

## Getting Started

### 1. Clone and install in editable mode

```bash
git clone https://github.com/mholmboe/atomipy.git
cd atomipy
pip install -e .
# Optional extras:
pip install -e ".[cif]"    # CIF/mmCIF support via gemmi
pip install -e ".[xrd]"    # XRD plotting via matplotlib + scipy
```

### 2. Verify your installation

```bash
python -c "import atomipy as ap; print('atomipy', ap.__version__)"
```

---

## Branch Strategy

We use a simple **feature-branch workflow** with two long-lived branches:

| Branch | Purpose |
|---|---|
| `main` | Stable, tagged releases only. Never push directly here. |
| `dev` | Active integration branch. All feature PRs merge here. |
| `feature/<short-name>` | Your working branch for a specific change. |
| `fix/<short-name>` | Bug fix branches. |
| `docs/<short-name>` | Documentation-only changes. |

### Typical workflow

```bash
# 1. Always start from the latest dev
git checkout dev
git pull origin dev

# 2. Create your feature branch
git checkout -b feature/improve-bond-angle-pbc

# 3. Do your work, commit often (see commit style below)
git add atomipy/bond_angle.py
git commit -m "feat: improve PBC handling for triclinic cells in bond_angle"

# 4. Push and open a Pull Request → dev
git push origin feature/improve-bond-angle-pbc
```

---

## Making Changes

- **One concern per branch** — keep branches focused. A PR should do one thing.
- **Write or update tests** for any new or changed functionality (see `tests/`).
- **Update the docstring** of any function you modify.
- **Do not commit large binary files** (`.gro`, `.pdb`, `.data` files > 1 MB) unless they are essential reference structures.

---

## Commit Message Style

Use the **conventional commits** format:

```
<type>: <short description in present tense>
```

| Type | When to use |
|---|---|
| `feat` | New function or feature |
| `fix` | Bug fix |
| `docs` | Documentation only |
| `refactor` | Code restructure, no behaviour change |
| `test` | Adding or updating tests |
| `perf` | Performance improvement |
| `chore` | Build system, CI, dependencies |

**Examples:**

```
feat: add OPC3 water model support to solvate
fix: wrap coordinates for triclinic boxes in bond_angle
docs: add MINFF link to README and ffparams README
perf: use numpy broadcasting in cell_list_dist_matrix
test: add CI smoke test for minff on Kaolinite PDB
```

---

## Pull Requests

1. Open a PR from your branch **into `dev`**, not `main`.
2. Fill in the PR description: what changed and why.
3. Request a review from at least one other collaborator.
4. All CI checks must pass before merging.
5. Use **Squash and Merge** to keep the `dev` history clean.

### Releasing to `main`

Only the repo owner merges `dev → main`. A new release PR is opened when `dev` is stable, a version bump is applied, and a tag is created:

```bash
git tag -a v0.96 -m "Release v0.96"
git push origin v0.96
```

---

## Running Tests

The CI runs automatically on every PR. To run locally before pushing:

```bash
# Core smoke test (import + bond_angle dispatch)
python scripts/test_threshold_logic.py

# MINFF typing on reference structures
python -c "
import atomipy as ap
atoms, box = ap.import_pdb('Kaolinite_GII_0.0487.pdb')
atoms = ap.minff(atoms, Box=box)
print(f'Typed {len(atoms)} atoms OK')
"

# Import smoke test (all submodules)
python -c "import atomipy as ap; print('Import OK:', ap.__version__)"
```

---

## Code Style

- Follow the conventions in [`CODING_STYLE_GUIDE.md`](./CODING_STYLE_GUIDE.md).
- Use **4-space indentation**, no tabs.
- All public functions must have a **NumPy-style docstring** with `Args`, `Returns`, and at least one `Examples` entry.
- Keep line length ≤ 120 characters.
- Run a quick lint check before committing:

```bash
pip install flake8
flake8 atomipy/ --max-line-length=120 --ignore=E203,W503
```

---

## Versioning

atomipy uses **semantic versioning** (`MAJOR.MINOR`):

- `MINOR` bump: new functions, backward-compatible changes.
- `MAJOR` bump: breaking API changes (rare).

The version is set in `setup.py`. Update it only in the release PR to `main`.
