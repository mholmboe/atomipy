#!/usr/bin/env python3
"""Structure → n2t generator for gmx x2top workflows.

This CLI loads a structure (PDB/GRO/XYZ), preserves the box/cell, and writes
an atom-name-to-type table via `atomipy.write_n2t`. Environments are merged so
nearly identical sites share the same type definitions.

Usage examples (run from repo root):
  python struct2n2t.py Kaolinite_GII_0.0487.pdb
  python struct2n2t.py conf/preem1.gro --output preem1_atomname2type.n2t
  python struct2n2t.py unitcell.xyz --output unitcell_atomname2type.n2t

What it does:
  - Autodetects format with `import_auto` and keeps any box/cell info.
  - Hands atoms and box to `write_n2t`, which infers bonds if needed.
  - Writes <input>_atomname2type.n2t unless --output is provided.
"""

import argparse
import os
import sys

import atomipy as ap


def main():
    parser = argparse.ArgumentParser(description="Generate a GROMACS .n2t file from a structure")
    parser.add_argument("structure", help="Input structure (.pdb, .gro, .xyz, etc.)")
    parser.add_argument("--output", "-o", help="Output .n2t path (default: <structure>_atomname2type.n2t)")
    args = parser.parse_args()

    if not os.path.isfile(args.structure):
        parser.error(f"Input file '{args.structure}' does not exist")

    atoms, Box = ap.import_auto(args.structure)

    # write_n2t internally ensures elements/masses and will infer bonds if required
    out_path = args.output
    if out_path is None:
        stem, _ = os.path.splitext(os.path.basename(args.structure))
        out_path = f"{stem}_atomname2type.n2t"

    try:
        n2t_path = ap.write_n2t(atoms, Box, n2t_file=out_path, verbose=True)
    except Exception as exc:
        print(f"Failed to generate .n2t: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
