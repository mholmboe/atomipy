#!/usr/bin/env python3
"""MINFF to n2t generator."""
# Example: python minff2n2t.py structure.gro --output minff_atomname2type.n2t

import argparse
from pathlib import Path

import atomipy as ap


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate MINFF n2t file from structure")
    parser.add_argument("input_file", help="Input structure file (gro, pdb, xyz)")
    parser.add_argument("--output", help="Output n2t file name")
    args = parser.parse_args()

    infile = Path(args.input_file)
    if not infile.is_file():
        raise SystemExit(f"Input file '{infile}' not found")

    atoms, Box = ap.import_auto(str(infile))
    atoms = ap.minff(atoms, Box)

    output = Path(args.output) if args.output else infile.with_name(f"{infile.stem}_minff_atomname2type.n2t")
    n2t_path = ap.write_n2t(atoms, Box=Box, n2t_file=str(output))
    print(f"Written {n2t_path}")


if __name__ == "__main__":
    main()
