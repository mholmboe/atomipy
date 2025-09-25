#!/usr/bin/env python3
"""Minimal MINFF n2t generator."""
# Example: python generate_n2t_example.py structure.gro --output structure_minff.n2t

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

    atoms, _, box = ap.import_auto(str(infile))
    atoms, _ = ap.minff(atoms, box)

    output = Path(args.output) if args.output else infile.with_name(f"{infile.stem}_minff.n2t")
    n2t_path = ap.write_n2t(atoms, n2t_file=str(output), box=box)
    print(f"Written {n2t_path}")


if __name__ == "__main__":
    main()
