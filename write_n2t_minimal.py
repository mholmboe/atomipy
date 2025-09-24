#!/usr/bin/env python3
"""Minimal CLI for converting a structure (.pdb/.gro) into an .n2t file.

The script loads coordinates, forwards the simulation cell so that
`atomipy.write_n2t` can honour periodic boundary conditions, and writes
the consolidated environment table to disk.
"""

import argparse
import os
import sys

import atomipy as ap


def main():
    parser = argparse.ArgumentParser(description="Generate a GROMACS .n2t file from a structure")
    parser.add_argument("structure", help="Input structure (.pdb, .gro, .xyz, etc.)")
    parser.add_argument("--output", "-o", help="Output .n2t path (default: <structure>.n2t)")
    args = parser.parse_args()

    if not os.path.isfile(args.structure):
        parser.error(f"Input file '{args.structure}' does not exist")

    atoms, cell, box = ap.import_auto(args.structure)

    # write_n2t internally ensures elements/masses and will infer bonds if required
    out_path = args.output
    if out_path is None:
        stem, _ = os.path.splitext(os.path.basename(args.structure))
        out_path = f"{stem}.n2t"

    try:
        n2t_path = ap.write_n2t(atoms, n2t_file=out_path, box=box, verbose=True)
    except Exception as exc:
        print(f"Failed to generate .n2t: {exc}", file=sys.stderr)
        return 1

    print(f"N2T file written to {n2t_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
