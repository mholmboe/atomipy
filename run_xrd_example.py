#!/usr/bin/env python3
"""Quick XRD runner using atomipy.diffraction.xrd.

This test script loads a structure (PDB/GRO/XYZ), prepares it for diffraction,
and plots the calculated XRD profile. Optional flags let you tweak the scan
range, wavelength, and whether to save the generated data files.

Usage examples (run from repo root):
  python run_xrd_example.py Kaolinite_GII_0.0487.pdb
  python run_xrd_example.py Kaolinite_GII_0.0487.pdb --two-theta 5 70 --angle-step 0.01
  python run_xrd_example.py Kaolinite_GII_0.0487.gro --save-output
  python run_xrd_example.py Kaolinite_GII_0.0487.xyz --box 52 52 60 --wavelength 0.71073 --neutral --no-plot

Tips:
  - Provide --box when the file lacks box/cell info (3, 6, or 9 numbers).
  - Use --save-output to emit xrd.dat/mat/txt/png alongside the plot.
  - Use --no-plot for headless or batch runs; data still computes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

import atomipy as ap


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an XRD pattern for a structure file using atomipy.diffraction.xrd.",
    )
    parser.add_argument("structure", help="Input structure file (.pdb, .gro, .xyz)")
    parser.add_argument(
        "--wavelength",
        type=float,
        default=1.54187,
        help="X-ray wavelength in Å (default: Cu Kα)",
    )
    parser.add_argument(
        "--two-theta",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(2.0, 90.0),
        help="2θ scan range in degrees (min max).",
    )
    parser.add_argument(
        "--angle-step",
        type=float,
        default=0.02,
        help="2θ step size in degrees.",
    )
    parser.add_argument(
        "--hkl-max",
        type=int,
        default=0,
        help="Optional override for max h,k,l indices (0 = auto from 2θ range).",
    )
    parser.add_argument(
        "--neutral",
        action="store_true",
        help="Force neutral atom scattering factors (skip ionic defaults).",
    )
    parser.add_argument(
        "--save-output",
        action="store_true",
        help="Write xrd.dat/mat/txt/png outputs in the working directory.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Compute pattern without showing the matplotlib window.",
    )
    parser.add_argument(
        "--box",
        nargs="+",
        type=float,
        help="Override box/cell values if missing in the structure (3, 6, or 9 numbers).",
    )
    return parser.parse_args()


def validate_box(box_values: Sequence[float] | None) -> List[float] | None:
    if box_values is None:
        return None
    if len(box_values) not in (3, 6, 9):
        raise ValueError("Box must contain 3, 6, or 9 values.")
    return [float(v) for v in box_values]


def prepare_atoms_for_xrd(atoms: Iterable[dict]) -> List[dict]:
    """Copy atoms and force the 'type' field to the element symbol for scattering lookups."""
    prepared: List[dict] = []
    for atom in atoms:
        atom_copy = atom.copy()
        element = atom_copy.get("element")
        if element:
            atom_copy["type"] = element
        prepared.append(atom_copy)
    return prepared


def main() -> int:
    args = parse_args()

    struct_path = Path(args.structure)
    if not struct_path.exists():
        print(f"Input structure not found: {struct_path}", file=sys.stderr)
        return 1

    try:
        box = validate_box(args.box)
    except ValueError as exc:
        print(f"Invalid box override: {exc}", file=sys.stderr)
        return 1

    try:
        atoms, file_box = ap.import_auto(str(struct_path))
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Failed to read structure: {exc}", file=sys.stderr)
        return 1

    file_box_validated: List[float] | None = None
    if file_box is not None:
        try:
            file_box_validated = validate_box(file_box)
        except ValueError as exc:
            print(f"Invalid box information in file: {exc}", file=sys.stderr)
            return 1

    # Prefer CLI-provided box/cell; otherwise use what the file contained.
    box = box if box is not None else file_box_validated
    if box is None:
        print("No box/cell information found. Provide --box with 3, 6, or 9 numbers.", file=sys.stderr)
        return 1

    if args.two_theta[0] >= args.two_theta[1]:
        print("--two-theta must be provided as 'min max' with min < max.", file=sys.stderr)
        return 1

    atoms = ap.element(atoms)  # Ensure elements are set
    atoms_for_xrd = prepare_atoms_for_xrd(atoms)

    print(f"Running XRD for {struct_path.name} ({len(atoms_for_xrd)} atoms)")
    print(f"2θ range: {args.two_theta[0]}-{args.two_theta[1]} deg, step {args.angle_step} deg")
    print(f"Wavelength: {args.wavelength} Å")

    try:
        result = ap.xrd(
            atoms=atoms_for_xrd,
            Box=box,
            wavelength=args.wavelength,
            angle_step=args.angle_step,
            two_theta_range=tuple(args.two_theta),
            neutral_atoms=args.neutral,
            hkl_max=args.hkl_max,
            plot=not args.no_plot,
            save_output=args.save_output,
        )
    except Exception as exc:  # pylint: disable=broad-except
        print(f"XRD calculation failed: {exc}", file=sys.stderr)
        return 1

    if isinstance(result, tuple) and len(result) == 3:
        two_theta, intensity, _ = result
    else:
        two_theta, intensity = result

    print(f"Computed {len(two_theta)} intensity points.")
    if args.save_output:
        print("Saved outputs: xrd.dat, xrd_results.mat, hkl_limits.txt, atomic_scattering_factors.txt, xrd_pattern.png")
    if args.no_plot:
        print("Plot display disabled (--no-plot).")

    return 0


if __name__ == "__main__":
    sys.exit(main())
