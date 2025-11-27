"""
Run Bond Valence Sum (BVS) and Shannon radii analysis on a structure file.

Usage:
  python run_BVS_analysis.py input.(pdb|gro|xyz) [--params bvparm2020.cif] [--csv output.csv]
       [--elements Al Si O H] [--oxidations 3 4 -2 1]

Outputs:
  - Console summary of GII, GII_noH, formal charge, top deviations
  - Optional CSV with per-atom BVS, deltas, bond contributions, and Shannon radii
"""

import argparse
import sys
from typing import List

import atomipy as ap


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run BVS + Shannon radii analysis on a structure file.")
    parser.add_argument("structure", help="Input structure file (.pdb, .gro, or .xyz)")
    parser.add_argument(
        "--params",
        default=None,
        help="Path to bond valence parameter table (default: bundled bvparm2020.cif).",
    )
    parser.add_argument(
        "--csv",
        dest="csv_out",
        default=None,
        help="Optional path to write per-atom results as CSV.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of worst offenders (by |delta|) to show.",
    )
    parser.add_argument(
        "--elements",
        nargs="+",
        help="Optional list of elements for oxidation hints, e.g., --elements Al Si O H",
    )
    parser.add_argument(
        "--oxidations",
        nargs="+",
        help="Optional list of oxidation states matching --elements, e.g., --oxidations 3 4 -2 1",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    elements: List[str] | None = args.elements
    ox_values = None
    if args.oxidations:
        ox_values = [int(float(x)) for x in args.oxidations]
    if elements and ox_values and len(elements) != len(ox_values):
        print("elements and oxidations must have the same length", file=sys.stderr)
        return 1

    try:
        analysis = ap.conf2bvs(
            args.structure,
            params_path=args.params,
            csv_path=args.csv_out,
            top_n=args.top_n,
            elements=elements,
            ox_values=ox_values,
        )
        results = analysis["results"]
        gii = analysis["gii"]
        gii_check = ap.global_instability_index(results)
        formal_charge = analysis["formal_charge"]
        gii_no_h = analysis["gii_no_h"]
        worst = analysis["top_worst"]
    except Exception as exc:  # pylint: disable=broad-except
        print(f"Error during BVS analysis: {exc}", file=sys.stderr)
        return 1

    print(f"BVS analysis for {args.structure}")
    print(f"Atoms: {len(results)}")
    print(f"Formal charge from assumed oxidation states: {formal_charge}")
    print(f"GII: {gii:.5f} (recomputed: {gii_check:.5f})")
    print(f"GII_noH (excluding H atoms): {gii_no_h:.5f}")
    print("Delta = actual BVS - |expected oxidation| (positive means overbonded)")

    if worst:
        print(f"Top {len(worst)} |delta| atoms:")
        for entry in worst:
            delta = entry.get("delta")
            sh_ionic = entry.get("shannon_ionic_radius")
            sh_crystal = entry.get("shannon_crystal_radius")
            cn = entry.get("cn")
            print(
                f"  #{entry.get('index'):>4} {entry.get('element') or 'X':>2} "
                f"BVS={entry.get('bvs', 0.0):6.3f} "
                f"exp={entry.get('expected_ox')} delta={delta:+6.3f} "
                f"cn={cn} sh_r_ionic={sh_ionic} sh_r_crystal={sh_crystal}"
            )

    if args.csv_out:
        print(f"Wrote CSV: {args.csv_out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
