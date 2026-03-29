#!/usr/bin/env python3
"""
Convert CIF/mmCIF to PDB with optional preprocessing and typing.

Workflow
--------
1. Import CIF/mmCIF (symmetry expansion enabled by default).
2. Optionally unfold layered structures in fractional z.
3. Fuse overlapping atom sites.
4. Optionally run BVS-based protonation-need diagnostics.
5. Optionally assign MINFF atom types and charges.
6. Write PDB.
"""

from __future__ import annotations

import argparse
import copy
import os
import sys
from typing import Iterable, Tuple

import numpy as np

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import atomipy as ap


DEFAULT_FUSE_RMAX = 0.85
DEFAULT_Z_SPLIT = 0.5
DEFAULT_AUTO_GAP_THRESHOLD = 0.25
DEFAULT_PROTONATION_DELTA_THRESHOLD = -0.50
DEFAULT_PROTONATION_TOP_N = 10


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Import CIF, expand full unit cell, fuse overlaps, optional BVS protonation check, "
            "optional MINFF typing, and write PDB."
        )
    )
    parser.add_argument("input_cif", help="Input .cif/.mmcif file.")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "Output PDB path (default: <input_stem>.pdb, or <input_stem>_minff.pdb "
            "when --assign-minff is enabled)."
        ),
    )
    parser.add_argument(
        "--no-expand-symmetry",
        action="store_true",
        help="Disable CIF symmetry expansion (default is expansion to full unit cell).",
    )
    parser.add_argument(
        "--fuse-rmax",
        type=float,
        default=DEFAULT_FUSE_RMAX,
        help=f"Overlap fusion cutoff in Angstrom (default: {DEFAULT_FUSE_RMAX}).",
    )
    parser.add_argument(
        "--fuse-criteria",
        choices=["average", "occupancy", "order"],
        default="average",
        help="How to set fused-site coordinates/properties (default: average).",
    )
    parser.add_argument(
        "--layer-z-unfold",
        action="store_true",
        help=(
            "Manual z-unfold for layered structures: shift atoms with zfrac >= z_split by -1.0 "
            "before Cartesian conversion."
        ),
    )
    parser.add_argument(
        "--no-auto-layer-z-unfold",
        action="store_true",
        help="Disable automatic layered z-unfold detection.",
    )
    parser.add_argument(
        "--z-split",
        type=float,
        default=DEFAULT_Z_SPLIT,
        help=f"Fractional z threshold for manual --layer-z-unfold (default: {DEFAULT_Z_SPLIT}).",
    )
    parser.add_argument(
        "--auto-gap-threshold",
        type=float,
        default=DEFAULT_AUTO_GAP_THRESHOLD,
        help=(
            "Minimum largest cyclic gap in fractional z required to auto-unfold "
            f"(default: {DEFAULT_AUTO_GAP_THRESHOLD})."
        ),
    )
    parser.add_argument(
        "--assign-minff",
        action="store_true",
        help="Assign MINFF atom types/charges. If omitted, keep CIF atom names/types.",
    )
    parser.add_argument(
        "--minff-log",
        default=None,
        help="Optional MINFF stats log file path (used only with --assign-minff).",
    )
    parser.add_argument(
        "--check-protonation",
        action="store_true",
        help="Run BVS analysis and report whether protonation is likely needed.",
    )
    parser.add_argument(
        "--protonation-delta-threshold",
        type=float,
        default=DEFAULT_PROTONATION_DELTA_THRESHOLD,
        help=(
            "BVS delta threshold for underbonded O protonation candidates "
            f"(default: {DEFAULT_PROTONATION_DELTA_THRESHOLD})."
        ),
    )
    parser.add_argument(
        "--protonation-top-n",
        type=int,
        default=DEFAULT_PROTONATION_TOP_N,
        help=(
            "Number of most underbonded oxygen sites to print when --check-protonation is enabled "
            f"(default: {DEFAULT_PROTONATION_TOP_N})."
        ),
    )
    return parser.parse_args()


def build_output_path(input_cif: str, output_path: str | None, assign_minff: bool) -> str:
    """Resolve output path from CLI args."""
    if output_path:
        return output_path
    stem, _ = os.path.splitext(input_cif)
    return f"{stem}_minff.pdb" if assign_minff else f"{stem}.pdb"


def unfold_layered_z_fractional(atoms: list[dict], split: float) -> int:
    """
    Shift atoms with zfrac >= split by -1.0 in fractional coordinates.

    Parameters
    ----------
    atoms : list of dict
        Atom list containing fractional z values in ``zfrac``.
    split : float
        Fractional threshold in [0, 1) used to decide which atoms are shifted.

    Returns
    -------
    int
        Number of shifted atoms.
    """
    moved = 0
    for atom in atoms:
        zfrac = atom.get("zfrac")
        if zfrac is None:
            continue
        zfrac_wrapped = zfrac % 1.0
        if zfrac_wrapped >= split:
            atom["zfrac"] = zfrac_wrapped - 1.0
            moved += 1
        else:
            atom["zfrac"] = zfrac_wrapped
    return moved


def detect_layered_z_unfold_split(
    atoms: Iterable[dict],
    min_gap: float,
) -> Tuple[float | None, float, int, str]:
    """
    Auto-detect whether z-unfold is likely needed using fractional-z distribution.

    Strategy
    --------
    1. Wrap zfrac values into [0, 1), sort, and compute cyclic gaps.
    2. Find the largest gap.
    3. If largest gap >= min_gap, set split at that gap midpoint.
    4. Shift atoms with zfrac >= split.

    Returns
    -------
    tuple
        (split, largest_gap, n_shift, reason)
        reason is one of:
        - ``ok``
        - ``insufficient_data``
        - ``below_gap_threshold``
        - ``degenerate_split``
    """
    z_values = [atom.get("zfrac") for atom in atoms if atom.get("zfrac") is not None]
    if len(z_values) < 2:
        return None, 0.0, 0, "insufficient_data"

    z_values = np.array(z_values, dtype=float) % 1.0
    z_values.sort()

    next_values = np.roll(z_values, -1)
    gaps = (next_values - z_values) % 1.0
    gap_index = int(np.argmax(gaps))
    largest_gap = float(gaps[gap_index])

    if largest_gap < min_gap:
        return None, largest_gap, 0, "below_gap_threshold"

    gap_start = float(z_values[gap_index])
    split = (gap_start + 0.5 * largest_gap) % 1.0
    n_shift = int(np.sum(z_values >= split))

    if n_shift <= 0 or n_shift >= len(z_values):
        return None, largest_gap, n_shift, "degenerate_split"

    return split, largest_gap, n_shift, "ok"


def maybe_unfold_layered_z(
    atoms: list[dict],
    Cell: list[float],
    args: argparse.Namespace,
) -> list[dict]:
    """
    Apply manual or auto z-unfold logic and rebuild Cartesian coordinates when needed.
    """
    did_unfold = False

    if args.layer_z_unfold:
        moved = unfold_layered_z_fractional(atoms, split=args.z_split)
        print(
            f"Manual layer z-unfold: shifted {moved} atoms with zfrac >= {args.z_split:.3f} by -1.0"
        )
        did_unfold = True
    elif not args.no_auto_layer_z_unfold:
        split, largest_gap, n_shift, reason = detect_layered_z_unfold_split(
            atoms,
            min_gap=args.auto_gap_threshold,
        )
        if split is not None:
            moved = unfold_layered_z_fractional(atoms, split=split)
            print(
                f"Auto layer z-unfold: largest_gap={largest_gap:.3f} >= "
                f"{args.auto_gap_threshold:.3f}, split={split:.3f}, shifted={moved} atoms"
            )
            did_unfold = True
        else:
            reason_map = {
                "insufficient_data": "insufficient fractional-z data",
                "below_gap_threshold": (
                    f"largest_gap={largest_gap:.3f} < {args.auto_gap_threshold:.3f}"
                ),
                "degenerate_split": (
                    f"degenerate split (largest_gap={largest_gap:.3f}, would shift {n_shift} atoms)"
                ),
            }
            print(f"Auto layer z-unfold: not applied ({reason_map.get(reason, reason)})")

    if did_unfold:
        ap.fractional_to_cartesian(atoms=atoms, Box=Cell, add_to_atoms=True)

    return atoms


def run_bvs_protonation_check(
    atoms: list[dict],
    Box: list[float],
    args: argparse.Namespace,
) -> None:
    """
    Run BVS analysis and print protonation diagnostics.

    Notes
    -----
    ``analyze_bvs`` calls ``element()``, which can overwrite atom ``type`` fields.
    To preserve CIF atom names when MINFF is disabled, this analysis runs on a deep copy.
    """
    bvs_atoms = copy.deepcopy(atoms)
    print("Running BVS analysis for protonation check ...")
    bvs_report = ap.analyze_bvs(bvs_atoms, Box, top_n=args.protonation_top_n)

    has_h = any((atom.get("element") or "").upper() == "H" for atom in atoms)
    underbonded_oxygen = [
        row
        for row in bvs_report["results"]
        if (row.get("element") or "").upper() == "O"
        and row.get("delta") is not None
        and row["delta"] < args.protonation_delta_threshold
    ]

    print(
        f"BVS protonation check: GII={bvs_report['gii']:.5f}, "
        f"GII_noH={bvs_report['gii_no_h']:.5f}, "
        f"underbonded_O={len(underbonded_oxygen)} "
        f"(delta < {args.protonation_delta_threshold:.3f}), has_H={has_h}"
    )

    if underbonded_oxygen and not has_h:
        print("Protonation likely needed: yes (underbonded O found and no H present).")
    elif underbonded_oxygen and has_h:
        print("Protonation may still be needed: underbonded O found despite existing H.")
    else:
        print("Protonation likely needed: no underbonded O below threshold.")

    if underbonded_oxygen:
        top_n = min(len(underbonded_oxygen), args.protonation_top_n)
        print(f"Top {top_n} underbonded O sites:")
        ranked = sorted(underbonded_oxygen, key=lambda row: row["delta"])[:top_n]
        for row in ranked:
            print(
                f"  #{row.get('index'):>4} O "
                f"BVS={row.get('bvs', 0.0):6.3f} "
                f"delta={row.get('delta', 0.0):+6.3f} "
                f"cn={row.get('cn')}"
            )


def main() -> int:
    """CLI entry point."""
    args = parse_args()

    if not os.path.exists(args.input_cif):
        print(f"Input file not found: {args.input_cif}", file=sys.stderr)
        return 1

    output_path = build_output_path(args.input_cif, args.output, args.assign_minff)

    print(f"Loading {args.input_cif} ...")
    atoms, Cell = ap.import_cif(
        args.input_cif,
        expand_symmetry=not args.no_expand_symmetry,
    )
    Box = ap.Cell2Box_dim(Cell)
    print(f"Imported {len(atoms)} atoms")

    atoms = maybe_unfold_layered_z(atoms, Cell, args)

    num_atoms_before_fuse = len(atoms)
    atoms = ap.fuse_atoms(atoms, Box, rmax=args.fuse_rmax, criteria=args.fuse_criteria)
    num_atoms_after_fuse = len(atoms)
    print(f"Fuse step: {num_atoms_before_fuse} -> {num_atoms_after_fuse} atoms")

    if args.check_protonation:
        run_bvs_protonation_check(atoms, Box, args)

    if args.assign_minff:
        print("Assigning MINFF atom types and charges ...")
        atoms = ap.minff(
            atoms,
            Box,
            log=bool(args.minff_log),
            log_file=args.minff_log,
        )
    elif args.minff_log:
        print("Warning: --minff-log ignored because --assign-minff was not set.")

    print(f"Writing PDB: {output_path}")
    ap.write_pdb(atoms, Box, output_path)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
