#!/usr/bin/env python3
"""
analyze_uc_conf_angles.py
=========================
Scan all .pdb files in atomipy/structures/minerals/UC_conf/ and report which
structures have bimodal angle type triplets — flagging where a single mean
angle value per triplet is insufficient and a bimodal treatment is warranted.

Usage:
    python scripts/analyze_uc_conf_angles.py

Output:
    Prints a Markdown table to stdout and saves to uc_conf_angle_analysis.md
"""

import sys
import os
from pathlib import Path
from collections import defaultdict
import numpy as np

# ---------------------------------------------------------------------------
# Path setup — allow running from repo root OR from scripts/ sub-dir
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import atomipy as ap
from atomipy.write_top import cluster_angles


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
UC_CONF_DIR = REPO_ROOT / "atomipy" / "structures" / "minerals" / "UC_conf"
BIMODAL_THRESHOLD = 30.0   # degrees — spread > this AND max-gap > half of it
MIN_COUNT_FOR_BIMODAL = 4  # need at least 4 values to call bimodal
OUTPUT_MD = REPO_ROOT / "uc_conf_angle_analysis.md"


# ---------------------------------------------------------------------------
# Helper: detect bimodal triplets from angle_type_triplets dict
# ---------------------------------------------------------------------------
def analyse_triplets(angle_type_triplets, threshold=BIMODAL_THRESHOLD):
    """Return (bimodal_list, unimodal_list) for a structure's triplets.

    Each element is a dict with keys:
        triplet, count, mean, std, min, max, spread, clusters
    """
    bimodal = []
    unimodal = []

    for triplet, values in sorted(angle_type_triplets.items()):
        n = len(values)
        if n == 0:
            continue

        arr = np.array(values, dtype=float)
        mean = float(np.mean(arr))
        std  = float(np.std(arr))
        vmin = float(arr.min())
        vmax = float(arr.max())
        spread = vmax - vmin

        clusters = cluster_angles(values, threshold=threshold)
        is_bimodal = (len(clusters) > 1) and (n >= MIN_COUNT_FOR_BIMODAL)

        entry = dict(
            triplet="-".join(triplet),
            count=n,
            mean=mean,
            std=std,
            vmin=vmin,
            vmax=vmax,
            spread=spread,
            clusters=clusters,
        )

        if is_bimodal:
            bimodal.append(entry)
        else:
            unimodal.append(entry)

    return bimodal, unimodal


# ---------------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------------
def main():
    pdb_files = sorted(UC_CONF_DIR.glob("*.pdb"))
    if not pdb_files:
        print(f"No .pdb files found in {UC_CONF_DIR}")
        sys.exit(1)

    print(f"Found {len(pdb_files)} .pdb files in {UC_CONF_DIR.relative_to(REPO_ROOT)}")
    print()

    rows = []  # one row per structure

    for pdb_path in pdb_files:
        name = pdb_path.stem
        print(f"  Processing: {name} ...", end=" ", flush=True)

        try:
            atoms, Box = ap.import_pdb(str(pdb_path))
            if not atoms or Box is None:
                print("SKIP (empty)")
                continue

            # Assign MINFF types (also calculates bonds/angles internally)
            atoms = ap.minff(atoms, Box, log=False)

            # get_structure_stats builds angle_type_triplets from atom['angles']
            # We need those triplets — extract them without calling get_structure_stats
            # (which also recalculates bonds redundantly). Instead, use the data
            # already on atoms after minff().
            angle_type_triplets = defaultdict(list)
            for i, atom in enumerate(atoms):
                if 'angles' in atom:
                    for angle_data in atom['angles']:
                        (n1_idx, n2_idx), angle_val = angle_data
                        t_center = atom.get('fftype', atom.get('type', 'X'))
                        t1 = atoms[n1_idx].get('fftype', atoms[n1_idx].get('type', 'X'))
                        t2 = atoms[n2_idx].get('fftype', atoms[n2_idx].get('type', 'X'))
                        # canonical ordering: min(t1,t2) as first, max as third
                        if t1 <= t2:
                            triplet = (t1, t_center, t2)
                        else:
                            triplet = (t2, t_center, t1)
                        angle_type_triplets[triplet].append(angle_val)

            bimodal, unimodal = analyse_triplets(angle_type_triplets)

            n_triplets = len(bimodal) + len(unimodal)
            verdict = "⚠️  BIMODAL" if bimodal else "✅ single-mean OK"
            print(f"{verdict}  ({n_triplets} triplet{'s' if n_triplets != 1 else ''})")

            rows.append(dict(
                name=name,
                n_triplets=n_triplets,
                bimodal=bimodal,
                unimodal=unimodal,
                verdict="BIMODAL" if bimodal else "OK",
            ))

        except Exception as e:
            print(f"ERROR: {e}")
            rows.append(dict(name=name, n_triplets=0, bimodal=[], unimodal=[], verdict="ERROR"))

    # -----------------------------------------------------------------------
    # Build Markdown output
    # -----------------------------------------------------------------------
    lines = []
    lines.append("# UC_conf Mineral Structure — Angle Triplet Bimodal Analysis")
    lines.append("")
    lines.append(
        "Threshold: spread > **{:.0f}°** AND max gap > **{:.0f}°** → flagged as bimodal.  ".format(
            BIMODAL_THRESHOLD, BIMODAL_THRESHOLD / 2
        )
    )
    lines.append(
        "Structures marked **⚠️ BIMODAL** should use `detect_bimodal=True` (LAMMPS/GROMACS) "
        "or two separate angle type entries. "
        "Structures marked **✅ OK** can safely use a single mean angle per triplet."
    )
    lines.append("")

    # Summary table
    lines.append("## Summary Table")
    lines.append("")
    lines.append("| Structure | Triplets | Bimodal triplets | Recommendation |")
    lines.append("|-----------|:--------:|:----------------:|----------------|")
    for r in rows:
        nm = r["name"]
        nb = len(r["bimodal"])
        nt = r["n_triplets"]
        if r["verdict"] == "ERROR":
            lines.append(f"| {nm} | — | — | ❌ Error during processing |")
        elif nb == 0:
            lines.append(f"| {nm} | {nt} | 0 | ✅ Single mean angle per triplet is fine |")
        else:
            lines.append(f"| {nm} | {nt} | **{nb}** | ⚠️ Use `detect_bimodal=True` |")
    lines.append("")

    # Detailed bimodal breakdown per structure
    lines.append("## Bimodal Triplet Details")
    lines.append("")
    any_bimodal = False
    for r in rows:
        if not r["bimodal"]:
            continue
        any_bimodal = True
        lines.append(f"### {r['name']}")
        lines.append("")
        lines.append("| Triplet (type1-center-type3) | n | mean (°) | std (°) | spread (°) | Cluster 1 (°) | Cluster 2 (°) |")
        lines.append("|------------------------------|:-:|:--------:|:-------:|:----------:|:-------------:|:-------------:|")
        for b in r["bimodal"]:
            c1_mean = f"{b['clusters'][0][0]:.1f}" if len(b['clusters']) > 0 else "—"
            c2_mean = f"{b['clusters'][1][0]:.1f}" if len(b['clusters']) > 1 else "—"
            lines.append(
                f"| `{b['triplet']}` | {b['count']} | {b['mean']:.1f} | {b['std']:.1f} | {b['spread']:.1f} | {c1_mean} | {c2_mean} |"
            )
        lines.append("")

    if not any_bimodal:
        lines.append("_No bimodal angle distributions detected in any UC_conf structure._")
        lines.append("")

    # Structures that are clean
    clean = [r for r in rows if r["verdict"] == "OK"]
    bimodal_structs = [r for r in rows if r["verdict"] == "BIMODAL"]
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- **{len(clean)}** structures → single mean angle OK")
    lines.append(f"- **{len(bimodal_structs)}** structures → bimodal treatment recommended")
    if bimodal_structs:
        lines.append("")
        lines.append("Structures requiring bimodal treatment:")
        for r in bimodal_structs:
            triplets_str = ", ".join(f"`{b['triplet']}`" for b in r["bimodal"])
            lines.append(f"  - **{r['name']}** — {triplets_str}")
    lines.append("")

    md_text = "\n".join(lines)

    # Save
    with open(OUTPUT_MD, "w", encoding="utf-8") as fh:
        fh.write(md_text)
    print()
    print(f"Saved Markdown report → {OUTPUT_MD.relative_to(REPO_ROOT)}")
    print()
    print(md_text)


if __name__ == "__main__":
    main()
