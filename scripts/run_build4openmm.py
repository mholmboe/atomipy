#!/usr/bin/env python3
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import atomipy as ap
import json
import math
import sys

# BVS-driven protonation controls
O_UNDERBONDED_DELTA = -0.5
BVS_GII_WARNING = 0.20

# Input parameters
input_cif = sys.argv[1] if len(sys.argv) > 1 else "Pyrophyllite_centered.cif"
output_name = "minff_system"
gii_output_file = f"{output_name}_gii.json"

# 1. Import CIF
print(f"Loading {input_cif}...")
atoms, cell = ap.import_cif(input_cif)
box = ap.Cell2Box_dim(cell)

# 2. BVS sanity check + protonation decision (avoid MINFF before H completion)
print("Running BVS sanity analysis...")
bvs_report = ap.analyze_bvs(atoms, box, top_n=10)
gii_before = bvs_report["gii"]
gii_no_h_before = bvs_report["gii_no_h"]
gii_after = None
gii_no_h_after = None
print(f"BVS GII: {gii_before:.4f} (no-H: {gii_no_h_before:.4f})")
if gii_before > BVS_GII_WARNING:
    print(f"Warning: GII {gii_before:.4f} > {BVS_GII_WARNING:.2f} (structure may be strained/disordered).")

# Protonate underbonded oxygen sites only when no H are present
has_h = any(a['element'] == 'H' for a in atoms)
if not has_h:
    underbonded_o = [
        row for row in bvs_report["results"]
        if row.get("element") == "O"
        and row.get("delta") is not None
        and row["delta"] < O_UNDERBONDED_DELTA
    ]

    if underbonded_o:
        print(f"Hydrogens missing. BVS found {len(underbonded_o)} underbonded O sites; protonating...")
        protonation_tag = "__BVS_O_PROTONATE__"
        original_types = [a.get("type") for a in atoms]

        # Tag only BVS-selected oxygen atoms as targets for add_H_atom
        for row in underbonded_o:
            atom_idx0 = row["index"] - 1  # BVS is 1-based
            atoms[atom_idx0]["type"] = protonation_tag

        n_before = len(atoms)
        atoms = ap.add_H_atom(
            atoms,
            box,
            target_type=protonation_tag,
            h_type='H',
            bond_length=1.0,
            coordination=3,
            max_h_per_atom=1
        )
        n_after = len(atoms)
        print(f"Added {n_after - n_before} H atoms from BVS-guided protonation.")

        # Restore original atom type labels for all pre-existing atoms
        for i, old_type in enumerate(original_types):
            atoms[i]["type"] = old_type

        print("Re-running BVS after protonation...")
        bvs_after = ap.analyze_bvs(atoms, box, top_n=10)
        gii_after = bvs_after["gii"]
        gii_no_h_after = bvs_after["gii_no_h"]
        print(f"BVS GII after protonation: {gii_after:.4f} (no-H: {gii_no_h_after:.4f})")
    else:
        print("Hydrogens missing, but BVS did not flag underbonded O sites. Skipping auto-protonation.")

# 3. Replication to > 2.5nm (X, Y) and > 3.0nm (Z)
# target_x = 25A, target_y = 25A, target_z = 30A
nx = math.ceil(25.0 / box[0])
ny = math.ceil(25.0 / box[1])
nz = 1 #math.ceil(30.0 / box[2])
print(f"Replicating {nx}x{ny}x{nz}...")
atoms, box, _ = ap.replicate_system(atoms, box, [nx, ny, nz])

# 4. Add vacuum space for solvation (e.g., 25A buffer in Z)
print("Extending box for solvation...")
box[2] += 25.0
# The solvate function will fill this new space automatically
limits = [0, 0, 0, box[0], box[1], box[2]]

# 5. Ionization (Add 10 Na+ and 10 Cl-)
print("Adding ions...")
na_ions = ap.ionize(ion_type='Na', resname='ION', limits=limits, num_ions=6, min_distance=3.0, solute_atoms=atoms)
atoms = ap.update(atoms, na_ions)
cl_ions = ap.ionize(ion_type='Cl', resname='ION', limits=limits, num_ions=6, min_distance=3.0, solute_atoms=atoms)
atoms = ap.update(atoms, cl_ions)

# 6. Solvation (Add TIP3P water to fill the box)
print("Solvating...")
# include_solute=True returns the merged system directly
atoms = ap.solvate(limits=limits, solute_atoms=atoms, solvent_type='tip3p', min_distance=2.5, include_solute=True)

# 7. Assign MINFF atom types and charges
print("Assigning MINFF parameters...")
atoms = ap.minff(atoms, box)

# 8. Write outputs for OpenMM
print(f"Writing {output_name}.pdb and {output_name}.psf...")
ap.write_pdb(atoms, box, f"{output_name}.pdb")
ap.write_psf(atoms, box, f"{output_name}.psf")

gii_payload = {
    "input_cif": input_cif,
    "gii": {
        "initial": gii_before,
        "initial_no_h": gii_no_h_before,
        "after_protonation": gii_after,
        "after_protonation_no_h": gii_no_h_after
    }
}
with open(gii_output_file, "w", encoding="utf-8") as fh:
    json.dump(gii_payload, fh, indent=2)
print(f"Wrote GII report to {gii_output_file}")

print("Done.")
