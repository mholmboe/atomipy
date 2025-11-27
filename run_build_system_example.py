#!/usr/bin/env python3
"""
Build a montmorillonite bilayer system with counter-ions and water.

This script creates a montmorillonite bilayer by:
1. Replicating a pyrophyllite unit cell
2. Creating two montmorillonite layers through isomorphic substitution
3. Positioning the layers with appropriate spacing
4. Adding sodium counter-ions
5. Solvating the interlayers with water
6. Writing output files
"""

import atomipy as ap

# ===== USER-DEFINED PARAMETERS =====
# Input structure
input_pdb = "Pyrophyllite_GII_0.071.pdb"
output_name = "MMT_bilayer"

# Replication factors (x, y, z)
replication = [6, 4, 1]

# Isomorphic substitutions (number of octahedral and tetrahedral substitutions)
n_subst = 16  # Number of octahedral Al3+ to Mg2+ substitutions
orig_o_type = 'Alo'
subst_o_type = 'Mgo'
cation = 'Na'

# Layer separation and spacing (in Angstroms)
basal_spacing = 22
nWaters = 360

# Step 1: Load and replicate the pyrophyllite unit cell
atoms, cell = ap.import_pdb(input_pdb)
box_dim = ap.Cell2Box_dim(cell)
box = box_dim  # Using cell parameters format for this script

# Center the structure along z and then translate to have z=0 at the center
atoms = ap.center(atoms, box_dim, dim="z")
atoms = ap.translate(atoms, [0, 0, -box_dim[2]/2])
# Write the centered structure to a PDB file
ap.write_pdb(atoms, box_dim, f"{output_name}_centered.pdb")

# Replicate the structure first
atoms_6x4x1, box_6x4x1, _ = ap.replicate_system(atoms, box, [6, 4, 1])

# Write the wrapped bilayer structure to a PDB file
ap.write_pdb(atoms_6x4x1, box_6x4x1, f"{output_name}_replicated.pdb")

# Step 2: Create first montmorillonite layer through isomorphic substitution
mmt1_atoms, _, _ = ap.substitute(
    atoms_6x4x1,
    box_6x4x1,
    num_oct_subst=n_subst,
    o1_type=orig_o_type,
    o2_type=subst_o_type,
    min_o2o2_dist=5.2
)
# Tag first layer residue name for clean separation
for atom in mmt1_atoms:
    atom['resname'] = 'MMT1'
    atom['molid'] = 1
z_coords_1 = [atom['z'] for atom in mmt1_atoms]

# Step 3: Create second montmorillonite layer with different substitution pattern
mmt2_atoms, _, _ = ap.substitute(
    atoms_6x4x1,
    box_6x4x1,
    num_oct_subst=n_subst,
    o1_type=orig_o_type,
    o2_type=subst_o_type,
    min_o2o2_dist=5.2
)

# Step 4: Compute the new box dimensions with basal spacing of 18.9 Å
full_box = box_6x4x1.copy()
# Set z-dimension to hold both layers
full_box[2] = 2 * basal_spacing

# Step 5: Position the two layers in the new cell
# Layer 1: centered at z=0 (split across bottom/top PBC)
# Layer 2: centered at z=basal_spacing (middle of box)
mmt2_atoms = ap.translate(mmt2_atoms, [0, 0, basal_spacing])

# Step 6: Combine the two layers
combined_atoms = ap.update(mmt1_atoms, mmt2_atoms)

# Write the combined bilayer structure to a PDB file
ap.write_pdb(combined_atoms, full_box, f"{output_name}_combined.pdb")

# Step 6.6: Wrap the two MMT layers into the full periodic box
combined_atoms = ap.wrap(combined_atoms, full_box, return_type='cartesian')

# Write the wrapped bilayer structure to a PDB file
ap.write_pdb(combined_atoms, full_box, f"{output_name}_wrapped.pdb")

# Step 7: Add 16 Na counter-ions in each interlayer
# First interlayer (lower)
# For the lower interlayer, use 0 to 1/4 of the basal spacing
lower_limits = [0, 0, 0, full_box[0], full_box[1], basal_spacing]
na_lower = ap.ionize(
    ion_type=cation,
    resname=cation,
    limits=lower_limits,
    num_ions=n_subst,
    min_distance=3.0,
    solute_atoms=combined_atoms,
    placement='bulk'
)

# Second interlayer (upper)
# For the upper interlayer, use 3/4 of the basal spacing to the full basal spacing
upper_limits = [0, 0, basal_spacing, full_box[0], full_box[1],  full_box[2]]
na_upper = ap.ionize(
    ion_type=cation,
    resname=cation,
    limits=upper_limits,
    num_ions=n_subst,
    min_distance=3.0,
    solute_atoms=combined_atoms,
    placement='bulk'
)

# Combine clay layers and ions
all_atoms = ap.update(combined_atoms, na_lower, na_upper)

# Write the wrapped bilayer structure to a PDB file
ap.write_pdb(all_atoms, full_box, f"{output_name}_ions.pdb")

# Wrap all atoms (clay + ions) into the periodic box
all_atoms = ap.wrap(all_atoms, full_box, return_type='cartesian')

# Write the wrapped bilayer structure to a PDB file
ap.write_pdb(all_atoms, full_box, f"{output_name}_wrapped_ions.pdb")

# Step 8: Solvate each interlayer with water
# Solvate lower interlayer - use the same limits as for ions
water_lower = ap.solvate(
    limits=lower_limits,
    max_solvent=nWaters,
    solute_atoms=all_atoms,
    min_distance=2.0,
    solvent_type='spce',
)

# Solvate upper interlayer - use the same limits as for ions
water_upper = ap.solvate(
    limits=upper_limits,
    max_solvent=nWaters,
    solute_atoms=all_atoms,
    min_distance=2.0,
    solvent_type='spce',
)

# Combine everything and update indices
final_atoms = ap.update(all_atoms, water_lower, water_upper)

# Step 8b: Write output files
# PDB file
ap.write_pdb(final_atoms, Box=full_box, file_path=f"{output_name}_final_atoms.pdb")

# Step 9: Apply MINFF force field

# Ensure the box is in the correct format
if len(full_box) > 6:
    # If new_box is in 9-element triclinic format, convert to 6-element
    cell = ap.Box_dim2Cell(full_box)
    minff_atoms = ap.minff(final_atoms, cell)
else:
    minff_atoms = ap.minff(final_atoms, full_box)

# Step 10: Write output files
# PDB file
ap.write_pdb(minff_atoms, Box=full_box, file_path=f"{output_name}.pdb")

# GROMACS topology file
if len(full_box) > 6:
    cell = ap.Box_dim2Cell(full_box)
    ap.write_itp(minff_atoms, cell, f"{output_name}.itp")
else:
    ap.write_itp(minff_atoms, full_box, f"{output_name}.itp")

# LAMMPS data file
if len(full_box) > 6:
    cell = ap.Box_dim2Cell(full_box)
    ap.write_lmp(minff_atoms, Box=cell, file_path=f"{output_name}.data")
else:
    ap.write_lmp(minff_atoms, Box=full_box, file_path=f"{output_name}.data")
