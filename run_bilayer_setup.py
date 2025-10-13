#!/usr/bin/env python
"""
Build script for creating a montmorillonite bilayer system with counter-ions and water.

Steps:
1. Replicate the pyrophyllite unit cell to 6x4x1
2. Create two montmorillonite layers through isomorphic substitution
3. Position the layers with appropriate basal spacing
4. Add sodium counter-ions in the interlayers
5. Solvate the interlayers with water
6. Apply MINFF force field and write output files
"""

import os
import sys
import numpy as np
import atomipy as ap

def main():
    # Input structure
    input_pdb = "Pyrophyllite_GII_0.071.pdb"
    output_name = "MMT_bilayer"
    
    # Step 1: Load and replicate the pyrophyllite unit cell
    print("Loading and replicating pyrophyllite unit cell...")
    atoms, cell, box_dim = ap.import_pdb(input_pdb)
    print(f"Loaded {len(atoms)} atoms with cell parameters: {cell}")
    box = box_dim  # Using cell parameters format for this script
    
    # Center the structure along z and then translate to have z=0 at the center
    atoms = ap.center(atoms, box_dim, dim="z")
    atoms = ap.translate(atoms, [0, 0, -box_dim[2]/2])
    # Write the centered structure to a PDB file
    ap.write_pdb(atoms, box_dim, f"{output_name}_centered.pdb")
    print(f"Wrote centered bilayer structure to {output_name}_centered.pdb")

    # Replicate the structure first
    atoms_6x4x1, box_6x4x1, _ = ap.replicate_system(atoms, box, [6, 4, 1])
    print(f"Created 6x4x1 pyrophyllite layer with {len(atoms_6x4x1)} atoms")

    # Write the wrapped bilayer structure to a PDB file
    ap.write_pdb(atoms_6x4x1, box_6x4x1, f"{output_name}_replicated.pdb")
    print(f"Wrote wrapped bilayer structure to {output_name}_replicated.pdb")
    
    # Step 2: Create first montmorillonite layer through isomorphic substitution
    print("Creating first montmorillonite layer...")
    mmt1_atoms, _, _ = ap.substitute(
        atoms_6x4x1, 
        box_6x4x1,
        num_oct_subst=16,
        o1_type='Alo',
        o2_type='Mgo',
        min_o2o2_dist=5.2
    )
    z_coords_1 = [atom['z'] for atom in mmt1_atoms]
    print(f"MMT1 z-range: {min(z_coords_1):.2f} to {max(z_coords_1):.2f} Å")
    
    # Step 3: Create second montmorillonite layer with different substitution pattern
    print("Creating second montmorillonite layer...")
    mmt2_atoms, _, _ = ap.substitute(
        atoms_6x4x1, 
        box_6x4x1,
        num_oct_subst=16,
        o1_type='Alo',
        o2_type='Mgo',
        min_o2o2_dist=5.2
    )
    z_coords_2 = [atom['z'] for atom in mmt2_atoms]
    print(f"MMT2 z-range before translation: {min(z_coords_2):.2f} to {max(z_coords_2):.2f} Å")
    
    # Step 4: Compute the new box dimensions with basal spacing of 18.9 Å
    basal_spacing = 18.9
    full_box = box_6x4x1.copy()
    # Set z-dimension to hold both layers
    full_box[2] = 2 * basal_spacing
    
    # Step 5: Position the two layers in the new cell
    # Layer 1: centered at z=0 (split across bottom/top PBC)
    # Layer 2: centered at z=basal_spacing (middle of box)
    mmt2_atoms = ap.translate(mmt2_atoms, [0, 0, basal_spacing])
    
    # Check if translation affected mmt1_atoms (should not!)
    z_coords_1_after = [atom['z'] for atom in mmt1_atoms]
    z_coords_2_after = [atom['z'] for atom in mmt2_atoms]
    print(f"After translation:")
    print(f"  MMT1 z-range: {min(z_coords_1_after):.2f} to {max(z_coords_1_after):.2f} Å")
    print(f"  MMT2 z-range: {min(z_coords_2_after):.2f} to {max(z_coords_2_after):.2f} Å")
    
    # Step 6: Combine the two layers
    combined_atoms = ap.update(mmt1_atoms, mmt2_atoms)
    print(f"Combined bilayer system has {len(combined_atoms)} atoms")
    print(f"Layer 1 centered at z = 0.00 Å (split across PBC)")
    print(f"Layer 2 centered at z = {basal_spacing:.2f} Å")
    print(f"Box height: {full_box[2]:.2f} Å")
    print(f"Two interlayers created at z ≈ 0-{basal_spacing:.1f} Å and z ≈ {basal_spacing:.1f}-{2*basal_spacing:.1f} Å")

    # Write the combined bilayer structure to a PDB file
    ap.write_pdb(combined_atoms, full_box, f"{output_name}_combined.pdb")
    print(f"Wrote combined bilayer structure to {output_name}_combined.pdb")
    
    # Step 6.6: Wrap the two MMT layers into the full periodic box
    print("Wrapping clay layers into periodic box...")
    combined_atoms = ap.wrap(combined_atoms, full_box, return_type='cartesian')
    print(f"Wrapped {len(combined_atoms)} atoms into box")
    
    # Write the wrapped bilayer structure to a PDB file
    ap.write_pdb(combined_atoms, full_box, f"{output_name}_wrapped.pdb")
    print(f"Wrote wrapped bilayer structure to {output_name}_wrapped.pdb")
    
    # Step 7: Add 16 Na counter-ions in each interlayer
    print("Adding sodium counter-ions...")
    
    # First interlayer (lower)
    # For the lower interlayer, use 0 to 1/4 of the basal spacing
    lower_limits = [0, 0, 0, full_box[0], full_box[1], basal_spacing]
    na_lower = ap.ionize(
        ion_type='Na',
        resname='NA',
        limits=lower_limits,
        num_ions=16,
        min_distance=3.0,
        solute_atoms=combined_atoms,
        placement='bulk'
    )
    print(f"Added {len(na_lower)} sodium ions to lower interlayer")
    
    # Second interlayer (upper)
    # For the upper interlayer, use 3/4 of the basal spacing to the full basal spacing
    upper_limits = [0, 0, basal_spacing, full_box[0], full_box[1],  full_box[2]]
    na_upper = ap.ionize(
        ion_type='Na',
        resname='NA',
        limits=upper_limits,
        num_ions=16,
        min_distance=3.0,
        solute_atoms=combined_atoms,
        placement='bulk'
    )
    print(f"Added {len(na_upper)} sodium ions to upper interlayer")
    
    # Combine clay layers and ions
    all_atoms = ap.update(combined_atoms, na_lower, na_upper)

        # Write the wrapped bilayer structure to a PDB file
    ap.write_pdb(all_atoms, full_box, f"{output_name}_ions.pdb")
    print(f"Wrote wrapped bilayer structure to {output_name}_ions.pdb")
    
    # Wrap all atoms (clay + ions) into the periodic box
    print("Wrapping all atoms (clay + ions) into periodic box...")
    all_atoms = ap.wrap(all_atoms, full_box, return_type='cartesian')
    print(f"Wrapped {len(all_atoms)} atoms into box")

        # Write the wrapped bilayer structure to a PDB file
    ap.write_pdb(all_atoms, full_box, f"{output_name}_wrapped_ions.pdb")
    print(f"Wrote wrapped bilayer structure to {output_name}_wrapped_ions.pdb")
    
    # Step 8: Solvate each interlayer with water
    print("Solvating interlayers...")
    
    # Solvate lower interlayer - use the same limits as for ions
    water_lower = ap.solvate(
        limits=lower_limits,
        max_solvent=360,
        solute_atoms=all_atoms,
        min_distance=2.0,
        solvent_type='spce'
    )
    print(f"Added {len(water_lower)} water atoms to lower interlayer")
    
    # Solvate upper interlayer - use the same limits as for ions
    water_upper = ap.solvate(
        limits=upper_limits,
        max_solvent=360,
        solute_atoms=all_atoms,
        min_distance=2.0,
        solvent_type='spce'
    )
    print(f"Added {len(water_upper)} water atoms to upper interlayer")
    
    # Combine everything and update indices
    final_atoms = ap.update(all_atoms, water_lower, water_upper)
    print(f"Final system has {len(final_atoms)} atoms")
    
    # Step 9: Apply MINFF force field
    print("Applying MINFF force field...")
    # Ensure the box is in the correct format
    if len(full_box) > 6:
        # If new_box is in 9-element triclinic format, convert to 6-element
        from atomipy.cell_utils import Box_dim2Cell
        cell_params = Box_dim2Cell(full_box)
        minff_atoms, itp_dir = ap.minff(final_atoms, cell_params)
    else:
        minff_atoms, itp_dir = ap.minff(final_atoms, full_box)
    
    # Step 10: Write output files
    print("Writing output files...")
    # PDB file
    ap.write_pdb(minff_atoms, Box=full_box, file_path=f"{output_name}.pdb")
    
    # GROMACS topology file
    print(f"itp_dir: {itp_dir}") # Debug print
    # Use cell_params which is already defined for the box parameter
    if len(full_box) > 6:
        # Use cell_params already defined
        from atomipy.cell_utils import Box_dim2Cell
        cell_params = Box_dim2Cell(full_box)
        ap.write_itp(minff_atoms, cell_params, f"{output_name}.itp")
    else:
        ap.write_itp(minff_atoms, full_box, f"{output_name}.itp")
    
    # LAMMPS data file
    if len(full_box) > 6:
        cell_params = Box_dim2Cell(full_box)  # Using Box_dim2Cell already imported above
        ap.write_lmp(minff_atoms, Box=cell_params, file_path=f"{output_name}.data")
    else:
        ap.write_lmp(minff_atoms, Box=full_box, file_path=f"{output_name}.data")
    
    print(f"\nSuccessfully created montmorillonite bilayer system '{output_name}'")
    print(f"Output file: {output_name}.pdb")

if __name__ == "__main__":
    main()
