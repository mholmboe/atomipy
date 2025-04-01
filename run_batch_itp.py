#!/usr/bin/env python3
"""
run_atomi.py - A beginner-friendly script to demonstrate atomipy functionality

This simple script shows the basic workflow for processing clay mineral structures
using the atomipy package. It demonstrates how to:
1. Import a structure file
2. Identify chemical elements
3. Create a supercell
4. Calculate bonds and angles
5. Assign specialized MINFF atom types
6. Generate a molecular topology file
"""

# Import the atomipy package and other required libraries
import atomipy as ap
import numpy as np
import os
from atomipy.cell_utils import Cell2Box_dim
from atomipy.minff import minff  # Directly import the minff function

def main():
    """Main function to demonstrate atomipy molecular modeling workflow"""
    print("AtomiPy Demo")
    print("======================================\n")

    # Step 1: Load a structure file in both PDB and GRO formats
    # ------------------------------------------
    pdb_file = "Kaolinite_GII_0.0487.pdb"
    
    print(f"Loading PDB structure from: {pdb_file}")
    # Read the PDB file - this returns a list of atom dictionaries and cell parameters
    pdb_atoms, cell = ap.import_conf.pdb(pdb_file)
    print(f"Successfully loaded {len(pdb_atoms)} atoms from PDB")
    # Convert cell parameters to box dimensions
    pdb_box_dim = Cell2Box_dim(cell)
    
    # Use the PDB structure for the rest of the script
    atoms = pdb_atoms
    box_dim = pdb_box_dim
    
    # Step 2: Assign chemical elements
    # -------------------------------
    # Each atom needs an element type (H, O, Si, Al, etc.)
    print("Assigning chemical elements to atoms...")
    ap.element(atoms)
    
    # Count elements in the structure
    elements = {}
    for atom in atoms:
        elem = atom.get('element', 'Unknown')
        elements[elem] = elements.get(elem, 0) + 1
    
    print("Element distribution:")
    for elem, count in elements.items():
        print(f"  {elem}: {count}")
    
    # Step 3: Create a supercell (replicate the structure)
    # ---------------------------------------------------
    # Replicate 5x3x4 times in x, y, z directions
    # The replicate function converts coordinates to fractional, replicates the unit cell,
    # then converts back to cartesian coordinates
    # Calculate replication based on desired size (approximately 25 Å in each dimension)
    target_size = 30  # Target size in Angstroms
    # Use the first 3 values from box_dim which represent x, y, z dimensions
    replicate_dims = np.ceil(np.array([target_size / box_dim[0], 
                                       target_size / box_dim[1], 
                                       target_size / box_dim[2]])).astype(int)
    print(f"Replicate dimensions: {replicate_dims} (target size: {target_size} Å)")
    print(f"\nCreating a {replicate_dims[0]}x{replicate_dims[1]}x{replicate_dims[2]} supercell...")
    
    replicated_atoms, replicated_box_dim, replicated_cell = ap.replicate.replicate_system(
        atoms, 
        box_dim, 
        replicate=replicate_dims,
        keep_molid=True, 
        renumber_index=True
    )
    
    print(f"Original atoms: {len(atoms)}, Replicated atoms: {len(replicated_atoms)}")
    
    # Step 4: Write the replicated structure to files
    # ----------------------------------------------
    print("\nSaving replicated structure...")
    ap.write_conf.gro(replicated_atoms, replicated_box_dim, "replicated_structure.gro")
    ap.write_conf.pdb(replicated_atoms, replicated_cell, "replicated_structure.pdb")
    
    # Step 5: Calculate bonds and angles in the structure
    # -------------------------------------------------
    print("\nCalculating bonds and angles...")
    # These parameters control how bonds are detected:
    # rmaxH = maximum bond length for hydrogen bonds (1.2 Å)
    # rmaxM = maximum bond length for metal-oxygen bonds (2.45 Å)
    replicated_atoms, Bond_index, Angle_index = ap.bond_angle(
        replicated_atoms, 
        replicated_box_dim, 
        rmaxH=1.2, 
        rmaxM=2.45
    )
    
    print(f"Found {len(Bond_index)} bonds and {len(Angle_index)} angles")
    
    # Write the system with bonds to a GRO file
    print("\nSaving structure with bond information...")
    ap.write_conf.gro(replicated_atoms, replicated_box_dim, "bonded_structure.gro")
    
    # Step 6: Assign MINFF atom types for force field
    # ---------------------------------------------
    print("\nAssigning specialized atom types using MINFF...")
    # MINFF classifies atoms based on their chemical environment
    # For example, oxygen atoms can be: Oh (hydroxyl), Op (bridging), Ow (water)
    # This function modifies atoms in-place (doesn't return anything)
    minff(replicated_atoms, replicated_box_dim)  # Use the directly imported minff function
    minff_atoms = replicated_atoms
    box_dim = replicated_box_dim
    
    # Count the different atom types
    atom_types = {}
    for atom in minff_atoms:
        atype = atom.get('type', 'Unknown')
        atom_types[atype] = atom_types.get(atype, 0) + 1
    
    print("Atom types after MINFF assignment:")
    for atype, count in atom_types.items():
        print(f"  {atype}: {count}")
    
    # Step 7: Generate a molecular topology file
    # ----------------------------------------
    print("\nGenerating molecular topology file...")
    # The topology (.itp) file contains all information needed for simulation:
    # - Atom definitions (type, charge, mass)
    # - Bond connections (only hydrogen bonds with current settings)
    # - Angle definitions
    ap.write_itp.write_itp(
        minff_atoms, 
        "molecular_topology.itp", 
        Box_dim=box_dim
    )
    
    # Step 8: Write final GRO file
    # ---------------------------
    print("Writing final structure to preem.gro...")
    ap.write_conf.gro(
        minff_atoms,
        box_dim,
        "preem.gro"
    )
    
    print("\nDone! Generated files:")
    print("1. replicated_structure.gro - Structure in GROMACS format")
    print("2. replicated_structure.pdb - Structure in PDB format")
    print("3. molecular_topology.itp - GROMACS topology file with H-bonds")
    print("4. preem.gro - Final structure with MINFF typing and charges")

if __name__ == "__main__":
    main()
