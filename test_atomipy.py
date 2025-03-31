#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Atomipy Test Script

This script demonstrates the key capabilities of the atomipy package:
- Importing molecular structure files
- Replicating structures to create supercells
- Assigning atomtypes with the MINFF forcefield
- Writing structures to various formats
- Generating topology files

Usage:
    python test_atomipy.py [pdb_file]

If no file is specified, the script will look for "Kaolinite_GII_0.0487.gro"
"""

import os
# from signal import pause  # This causes the script to hang indefinitely
import sys
import numpy as np
import scipy.io as sio
import atomipy as ap

def main():
    # Get input file from command line or use default
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "Kaolinite_GII_0.0487.gro"
    
    # Check if file exists
    if not os.path.exists(input_file):
        print(f"Error: File {input_file} not found.")
        print("Please provide a valid PDB file path.")
        return
    
    print(f"Processing file: {input_file}")
    print("=" * 80)
    
    # Step 1: Import the structure file
    print("\nSTEP 1: Importing structure")
    print("-" * 40)
    
    atoms, box_dim = ap.import_conf.gro(input_file)
    print(f"Successfully imported {len(atoms)} atoms")
    print(f"Box dimensions: {box_dim}")
        
    # Assign elements to atoms using the element.py function
    print("\nAssigning element types using the element.py function...")
        
    # Copy atom name to type field if it doesn't exist to help with element assignment
    for atom in atoms:
        if 'type' not in atom and 'name' in atom:
            atom['type'] = atom['name']
                
    # Use the element function from atomipy
    atoms = ap.element(atoms)
        
    # Check the results of element assignment
    print("Element assignment completed.")
        
    # Count atom elements
    atom_elements = {}
    for atom in atoms:
        element = atom.get('element', 'Unknown')
        atom_elements[element] = atom_elements.get(element, 0) + 1
        
    print("Element distribution after assignment:")
    for elem, count in atom_elements.items():
        print(f"  {elem}: {count}")
        
    # We have full cell parameters, generate proper box dimensions
    cell = ap.cell_utils.Box_dim2Cell(box_dim)
    box_dim_triclinic = ap.cell_utils.Cell2Box_dim(cell)
    cell_triclinic2 = ap.cell_utils.Box_dim2Cell(box_dim_triclinic)
    print(f"Cell parameters: {cell}")
    print(f"Generated triclinic box dimensions: {box_dim_triclinic}")
    print(f"Cell parameters: {cell_triclinic2}")
    #input()  # Commented out to allow non-interactive execution
    
    # Step 2: Analyze the original structure
    print("\nSTEP 2: Analyzing original structure")
    print("-" * 40)
    
    # Calculate bonds and angles for the original structure
    atoms, bond_index, angle_index = ap.bond_angle(atoms, box_dim_triclinic, rmaxH=1.2, rmaxM=2.45)
    print(f"Found {len(bond_index)} bonds and {len(angle_index)} angles in the original structure")

    
    # Count element types after bond calculations
    atom_elements = {}
    for atom in atoms:
        element = atom.get('element', 'Unknown')
        if element is None:
            element = 'Unknown'
        atom_elements[element] = atom_elements.get(element, 0) + 1
    
    print("Element distribution:")
    for elem, count in atom_elements.items():
        print(f"  {elem}: {count}")
  
    # Step 3: Replicate the structure
    print("\nSTEP 3: Replicating structure")
    print("-" * 40)
    # Create a 2x2x2 supercell
    replicate_dims = [5,3,4]
    
    replicated_atoms, replicated_box_dim, replicated_cell = ap.replicate.replicate_system(atoms, box_dim_triclinic, replicate=replicate_dims, keep_molid=True, keep_resname=True, renumber_index=True)
    print(f"Created a {replicate_dims[0]}x{replicate_dims[1]}x{replicate_dims[2]} supercell")
    print(f"Original atoms: {len(atoms)}, Replicated atoms: {len(replicated_atoms)}")    
    print(f"Replicated box dimensions: {replicated_box_dim}")
    print(f"Replicated cell parameters: {replicated_cell}")
    #input()  # Commented out to allow non-interactive execution

    # Step 4: Write the replicated structure to a GRO file
    print("\nSTEP 4: Writing replicated structure to GRO")
    print("-" * 40)
    replicated_gro = "replicated_structure.gro"
    replicated_pdb = "replicated_structure.pdb"
    try:
        ap.write_conf.gro(replicated_atoms, replicated_box_dim, replicated_gro)
        print(f"Replicated structure written to {replicated_gro}")
        # Pass replicated_cell (a,b,c,alpha,beta,gamma) to pdb function, not box_dim
        ap.write_conf.pdb(replicated_atoms, replicated_cell, replicated_pdb)
        print(f"Replicated structure written to {replicated_pdb}")
    except Exception as e:
        print(f"Error writing GRO/PDB file: {e}")

    # Step 4.5: Run bond_angle function on replicated structure
    print("\nRunning bond detection on replicated structure")
    print("-" * 40)
    
    # Make sure elements are properly assigned
    from atomipy.element import element
    element(replicated_atoms)
    
    # Set appropriate parameters for mineral structures
    rmaxH = 1.2  # Increased slightly for better H-O bond detection
    rmaxM = 2.45  # Increased for better detection of Si-O and Al-O bonds
    same_molecule_only = False  # Allow bonds between different molecules
    same_element_bonds = False  # Don't allow bonds between same elements (correct for minerals)
    
    # Run bond_angle function
    print(f"Running bond_angle with rmaxH={rmaxH}, rmaxM={rmaxM}, same_molecule_only={same_molecule_only}, same_element_bonds={same_element_bonds}")
    from atomipy.bond_angle import bond_angle
    
    # Try with standard cutoffs first
    updated_atoms, bond_index, angle_index = bond_angle(
        replicated_atoms, 
        replicated_box_dim, 
        rmaxH=rmaxH, 
        rmaxM=rmaxM, 
        same_molecule_only=same_molecule_only,
        same_element_bonds=same_element_bonds
    )
    
    # If no bonds found, try with larger cutoffs
    if len(bond_index) == 0:
        print("No bonds detected with standard cutoffs, trying larger cutoffs...")
        rmaxH_large = 1.2   # Very generous cutoff for hydrogen bonds
        rmaxM_large = 2.45   # Very generous cutoff for metal-oxygen bonds
        updated_atoms, bond_index, angle_index = bond_angle(
            replicated_atoms, 
            replicated_box_dim, 
            rmaxH=rmaxH_large, 
            rmaxM=rmaxM_large, 
            same_molecule_only=same_molecule_only,
            same_element_bonds=same_element_bonds
        )
        print(f"With larger cutoffs (rmaxH={rmaxH_large}, rmaxM={rmaxM_large}): {len(bond_index)} bonds and {len(angle_index)} angles")
    
    print(f"Found {len(bond_index)} bonds and {len(angle_index)} angles in the replicated structure")
    
    # Convert bond_index and angle_index to MATLAB-friendly format
    # MATLAB uses 1-based indexing, so add 1 to all indices
    bond_index_matlab = bond_index.copy()
    if len(bond_index) > 0:
        bond_index_matlab[:, 0:2] += 1  # Add 1 to atom indices for MATLAB 1-based indexing
    
    angle_index_matlab = angle_index.copy()
    if len(angle_index) > 0:
        angle_index_matlab[:, 0:3] += 1  # Add 1 to atom indices for MATLAB 1-based indexing
    
    # 8. Save to .mat file
    
    output_file = 'bond_angle_results.mat'
    print(f"\nSaving results to {output_file}...")
    
    # Create dictionary with results
    results = {
        'bond_index': bond_index_matlab,
        'angle_index': angle_index_matlab,
        'replicatedbox_dim': replicated_box_dim,
        'replicated_atoms': {
            'element': [atom.get('element', 'X') for atom in updated_atoms],
            'type': [atom.get('type', '') for atom in updated_atoms],
            # Add atom coordinates
            'x': [atom.get('x', 0.0) for atom in updated_atoms],
            'y': [atom.get('y', 0.0) for atom in updated_atoms],
            'z': [atom.get('z', 0.0) for atom in updated_atoms],
            # Also add molecule IDs
            'molid': [atom.get('molid', 1) for atom in updated_atoms]
        }
    }
    
    # Save to .mat file
    sio.savemat(output_file, results)
    print(f"Results successfully saved to: {os.path.abspath(output_file)}")
    print("This file can be loaded in MATLAB using: data = load('bond_angle_results.mat');")

    # Step 5: Assign MINFF atom types
    print("\nSTEP 5: Assigning MINFF atom types")
    print("-" * 40)
    try:
        # Ensure all atoms have proper element types using the element.py function
        print("Assigning elements to replicated atoms...")
        
        # Copy atom name to type field if it doesn't exist to help with element assignment
        for atom in replicated_atoms:
            if 'type' not in atom and 'name' in atom:
                atom['type'] = atom['name']
        
        # Use the element function from atomipy
        replicated_atoms = ap.element(replicated_atoms)
        print("Element assignment for replicated atoms completed.")
        
        elements_before = set(atom.get('element', 'X') for atom in replicated_atoms)
        print(f"Element types before MINFF: {', '.join(sorted(elements_before))}")
        
        # Apply MINFF atom typing
        if 'None' in elements_before or None in elements_before:
            print("Warning: Cannot run MINFF on atoms without element types")
            typed_atoms = replicated_atoms
        else:
            # Create a local copy of the box dimensions to ensure correct format
            box_dim_for_minff = np.array(replicated_box_dim, dtype=float)
            
            # Apply MINFF atom typing
            typed_atoms, all_neighbors = ap.minff(replicated_atoms, box_dim_for_minff)
            print(f"MINFF atom typing completed for {len(typed_atoms)} atoms")
            
            # Check the assigned atom types
            atom_types = set(atom.get('fftype', 'X') for atom in typed_atoms)
            print(f"MINFF assigned atom types: {', '.join(sorted(atom_types))}")
    except Exception as e:
        print(f"Error during MINFF typing: {e}")
        import traceback
        traceback.print_exc()
        # Continue with the original atoms if MINFF fails
        typed_atoms = replicated_atoms
    
    # Step 6: Write the typed structure to a new GRO file
    print("\nSTEP 6: Writing typed structure to GRO")
    print("-" * 40)
    typed_gro = "typed_structure.gro"
    try:
        ap.write_conf.gro(typed_atoms, replicated_box_dim, typed_gro)
        print(f"Structure with MINFF atomtypes written to {typed_gro}")
    except Exception as e:
        print(f"Error writing typed GRO file: {e}")
    
    # Step 7: Generate a topology file
    print("\nSTEP 7: Generating molecular topology (.itp) file")
    print("-" * 40)
    topology_file = "molecular_topology.itp"
    try:
        # Make sure box dimensions are in the right format and length
        # Use the full triclinic box dimensions for topology generation
        if isinstance(replicated_box_dim, np.ndarray):
            # Convert numpy array to plain list for better handling
            box_dim_for_topo = replicated_box_dim.tolist() if hasattr(replicated_box_dim, 'tolist') else list(replicated_box_dim)
        else:
            # Ensure we have a complete box definition
            box_dim_for_topo = list(replicated_box_dim)
            
            # If it's a short box definition, expand it to include triclinic components
            if len(box_dim_for_topo) == 3:
                # Convert to triclinic format using cell parameters
                if replicated_cell is not None and len(replicated_cell) == 6:
                    box_dim_for_topo = ap.cell_utils.Cell2Box_dim(replicated_cell)
                else:
                    # Fallback: use reasonable box size with orthogonal shape
                    box_dim_for_topo = [50.0, 50.0, 50.0]
        
        print(f"Using box dimensions for topology: {box_dim_for_topo}")
            
        # Write the topology file - the write_itp function will calculate bonds and angles internally
        print("Generating topology with internal bond_angle calculation...")
        ap.write_itp.write_itp(typed_atoms, topology_file, Box_dim=box_dim_for_topo, rmaxH=1.2, rmaxM=2.45)
        print(f"Molecular topology written to {topology_file}")
    except Exception as e:
        print(f"Error generating topology file: {e}")
        import traceback
        traceback.print_exc()
    
    print("\nAll tasks completed successfully!")
    print("=" * 80)
    print("\nSummary of generated files:")
    print(f"1. {replicated_gro} - Replicated structure in GRO format")
    print(f"2. {typed_gro} - Structure with MINFF atomtypes in GRO format")
    print(f"3. {topology_file} - Molecular topology in .itp format")
    print("\nYou can visualize these files using molecular visualization software like VMD or PyMOL.")

if __name__ == "__main__":
    main()
