#!/usr/bin/env python3
"""
run_batch_process.py - Script to batch process multiple GRO files


This script processes multiple GRO files from a specified folder by:
1. Loading each GRO file (preem1.gro through preem45.gro)
2. Assigning MINFF atom types to each structure
3. Generating new GRO files with the processed structures
4. Creating molecular topology (.itp) files
5. Tracking total charge for each system and saving to a summary file
"""

# Import the atomipy package and other required libraries
import atomipy as ap
import numpy as np
import os
import glob

def main():
    """Main function to batch process GRO files"""
    print("AtomiPy Batch Processing")
    print("======================================\n")

    # Create output directory if it doesn't exist
    os.makedirs("output", exist_ok=True)
    
    # Open the charges summary file
    with open("output/charge_summary.txt", "w") as summary_file, open("output/atom_type_changes.txt", "w") as type_changes_file:
        summary_file.write("Filename\tTotal Charge\tAtom Count\tAtom Types\n")
        type_changes_file.write("Filename\tChanged Atom Count\tAtom Type Changes\n")
        
        # Loop through all preem{1-45}.gro files in the conf directory
        for i in range(1,46): #(1, 46):
            input_file = f"conf/preem{i}.gro"
            
            # Skip if the file doesn't exist
            if not os.path.exists(input_file):
                print(f"Warning: {input_file} does not exist. Skipping...")
                continue
                
            print(f"\nProcessing file {i}/45: {input_file}")
            
            # Step 1: Load the GRO file
            # ------------------------------
            print(f"Loading structure from: {input_file}")
            try:
                atoms, Box_dim = ap.import_conf.gro(input_file)
                print(f"Successfully loaded {len(atoms)} atoms")
                
                # Store original atom types before processing
                original_types = []
                for atom in atoms:
                    original_types.append(atom.get('type', ''))
            except Exception as e:
                print(f"Error loading {input_file}: {e}")
                continue
            
            # Step 2: Assign chemical elements to atoms
            # --------------------------------------
            print("Assigning chemical elements to atoms...")
            atoms = ap.element(atoms)
            
            # Count elements for reporting
            element_counts = {}
            for atom in atoms:
                if 'element' in atom:
                    elem = atom['element']
                    element_counts[elem] = element_counts.get(elem, 0) + 1
            
            print("Element distribution:")
            for elem, count in element_counts.items():
                print(f"  {elem}: {count}")
            
            # Step 3: Calculate bonds and angles in the structure
            # -------------------------------------------------
            print("\nCalculating bonds and angles...")
            # These parameters control how bonds are detected:
            # rmaxH = maximum bond length for hydrogen bonds (1.2 Å)
            # rmaxM = maximum bond length for metal-oxygen bonds (2.45 Å)
            atoms, Bond_index, Angle_index = ap.bond_angle(
                atoms, 
                Box=Box_dim, 
                rmaxH=1.2, 
                rmaxM=2.45
            )
            
            print(f"Found {len(Bond_index)} bonds and {len(Angle_index)} angles")
            
            # Step 4: Assign MINFF atom types for force field
            # ---------------------------------------------
            print("\nAssigning specialized atom types using MINFF...")
            # MINFF classifies atoms based on their chemical environment
            # This function modifies atoms in-place (doesn't return anything)
            ap.minff(atoms, Box=Box_dim)  # Use the minff function from the ap package
            
            # Calculate total charge
            total_charge = sum(atom.get('charge', 0.0) for atom in atoms)
            print(f"Total charge: {total_charge}")
            
            # Count atom types
            type_counts = {}
            for atom in atoms:
                if 'fftype' in atom:
                    atom_type = atom['fftype']
                    type_counts[atom_type] = type_counts.get(atom_type, 0) + 1
            
            print("Atom types after MINFF assignment:")
            type_summary = ""
            for atom_type, count in type_counts.items():
                print(f"  {atom_type}: {count}")
                type_summary += f"{atom_type}:{count},"
            
            # Write to the summary file
            summary_file.write(f"{input_file}\t{total_charge:.6f}\t{len(atoms)}\t{type_summary[:-1]}\n")
            
            # Compare original atom types with new atom types after MINFF assignment
            new_types = [atom.get('type', '') for atom in atoms]
            type_changes = []
            for idx, (orig, new) in enumerate(zip(original_types, new_types)):
                if orig != new:
                    type_changes.append(f"Atom {idx+1}: {orig} → {new}")
            
            # Write type changes to file
            if type_changes:
                changes_summary = "; ".join(type_changes[:10])
                if len(type_changes) > 10:
                    changes_summary += f"; ...and {len(type_changes)-10} more"
                type_changes_file.write(f"{input_file}\t{len(type_changes)}\t{changes_summary}\n")
            else:
                type_changes_file.write(f"{input_file}\t0\tNo changes\n")
            
            # Step 5: Generate output files
            # ---------------------------
            # Output GRO file
            output_gro = f"output/pypreem{i}.gro"
            print(f"\nWriting processed structure to {output_gro}...")
            ap.write_gro(atoms, Box=Box_dim, file_path=output_gro)
            
            # Output ITP file
            output_itp = f"output/pymin{i}.itp"
            print(f"Generating molecular topology file {output_itp}...")
            ap.write_itp(atoms, Box=Box_dim, file_path=output_itp)
            
            print(f"Completed processing {input_file}")
    
    print("\nAll processing complete!")
    print("Output files have been saved to the 'output' directory.")
    print("Charge summary has been saved to 'output/charge_summary.txt'.")
    print("Atom type changes have been saved to 'output/atom_type_changes.txt'.")

if __name__ == "__main__":
    main()
