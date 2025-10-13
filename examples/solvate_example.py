#!/usr/bin/env python
"""
Example script demonstrating how to use the atomipy solvation functions.

This script:
1. Loads a molecular structure
2. Creates a solvent box
3. Combines them using the solvate function
4. Updates the atom indices and molecule IDs
5. Writes the solvated structure to a PDB file
"""

import os
import sys
import numpy as np

# Add the parent directory to the path to import atomipy
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import atomipy as ap

def main():
    """Run the example script."""
    # Create a simple molecular structure (a small cube of atoms)
    solute_atoms = []
    
    # Create a 3x3x3 cube of carbon atoms
    for x in range(3):
        for y in range(3):
            for z in range(3):
                atom = {
                    'index': len(solute_atoms) + 1,
                    'molid': 1,
                    'resname': 'CUB',
                    'type': 'C',
                    'element': 'C',
                    'x': x * 3.0,
                    'y': y * 3.0,
                    'z': z * 3.0
                }
                solute_atoms.append(atom)
    
    # Set box dimensions (add padding around the cube)
    box = [10.0, 10.0, 10.0]
    
    print(f"Created solute with {len(solute_atoms)} atoms")
    
    # Solvate the structure with a shell of water molecules
    print("Solvating the structure...")
    solvated_atoms = ap.solvate(
        limits=box,
        min_distance=2.0,
        max_solvent='shell10',  # Create a 10 Å shell
        solute_atoms=solute_atoms,
        solvent_type='spce'
    )
    
    print(f"Created solvated structure with {len(solvated_atoms)} atoms")
    
    # Update atom indices and molecule IDs
    print("Updating atom indices and molecule IDs...")
    updated_atoms = ap.update(solvated_atoms)
    
    # Count molecules by residue name
    residue_counts = {}
    molids = set()
    for atom in updated_atoms:
        molids.add(atom['molid'])
        residue = atom.get('resname', 'Unknown')
        if residue not in residue_counts:
            residue_counts[residue] = 0
        
        # Only count each molecule once (by checking the first atom in each molecule)
        if any(a['molid'] == atom['molid'] and a['index'] < atom['index'] for a in updated_atoms):
            continue
        
        residue_counts[residue] += 1
    
    print("\nMolecule counts:")
    for resname, count in residue_counts.items():
        print(f"  {resname}: {count}")
    
    print(f"Total molecules: {len(molids)}")
    
    # Write the solvated structure to a PDB file
    output_file = "solvated_cube.pdb"
    ap.write_pdb(updated_atoms, Box=box, file_path=output_file)
    print(f"\nWrote solvated structure to {output_file}")

if __name__ == "__main__":
    main()
