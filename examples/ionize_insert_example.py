#!/usr/bin/env python
"""
Example script demonstrating how to use the atomipy ionize and insert functions.

This script:
1. Creates a simple clay mineral model
2. Adds ions to neutralize the system
3. Inserts molecules into the simulation box
4. Updates atom indices and writes the final structure
"""

import os
import sys
import numpy as np

# Add the parent directory to the path to import atomipy
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import atomipy as ap

def create_clay_sheet():
    """Create a simple clay sheet model."""
    # Create a 4x4 sheet of atoms
    atoms = []
    
    # Create clay framework with Al and Si atoms
    for x in range(4):
        for y in range(4):
            # Silicon layer (tetrahedral)
            si = {
                'index': len(atoms) + 1,
                'molid': 1,
                'resname': 'CLY',
                'type': 'Si',
                'element': 'Si',
                'x': x * 5.0,
                'y': y * 5.0,
                'z': 5.0
            }
            atoms.append(si)
            
            # Aluminum layer (octahedral) with 25% Mg substitution
            if (x + y) % 4 == 0:  # 25% substitution
                mg = {
                    'index': len(atoms) + 1,
                    'molid': 1,
                    'resname': 'CLY',
                    'type': 'Mg',
                    'element': 'Mg',
                    'x': x * 5.0,
                    'y': y * 5.0,
                    'z': 2.5
                }
                atoms.append(mg)
            else:
                al = {
                    'index': len(atoms) + 1,
                    'molid': 1,
                    'resname': 'CLY',
                    'type': 'Al',
                    'element': 'Al',
                    'x': x * 5.0,
                    'y': y * 5.0,
                    'z': 2.5
                }
                atoms.append(al)
            
            # Oxygen atoms surrounding Si and Al
            for dx, dy in [(0.8, 0), (-0.8, 0), (0, 0.8), (0, -0.8)]:
                o = {
                    'index': len(atoms) + 1,
                    'molid': 1,
                    'resname': 'CLY',
                    'type': 'O',
                    'element': 'O',
                    'x': x * 5.0 + dx,
                    'y': y * 5.0 + dy,
                    'z': 3.75  # Between Si and Al layers
                }
                atoms.append(o)
    
    # Set box dimensions (with some padding)
    box = [22.0, 22.0, 22.0]
    
    return atoms, box

def create_water_molecule():
    """Create a single water molecule template."""
    water = [
        {
            'index': 1,
            'molid': 1,
            'resname': 'SOL',
            'type': 'OW',
            'element': 'O',
            'x': 0.0,
            'y': 0.0,
            'z': 0.0
        },
        {
            'index': 2,
            'molid': 1,
            'resname': 'SOL',
            'type': 'HW1',
            'element': 'H',
            'x': 0.0,
            'y': 0.1,
            'z': 0.9
        },
        {
            'index': 3,
            'molid': 1,
            'resname': 'SOL',
            'type': 'HW2',
            'element': 'H',
            'x': 0.9,
            'y': 0.0,
            'z': -0.3
        }
    ]
    return water

def main():
    """Run the example script."""
    # Create a clay sheet model
    clay_atoms, box = create_clay_sheet()
    print(f"Created clay sheet with {len(clay_atoms)} atoms")
    
    # Add sodium ions to neutralize the Mg substitutions (4 Mg atoms)
    print("Adding sodium ions...")
    ions = ap.ionize(
        ion_type='Na',
        resname='NA',
        limits=box,
        num_ions=4,  # 4 Mg substitutions require 4 Na+ ions
        min_distance=3.0,
        solute_atoms=clay_atoms,
        placement='surface'  # Place ions near the clay surface
    )
    print(f"Added {len(ions)} sodium ions")
    
    # Combine clay and ions
    combined = ap.update(clay_atoms, ions)
    
    # Create water molecule template
    water_template = create_water_molecule()
    
    # Insert water molecules into the system
    print("Inserting water molecules...")
    waters = ap.insert(
        molecule_atoms=water_template,
        limits=box,
        rotate='random',
        min_distance=2.0,
        num_molecules=20,  # Insert 20 water molecules
        solute_atoms=combined
    )
    print(f"Inserted {len(waters) // 3} water molecules")
    
    # Combine all atoms and update indices
    final_atoms = ap.update(combined, waters)
    
    # Count molecules by residue name
    residue_counts = {}
    for atom in final_atoms:
        residue = atom.get('resname', 'Unknown')
        if residue not in residue_counts:
            residue_counts[residue] = 0
        
        # Only count each molecule once
        if atom['index'] == 1 or any(a['index'] == atom['index']-1 and a['molid'] != atom['molid'] for a in final_atoms):
            residue_counts[residue] += 1
    
    print("\nSystem composition:")
    for resname, count in residue_counts.items():
        print(f"  {resname}: {count} molecules")
    
    # Write the final structure to a PDB file
    output_file = "clay_ions_water.pdb"
    ap.write_pdb(final_atoms, Box=box, file_path=output_file)
    print(f"\nWrote final structure to {output_file}")

if __name__ == "__main__":
    main()
