#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple script to run bond_angle.py on preem.gro and save results to a MATLAB .mat file
Using the improved atomipy.element to assign chemical elements and adjusting indices for MATLAB
"""

import os
import sys
import numpy as np

# Add main project directory to path to import atomipy
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import atomipy as ap

# Import cell_utils for format conversions
from atomipy import cell_utils

# Try to import scipy.io for saving .mat files
try:
    import scipy.io as sio
except ImportError:
    print("scipy is not installed. Installing it now...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scipy"])
    import scipy.io as sio

def main():
    # 1. Import preem.gro
    print("Importing preem.gro file...")
    file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'preem.gro')
    atoms, box_dim = ap.import_gro(file_path)
    print(f"Successfully imported {len(atoms)} atoms")
    print(f"Box dimensions: {box_dim}")
    
    # 2. Transfer atname to type field for better element detection
    print("\nPreparing atoms for element assignment...")
    for i in range(len(atoms)):
        # In GRO files, atom type information is typically in the atname field
        if 'atname' in atoms[i] and not atoms[i].get('type', ''):
            atoms[i]['type'] = atoms[i]['atname']
    
    # 3. Assign elements using the improved element.py function
    print("Assigning elements to atoms...")
    atoms = ap.element(atoms)  # Process all atoms at once
    
    # 3a. Write full structure to output file (before subsetting)
    print("\nWriting full structure to output file...")
    out_gro_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out.gro')
    
    # Write GRO file
    from atomipy import write_module
    write_module.gro(atoms, box_dim, out_gro_file)
    print(f"Full structure written to {out_gro_file}")
    
    # Process only a subset of atoms for faster testing
    # Use a smaller sample to reduce computational burden
    # atoms = atoms[:960]  # Process only first 960 atoms
    
    # Check if molid exists in the atoms
    has_molid = any('molid' in atom for atom in atoms)
    print(f"\nMolecule IDs found in data: {'Yes' if has_molid else 'No'}")
    
    # 4. Run bond_angle
    print("\nCalculating bonds and angles...")
    print("Excluding bonds between atoms of the same element...")
    print("Restricting bonds to atoms in the same molecule (same molid)...")
    
    # Use original cutoff values: rmaxH=1.2 and rmaxM=2.45 Angstrom
    updated_atoms, bond_index, angle_index = ap.bond_angle(atoms, box_dim, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True)
    
    # 5. Report results
    print(f"Found {len(bond_index)} bonds and {len(angle_index)} angles")
    
    # 6. Adjust indices for MATLAB (add 1 to all atom indices)
    print("\nAdjusting atom indices for MATLAB compatibility (0-based → 1-based)...")
    
    # Create copies to avoid modifying the original arrays
    bond_index_matlab = bond_index.copy()
    angle_index_matlab = angle_index.copy()
    
    # For bond_index: first two columns are atom indices
    bond_index_matlab[:, 0] += 1  # atom1_idx
    bond_index_matlab[:, 1] += 1  # atom2_idx
    
    # For angle_index: first three columns are atom indices
    angle_index_matlab[:, 0] += 1  # atom1_idx
    angle_index_matlab[:, 1] += 1  # atom2_idx (center atom)
    angle_index_matlab[:, 2] += 1  # atom3_idx
    
    # 7. Save to .mat file
    output_file = 'bond_angle_results.mat'
    print(f"\nSaving results to {output_file}...")
    
    # Create dictionary with results
    results = {
        'bond_index': bond_index_matlab,
        'angle_index': angle_index_matlab,
        'box_dim': box_dim,
        'atoms': {
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
    
    # 8. For verification, show index ranges
    print("\nVerification of index adjustment:")
    print(f"  Original bond_index indices: {int(bond_index[:,0].min())}-{int(bond_index[:,0].max())}, {int(bond_index[:,1].min())}-{int(bond_index[:,1].max())}")
    print(f"  MATLAB bond_index indices: {int(bond_index_matlab[:,0].min())}-{int(bond_index_matlab[:,0].max())}, {int(bond_index_matlab[:,1].min())}-{int(bond_index_matlab[:,1].max())}")
    print(f"  Original angle_index indices: {int(angle_index[:,0].min())}-{int(angle_index[:,0].max())}, {int(angle_index[:,1].min())}-{int(angle_index[:,1].max())}, {int(angle_index[:,2].min())}-{int(angle_index[:,2].max())}")
    print(f"  MATLAB angle_index indices: {int(angle_index_matlab[:,0].min())}-{int(angle_index_matlab[:,0].max())}, {int(angle_index_matlab[:,1].min())}-{int(angle_index_matlab[:,1].max())}, {int(angle_index_matlab[:,2].min())}-{int(angle_index_matlab[:,2].max())}")

if __name__ == '__main__':
    main()
