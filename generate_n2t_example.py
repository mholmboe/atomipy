#!/usr/bin/env python3
"""
generate_n2t_example.py - Script to demonstrate n2t file generation.

The script loads a structure file, optionally applies MINFF or CLAYFF
typing, and generates a GROMACS-compatible n2t (atom name to type) file.
Box dimensions are forwarded so neighbour distances respect periodic
boundary conditions.
"""

import os
import argparse
import atomipy as ap

def main():
    """Main function to demonstrate n2t file generation"""
    parser = argparse.ArgumentParser(description='Generate n2t file from structure')
    parser.add_argument('input_file', help='Input structure file (gro, pdb, xyz)')
    parser.add_argument('--forcefield', choices=['minff', 'clayff', 'auto'], default='auto',
                        help='Forcefield to use (minff, clayff, or auto-detect)')
    parser.add_argument('--output', help='Output n2t file name')
    args = parser.parse_args()
    
    print("AtomiPy N2T Generator")
    print("====================\n")
    
    # Check if input file exists
    if not os.path.isfile(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        return
    
    # Auto-detect file format and import
    print(f"Loading structure from {args.input_file}...")
    atoms, cell, box = ap.import_auto(args.input_file)
    print(f"Loaded {len(atoms)} atoms")
    
    # Ensure elements are assigned
    print("Assigning elements...")
    atoms = ap.element(atoms)
    
    # Set atomic masses
    print("Assigning masses...")
    from atomipy.mass import set_atomic_masses
    atoms = set_atomic_masses(atoms)
    
    # Process with forcefield if specified
    if args.forcefield != 'auto':
        print(f"Processing structure with {args.forcefield.upper()} forcefield...")
        if args.forcefield == 'minff':
            atoms, _ = ap.minff(atoms, box)
        else:  # clayff
            atoms, _ = ap.clayff(atoms, box)
    else:
        # If we're not processing with a forcefield, at least assign formal charges
        print("Auto-detecting structure types and assigning formal charges...")
        atoms = ap.assign_formal_charges(atoms)
    
    # Generate output file name if not provided
    output_file = args.output
    if not output_file:
        # Use input file name with .n2t extension
        base_name = os.path.splitext(os.path.basename(args.input_file))[0]
        output_file = f"{base_name}_{args.forcefield}.n2t"
    
    # Generate n2t file
    print(f"Generating n2t file with {args.forcefield} parameters...")
    n2t_path = ap.write_n2t(atoms, n2t_file=output_file, box=box)
    print(f"N2T file saved to: {n2t_path}")
    
    print("\nDone! You can use this file with GROMACS utilities like gmx x2top.")

if __name__ == "__main__":
    main()
