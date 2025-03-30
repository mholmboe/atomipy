import atomipy as ap
import numpy as np
import time
import argparse

# Define a mapping for the custom atom types to elements
atom_type_to_element = {
    'Mgo': 'Mg',
    'Si': 'Si',
    'Al': 'Al',
    'Omg': 'O',
    'Op': 'O',
    'Ohmg': 'O',
    'Ob': 'O',
    'H': 'H'
}

def assign_element_from_type(atom):
    """Assign element based on atom type"""
    atom_type = atom.get('type', '')
    if not atom_type and 'resname' in atom:
        # In this GRO file, the type is in column 3 (atom name)
        atom_name = atom.get('atname', '')
        if atom_name:
            atom['type'] = atom_name
            atom_type = atom_name
    
    # Get the element from our mapping
    if atom_type in atom_type_to_element:
        atom['element'] = atom_type_to_element[atom_type]
    else:
        # Try standard element guessing as fallback
        atom = ap.element(atom)
    
    return atom

# Import the GRO file
print("Importing preem.gro...")
atoms, box_dim = ap.import_gro('preem.gro')

# Print some information about the atoms
print(f"Successfully imported {len(atoms)} atoms")
print(f"Box dimensions: {box_dim}")

# Print the first few atoms before element assignment
print("\nFirst 5 atoms before element assignment:")
for i, atom in enumerate(atoms[:5]):
    print(f"Atom {i+1}: {atom['resname']} {atom.get('atname', '?')} at ({atom['x']:.3f}, {atom['y']:.3f}, {atom['z']:.3f})")

# Assign elements based on atom types
print("\nAssigning elements...")
for i in range(len(atoms)):
    atoms[i] = assign_element_from_type(atoms[i])

# Print the first few atoms again with elements
print("\nFirst 5 atoms with elements:")
for i, atom in enumerate(atoms[:5]):
    print(f"Atom {i+1}: {atom['resname']} {atom.get('element', '?')} at ({atom['x']:.3f}, {atom['y']:.3f}, {atom['z']:.3f})")

# Count elements in the system
element_counts = {}
for atom in atoms:
    element = atom.get('element', 'Unknown')
    if element in element_counts:
        element_counts[element] += 1
    else:
        element_counts[element] = 1

print("\nElement distribution:")
for element, count in element_counts.items():
    print(f"{element}: {count} atoms")

# Add command line argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Import and analyze a GRO file')
    parser.add_argument('--calc-bonds', action='store_true', help='Calculate bonds (slow for large systems)')
    parser.add_argument('--sample', type=int, default=0, help='Only calculate bonds for a sample of N atoms')
    args = parser.parse_args()
    
    # Calculate bonds if requested
    if args.calc_bonds:
        print("\nCalculating bonds (this might take a while)...")
        start_time = time.time()
        
        if args.sample > 0 and args.sample < len(atoms):
            # Only calculate bonds for a subset of atoms
            print(f"Calculating bonds for a sample of {args.sample} atoms")
            sample_atoms = atoms[:args.sample]
            sample_atoms = ap.bond_angle(sample_atoms, box_dim)
            
            # Print coordination numbers for the sample atoms
            print(f"\nCoordination numbers for {min(5, args.sample)} sample atoms:")
            for i, atom in enumerate(sample_atoms[:5]):
                neighbors = len(atom.get('neigh', []))
                element = atom.get('element', '?')
                print(f"Atom {i+1} ({element}): {neighbors} neighbors")
        else:
            # Calculate bonds for all atoms
            atoms = ap.bond_angle(atoms, box_dim)
            
            # Print coordination numbers for the first few atoms
            print("\nCoordination numbers for first 5 atoms:")
            for i, atom in enumerate(atoms[:5]):
                neighbors = len(atom.get('neigh', []))
                element = atom.get('element', '?')
                print(f"Atom {i+1} ({element}): {neighbors} neighbors")
        
        end_time = time.time()
        print(f"\nTime taken for bond calculation: {end_time - start_time:.2f} seconds")
    else:
        print("\nSkipping bond calculation. Use --calc-bonds to calculate bonds.")
        print("For large systems, you can use --sample N to only calculate bonds for N atoms.")
