import atomipy as ap
import numpy as np
import time
import os

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

def analyze_bond_angles(atoms):
    """Analyze bond angles and print statistics"""
    # Count total bonds
    total_bonds = 0
    bond_lengths = []
    
    # Count total angles
    total_angles = 0
    angle_values = []
    
    # Coordination numbers
    coordination = {}
    
    for atom in atoms:
        # Count bonds
        neighbors = len(atom.get('neigh', []))
        total_bonds += neighbors
        
        # Record bond lengths
        for bond in atom.get('bonds', []):
            bond_lengths.append(bond[1])
        
        # Record angle values
        for angle in atom.get('angles', []):
            angle_values.append(angle[1])
            total_angles += 1
        
        # Record coordination numbers by element
        element = atom.get('element', 'Unknown')
        if element not in coordination:
            coordination[element] = []
        coordination[element].append(neighbors)
    
    # Calculate statistics
    print(f"\nTotal bonds identified: {total_bonds // 2}")  # Divide by 2 because each bond is counted twice
    
    if bond_lengths:
        print(f"Average bond length: {np.mean(bond_lengths):.3f} Å")
        print(f"Min bond length: {np.min(bond_lengths):.3f} Å")
        print(f"Max bond length: {np.max(bond_lengths):.3f} Å")
    
    print(f"\nTotal angles calculated: {total_angles}")
    
    if angle_values:
        print(f"Average angle: {np.mean(angle_values):.2f}°")
        print(f"Min angle: {np.min(angle_values):.2f}°")
        print(f"Max angle: {np.max(angle_values):.2f}°")
    
    print("\nAverage coordination numbers by element:")
    for element, coords in coordination.items():
        print(f"{element}: {np.mean(coords):.2f}")

# Main execution
print("Importing preem.gro...")
file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'preem.gro')
atoms, box_dim = ap.import_gro(file_path)

print(f"Successfully imported {len(atoms)} atoms")
print(f"Box dimensions: {box_dim}")

# Assign elements based on atom types
print("\nAssigning elements...")
for i in range(len(atoms)):
    atoms[i] = assign_element_from_type(atoms[i])

# Count elements in the system
element_counts = {}
for atom in atoms:
    element = atom.get('element', 'Unknown')
    if element in element_counts:
        element_counts[element] += 1
    else:
        element_counts[element] = 1

print("\nElement distribution:")
for element, count in sorted(element_counts.items()):
    print(f"{element}: {count} atoms")

# Calculate bonds and angles
print("\nCalculating bonds and angles (this might take a while)...")
start_time = time.time()

# You can use a smaller sample for testing or large systems
sample_size = len(atoms)  # Use all atoms - set to a smaller number if it's too slow
sample_atoms = atoms[:sample_size]

# Run the bond_angle function
sample_atoms, Bond_index, Angle_index = ap.bond_angle(sample_atoms, box_dim)

# Print information about the Bond_index and Angle_index
print(f"Bond_index shape: {Bond_index.shape}")
print(f"Number of unique bonds: {len(Bond_index)}")
print(f"Angle_index shape: {Angle_index.shape}")
print(f"Number of angles: {len(Angle_index)}")

# Print a sample of the Angle_index array to verify format
if len(Angle_index) > 0:
    print("\nSample angle data (first row of Angle_index):")
    print(f"  Atoms: {int(Angle_index[0,0])}-{int(Angle_index[0,1])}-{int(Angle_index[0,2])}")
    print(f"  Angle: {Angle_index[0,3]:.2f}°")
    print(f"  Vector 1->2: [{Angle_index[0,4]:.3f}, {Angle_index[0,5]:.3f}, {Angle_index[0,6]:.3f}]")
    print(f"  Vector 2->3: [{Angle_index[0,7]:.3f}, {Angle_index[0,8]:.3f}, {Angle_index[0,9]:.3f}]")

end_time = time.time()
print(f"\nBond calculation completed in {end_time - start_time:.2f} seconds for {sample_size} atoms")

# Analyze the results
analyze_bond_angles(sample_atoms)

# Print some specific examples
print("\nExamples of atoms with their bonds and angles:")
for i, atom in enumerate(sample_atoms[:5]):
    print(f"\nAtom {i}: {atom.get('element', '?')}")
    print(f"  Neighbors: {len(atom.get('neigh', []))}")
    print(f"  Bonds: {len(atom.get('bonds', []))}")
    if atom.get('bonds', []):
        print(f"    First bond: to atom {atom['bonds'][0][0]} with length {atom['bonds'][0][1]:.3f} Å")
    print(f"  Angles: {len(atom.get('angles', []))}")
    if atom.get('angles', []):
        print(f"    First angle: between atoms {atom['angles'][0][0][0]} and {atom['angles'][0][0][1]} with value {atom['angles'][0][1]:.2f}°")
