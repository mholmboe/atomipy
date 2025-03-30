import atomipy as ap
import numpy as np
import os

# Import the preem.gro file
print("Importing preem.gro...")
file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'preem.gro')
atoms, box_dim = ap.import_gro(file_path)

# Define element mapping
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

# Assign elements based on atom types
for i in range(len(atoms)):
    atom_type = atoms[i].get('type', '')
    if not atom_type and 'resname' in atoms[i]:
        atom_name = atoms[i].get('atname', '')
        if atom_name:
            atoms[i]['type'] = atom_name
            atom_type = atom_name
    
    if atom_type in atom_type_to_element:
        atoms[i]['element'] = atom_type_to_element[atom_type]
    else:
        # Try standard element guessing as fallback
        atoms[i] = ap.element(atoms[i])

print(f"Successfully imported {len(atoms)} atoms")

# Run the bond_angle function
print("\nRunning bond_angle function...")
atoms, Bond_index, Angle_index = ap.bond_angle(atoms, box_dim)

# Print first 10 lines of Bond_index
print("\nFirst 10 lines of Bond_index:")
print("----------------------------")
print("atom1  atom2  distance")
print("----------------------------")
for i in range(min(10, len(Bond_index))):
    print(f"{int(Bond_index[i,0]):5d} {int(Bond_index[i,1]):5d} {Bond_index[i,2]:.4f}")

# Print first 10 lines of Angle_index
print("\nFirst 10 lines of Angle_index:")
print("------------------------------------------------------------------")
print("atom1  atom2  atom3   angle     dx12     dy12     dz12     dx23     dy23     dz23")
print("------------------------------------------------------------------")
for i in range(min(10, len(Angle_index))):
    print(f"{int(Angle_index[i,0]):5d} {int(Angle_index[i,1]):5d} {int(Angle_index[i,2]):5d} {Angle_index[i,3]:7.2f} " +
          f"{Angle_index[i,4]:8.3f} {Angle_index[i,5]:8.3f} {Angle_index[i,6]:8.3f} " +
          f"{Angle_index[i,7]:8.3f} {Angle_index[i,8]:8.3f} {Angle_index[i,9]:8.3f}")
