import os
import numpy as np
from datetime import datetime
from .bond_angle import bond_angle


def write_itp(atoms, file_path, Box_dim=None, molecule_name=None, nrexcl=1, comment=None, 
          rmaxH=1.2, rmaxM=2.45, explicit_bonds=0, explicit_angles=1, KANGLE=500):
    """
    Write atoms to a Gromacs molecular topology (.itp) file.
    
    This function takes a list of atom dictionaries from atomipy and outputs a formatted
    Gromacs topology file containing atom, bond, and angle definitions similar to the 
    MATLAB write_minff_itp.m script.
    
    Args:
        atoms: List of atom dictionaries.
        file_path: Output file path for the .itp file.
        Box_dim: Box dimensions (required for bond/angle calculations if not already in atoms).
        molecule_name: Name of the molecule (default: derived from atoms[0].resname).
        nrexcl: Number of exclusions (default: 3).
        comment: Optional comment to include in the header.
        rmaxH: Maximum bond distance for hydrogen bonds (default: 1.25 Å).
        rmaxM: Maximum bond distance for non-hydrogen bonds (default: 2.45 Å).
        explicit_bonds: If 1, include bond parameters in the file (default: 0).
        explicit_angles: If 1, include angle parameters in the file (default: 1).
        KANGLE: Force constant for generic angles in kJ/(mol·rad²) (default: 500).
        
    Returns:
        None
    
    Example:
        write_itp(atoms, "molecule.itp", Box_dim=[50, 50, 50], molecule_name="MMT")
    """
    # If file doesn't have .itp extension, add it
    if not file_path.endswith('.itp'):
        file_path = file_path + '.itp'
    
    nAtoms = len(atoms)
    
    # Get unique atom types
    atom_types = set(atom.get('type', '') for atom in atoms)
    atom_types = [a for a in atom_types if a]  # Filter out empty strings
    
    # Check if atoms have masses, if not add default values
    for atom in atoms:
        if 'mass' not in atom:
            # Add default mass based on atom type (in reality, you'd want a proper lookup table)
            atom['mass'] = 0.0
    
    # Check if atoms have charge values
    for atom in atoms:
        if 'charge' not in atom:
            atom['charge'] = 0.0
            
    # Extract molecule name if not provided
    if molecule_name is None and atoms and 'resname' in atoms[0]:
        molecule_name = atoms[0]['resname']
    elif molecule_name is None:
        molecule_name = "MOL"
    
    # Use bond_angle function to calculate bonds and angles
    if Box_dim is None:
        raise ValueError("Box_dim is required to calculate bonds and angles using bond_angle function")
    
    # Add debug output for atom coordinates and box dimensions
    print(f"write_itp: Using box dimensions: {Box_dim}")
    print(f"write_itp: Sample atoms (first 3 with coordinates):")
    for i, atom in enumerate(atoms[:3]):
        print(f"  Atom {i+1}: ({atom.get('x', 'None')}, {atom.get('y', 'None')}, {atom.get('z', 'None')}), Element: {atom.get('element', 'None')}, Type: {atom.get('type', 'None')}")
    
    # Make sure Box_dim is correctly formatted
    if isinstance(Box_dim, list):
        Box_dim = np.array(Box_dim, dtype=float)
    
    # Call the bond_angle function with the provided rmaxH and rmaxM parameters
    # Note: bond_angle function expects coordinates in Angstroms
    # Important: Use same_molecule_only=False to allow bonds between different molecules
    # Keep same_element_bonds=False (default) which is correct for mineral structures
    print(f"write_itp: Calling bond_angle with rmaxH={rmaxH}, rmaxM={rmaxM}, same_molecule_only=False")
    updated_atoms, Bond_index, Angle_index = bond_angle(atoms, Box_dim, rmaxH=rmaxH, rmaxM=rmaxM, same_molecule_only=False)
    print(f"write_itp: bond_angle found {len(Bond_index)} bonds and {len(Angle_index)} angles")
    
    # Convert bond and angle indices to 1-based if they're not already
    # Each bond is [atom1_idx, atom2_idx, distance]
    # Each angle is [atom1_idx, atom2_idx, atom3_idx, angle_value]
    
    # Check if Bond_index is not empty and contains 0-based indices
    if isinstance(Bond_index, np.ndarray) and Bond_index.size > 0:
        # For numpy arrays, check the first element's first value
        if np.min(Bond_index[:, 0]) == 0:  # Check if minimum index is 0 (0-based)
            Bond_index = np.array([
                [int(i)+1, int(j)+1, dist] for i, j, dist in Bond_index
            ])
            print("Converted Bond_index from numpy array to 1-based indexing")
    elif Bond_index and len(Bond_index) > 0:
        # For Python lists
        if min(int(bond[0]) for bond in Bond_index) == 0:
            Bond_index = [[int(i)+1, int(j)+1, dist] for i, j, dist in Bond_index]
            print("Converted Bond_index from list to 1-based indexing")
    
    # Similar check for Angle_index
    if isinstance(Angle_index, np.ndarray) and Angle_index.size > 0:
        # For numpy arrays, check the first element's first value
        if np.min(Angle_index[:, 0]) == 0:  # Check if minimum index is 0 (0-based)
            # Convert all indices to 1-based, preserving other data
            # First determine shape to handle correctly
            cols = Angle_index.shape[1]
            new_angles = []
            for angle in Angle_index:
                # First 3 elements are always atom indices
                new_angle = [int(angle[0])+1, int(angle[1])+1, int(angle[2])+1]
                # Add any remaining values unchanged (e.g., angle value)
                if cols > 3:
                    new_angle.extend(angle[3:])
                new_angles.append(new_angle)
            Angle_index = np.array(new_angles)
            print(f"Converted Angle_index from numpy array to 1-based indexing (shape: {Angle_index.shape})")
    elif Angle_index and len(Angle_index) > 0:
        # For Python lists
        if min(int(angle[0]) for angle in Angle_index) == 0:
            new_angles = []
            for angle in Angle_index:
                # First 3 elements are always atom indices
                new_angle = [int(angle[0])+1, int(angle[1])+1, int(angle[2])+1]
                # Add any remaining values unchanged (e.g., angle value)
                if len(angle) > 3:
                    new_angle.extend(angle[3:])
                new_angles.append(new_angle)
            Angle_index = new_angles
            print("Converted Angle_index from list to 1-based indexing")
    
    # Find atom indices for special types (similar to MATLAB script)
    ind_H = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('H')]
    ind_O = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('O')]
    ind_Al = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Al')]
    ind_Si = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Si')]
    ind_Mgo = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Mg')]
    
    # Filter Bond_index to only include bonds with at least one hydrogen atom
    if Bond_index is not None and len(Bond_index) > 0:
        total_bonds = len(Bond_index)
        Bond_index = [bond for bond in Bond_index if int(bond[0]) in ind_H or int(bond[1]) in ind_H]
        print(f"write_itp: Filtered to {len(Bond_index)} hydrogen bonds (from {total_bonds} total bonds)")
    
    # Calculate total charge
    total_charge = sum(atom.get('charge', 0.0) for atom in atoms)
    total_charge = round(total_charge, 6)
        
    # Open the file for writing
    with open(file_path, 'w') as f:
        # Write header
        f.write("; Gromacs topology file\n")
        f.write(f"; File generated by atomipy on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        # Add custom comment if provided
        if comment:
            f.write(f"; {comment}\n")
        
        # Add info about bonds and angles if available
        if Bond_index is not None and Angle_index is not None:
            f.write(f"; Structure with {nAtoms} atoms, {len(Bond_index)} bonds, {len(Angle_index)} angles\n")
            print(f"write_itp: Writing {len(Bond_index)} bonds and {len(Angle_index)} angles to file")
        
        f.write("\n")
        
        # Write moleculetype section
        f.write("[ moleculetype ]\n")
        f.write("; molname   nrexcl\n")
        
        # Use first 3 chars of molecule name if possible
        mol_name_short = molecule_name[:3] if len(molecule_name) >= 3 else molecule_name
        f.write(f"{mol_name_short.upper()}         {nrexcl}\n\n")
        
        # Write atoms section
        f.write("[ atoms ]\n")
        f.write("; id   attype  resnr resname  atname   cgnr        charge      mass\n")
        
        for i, atom in enumerate(atoms, 1):
            # Get values with defaults for missing fields
            at_type = atom.get('fftype', atom.get('type', ''))
            if at_type is None:
                at_type = 'X'  # Default type if none exists
                
            res_nr = atom.get('molid', atom.get('resid', 1))
            if res_nr is None:
                res_nr = 1
                
            res_name = atom.get('resname', molecule_name)
            if res_name is None:
                res_name = 'UNK'
            else:
                res_name = res_name[:3].upper()
                
            at_name = atom.get('type', '')
            if at_name is None:
                at_name = at_type
                
            charge = round(atom.get('charge', 0.0), 6)
            mass = round(atom.get('mass', 0.0), 6)
            
            # Write the atom line
            f.write(f"{i:<7} {at_type:<7} {res_nr:<7} {res_name:<7} {at_name:<7} {i:<7}  {charge:>10.6f}    {mass:>7.4f}\n")
        
        # Write bonds section if we have bonds
        if Bond_index is not None and len(Bond_index) > 0:
            f.write("\n[ bonds ]\n")
            if explicit_bonds == 1:
                f.write("; i     j       funct   length  force.c.\n")
            else:
                f.write("; i     j       funct\n")
            
            # Sort bonds by first atom index (they're already filtered to H-bonds)
            Bond_index = sorted(Bond_index, key=lambda x: x[0])
            print(f"write_itp: First 3 H-bonds: {Bond_index[:3] if len(Bond_index) >= 3 else Bond_index}")
            
            for bond in Bond_index:
                a1, a2, dist = int(bond[0]), int(bond[1]), float(bond[2])
                
                if explicit_bonds == 1:
                    # H-O bonds use specific parameters
                    r = 0.09572  # H bond length in nm (standard value for OH bonds)
                    kb = 441050  # Force constant
                    
                    # Write bond with parameters
                    at1_type = atoms[a1-1].get('fftype', atoms[a1-1].get('type', ''))
                    at2_type = atoms[a2-1].get('fftype', atoms[a2-1].get('type', ''))
                    f.write(f"{a1:<5} {a2:<5} {1:<5} {r:<8.4f} {kb:<8.4f} ; {at1_type}-{at2_type}\n")
                else:
                    # Write bond without parameters
                    at1_type = atoms[a1-1].get('fftype', atoms[a1-1].get('type', ''))
                    at2_type = atoms[a2-1].get('fftype', atoms[a2-1].get('type', ''))
                    f.write(f"{a1:<5} {a2:<5} {1:<5} ; {dist/10:<8.4f} {at1_type}-{at2_type}\n")
        
        # Write angles section if we have angles
        if Angle_index is not None and len(Angle_index) > 0:
            f.write("\n[ angles ]\n")
            if explicit_angles == 1:
                f.write("; i    j   k   type   theta   force.c.\n")
            else:
                f.write("; i    j   k   type\n")
            
            # Sort angles by middle atom index
            Angle_index = sorted(Angle_index, key=lambda x: x[1])
            
            for angle in Angle_index:
                a1, a2, a3 = int(angle[0]), int(angle[1]), int(angle[2])
                angle_val = float(angle[3]) if len(angle) > 3 else 0.0
                
                if explicit_angles == 1:
                    # Determine angle parameters based on atom types
                    h_count = sum(1 for a in [a1, a2, a3] if a in ind_H)
                    
                    if h_count == 1:
                        if any(a in ind_Mgo for a in [a1, a2, a3]):
                            adeg = 110.0
                            ktheta = 50.208
                        elif any(a in ind_Al for a in [a1, a2, a3]):
                            adeg = 110.0
                            ktheta = 125.52
                        else:
                            adeg = 110.0
                            ktheta = 125.52
                    elif h_count == 2:
                        adeg = 109.47  # SPC water
                        ktheta = 383.0
                    else:
                        adeg = angle_val
                        ktheta = KANGLE
                    
                    # Write angle with parameters
                    at1_type = atoms[a1-1].get('type', '')
                    at2_type = atoms[a2-1].get('type', '')
                    at3_type = atoms[a3-1].get('type', '')
                    f.write(f"{a1:<5} {a2:<5} {a3:<5} {1:<5} {adeg:<6.2f}   {ktheta:<8.2f} ; {at1_type}-{at2_type}-{at3_type}\n")
                else:
                    # Write angle without parameters
                    at1_type = atoms[a1-1].get('fftype', atoms[a1-1].get('type', ''))
                    at2_type = atoms[a2-1].get('fftype', atoms[a2-1].get('type', ''))
                    at3_type = atoms[a3-1].get('fftype', atoms[a3-1].get('type', ''))
                    f.write(f"{a1:<5} {a2:<5} {a3:<5} {1:<5} ; {angle_val:<6.2f} {at1_type}-{at2_type}-{at3_type}\n")
        
        # Write position restraints section
        f.write("\n#ifdef POSRES  \n")
        f.write("[ position_restraints ]\n")
        f.write("; atom  type      fx      fy      fz\n")
        
        for i, atom in enumerate(atoms, 1):
            if i in ind_Al or i in ind_Si or i in ind_Mgo:  # Equivalent to ind_Oct in MATLAB
                f.write(f"{i:<6}\t{1:<6}\t{100:<6}\t{100:<6}\t{10000:<6}\n")
        
        f.write("#endif\n\n")
        
        # Add POSRES for specific force field if needed
        f.write("#ifdef POSRES_MINFF \n")
        f.write("[ position_restraints ]\n")
        f.write("; atom  type      fx      fy      fz\n")
        
        for i, atom in enumerate(atoms, 1):
            at_type = atom.get('type', '')
            # Include all atoms except hydrogen
            if not at_type.startswith('H'):
                f.write(f"{i:<6}\t{1:<6}\t{1000:<6}\t{1000:<6}\t{1000:<6}\n")
        
        f.write("#endif\n")
        
        # All sections (moleculetype, atoms, bonds, angles, and position restraints) are complete
        
        f.write("\n")
        
        # Check if bonds are defined in the atoms
        has_bonds = any('bonds' in atom and atom['bonds'] for atom in atoms)
        
        if has_bonds:
            # Write bonds section if bonds are defined
            f.write("[ bonds ]\n")
            f.write("; i     j       funct   length  force.c.\n")
            
            # This section is handled in the main bond processing code above
            
            f.write("\n")
        
        # We'll skip processing of 'angles' attribute on atoms - this is handled by Angle_index
        # Getting angles from the Angle_index is more reliable and consistent


def from_atom_types(atom_types, charges, masses, file_path, molecule_name='MOL', nrexcl=1, comment=None):
    """
    Write a Gromacs topology file directly from atom types, charges, and masses.
    
    This is a convenience function that doesn't require full atom dictionaries.
    
    Args:
        atom_types: List of atom type strings.
        charges: List of charges corresponding to atom types.
        masses: List of masses corresponding to atom types.
        file_path: Output file path for the .itp file.
        molecule_name: Name of the molecule (default: 'MOL').
        nrexcl: Number of exclusions (default: 1).
        comment: Optional comment to include in the header.
        
    Returns:
        None
        
    Example:
        ap.write_itp.from_atom_types(
            ['Al', 'Si', 'O', 'H'],
            [1.782, 1.884, -1.065, 0.4],
            [26.98, 28.09, 16.0, 1.01],
            "simple.itp",
            molecule_name="MIN"
        )
    """
    # Create simple atom dictionaries
    atoms = []
    for i, (atom_type, charge, mass) in enumerate(zip(atom_types, charges, masses), start=1):
        atoms.append({
            'index': i,
            'fftype': atom_type,
            'atname': atom_type,
            'resname': molecule_name,
            'resnr': 1,
            'charge': charge,
            'mass': mass,
            'cgnr': i
        })
    
    # Call the main write function
    write(atoms, file_path, molecule_name, nrexcl, comment)
