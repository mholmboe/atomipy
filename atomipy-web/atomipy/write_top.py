import os
import numpy as np
from datetime import datetime
from .bond_angle import bond_angle


def itp(atoms, Box_dim, file_path=None, molecule_name=None, nrexcl=1, comment=None, 
          rmaxH=1.2, rmaxM=2.45, explicit_bonds=0, explicit_angles=1, KANGLE=500):
    """
    Write atoms to a Gromacs molecular topology (.itp) file.
    
    This function takes a list of atom dictionaries from atomipy and outputs a formatted
    Gromacs topology file containing atom, bond, and angle definitions similar to the 
    MATLAB write_minff_itp.m script.
    
    Args:
        atoms: List of atom dictionaries.
        Box_dim: Box dimensions (required for bond/angle calculations if not already in atoms).
        file_path: Output file path for the .itp file.
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
        itp(atoms, Box_dim=[50, 50, 50], file_path="molecule.itp", molecule_name="MMT")
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
 
    # Make sure Box_dim is correctly formatted
    if isinstance(Box_dim, list):
        Box_dim = np.array(Box_dim, dtype=float)
    
    # Call the bond_angle function with the provided rmaxH and rmaxM parameters
    # Note: bond_angle function expects coordinates in Angstroms
    # Important: Use same_molecule_only=False to allow bonds between different molecules
    # Keep same_element_bonds=False (default) which is correct for mineral structures
    print(f"write_itp: Calling bond_angle with rmaxH={rmaxH}, rmaxM={rmaxM}, same_molecule_only=False")
    updated_atoms, Bond_index, Angle_index = bond_angle(atoms, Box_dim, rmaxH=rmaxH, rmaxM=rmaxM, same_element_bonds=False, same_molecule_only=True)
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
        f.write(f"; Total charge of the system is {total_charge:.5f}\n")
        
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
        
        # We'll skip processing of 'angles' attribute on atoms - this is handled by Angle_index
        # Getting angles from the Angle_index is more reliable and consistent


def psf(atoms, Box_dim, file_path,segid=None, rmaxH=1.2, rmaxM=2.45, comment=None):
    """
    Write atoms to a NAMD/CHARMM PSF topology file.
    
    This function takes a list of atom dictionaries from atomipy and outputs a formatted
    PSF topology file containing atom, bond, and angle definitions for use with 
    NAMD or CHARMM molecular dynamics software.
    
    Args:
        atoms: List of atom dictionaries.
        Box_dim: Box dimensions (required for bond/angle calculations if not already in atoms).
        file_path: Output file path for the .psf file.
        segid: Segment ID to use in the PSF file (default: derived from atoms[0].resname).
        rmaxH: Maximum bond distance for hydrogen bonds (default: 1.2 Å).
        rmaxM: Maximum bond distance for non-hydrogen bonds (default: 2.45 Å).
        comment: Optional comment to include in the header.
        
    Returns:
        None
    
    Example:
        psf(atoms, Box_dim=[50, 50, 50], file_path="molecule.psf", segid="CLAY")
    """
    # If file doesn't have .psf extension, add it
    if not file_path.endswith('.psf'):
        file_path = file_path + '.psf'
    
    nAtoms = len(atoms)
    
    # Extract segment ID if not provided
    if segid is None and atoms and 'resname' in atoms[0]:
        segid = atoms[0]['resname'][:4].upper()  # PSF typically uses 4-char segids
    elif segid is None:
        segid = "CLAY"
    
    # Use bond_angle function to calculate bonds and angles
    if Box_dim is None:
        raise ValueError("Box_dim is required to calculate bonds and angles using bond_angle function")
    
    # Make sure Box_dim is correctly formatted
    if isinstance(Box_dim, list):
        Box_dim = np.array(Box_dim, dtype=float)
    
    # Call the bond_angle function with the provided rmaxH and rmaxM parameters
    print(f"write_psf: Calling bond_angle with rmaxH={rmaxH}, rmaxM={rmaxM}")
    updated_atoms, Bond_index, Angle_index = bond_angle(atoms, Box_dim, rmaxH=rmaxH, rmaxM=rmaxM, same_element_bonds=False, same_molecule_only=True)
    print(f"write_psf: bond_angle found {len(Bond_index)} bonds and {len(Angle_index)} angles")
    
    # PSF uses 0-based indexing for its pointers internally but displays as 1-based
    # Ensure we're working with 0-based indices for internal calculations
    if isinstance(Bond_index, np.ndarray) and Bond_index.size > 0:
        if np.min(Bond_index[:, 0]) == 1:  # Check if minimum index is 1 (1-based)
            Bond_index = np.array([
                [int(i)-1, int(j)-1, dist] for i, j, dist in Bond_index
            ])
    elif Bond_index and len(Bond_index) > 0:
        if min(int(bond[0]) for bond in Bond_index) == 1:
            Bond_index = [[int(i)-1, int(j)-1, dist] for i, j, dist in Bond_index]
    
    # Similar check for Angle_index
    if isinstance(Angle_index, np.ndarray) and Angle_index.size > 0:
        if np.min(Angle_index[:, 0]) == 1:  # Check if minimum index is 1 (1-based)
            cols = Angle_index.shape[1]
            new_angles = []
            for angle in Angle_index:
                new_angle = [int(angle[0])-1, int(angle[1])-1, int(angle[2])-1]
                if cols > 3:
                    new_angle.extend(angle[3:])
                new_angles.append(new_angle)
            Angle_index = np.array(new_angles)
    elif Angle_index and len(Angle_index) > 0:
        if min(int(angle[0]) for angle in Angle_index) == 1:
            new_angles = []
            for angle in Angle_index:
                new_angle = [int(angle[0])-1, int(angle[1])-1, int(angle[2])-1]
                if len(angle) > 3:
                    new_angle.extend(angle[3:])
                new_angles.append(new_angle)
            Angle_index = new_angles
    
        # Find atom indices for special types (similar to MATLAB script)
    # For PSF, since Bond_index was already converted to 0-based, we need 0-based indices for filtering
    ind_H = [i-1 for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('H')]
    ind_O = [i-1 for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('O')]
    ind_Al = [i-1 for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Al')]
    ind_Si = [i-1 for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Si')]
    ind_Mgo = [i-1 for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Mg')]
    
    # Filter Bond_index to only include bonds with at least one hydrogen atom
    if Bond_index is not None and len(Bond_index) > 0:
        total_bonds = len(Bond_index)
        Bond_index = [bond for bond in Bond_index if int(bond[0]) in ind_H or int(bond[1]) in ind_H]
        print(f"write_psf: Filtered to {len(Bond_index)} hydrogen bonds (from {total_bonds} total bonds)")
    
    # Calculate total charge
    total_charge = sum(atom.get('charge', 0.0) for atom in atoms)
    total_charge = round(total_charge, 6)
    
    # Open the file for writing
    with open(file_path, 'w') as f:
        # Write PSF header
        f.write("PSF\n\n")
        f.write(f"       2 !NTITLE\n")
        header_comment = f"REMARKS Generated by atomipy on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        f.write(f" {header_comment}\n")
        f.write(f" REMARKS Total charge of the system is {total_charge:.5f}\n")
        if comment:
            f.write(f" REMARKS {comment}\n")
        f.write("\n")
        
        # Write atoms section
        f.write(f"{nAtoms:8d} !NATOM\n")
        
        for i, atom in enumerate(atoms, 1):
            # Get values with defaults for missing fields
            res_nr = atom.get('molid', atom.get('resid', 1))
            if res_nr is None:
                res_nr = 1
                
            res_name = atom.get('resname', 'UNK')
            if res_name is None:
                res_name = 'UNK'
            else:
                res_name = res_name[:4].upper()  # PSF typically uses 4-char residue names
                
            at_name = atom.get('type', '')
            if at_name is None:
                at_name = 'X'
            
            at_type = atom.get('fftype', atom.get('type', ''))
            if at_type is None:
                at_type = 'X'
                
            charge = atom.get('charge', 0.0)
            mass = atom.get('mass', 0.0)
            
            # PSF format:
            # atom-ID  segment-name  residue-ID  residue-name  atom-name  atom-type  charge  mass  0
            f.write(f"{i:8d} {segid:4s} {res_nr:4d} {res_name:4s} {at_name:4s} {at_type:4s} {charge:10.6f} {mass:10.4f} 0\n")
        
        # Write bonds section
        if Bond_index is not None:
            nbonds = len(Bond_index)
            f.write(f"\n{nbonds:8d} !NBOND: bonds\n")
            
            # PSF format - 8 columns per line, each with 8 characters
            # Bond indices are 1-based in the output file
            line_items = 0
            for bond in Bond_index:
                # Since we're already working with 0-based indices after conversion
                # And PSF requires 1-based indices in the file, add 1 here
                atom1 = int(bond[0]) + 1  # Convert to 1-based for output
                atom2 = int(bond[1]) + 1
                f.write(f"{atom1:8d}{atom2:8d}")
                line_items += 2
                if line_items >= 8:  # 4 bonds per line (8 indices)
                    f.write("\n")
                    line_items = 0
            
            # Add a newline if the last line was incomplete
            if line_items > 0:
                f.write("\n")
        else:
            f.write("\n0 !NBOND: bonds\n")
        
        # Write angles section
        if Angle_index is not None:
            nangles = len(Angle_index)
            f.write(f"\n{nangles:8d} !NTHETA: angles\n")
            
            # PSF format - 9 columns per line, each with 8 characters
            # Angle indices are 1-based in the output file
            line_items = 0
            for angle in Angle_index:
                atom1 = int(angle[0]) + 1  # Convert back to 1-based for output
                atom2 = int(angle[1]) + 1
                atom3 = int(angle[2]) + 1
                f.write(f"{atom1:8d}{atom2:8d}{atom3:8d}")
                line_items += 3
                if line_items >= 9:  # 3 angles per line (9 indices)
                    f.write("\n")
                    line_items = 0
            
            # Add a newline if the last line was incomplete
            if line_items > 0:
                f.write("\n")
        else:
            f.write("\n0 !NTHETA: angles\n")
        
        # Write dihedrals section (empty as we don't calculate them)
        f.write("\n0 !NPHI: dihedrals\n")
        
        # Write impropers section (empty as we don't calculate them)
        f.write("\n0 !NIMPHI: impropers\n")
        
        # Write donors section (empty)
        f.write("\n0 !NDON: donors\n")
        
        # Write acceptors section (empty)
        f.write("\n0 !NACC: acceptors\n")
        
        # Write non-bonded exclusion section (empty)
        f.write("\n0 !NNB\n\n")
        
        # Write zeros for the excluded atoms list (required by PSF format)
        f.write(" ".join("0" for _ in range(min(10, nAtoms))) + "\n\n")
        
        # Write zeros for the 1-4 interaction list (required by PSF format)
        f.write(" ".join("0" for _ in range(min(10, nAtoms))) + "\n\n")
        
        # Write group section (empty)
        f.write("0 !NGRP\n")
        f.write("0 0 0\n\n")
        
    print(f"write_psf: Wrote {nAtoms} atoms, {0 if Bond_index is None else len(Bond_index)} bonds, {0 if Angle_index is None else len(Angle_index)} angles to {file_path}")
    return file_path


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
    write_itp(atoms, file_path, molecule_name, nrexcl, comment)


def lmp(atoms, Box_dim, file_path, forcefield=None, rmaxH=1.2, rmaxM=2.45, comment=None):
    """
    Write a LAMMPS data file (.data) from the atom data.
    
    This function takes a list of atom dictionaries from atomipy and outputs a formatted
    LAMMPS data file containing atom, bond, and angle definitions for use with
    LAMMPS molecular dynamics software using the 'real' units system.
    
    Args:
        atoms: List of atom dictionaries.
        Box_dim: Box dimensions (required for bond/angle calculations if not already in atoms).
        file_path: Output file path for the .data file.
        forcefield: Dictionary of forcefield parameters or None (default: None).
        rmaxH: Maximum bond distance for hydrogen bonds (default: 1.2 Å).
        rmaxM: Maximum bond distance for non-hydrogen bonds (default: 2.45 Å).
        comment: Optional comment to include in the header.
        
    Returns:
        The file path to which the data was written.
    
    Example:
        lmp(atoms, Box_dim=[50, 50, 50], file_path="molecule.data")
    """
    import numpy as np
    from datetime import datetime
    from .bond_angle import bond_angle
    
    # Constants for bond and angle parameters in LAMMPS real units

    bHdist=0.97888                      # OPC3 water model, in Å
    KANGLE_WAT=383 /2/4.184             # Dummy value for the rigid OPC3
    ANGLE_WAT=109.47                    # Angle value for the rigid OPC3
    # Convert Gromacs units [kJ/mol] to LAMMPS real units [kcal/mol]
    kbH = 441050 / (2 * 4.184 * 10**2)  # For H bonds
    kbM = 0 / (2 * 4.184 * 10**2)       # For other bonds
    KANGLE_H = 125.52 / (2 * 4.184)      # For angles with H
    KANGLE_M = 500.00 / (2 * 4.184)       # For other angles
    ANGLE_H = 110                        # Angle value for H-containing angles

    
    # Check if we need to convert file_path to include .data extension
    if not file_path.endswith('.data') and not file_path.endswith('.data'):
        file_path = file_path + '.data'
    
    # Count number of atoms
    nAtoms = len(atoms)
    print(f"write_lmp: Processing {nAtoms} atoms")
    
    # Call the bond_angle function with the provided rmaxH and rmaxM parameters
    print(f"write_lmp: Calling bond_angle with rmaxH={rmaxH}, rmaxM={rmaxM}")
    updated_atoms, Bond_index, Angle_index = bond_angle(atoms, Box_dim, rmaxH=rmaxH, rmaxM=rmaxM, same_element_bonds=False, same_molecule_only=True)
    print(f"write_lmp: bond_angle found {len(Bond_index)} bonds and {len(Angle_index)} angles")
    
    # Ensure Bond_index is using 1-based indexing for consistent processing
    if isinstance(Bond_index, np.ndarray) and Bond_index.size > 0:
        if np.min(Bond_index[:, 0]) == 0:  # Check if minimum index is 0 (0-based)
            Bond_index = np.array([
                [int(i)+1, int(j)+1, dist] for i, j, dist in Bond_index
            ])
    elif Bond_index and len(Bond_index) > 0:
        if min(int(bond[0]) for bond in Bond_index) == 0:
            Bond_index = [[int(i)+1, int(j)+1, dist] for i, j, dist in Bond_index]
    
    # Similar check for Angle_index
    if isinstance(Angle_index, np.ndarray) and Angle_index.size > 0:
        if np.min(Angle_index[:, 0]) == 0:  # Check if minimum index is 0 (0-based)
            cols = Angle_index.shape[1]
            new_angles = []
            for angle in Angle_index:
                new_angle = [int(angle[0])+1, int(angle[1])+1, int(angle[2])+1]
                if cols > 3:
                    new_angle.extend(angle[3:])
                new_angles.append(new_angle)
            Angle_index = np.array(new_angles)
    elif Angle_index and len(Angle_index) > 0:
        if min(int(angle[0]) for angle in Angle_index) == 0:
            new_angles = []
            for angle in Angle_index:
                new_angle = [int(angle[0])+1, int(angle[1])+1, int(angle[2])+1]
                if len(angle) > 3:
                    new_angle.extend(angle[3:])
                new_angles.append(new_angle)
            Angle_index = new_angles
            
    # Find atom indices for special types
    ind_H = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('H')]
    ind_O = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('O')]
    ind_Al = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Al')]
    ind_Si = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Si')]
    ind_Mgo = [i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('Mg')]
    ind_Edge = list(set(ind_H + [i for i in ind_O if atoms[i-1].get('type', '').startswith('Osih')]))
    
    # Filter Bond_index to only include bonds with at least one hydrogen atom
    if Bond_index is not None and len(Bond_index) > 0:
        total_bonds = len(Bond_index)
        Bond_index = [bond for bond in Bond_index if int(bond[0]) in ind_H or int(bond[1]) in ind_H]
        print(f"write_lmp: Filtered to {len(Bond_index)} hydrogen bonds (from {total_bonds} total bonds)")
    
    # Sort bonds by first atom index
    if Bond_index is not None and len(Bond_index) > 0:
        Bond_index = sorted(Bond_index, key=lambda x: x[0])
        
    # Calculate total charge
    total_charge = sum(atom.get('charge', 0.0) for atom in atoms)
    total_charge = round(total_charge, 6)
    print(f"write_lmp: Total charge: {total_charge}")
    
    # Get unique atom types
    atom_types = sorted(list(set(atom.get('type', '') for atom in atoms)))
    natom_types = len(atom_types)
    print(f"write_lmp: Found {natom_types} unique atom types")
    
    # Create mapping from atom type to type index
    atom_type_map = {atype: i+1 for i, atype in enumerate(atom_types)}
    
    # Get masses for each atom type
    masses = {}
    for atype in atom_types:
        # Find first atom with this type to get mass
        for atom in atoms:
            if atom.get('type', '') == atype:
                masses[atype] = atom.get('mass', 0.0)
                break
    
    # Process bonds to get unique bond types
    bond_types = []
    unique_bond_info = []
    bond_distances = {}  # Store distances for each bond type
    if Bond_index is not None and len(Bond_index) > 0:
        # Create string representation of each bond type (only atom types, not distances)
        bond_info = []
        bond_data = []  # Store bond data including distances
        for bond in Bond_index:
            a1, a2 = int(bond[0]), int(bond[1])
            type1 = atoms[a1-1].get('type', '')
            type2 = atoms[a2-1].get('type', '')
            dist = float(bond[2])
            # Order types alphabetically for consistency
            if type1 > type2:
                type1, type2 = type2, type1
            # Store type pair as unique identifier
            type_pair = f"{type1} {type2}"
            bond_info.append(type_pair)
            # Store all information for later use
            bond_data.append((type_pair, dist))
        
        # Get unique bond types (based only on atom type pairs)
        seen = set()
        for type_pair in bond_info:
            if type_pair not in seen:
                seen.add(type_pair)
                unique_bond_info.append(type_pair)
                bond_distances[type_pair] = []  # Initialize list for distances
        
        # Collect distances for each bond type
        for type_pair, dist in bond_data:
            bond_distances[type_pair].append(dist)
        
        # Map each bond to its type
        bond_type_map = {info: i+1 for i, info in enumerate(unique_bond_info)}
        for i, bond in enumerate(Bond_index):
            a1, a2 = int(bond[0]), int(bond[1])
            type1 = atoms[a1-1].get('type', '')
            type2 = atoms[a2-1].get('type', '')
            if type1 > type2:
                type1, type2 = type2, type1
            type_pair = f"{type1} {type2}"
            bond_types.append(bond_type_map[type_pair])
    
    # Calculate bond parameters
    bond_params = []
    for i, type_pair in enumerate(unique_bond_info):
        type1, type2 = type_pair.split()
        # Calculate average distance for this bond type
        avg_dist = sum(bond_distances[type_pair]) / len(bond_distances[type_pair])
        # Assign bond constants
        if avg_dist < 1.25:  # H-bond
            k = kbH
            r0 = bHdist
        else:  # Non-H bond
            k = kbM
            r0 = avg_dist
        bond_params.append((i+1, k, r0, f"{type1}-{type2}"))
    
    # Process angles to get unique angle types
    angle_types = []
    unique_angle_info = []
    angle_values = {}  # Store angle values for each angle type
    if Angle_index is not None and len(Angle_index) > 0:
        # Create string representation of each angle type (only atom types, not angles)
        angle_info = []
        angle_data = []  # Store angle data including values
        for angle in Angle_index:
            a1, a2, a3 = int(angle[0]), int(angle[1]), int(angle[2])
            type1 = atoms[a1-1].get('type', '')
            type2 = atoms[a2-1].get('type', '')
            type3 = atoms[a3-1].get('type', '')
            angle_val = float(angle[3]) if len(angle) > 3 else 0.0
            
            # Order first and third atom types alphabetically for consistency
            # (keeping the middle/central atom as the second one)
            if type1 > type3:
                type1, type3 = type3, type1
            
            # Store type triplet as unique identifier
            type_triplet = f"{type1} {type2} {type3}"
            angle_info.append(type_triplet)
            # Store all information for later use
            angle_data.append((type_triplet, angle_val))
        
        # Get unique angle types (based only on atom type triplets)
        seen = set()
        for type_triplet in angle_info:
            if type_triplet not in seen:
                seen.add(type_triplet)
                unique_angle_info.append(type_triplet)
                angle_values[type_triplet] = []  # Initialize list for angle values
        
        # Collect angle values for each angle type
        for type_triplet, angle_val in angle_data:
            angle_values[type_triplet].append(angle_val)
        
        # Map each angle to its type
        angle_type_map = {info: i+1 for i, info in enumerate(unique_angle_info)}
        for i, angle in enumerate(Angle_index):
            a1, a2, a3 = int(angle[0]), int(angle[1]), int(angle[2])
            type1 = atoms[a1-1].get('type', '')
            type2 = atoms[a2-1].get('type', '')
            type3 = atoms[a3-1].get('type', '')
            if type1 > type3:
                type1, type3 = type3, type1
            type_triplet = f"{type1} {type2} {type3}"
            angle_types.append(angle_type_map[type_triplet])
    
    # Calculate angle parameters
    angle_params = []
    for i, type_triplet in enumerate(unique_angle_info):
        type1, type2, type3 = type_triplet.split()
        # Check if any of the atoms is hydrogen
        has_h = 'H' in type1 or 'H' in type2 or 'H' in type3
        
        # Calculate average angle for this angle type if we have values
        if angle_values[type_triplet] and any(val > 0 for val in angle_values[type_triplet]):
            # Filter out zero values which might indicate missing angle data
            valid_angles = [val for val in angle_values[type_triplet] if val > 0]
            avg_angle = sum(valid_angles) / len(valid_angles) if valid_angles else 109.5
        else:
            avg_angle = 109.5  # Default value if no valid angles are found
        
        # Assign angle constants
        # Check if this is a water angle (atom types start with 'ow' and 'hw')
        type1_lower = type1.lower()
        type2_lower = type2.lower()
        type3_lower = type3.lower()
        
        is_water_angle = (type2_lower.startswith('ow') and 
                         (type1_lower.startswith('hw') or type3_lower.startswith('hw')))
        
        if is_water_angle:
            # Special parameters for water angles
            k = KANGLE_WAT
            theta0 = ANGLE_WAT
        elif has_h:
            # Parameters for hydrogen-containing angles
            k = KANGLE_H
            theta0 = ANGLE_H
        else:
            # Parameters for other angles
            k = KANGLE_M
            theta0 = avg_angle  # Use calculated average
        
        angle_params.append((i+1, k, theta0, f"{type1}-{type2}-{type3}"))
    
    # Open the file for writing
    with open(file_path, 'w') as f:
        # Write header
        f.write(f"LAMMPS data file generated by atomipy on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Total charge of the system is {total_charge:.5f}\n")
        
        # Add custom comment if provided
        if comment:
            f.write(f"# {comment}\n")
        f.write("\n")
        
        # Write counts section
        f.write(f"{nAtoms} atoms\n")
        # Safely check if Bond_index exists and has elements
        if Bond_index is not None and (isinstance(Bond_index, np.ndarray) and Bond_index.size > 0 or 
                                    not isinstance(Bond_index, np.ndarray) and len(Bond_index) > 0):
            f.write(f"{len(Bond_index)} bonds\n")
        else:
            f.write("0 bonds\n")
            
        # Safely check if Angle_index exists and has elements
        if Angle_index is not None and (isinstance(Angle_index, np.ndarray) and Angle_index.size > 0 or 
                                      not isinstance(Angle_index, np.ndarray) and len(Angle_index) > 0):
            f.write(f"{len(Angle_index)} angles\n")
        else:
            f.write("0 angles\n")
            
        f.write("0 dihedrals\n")  # Not implemented
        f.write("0 impropers\n\n")  # Not implemented
        
        # Write types section
        f.write(f"{natom_types} atom types\n")
        f.write(f"{len(bond_params)} bond types\n")
        f.write(f"{len(angle_params)} angle types\n")
        f.write("0 dihedral types\n")  # Not implemented
        f.write("0 improper types\n\n")  # Not implemented
        
        # Write box dimensions
        f.write(f"0.0 {Box_dim[0]:.6f} xlo xhi\n")
        f.write(f"0.0 {Box_dim[1]:.6f} ylo yhi\n")
        f.write(f"0.0 {Box_dim[2]:.6f} zlo zhi\n")
        
        # Add tilt factors if box is triclinic (more than 3 dimensions provided)
        # Check for various formats of box dimensions:
        # Standard triclinic box might have 9 elements: [Lx, Ly, Lz, xy, xz, yz]
        # Some formats might include [Lx, Ly, Lz, 0, 0, xy, 0, xz, yz]
        if len(Box_dim) > 5:  # At least has tilt factors
            if len(Box_dim) >= 9:
                # Format with 9 or more elements
                xy = Box_dim[5] if len(Box_dim) == 6 else Box_dim[6]
                xz = Box_dim[6] if len(Box_dim) == 6 else Box_dim[8]
                yz = Box_dim[7] if len(Box_dim) == 6 else Box_dim[9] if len(Box_dim) > 9 else 0.0
            else:
                # Format with 6 elements
                xy = Box_dim[3]
                xz = Box_dim[4]
                yz = Box_dim[5]
                
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n")
        f.write("\n")
        
        # Write masses section
        f.write("Masses\n\n")
        for i, atype in enumerate(atom_types, 1):
            f.write(f"{i} {masses[atype]:.6f}  # {atype}\n")
        f.write("\n")
        
        # Write pair coefficients section if forcefield is provided
        if forcefield is not None:
            f.write("Pair Coeffs\n\n")
            for i, atype in enumerate(atom_types, 1):
                # Get epsilon and sigma from forcefield
                if atype in forcefield:
                    epsilon = forcefield[atype].get('epsilon', 0.0)
                    sigma = forcefield[atype].get('sigma', 0.0)
                else:
                    epsilon = 0.0
                    sigma = 0.0
                f.write(f"{i} {epsilon:.6f} {sigma:.6f}  # {atype}\n")
            f.write("\n")
        
        # Write atoms section
        f.write("Atoms\n\n")
        for i, atom in enumerate(atoms, 1):
            atom_type = atom_type_map[atom.get('type', '')]
            molid = atom.get('molid', atom.get('resid', 1))
            charge = atom.get('charge', 0.0)
            x, y, z = atom.get('x', 0.0), atom.get('y', 0.0), atom.get('z', 0.0)
            f.write(f"{i}\t{molid}\t{atom_type}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}  # {atom.get('type', '')}\n")
        f.write("\n")
        
        # Write bond coefficients if we have bonds
        if Bond_index is not None and len(bond_params) > 0:
            f.write("Bond Coeffs\n\n")
            for bond_id, k, r0, comment in bond_params:
                f.write(f"{bond_id} {k:.6f} {r0:.6f}  # {comment}\n")
            f.write("\n")
            
            # Write bonds section
            f.write("Bonds\n\n")
            for i, bond in enumerate(Bond_index, 1):
                a1, a2 = int(bond[0]), int(bond[1])
                bond_type = bond_types[i-1]
                type1 = atoms[a1-1].get('type', '')
                type2 = atoms[a2-1].get('type', '')
                f.write(f"{i}\t{bond_type}\t{a1}\t{a2}  # {type1}-{type2} {bond[2]:.3f}\n")
            f.write("\n")
        
        # Write angle coefficients if we have angles
        if Angle_index is not None and len(angle_params) > 0:
            f.write("Angle Coeffs\n\n")
            for angle_id, k, theta0, comment in angle_params:
                f.write(f"{angle_id} {k:.6f} {theta0:.6f}  # {comment}\n")
            f.write("\n")
            
            # Write angles section
            f.write("Angles\n\n")
            for i, angle in enumerate(Angle_index, 1):
                a1, a2, a3 = int(angle[0]), int(angle[1]), int(angle[2])
                angle_type = angle_types[i-1]
                type1 = atoms[a1-1].get('type', '')
                type2 = atoms[a2-1].get('type', '')
                type3 = atoms[a3-1].get('type', '')
                angle_val = float(angle[3]) if len(angle) > 3 else 0.0
                f.write(f"{i}\t{angle_type}\t{a1}\t{a2}\t{a3}  # {type1}-{type2}-{type3} {angle_val:.2f}\n")
            f.write("\n")
    
    # Print success message
    print(f"write_lmp: Successfully wrote LAMMPS data file to {file_path}")
    return file_path
