import numpy as np
import os

def pdb(file_path):
    """Import atoms from a PDB file.
    
    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, neigh, bonds, angles, element, type, fftype.
       cell: a 1x6 list [a, b, c, alpha, beta, gamma] if available from CRYST1 record.
    """
    atoms = []
    cell = None
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("CRYST1"):
                # Parse cell parameters from CRYST1 line
                a = float(line[6:15].strip())
                b = float(line[15:24].strip())
                c = float(line[24:33].strip())
                alpha = float(line[33:40].strip())
                beta = float(line[40:47].strip())
                gamma = float(line[47:54].strip())
                cell = [a, b, c, alpha, beta, gamma]
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                index = int(line[6:11].strip())
                # Use atom name from columns 12-16 as a preliminary guess for element
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                # Extract residue sequence number (columns 23-26) to use as molecule ID
                try:
                    # PDB format has residue sequence number in columns 23-26
                    molid = int(line[22:26].strip())
                except (ValueError, IndexError):
                    # Default to 1 if conversion fails
                    molid = 1
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom = {
                    "molid": molid,
                    "index": index,
                    "resname": resname,
                    "x": x,
                    "y": y,
                    "z": z,
                    "neigh": [],
                    "bonds": [],
                    "angles": [],
                    "element": None,
                    "type": None,
                    "fftype": None
                }
                atoms.append(atom)
    return atoms, cell


def gro(file_path):
    """Import atoms from a Gromacs .gro file.

    Gromacs .gro files store coordinates in nanometers (nm), but atomipy uses Angstroms (Å).
    This function automatically converts the coordinates and box dimensions from nm to Å.

    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, vx, vy, vz, neigh, bonds, angles, element, type, fftype.
                Coordinates (x, y, z) are converted to Angstroms.
       Box_dim: a 1x9 list representing the triclinic cell dimensions in Angstroms.
    """
    atoms = []
    Box_dim = None
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # First line is title, second line is number of atoms, last line is box dimensions
    num_atoms = int(lines[1].strip())
    atom_lines = lines[2:2+num_atoms]
    box_line = lines[2+num_atoms].strip()

    # Conversion factor from nm to Angstroms
    nm_to_angstrom = 10.0

    # Parse atom lines
    for line in atom_lines:
        # GRO format: residue number (0-5), residue name (5-10), atom name (10-15), atom number (15-20), 
        # x (20-28), y (28-36), z (36-44), and optionally vx (44-52), vy (52-60), vz (60-68)
        try:
            # Extract residue number (columns 0-5) to use as molecule ID
            molid = int(line[0:5].strip())
            index = int(line[15:20].strip())
            resname = line[5:10].strip()
            # Extract atom name (columns 10-15)
            atname = line[10:15].strip()
            # Convert coordinates from nm to Angstroms
            x = float(line[20:28].strip()) * nm_to_angstrom
            y = float(line[28:36].strip()) * nm_to_angstrom
            z = float(line[36:44].strip()) * nm_to_angstrom
        except Exception as e:
            continue
        # Check if velocities are present (line length >= 68 characters)
        if len(line) >= 68:
            try:
                # Also convert velocities from nm/ps to Å/ps
                vx = float(line[44:52].strip()) * nm_to_angstrom
                vy = float(line[52:60].strip()) * nm_to_angstrom
                vz = float(line[60:68].strip()) * nm_to_angstrom
            except Exception as e:
                vx, vy, vz = None, None, None
        else:
            vx, vy, vz = None, None, None

        atom = {
            "molid": molid,
            "index": index,
            "resname": resname,
            "atname": atname,  # Add atom name
            "x": x,
            "y": y,
            "z": z,
            "vx": vx,
            "vy": vy,
            "vz": vz,
            "neigh": [],
            "bonds": [],
            "angles": [],
            "element": None,
            "type": atname,  # Use atom name as initial type
            "fftype": None,
            "is_nm": False  # Mark that coordinates are now in Angstroms
        }
        atoms.append(atom)

    # Parse box dimensions
    try:
        values = [float(val) for val in box_line.split()]
        # Convert box dimensions from nm to Angstroms
        values = [val * nm_to_angstrom for val in values]
    except Exception as e:
        values = []

    if len(values) == 3:
        # Only box lengths are provided, assume orthorhombic
        Box_dim = [values[0], 0, 0, 0, values[1], 0, 0, 0, values[2]]
    elif len(values) == 9:
        Box_dim = values
    else:
        Box_dim = values  # arbitrary

    return atoms, Box_dim


def auto(file_path):
    """Automatically detect file format and import atoms.
    
    This function will try to detect whether the file is a PDB or GRO file based on the file extension
    and call the appropriate import function.
    
    Args:
        file_path: Path to the input file (PDB or GRO)
        
    Returns:
        atoms: List of atom dictionaries
        box: Either cell (for PDB) or Box_dim (for GRO) depending on the file format
    """
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == '.pdb':
        return pdb(file_path)
    elif ext == '.gro':
        return gro(file_path)
    else:
        # Try to detect the format by checking file contents
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('REMARK') or 'PDB' in first_line:
                return pdb(file_path)
            else:
                # Default to GRO if can't determine
                return gro(file_path)
