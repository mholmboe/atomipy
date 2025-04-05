import numpy as np
import os
from . import element as element_module
from .cell_utils import Cell2Box_dim, Box_dim2Cell  # Will use lowercase box_dim in our code

def pdb(file_path):
    """Import atoms from a PDB file.
    
    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, neigh, bonds, angles, element, type, fftype.
       cell: a 1x6 list [a, b, c, alpha, beta, gamma] if available from CRYST1 record.
       box_dim: a 1x3 list for orthogonal cells or a 1x9 list representing the triclinic cell dimensions in Angstroms.
    """
    atoms = []
    cell = None
    box_dim = None
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
                # Extract atom name and type information from columns 12-16
                atom_name = line[12:16].strip()
                atom_type = atom_name.strip() # Store the full atom name as type
                
                # Extract residue name from columns 17-20
                resname = line[17:20].strip()
                
                # Extract element from columns 76-78 if available (standard PDB format)
                element_from_pdb = None
                if len(line) >= 78:
                    element_from_pdb = line[76:78].strip()
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
                    "element": element_from_pdb,  # Will be properly set by element function later
                    "type": atom_type,   # Set the atom type based on name
                    "fftype": None
                }
                atoms.append(atom)
    
    # Now use the element.py function to properly determine elements
    element_module.element(atoms)
    
    # Convert cell to box_dim if cell is available
    box_dim = None
    if cell is not None:
        box_dim = Cell2Box_dim(cell)
    
    return atoms, cell, box_dim


def gro(file_path):
    """Import atoms from a Gromacs .gro file.

    Gromacs .gro files store coordinates in nanometers (nm), but atomipy uses Angstroms (Å).
    This function automatically converts the coordinates and box dimensions from nm to Å.

    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, vx, vy, vz, neigh, bonds, angles, element, type, fftype.
                Coordinates (x, y, z) are converted to Angstroms.
       cell: a 1x6 list [a, b, c, alpha, beta, gamma] derived from box_dim.
       box_dim: a 1x3 list for orthogonal cells or a 1x9 list representing the triclinic cell dimensions in Angstroms.
    """
    atoms = []
    box_dim = None
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
            "type": atname,  # Store atom name in type field
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
        box_dim = [values[0],values[1], values[2]]
    elif len(values) == 9:
        box_dim = values
    else:
        box_dim = values  # arbitrary

    # Convert box_dim to cell
    cell = None
    if box_dim is not None:
        cell = Box_dim2Cell(box_dim)  # Using function name as imported
    
    return atoms, cell, box_dim


def xyz(file_path):
    """Import atoms from an XYZ file.
    
    XYZ format has the following structure:
    - First line: number of atoms
    - Second line: comment line, may contain box_dim or Cell info starting with #
    - Remaining lines: atom entries in format: Element X Y Z
    
    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, neigh, bonds, angles, element, type, fftype.
       cell: a 1x6 list [a, b, c, alpha, beta, gamma] derived from box_dim or directly from comment.
       box_dim: a 1x3 or 1x9 list representing the box dimensions in Angstroms if available from comment line.
    """
    atoms = []
    box_dim = None
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # First line is the number of atoms
    try:
        num_atoms = int(lines[0].strip())
    except (ValueError, IndexError):
        raise ValueError("Invalid XYZ file: First line must contain the number of atoms")
    
    # Second line is a comment line, which may contain box dimensions
    comment_line = lines[1].strip()
    if comment_line.startswith('#'):
        # Try to extract box dimensions from comment line
        try:
            # Parse out the values after the # symbol
            values = comment_line.split('#')[1].strip().split()
            if len(values) == 9:
                # Assuming format is a 1x9 box_dim array
                box_dim = [float(val) for val in values]
            elif len(values) == 6:
                # Assuming format is a 1x6 Cell array [a, b, c, alpha, beta, gamma]
                box_dim = [float(val) for val in values]
            elif len(values) == 3:
                # Assuming format is a 1x3 box_dim array for orthogonal box
                box_dim = [float(val) for val in values]
        except (ValueError, IndexError):
            # If parsing fails, ignore and continue without box dimensions
            pass
    
    # Parse atom lines (starting from line 3)
    atom_lines = lines[2:2+num_atoms]
    for i, line in enumerate(atom_lines, start=1):
        parts = line.strip().split()
        if len(parts) < 4:
            continue  # Skip invalid lines
        
        element = parts[0].strip()
        try:
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
        except (ValueError, IndexError):
            continue  # Skip invalid lines
        
        atom = {
            "molid": 1,  # Default molecule ID
            "index": i,  # Use sequential numbering
            "resname": "UNK",  # Default residue name
            "x": x,
            "y": y,
            "z": z,
            "neigh": [],
            "bonds": [],
            "angles": [],
            "element": element,  # Set element directly from XYZ
            "type": element,     # Use element as type by default
            "fftype": None
        }
        atoms.append(atom)
    
    # Ensure the element is properly set for all atoms
    element_module.element(atoms)
    
    # Convert box_dim to cell
    cell = None
    if box_dim is not None:
        cell = Box_dim2Cell(box_dim)  # Using function name as imported
    
    return atoms, cell, box_dim


def auto(file_path):
    """Automatically detect file format and import atoms.
    
    This function will try to detect whether the file is a PDB, GRO, or XYZ file based on the file extension
    and call the appropriate import function.
    
    Args:
        file_path: Path to the input file (PDB, GRO, or XYZ)
        
    Returns:
        atoms: List of atom dictionaries
        cell: A 1x6 list [a, b, c, alpha, beta, gamma] representing cell parameters
        box_dim: A 1x3 list for orthogonal cells or a 1x9 list for triclinic cells
    """
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == '.pdb':
        return pdb(file_path)
    elif ext == '.gro':
        return gro(file_path)
    elif ext == '.xyz':
        return xyz(file_path)
    else:
        # Try to detect the format by checking file contents
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('REMARK') or 'PDB' in first_line:
                return pdb(file_path)
            else:
                try:
                    # If first line is a number, it's likely an XYZ file
                    int(first_line)
                    return xyz(file_path)
                except ValueError:
                    # Default to GRO if can't determine
                    return gro(file_path)
