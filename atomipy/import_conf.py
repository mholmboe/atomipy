import numpy as np
import os
from . import element as element_module
from .cell_utils import Cell2Box_dim, Box_dim2Cell  # Will use lowercase Box_dim in our code

def pdb(file_path):
    """Import atoms from a PDB file.
    
    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, neigh, bonds, angles, element, type, fftype.
       Cell: a 1x6 list [a, b, c, alpha, beta, gamma] if available from CRYST1 record.
       Box_dim: a 1x3 list for orthogonal cells or a 1x9 list representing the triclinic Cell dimensions in Angstroms.
    """
    atoms = []
    Cell = None
    Box_dim = None
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("CRYST1"):
                # Parse Cell parameters from CRYST1 line
                a = float(line[6:15].strip())
                b = float(line[15:24].strip())
                c = float(line[24:33].strip())
                alpha = float(line[33:40].strip())
                beta = float(line[40:47].strip())
                gamma = float(line[47:54].strip())
                Cell = [a, b, c, alpha, beta, gamma]
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                # PDB format column specifications (1-indexed based on documentation)
                # Serial:       7-11
                # Atom name:   13-16 (Atom type in user's description)
                # AltLoc:      17
                # ResName:     18-20
                # ChainID:     22
                # ResSeq:      23-26
                # X:           31-38
                # Y:           39-46
                # Z:           47-54
                # Occupancy:   55-60
                # TempFactor:  61-66
                # Element:     77-78 (right-justified)
                # Charge:      79-80

                try:
                    index = int(line[6:11].strip())          # Cols 7-11
                    atom_name = line[12:16].strip()        # Cols 13-16
                    resname = line[17:20].strip()          # Cols 18-20
                    
                    try:
                        molid = int(line[22:26].strip())     # Cols 23-26 (Residue sequence number as molid)
                    except (ValueError, IndexError):
                        molid = 1 # Default if not present or invalid

                    x = float(line[30:38].strip())           # Cols 31-38
                    y = float(line[38:46].strip())           # Cols 39-46
                    z = float(line[46:54].strip())           # Cols 47-54

                    occupancy = 1.0 # Default occupancy
                    if len(line) >= 60:                      # Check if line is long enough for occupancy
                        try:
                            occupancy_str = line[54:60].strip() # Cols 55-60
                            if occupancy_str: # Ensure not empty before float conversion
                                occupancy = float(occupancy_str)
                        except ValueError:
                            pass # Keep default if parsing fails

                    temp_factor = 0.0 # Default temperature factor
                    if len(line) >= 66:                      # Check if line is long enough for temp_factor
                        try:
                            temp_factor_str = line[60:66].strip() # Cols 61-66
                            if temp_factor_str:
                                temp_factor = float(temp_factor_str)
                        except ValueError:
                            pass # Keep default

                    element_symbol = "" # Default element symbol
                    if len(line) >= 78:                      # Check for element symbol
                        element_symbol = line[76:78].strip().upper() # Cols 77-78, ensure upper for consistency
                    
                    charge_str = "" # Default charge
                    if len(line) >= 80:                      # Check for charge
                        charge_str = line[78:80].strip()     # Cols 79-80

                    atom = {
                        "molid": molid,
                        "index": index,
                        "resname": resname,
                        "x": x,
                        "y": y,
                        "z": z,
                        "occupancy": occupancy,
                        "temp_factor": temp_factor,
                        "element": element_symbol, # Explicitly from PDB cols 77-78
                        "charge": charge_str,
                        "type": atom_name,   # Original atom name from PDB, often used as type
                        "neigh": [],
                        "bonds": [],
                        "angles": [],
                        "fftype": None 
                    }
                    atoms.append(atom)
                except Exception as e:
                    # print(f"Warning: Could not parse ATOM/HETATM line: {line.strip()} - Error: {e}")
                    continue # Skip malformed ATOM/HETATM lines
    
    # Now use the element.py function to properly determine elements
    element_module.element(atoms)
    
    # Convert Cell to Box_dim if Cell is available
    Box_dim = None
    if Cell is not None:
        Box_dim = Cell2Box_dim(Cell)
    
    return atoms, Cell, Box_dim


def gro(file_path):
    """Import atoms from a Gromacs .gro file.

    Gromacs .gro files store coordinates in nanometers (nm), but atomipy uses Angstroms (Å).
    This function automatically converts the coordinates and Box dimensions from nm to Å.

    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, vx, vy, vz, neigh, bonds, angles, element, type, fftype.
                Coordinates (x, y, z) are converted to Angstroms.
       Cell: a 1x6 list [a, b, c, alpha, beta, gamma] derived from Box_dim.
       Box_dim: a 1x3 list for orthogonal cells or a 1x9 list representing the triclinic Cell dimensions in Angstroms.
    """
    atoms = []
    Box_dim = None
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # First line is title, second line is number of atoms, last line is Box dimensions
    num_atoms = int(lines[1].strip())
    atom_lines = lines[2:2+num_atoms]
    Box_line = lines[2+num_atoms].strip()

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

    # Parse Box dimensions
    try:
        values = [float(val) for val in Box_line.split()]
        # Convert Box dimensions from nm to Angstroms
        values = [val * nm_to_angstrom for val in values]
    except Exception as e:
        values = []

    if len(values) == 3:
        # Only Box lengths are provided, assume orthorhombic
        Box_dim = [values[0],values[1], values[2]]
    elif len(values) == 9:
        Box_dim = values
    else:
        Box_dim = values  # arbitrary

    # Convert Box_dim to Cell
    Cell = None
    if Box_dim is not None:
        Cell = Box_dim2Cell(Box_dim)  # Using function name as imported
    
    return atoms, Cell, Box_dim


def xyz(file_path):
    """Import atoms from an XYZ file.
    
    XYZ format has the following structure:
    - First line: number of atoms
    - Second line: comment line, may contain Box_dim or Cell info starting with #
    - Remaining lines: atom entries in format: Element X Y Z
    
    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, neigh, bonds, angles, element, type, fftype.
       Cell: a 1x6 list [a, b, c, alpha, beta, gamma] derived from Box_dim or directly from comment.
       Box_dim: a 1x3 or 1x9 list representing the Box dimensions in Angstroms if available from comment line.
    """
    atoms = []
    Box_dim = None
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # First line is the number of atoms
    try:
        num_atoms = int(lines[0].strip())
    except (ValueError, IndexError):
        raise ValueError("Invalid XYZ file: First line must contain the number of atoms")
    
    # Second line is a comment line, which may contain Box dimensions
    comment_line = lines[1].strip()
    if comment_line.startswith('#'):
        # Try to extract Box dimensions from comment line
        try:
            # Parse out the values after the # symbol
            values = comment_line.split('#')[1].strip().split()
            if len(values) == 9:
                # Assuming format is a 1x9 Box_dim array
                Box_dim = [float(val) for val in values]
            elif len(values) == 6:
                # Assuming format is a 1x6 Cell array [a, b, c, alpha, beta, gamma]
                Box_dim = [float(val) for val in values]
            elif len(values) == 3:
                # Assuming format is a 1x3 Box_dim array for orthogonal Box
                Box_dim = [float(val) for val in values]
        except (ValueError, IndexError):
            # If parsing fails, ignore and continue without Box dimensions
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
    
    # Convert Box_dim to Cell
    Cell = None
    if Box_dim is not None:
        Cell = Box_dim2Cell(Box_dim)  # Using function name as imported
    
    return atoms, Cell, Box_dim


def auto(file_path):
    """Automatically detect file format and import atoms.
    
    This function will try to detect whether the file is a PDB, GRO, or XYZ file based on the file extension
    and call the appropriate import function.
    
    Args:
        file_path: Path to the input file (PDB, GRO, or XYZ)
        
    Returns:
        atoms: List of atom dictionaries
        Cell: A 1x6 list [a, b, c, alpha, beta, gamma] representing Cell parameters
        Box_dim: A 1x3 list for orthogonal cells or a 1x9 list for triclinic cells
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
