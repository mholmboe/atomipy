import numpy as np


def import_atom_gro(file_path):
    """Import atoms from a Gromacs .gro file.

    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, vx, vy, vz, neigh, bonds, angles, element, type, fftype.
       Box_dim: a 1x9 list representing the triclinic cell dimensions.
    """
    atoms = []
    Box_dim = None
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # First line is title, second line is number of atoms, last line is box dimensions
    num_atoms = int(lines[1].strip())
    atom_lines = lines[2:2+num_atoms]
    box_line = lines[2+num_atoms].strip()

    # Parse atom lines
    for line in atom_lines:
        # GRO format: residue number (0-5), residue name (5-10), atom name (10-15), atom number (15-20), 
        # x (20-28), y (28-36), z (36-44), and optionally vx (44-52), vy (52-60), vz (60-68)
        try:
            index = int(line[15:20].strip())
            resname = line[5:10].strip()
            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())
        except Exception as e:
            continue
        # Check if velocities are present (line length >= 68 characters)
        if len(line) >= 68:
            try:
                vx = float(line[44:52].strip())
                vy = float(line[52:60].strip())
                vz = float(line[60:68].strip())
            except Exception as e:
                vx, vy, vz = None, None, None
        else:
            vx, vy, vz = None, None, None

        atom = {
            "molid": 1,
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
            "type": None,
            "fftype": None
        }
        atoms.append(atom)

    # Parse box dimensions
    try:
        values = [float(val) for val in box_line.split()]
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
