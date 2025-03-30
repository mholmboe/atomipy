import numpy as np


def import_atom_pdb(file_path):
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
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom = {
                    "molid": 1,
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
