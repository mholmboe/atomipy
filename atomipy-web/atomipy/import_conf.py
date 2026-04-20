import numpy as np
import os
from . import element as element_module
from .cell_utils import Cell2Box_dim, Box_dim2Cell  # Will use lowercase Box in our code


def pdb(file_path):
    """Import atoms from a PDB file.
    
    
    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, neigh, bonds, angles, element, type, fftype.
       Cell: a 1x6 list [a, b, c, alpha, beta, gamma] if available from CRYST1 record.

    Examples
    --------
    atoms, Cell = pdb("structure.pdb")
    atoms, Cell = pdb(file_path="structure.pdb")
    """
    atoms = []
    Cell = None
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
                    index = int(line[6:11])        # Cols 7-11
                    atname = line[12:16].strip()           # Cols 13-16
                    resname = line[17:20].strip()          # Cols 18-20
                    
                    try:
                        molid = int(line[22:26])     # Cols 23-26 (Residue sequence number as molid)
                    except (ValueError, IndexError):
                        molid = 1 # Default if not present or invalid

                    x = float(line[30:38])           # Cols 31-38
                    y = float(line[38:46])           # Cols 39-46
                    z = float(line[46:54])           # Cols 47-54

                    occupancy = 1.0 # Default occupancy
                    if len(line) >= 60:                      # Check if line is long enough for occupancy
                        try:
                            occupancy = float(line[54:60]) # Cols 55-60
                        except ValueError:
                            pass # Keep default if parsing fails

                    temp_factor = 0.0 # Default temperature factor
                    if len(line) >= 66:                      # Check if line is long enough for temp_factor
                        try:
                            temp_factor = float(line[60:66]) # Cols 61-66
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
                        "type": atname,   # Original atom name from PDB, often used as type
                        "neigh": [],
                        "bonds": [],
                        "angles": [],
                        "fftype": atname   # Original atom name from PDB, often used as type 
                    }
                    atoms.append(atom)
                except Exception as e:
                    # print(f"Warning: Could not parse ATOM/HETATM line: {line.strip()} - Error: {e}")
                    continue # Skip malformed ATOM/HETATM lines
    
    # Save original atom types before element() overwrites them
    original_types = [atom.get('type') for atom in atoms]

    # Now use the element.py function to properly determine elements
    # logic: element() populates the 'element' field but currently overwrites 'type'
    element_module.element(atoms)
    
    # Restore original atom types so that 'type' reflects 'atname' (e.g. 'Alo')
    # while 'element' reflects the chemical element (e.g. 'Al')
    for atom, orig_type in zip(atoms, original_types):
        if orig_type:
            atom['type'] = orig_type
    
    return atoms, Cell


def gro(file_path):
    """Import atoms from a Gromacs .gro file.

    Gromacs .gro files store coordinates in nanometers (nm), but atomipy uses Angstroms (Å).
    This function automatically converts the coordinates and Box dimensions from nm to Å.

    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, vx, vy, vz, neigh, bonds, angles, element, type, fftype.
                Coordinates (x, y, z) are converted to Angstroms.
       Box_dim: 1x3 (orthogonal) or 1x9 (triclinic) Box dimensions in Angstroms.

    Examples
    --------
    atoms, Box_dim = gro("structure.gro")
    atoms, Box_dim = gro(file_path="unitcell.gro")
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
            "fftype": atname,  # Store atom name in type field
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

    if len(values) in (3, 9):
        Box_dim = values
    elif values:
        Box_dim = values  # fall back to whatever was provided

    # Save original atom types before element() overwrites them
    original_types = [atom.get('type') for atom in atoms]

    # Now use the element.py function to properly determine elements
    # logic: element() populates the 'element' field but currently overwrites 'type'
    element_module.element(atoms)
    
    # Restore original atom types so that 'type' reflects 'atname' (e.g. 'Alo')
    # while 'element' reflects the chemical element (e.g. 'Al')
    for atom, orig_type in zip(atoms, original_types):
        if orig_type:
            atom['type'] = orig_type

    # Convert Box to Cell
    return atoms, Box_dim


def xyz(file_path):
    """Import atoms from an XYZ file.
    
    XYZ format has the following structure:
    - First line: number of atoms
    - Second line: comment line, may contain Box_dim or Cell info starting with #
    - Remaining lines: atom entries in format: Element X Y Z
    
    Returns:
       atoms: list of dictionaries, each with keys: molid, index, resname, x, y, z, neigh, bonds, angles, element, type, fftype.
       Cell: a 1x6 list [a, b, c, alpha, beta, gamma] derived from Box_dim or directly from comment.

    Examples
    --------
    atoms, Cell = xyz("structure.xyz")
    atoms, Cell = xyz(file_path="with_box.xyz")
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
            if len(values) in (3, 9):
                Box_dim = [float(val) for val in values]
            elif len(values) == 6:
                # If Cell provided directly, convert to Box_dim for consistency
                Box_dim = Cell2Box_dim([float(val) for val in values])
        except (ValueError, IndexError):
            # If parsing fails, ignore and continue without Box dimensions
            pass
    
    # Parse atom lines (starting from line 3)
    atom_lines = lines[2:2+num_atoms]
    for i, line in enumerate(atom_lines, start=1):
        parts = line.strip().split()
        if len(parts) < 4:
            continue  # Skip invalid lines
        
        atname = parts[0].strip()
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
            "element": None,    # Set element directly from XYZ
            "type": atname,     # Use element as type by default
            "fftype": atname    # Use element as fftype by default
        }
        atoms.append(atom)
    
    # Ensure the element is properly set for all atoms
    element_module.element(atoms)
    
    # Convert Box to Cell
    Cell = None
    if Box_dim is not None:
        Cell = Box_dim2Cell(Box_dim)
    
    return atoms, Cell


def cif(file_path, expand_symmetry=True):
    """Import atoms from a CIF or mmCIF file using GEMMI.

    Handles both small-molecule CIF (minerals, crystals from e.g. COD/ICSD)
    and macromolecular mmCIF. Fractional coordinates are converted to
    Cartesian (Angstroms). By default, symmetry operations are applied to
    generate the full unit cell from the asymmetric unit.

    Parameters
    ----------
    file_path : str
        Path to the .cif or .mmcif file.
    expand_symmetry : bool, optional
        If True (default), apply space-group symmetry operations to expand
        the asymmetric unit into the full unit cell.  Set to False if the
        CIF already contains all atoms in the unit cell.

    Returns
    -------
    atoms : list of dict
        Standard atomipy atom list with keys: molid, index, resname,
        x, y, z, element, type, fftype, neigh, bonds, angles, xfrac,
        yfrac, zfrac, occupancy.
    Cell : list
        [a, b, c, alpha, beta, gamma] in Angstroms / degrees.

    Raises
    ------
    ImportError
        If GEMMI is not installed (``pip install gemmi``).

    Notes
    -----
    - Uncertainty values in parentheses (e.g. ``5.1793(4)``) are
      automatically stripped.
    - For small-molecule CIF the space group is read from
      ``_symmetry_space_group_name_H-M`` or
      ``_space_group_name_H-M_alt``.  If neither is found, ``P 1``
      is assumed (no expansion).
    - For macromolecular mmCIF, GEMMI's ``read_structure`` is used
      which already provides Cartesian coordinates.

    Examples
    --------
    atoms, Cell = cif("Kaolinite.cif")
    atoms, Cell = cif("quartz.cif", expand_symmetry=False)
    """
    try:
        import gemmi
    except ImportError:
        raise ImportError(
            "GEMMI is required for CIF import.  Install it with: "
            "pip install gemmi"
        )

    from .transform import fractional_to_cartesian

    # Helper: strip CIF uncertainty parentheses, e.g. '5.1793(4)' -> 5.1793
    def _cifval(s):
        if s is None or s == '?' or s == '.':
            return None
        s = str(s).split('(')[0].strip()
        return float(s)

    ext = os.path.splitext(file_path)[1].lower()

    # ------------------------------------------------------------------
    # Try macromolecular mmCIF first (GEMMI auto-detects the format)
    # ------------------------------------------------------------------
    try:
        st = gemmi.read_structure(file_path)
        if st and len(st) > 0 and len(st[0]) > 0:
            # Successfully parsed as a macromolecular structure
            model = st[0]

            Cell = [
                st.cell.a, st.cell.b, st.cell.c,
                st.cell.alpha, st.cell.beta, st.cell.gamma
            ]

            atoms = []
            idx = 1
            for chain in model:
                for residue in chain:
                    for ga in residue:
                        atom = {
                            'index':   idx,
                            'molid':   int(residue.seqid.num) if residue.seqid.num else 1,
                            'resname': residue.name.strip(),
                            'type':    ga.name.strip(),
                            'fftype':  ga.name.strip(),
                            'element': ga.element.name if ga.element else '',
                            'x':       ga.pos.x,
                            'y':       ga.pos.y,
                            'z':       ga.pos.z,
                            'occupancy': ga.occ,
                            'neigh': [], 'bonds': [], 'angles': [],
                        }
                        atoms.append(atom)
                        idx += 1

            # Fill in element field if missing
            element_module.element(atoms)
            return atoms, Cell
    except Exception:
        pass  # Not a macromolecular format — fall through to small-molecule

    # ------------------------------------------------------------------
    # Small-molecule / mineral CIF
    # ------------------------------------------------------------------
    doc = gemmi.cif.read_file(file_path)
    block = doc.sole_block()

    # ---- Cell parameters ----
    a     = _cifval(block.find_value('_cell_length_a'))
    b     = _cifval(block.find_value('_cell_length_b'))
    c     = _cifval(block.find_value('_cell_length_c'))
    alpha = _cifval(block.find_value('_cell_angle_alpha'))
    beta  = _cifval(block.find_value('_cell_angle_beta'))
    gamma = _cifval(block.find_value('_cell_angle_gamma'))

    if any(v is None for v in [a, b, c, alpha, beta, gamma]):
        raise ValueError(
            f"CIF file {file_path} is missing one or more _cell_* fields"
        )
    Cell = [a, b, c, alpha, beta, gamma]

    # ---- Atom sites (fractional coordinates) ----
    # Try the standard CIF tags; fall back to alternatives
    label_tag   = '_atom_site_label'
    symbol_tag  = '_atom_site_type_symbol'
    fx_tag      = '_atom_site_fract_x'
    fy_tag      = '_atom_site_fract_y'
    fz_tag      = '_atom_site_fract_z'
    occ_tag     = '_atom_site_occupancy'

    # Build a Table from the block for the tags we need
    tags_wanted = [label_tag, symbol_tag, fx_tag, fy_tag, fz_tag, occ_tag]
    # Not all CIF files have all fields; figure out which exist
    tags_present = [t for t in tags_wanted if block.find_value(t) is not None
                    or block.find([t]) is not None]

    # Use GEMMI's loop search
    loop = block.find(tags_wanted[:5])  # At least label + coords

    raw_atoms = []
    for row in loop:
        label  = str(row[0]).strip() if row[0] else 'X'
        symbol = str(row[1]).strip() if len(row) > 1 and row[1] else label
        # Strip charge suffixes like '2+', '3-', digits from symbol
        clean_symbol = ''.join(c for c in symbol if c.isalpha())
        if not clean_symbol:
            clean_symbol = ''.join(c for c in label if c.isalpha())

        fx = _cifval(row[2]) if len(row) > 2 else 0.0
        fy = _cifval(row[3]) if len(row) > 3 else 0.0
        fz = _cifval(row[4]) if len(row) > 4 else 0.0

        # Try to get occupancy from the loop or as a separate column
        occ = 1.0
        if loop.width() > 5:
            occ_val = _cifval(row[5])
            if occ_val is not None:
                occ = occ_val

        raw_atoms.append({
            'label':   label,
            'element': clean_symbol,
            'xfrac':   fx,
            'yfrac':   fy,
            'zfrac':   fz,
            'occupancy': occ,
        })

    n_asym = len(raw_atoms)

    # ---- Symmetry expansion (asymmetric unit -> full unit cell) ----
    if expand_symmetry:
        # Find space group
        sg_str = (block.find_value('_symmetry_space_group_name_H-M')
                  or block.find_value('_space_group_name_H-M_alt')
                  or block.find_value('_symmetry_space_group_name_Hall'))

        if sg_str and sg_str not in ('?', '.'):
            sg_str = sg_str.strip("'\" ")
            try:
                sg = gemmi.SpaceGroup(sg_str)
            except ValueError:
                print(f"Warning: Could not recognise space group '{sg_str}', "
                      f"assuming P 1 (no symmetry expansion)")
                sg = gemmi.SpaceGroup('P 1')
        else:
            sg = gemmi.SpaceGroup('P 1')

        ops = sg.operations()
        expanded = []
        seen = set()
        for atom in raw_atoms:
            for op in ops:
                fxyz = op.apply_to_xyz(
                    [atom['xfrac'], atom['yfrac'], atom['zfrac']]
                )
                # Wrap into [0, 1)
                fxyz = [f % 1.0 for f in fxyz]
                # De-duplicate within tolerance
                key = (
                    atom['element'],
                    round(fxyz[0], 4),
                    round(fxyz[1], 4),
                    round(fxyz[2], 4),
                )
                if key not in seen:
                    seen.add(key)
                    new = dict(atom)
                    new['xfrac'], new['yfrac'], new['zfrac'] = fxyz
                    expanded.append(new)
        raw_atoms = expanded

    # ---- Build atomipy atom list ----
    atoms = []
    for i, ra in enumerate(raw_atoms, 1):
        atom = {
            'index':     i,
            'molid':     1,
            'resname':   'MIN',
            'type':      ra['label'],
            'fftype':    ra['label'],
            'element':   ra['element'],
            'xfrac':     ra['xfrac'],
            'yfrac':     ra['yfrac'],
            'zfrac':     ra['zfrac'],
            'x':         0.0,
            'y':         0.0,
            'z':         0.0,
            'occupancy': ra.get('occupancy', 1.0),
            'neigh': [], 'bonds': [], 'angles': [],
        }
        atoms.append(atom)

    # Convert fractional -> Cartesian using existing transform module
    fractional_to_cartesian(atoms=atoms, Box=Cell, add_to_atoms=True)

    # Ensure element field is properly set
    original_types = [atom.get('type') for atom in atoms]
    element_module.element(atoms)
    for atom, orig_type in zip(atoms, original_types):
        if orig_type:
            atom['type'] = orig_type

    print(f"Imported {len(atoms)} atoms from CIF "
          f"(asymmetric unit: {n_asym}, expand_symmetry={expand_symmetry}, "
          f"Cell: {a:.3f} {b:.3f} {c:.3f} "
          f"{alpha:.1f} {beta:.1f} {gamma:.1f})")

    return atoms, Cell


def auto(file_path):
    """Automatically detect file format and import atoms.
    
    This function will try to detect whether the file is a PDB, GRO, or XYZ file based on the file extension
    and call the appropriate import function.
    
    Args:
        file_path: Path to the input file (PDB, GRO, or XYZ)
        
    Returns:
        atoms: List of atom dictionaries
        Cell or Box_dim depending on format:
            - pdb/xyz: Cell
            - gro: Box_dim
    """
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == '.pdb':
        return pdb(file_path)
    elif ext == '.gro':
        return gro(file_path)
    elif ext == '.xyz':
        return xyz(file_path)
    elif ext in ('.cif', '.mmcif', '.mcif'):
        return cif(file_path)
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
