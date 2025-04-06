import os

def pdb(atoms, box, file_path):
    """Write atoms and cell dimensions to a PDB file.

    Args:
       atoms: list of atom dictionaries.
       box: a 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
            a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
            a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
            Note that for orthogonal boxes Cell = Box_dim.
       file_path: output filepath.
    """
    # Possibly convert Box_dim into [a,b,c,alpha,beta,gamma] form
    if box is not None:
        if len(box) == 9:
            Box_dim = box 
            # Convert from Box_dim format to Cell format
            Cell = Box_dim2Cell(Box_dim)
        elif len(box) == 6:
            Cell = box
        elif len(box) == 3:
            # Orthogonal box
            Cell = list(box) + [90.0, 90.0, 90.0]
    
    with open(file_path, 'w') as f:
        # Add REMARK lines
        f.write("REMARK    GENERATED BY ATOMIPY\n")
        f.write("REMARK    THIS IS A SIMULATION BOX\n")
        
        if Cell is not None:
            # Write CRYST1 record; formatted width according to PDB spec
            a, b, c, alpha, beta, gamma = Cell
            f.write(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1\n")
        
        # Add MODEL record
        f.write("MODEL        1\n")
        
        for i, atom in enumerate(atoms, start=1):
            # Format the ATOM record according to PDB specification
            index = atom.get('index', i)
            
            # Get atom name from 'type' field, fall back to element
            # PDB atom names (cols 13-16). Rules: up to 4 chars.
            # - Alignment is tricky: usually left-aligned for elements <= 2 chars,
            #   right-aligned (?) or starting col 13 for longer names.
            #   Let's try space-padding and left-alignment for simplicity first.
            raw_atomname = atom.get('type', atom.get('element', 'X'))
            if len(raw_atomname) == 1: # Single char element, place in col 14
                pdb_atomname = f" {raw_atomname}  "
            elif len(raw_atomname) > 4: # Truncate
                pdb_atomname = raw_atomname[:4]
            else: # Left-align others
                pdb_atomname = f"{raw_atomname:<4}"
            
            alt_loc = ' ' # Column 17: Alternate location indicator
            
            # Get residue name: use resname if available, otherwise UNK
            # PDB residue names (cols 18-20). 3 chars, right-aligned.
            raw_resname = atom.get('resname', 'UNK')
            pdb_resname = f"{raw_resname[:3]:>3}" 
            
            chain_id = 'A' # Column 22: Chain identifier (default to 'A')
            res_seq = atom.get('resid', 1) # Column 23-26: Residue sequence number
            icode = ' ' # Column 27: Code for insertion of residues
            
            x = atom.get('x', 0.0)
            y = atom.get('y', 0.0)
            z = atom.get('z', 0.0)
            
            occupancy = 1.00 # Columns 55-60
            temp_factor = 0.00 # Columns 61-66
            
            # Get element symbol (cols 77-78), right-justified
            element_symbol = atom.get('element', 'X')[:2] # Max 2 chars
            element_symbol_pdb = f"{element_symbol:>2}"
            
            charge_str = '  ' # Columns 79-80 (Not handled yet)
            
            # Construct the ATOM line using f-string formatting for precise columns
            f.write(f"ATOM  {index:5d} {pdb_atomname}{alt_loc}{pdb_resname} {chain_id}{res_seq:4d}{icode}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp_factor:6.2f}          {element_symbol_pdb}{charge_str}\n")
        
        # Add ENDMDL record
        f.write("ENDMDL\n")
        
        f.write("END\n")


def gro(atoms, box, file_path):
    """Write atoms and box dimensions to a Gromacs .gro file.

    Gromacs .gro files store coordinates in nanometers (nm), but atomipy uses Angstroms (Å).
    This function automatically converts the coordinates and box dimensions from Å to nm.

    Args:
       atoms: list of atom dictionaries with coordinates in Angstroms.
       box: a 1x6 or 1x9 list representing cell dimensions in Angstroms, either as 
            a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
            a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
            Note that for orthogonal boxes Cell = Box_dim.
       file_path: output filepath.
    """
    # Conversion factor from Angstroms to nm
    angstrom_to_nm = 0.1

    # Possibly convert Box_dim into [a,b,c,alpha,beta,gamma] form
    if box is not None:
        if len(box) == 9:
            Box_dim = box 
        elif len(box) == 6:
            Cell = box
            # Convert from Box_dim format to Cell format
            Box_dim = Cell2Box_dim(Cell)
        elif len(box) == 3:
            # Orthogonal box
            Box_dim = box
    
    with open(file_path, 'w') as f:
        # Write title
        f.write("Generated by atomipy\n")
        # Write number of atoms
        f.write(f"{len(atoms)}\n")
        # Write each atom line with GRO format: residue number (5 chars), residue name (5 chars), atom name (5 chars), atom number (5 chars), x (8.3f), y (8.3f), z (8.3f) and optionally velocities vx, vy, vz (8.4f each)
        for i, atom in enumerate(atoms, start=1):
            # Use molid if available, otherwise default to 1
            resnum = atom.get('molid', 1)  # Use molecule ID as residue number in GRO format
            resname = atom.get('resname', 'UNK')
            
            # Use 'type' field for atom name, fall back to element or first character of resname
            atomname = atom.get('type')
            if atomname is None:
                atomname = atom.get('element') if atom.get('element') is not None else resname[0]
                
            index = atom.get('index', i)
            
            # Convert coordinates from Angstroms to nm for .gro format
            x = atom.get('x', 0.0) * angstrom_to_nm
            y = atom.get('y', 0.0) * angstrom_to_nm
            z = atom.get('z', 0.0) * angstrom_to_nm
            
            # Check if velocities are present and convert them too
            vx = atom.get('vx', None)
            vy = atom.get('vy', None)
            vz = atom.get('vz', None)
            
            if vx is not None and vy is not None and vz is not None:
                # Convert velocities from Å/ps to nm/ps
                vx *= angstrom_to_nm
                vy *= angstrom_to_nm
                vz *= angstrom_to_nm
                line = f"{resnum:5d}{resname:<5s}{atomname:>5s}{index:5d}{x:8.3f}{y:8.3f}{z:8.3f}{vx:8.4f}{vy:8.4f}{vz:8.4f}\n"
            else:
                line = f"{resnum:5d}{resname:<5s}{atomname:>5s}{index:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
                
            f.write(line)
            
        if Box_dim is not None:
            # Convert box dimensions from Angstroms to nm
            box_dim_nm = [val * angstrom_to_nm for val in Box_dim]
            
            # Handle different Box_dim formats
            if len(Box_dim) == 3:  # Orthogonal box with just 3 dimensions
                # For orthogonal boxes, GROMACS expects just 3 values
                box_str = '   '.join(f"{val:.5f}" for val in box_dim_nm[:3])
            elif len(Box_dim) == 9:  # Triclinic box with 9 dimensions
                # For triclinic boxes, GROMACS expects all 9 values
                box_str = '   '.join(f"{val:.5f}" for val in box_dim_nm)
            else:
                # Handle unexpected Box_dim length
                print(f"Warning: Box_dim has unexpected length {len(Box_dim)}. Expected 3 or 9.")
                # Default to using whatever was provided
                box_str = '   '.join(f"{val:.5f}" for val in box_dim_nm)
                
            f.write(box_str + "\n")
        else:
            f.write("\n")


def xyz(atoms, box=None, file_path=None):
    """Write atoms and cell dimensions to an XYZ file.

    XYZ format has the following structure:
    - First line: number of atoms
    - Second line: comment line with box dimensions (if provided). Box dimensions are a 1x6 or 1x9 list representing cell dimensions (in Angstroms), either as 
            a Cell variable having cell parameters array [a, b, c, alpha, beta, gamma], or as 
            a Box_dim variable having box dimensions [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells.
            Note that for orthogonal boxes Cell = Box_dim.
    - Remaining lines: atom entries in format: Element X Y Z

    Args:
       atoms: list of atom dictionaries.
       box: Optional 1x6 list [a, b, c, alpha, beta, gamma] or Box_dim (1x3 or 1x9 list).
           Default is to write Cell parameters if provided.
       file_path: output filepath.
    """
    # Initialize Cell to None in case box is None
    Cell = None
    
    # Possibly convert Box_dim into [a,b,c,alpha,beta,gamma] form
    if box is not None:
        if len(box) == 9:
            Box_dim = box 
            # Convert from Box_dim format to Cell format
            Cell = Box_dim2Cell(Box_dim)
        elif len(box) == 6:
            Cell = box
        elif len(box) == 3:
            # Orthogonal box
            Cell = list(box) + [90.0, 90.0, 90.0]


    with open(file_path, 'w') as f:
        # Write number of atoms on first line
        f.write(f"{len(atoms)}\n")
        
        # Write comment line with box dimensions if available
        if Cell is not None:
            if len(Cell) == 6:  # Cell parameters format
                f.write(f"#    {Cell[0]:.5f}   {Cell[1]:.5f}   {Cell[2]:.5f}   {Cell[3]:.5f}   {Cell[4]:.5f}   {Cell[5]:.5f}\n")
            elif len(Cell) == 3:  # Orthogonal box
                f.write(f"#    {Cell[0]:.5f}   {Cell[1]:.5f}   {Cell[2]:.5f}\n")
            elif len(Cell) == 9:  # Triclinic box
                # Write all 9 values
                box_str = '   '.join(f"{val:.5f}" for val in Cell)
                f.write(f"#    {box_str}\n")
            else:
                # No recognized format, just write a placeholder comment
                f.write("# Generated by atomipy\n")
        else:
            # No cell information, just write a placeholder comment
            f.write("# Generated by atomipy\n")
            
        # Write atom entries
        for atom in atoms:
            # Get element symbol from 'element' field
            element = atom.get('element', 'X')  # Default to 'X' if no element
            
            # Get coordinates
            x = atom.get('x', 0.0)
            y = atom.get('y', 0.0)
            z = atom.get('z', 0.0)
            
            # Write atom entry: Element X Y Z with right-aligned coordinates and extra spacing
            f.write(f"{element:<2}         {x:>12.5f}    {y:>12.5f}    {z:>12.5f}\n")


def auto(atoms, box, file_path):
    """Automatically choose the appropriate write function based on file extension.
    
    This function will analyze the file extension and call either write_pdb, write_gro, or write_xyz
    based on the detected format.
    
    Args:
        atoms: List of atom dictionaries
        box: Either cell parameters for PDB (1x6 list) or Box_dim for GRO/XYZ (1x3 or 1x9 list)
        file_path: Path for the output file
    """
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()
    
    if ext == '.pdb':
        # For PDB, we need 1x6 cell format
        if box is not None and len(box) == 9:
            # Convert Box_dim to cell parameters if needed
            # This is a simplistic conversion and may not be accurate for all cases
            a = box[0]
            b = box[4]
            c = box[8]
            alpha = beta = gamma = 90.0  # Default to orthogonal
            cell = [a, b, c, alpha, beta, gamma]
        else:
            cell = box
        pdb(atoms, cell, file_path)
    elif ext == '.gro':
        # For GRO, we need 1x9 Box_dim format
        if box is not None and len(box) == 6:
            # Convert cell parameters to Box_dim if needed
            # This is a simplistic conversion assuming orthogonal cell
            a, b, c = box[0], box[1], box[2]
            Box_dim = [a, 0, 0, 0, b, 0, 0, 0, c]
        else:
            Box_dim = box
        gro(atoms, Box_dim, file_path)
    elif ext == '.xyz':
        # For XYZ, we'll use the box parameters as provided
        # The xyz function handles different formats internally
        xyz(atoms, box, file_path)
    else:
        # Default to PDB if the extension is unrecognized
        pdb(atoms, box, file_path)
