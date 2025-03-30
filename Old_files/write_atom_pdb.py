import os


def write_atom_pdb(atoms, cell, file_path):
    """Write atoms and cell dimensions to a PDB file.

    Args:
       atoms: list of atom dictionaries.
       cell: 1x6 list [a, b, c, alpha, beta, gamma] or None.
       file_path: output filepath.
    """
    with open(file_path, 'w') as f:
        if cell is not None:
            # Write CRYST1 record; formatted width according to PDB spec
            a, b, c, alpha, beta, gamma = cell
            f.write(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1\n")
        for atom in atoms:
            # Format the ATOM record: using index, resname and coordinates
            index = atom.get('index', 0)
            resname = atom.get('resname', 'UNK')
            x = atom.get('x', 0.0)
            y = atom.get('y', 0.0)
            z = atom.get('z', 0.0)
            # Use element if available, else guess from resname
            element = atom.get('element') if atom.get('element') is not None else resname[0]
            # Standard PDB format for ATOM
            line = (f"ATOM  {index:5d} {resname:^4s} MOL A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n")
            f.write(line)
        f.write("END\n")
