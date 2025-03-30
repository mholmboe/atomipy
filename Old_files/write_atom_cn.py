def write_atom_cn(atoms, file_path):
    """Write the coordination number (CN) for each atom to a text file.
    
    For each atom, the CN is defined as the number of bonded nearest neighbours as stored in the 'neigh' field.
    The output file will have one line per atom with the atom's index, residue name, chemical element, and CN.
    
    Args:
       atoms: list of atom dictionaries. Each atom dictionary should have a 'neigh' field.
       file_path: output file path.
    """
    # Update each atom's dictionary with its coordination number (CN)
    for atom in atoms:
        atom['cn'] = len(atom.get('neigh', []))

    with open(file_path, 'w') as f:
        # Write header
        f.write("Index Resname Element CN\n")
        for atom in atoms:
            index = atom.get('index', 0)
            resname = atom.get('resname', 'UNK')
            element = atom.get('element', 'X')
            cn = atom.get('cn', 0)
            f.write(f"{index:5d} {resname:6s} {element:>2s} {cn:3d}\n")
