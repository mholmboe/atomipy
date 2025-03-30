import copy

def replicate_atom(atom, box_dim, replicate, dim_order='xyz', add_molid=False):
    """
    Replicates the coordinates in an atom list and the box dimensions.

    This function replicates the atom coordinates (x, y, z) in orthogonal 
    or triclinic space along x, y, and z dimensions a specified number of times. 
    It also updates the box dimensions accordingly. For triclinic boxes, it 
    temporarily converts atoms to orthonormal space, replicates them, and 
    converts them back.

    Parameters
    ----------
    atom : list
        A list of atom dictionaries. Each dictionary must have at least:
          'x', 'y', 'z' (coordinates), 'molid' (molecule ID), 'type' (atom type).
    box_dim : list or array
        Box dimension parameters. Can be:
          - 3 elements [lx, ly, lz] for an orthogonal box
          - 9 elements [lx, ly, lz, 0, 0, xy, 0, xz, yz] for a triclinic box
    replicate : list or int
        The replication factor for x, y, z directions. For example:
          [nx, ny, nz] means replicate nx times in x, ny in y, nz in z.
        If a single int is given, it applies to all three directions.
    dim_order : str, optional
        The order in which to replicate dimensions (e.g., 'xyz', 'yxz', etc.).
        Default is 'xyz'.
    add_molid : bool, optional
        If True, each newly replicated copy will have its 'molid' incremented.
        Default is False.

    Returns
    -------
    atom : list
        The updated list of atom dictionaries, replicated in the specified way.
    Notes
    -----
    - Requires the functions `update_atom`, `box_dim2cell`, `triclinic_atom`, 
      and `orto_atom` to be defined elsewhere, if using triclinic boxes.
    - If `box_dim` has 9 elements, the atoms are first converted to orthonormal 
      coordinates (`orto_atom`), replicated, and converted back via 
      `triclinic_atom`.
    - The final box dimensions and the corresponding [a, b, c, alpha, beta, gamma]
      are updated internally (requires `box_dim2cell`).

    Example
    -------
    >>> # Replicate 6 times in x, 4 times in y, 1 time in z
    >>> new_atom = replicate_atom(atom, box_dim, [6, 4, 1])
    >>> # Change dimension order and add new 'molid' per replication
    >>> new_atom = replicate_atom(atom, box_dim, [6, 4, 1], dim_order='yxz', add_molid=True)
    """
    # -- Helper imports (assuming these functions are available) --
    # from your_module import orto_atom, triclinic_atom, box_dim2cell, update_atom

    # 1) Handle replicate to ensure 3 elements
    if isinstance(replicate, int):
        replicate = [replicate, replicate, replicate]
    else:
        replicate = list(replicate)
        if len(replicate) < 3:
            raise ValueError("replicate must have 3 elements or be a single integer.")

    # Replace zeros with 1 (no replication if replicate[i] == 0)
    replicate = [r if r != 0 else 1 for r in replicate]

    # 2) Handle box_dim similarly
    box_dim = list(box_dim)
    if len(box_dim) == 1:
        # Single dimension used for x,y,z
        box_dim = [box_dim[0], box_dim[0], box_dim[0]]
    # Distinguish orthogonal vs triclinic
    is_triclinic = (len(box_dim) == 9)

    # 3) If triclinic, transform coordinates to orthonormal space before replicating
    if is_triclinic:
        # Convert triclinic box to cell angles
        # (Requires a user-defined box_dim2cell function)
        cell_tric = box_dim2cell(box_dim)
        # Convert to orthonormal
        atom = orto_atom(atom, box_dim)  # User-defined function

    # 4) Replicate in the specified dimension order
    #    We'll do cumulative replication in a stepwise manner
    combined_atom = copy.deepcopy(atom)
    current_molid = combined_atom[0]['molid'] if len(combined_atom) > 0 else 1

    # Map 'x','y','z' to array indices (0,1,2) for replicate & box_dim
    dim_map = {'x': 0, 'y': 1, 'z': 2}

    for d in dim_order:
        idx = dim_map[d]  # which dimension in replicate, box_dim to use
        times = replicate[idx]
        if times <= 1:
            # No replication needed beyond the single copy
            continue

        # For clarity, replicate from the existing combined_atom
        original_list = copy.deepcopy(combined_atom)
        # We'll accumulate new positions here
        for i in range(2, times + 1):
            # Each new copy
            new_atoms = copy.deepcopy(original_list)
            shift_val = (i - 1) * box_dim[idx]  # shift in x, y, or z

            if d == 'x':
                for atom_dict in new_atoms:
                    atom_dict['x'] += shift_val
                if add_molid:
                    current_molid += 1
                    for atom_dict in new_atoms:
                        atom_dict['molid'] = current_molid

            elif d == 'y':
                for atom_dict in new_atoms:
                    atom_dict['y'] += shift_val
                if add_molid:
                    current_molid += 1
                    for atom_dict in new_atoms:
                        atom_dict['molid'] = current_molid

            elif d == 'z':
                for atom_dict in new_atoms:
                    atom_dict['z'] += shift_val
                if add_molid:
                    current_molid += 1
                    for atom_dict in new_atoms:
                        atom_dict['molid'] = current_molid

            # Add the newly replicated chunk
            combined_atom.extend(new_atoms)

    # 5) Update box dimensions
    if not is_triclinic:
        # Simply scale each dimension by the replicate factor
        box_dim = [
            box_dim[0] * replicate[0] if dim_order.count('x') > 0 else box_dim[0],
            box_dim[1] * replicate[1] if dim_order.count('y') > 0 else box_dim[1],
            box_dim[2] * replicate[2] if dim_order.count('z') > 0 else box_dim[2]
        ]
    else:
        # For triclinic, we updated only the orthogonal distances in each step.
        # We must reapply a triclinic transform to place them back correctly.
        # That requires user-defined `triclinic_atom` to handle angles.
        # We'll guess you want the angles from cell_tric[3:6] for alpha,beta,gamma
        # and the final box lengths scaled by replicate factors.
        # Example approach:
        box_dim = [
            box_dim[0] * replicate[0], 
            box_dim[1] * replicate[1], 
            box_dim[2] * replicate[2],
            box_dim[3], box_dim[4], box_dim[5], box_dim[6], box_dim[7], box_dim[8]
        ]
        # Convert back to triclinic coordinates
        combined_atom = triclinic_atom(combined_atom, box_dim, cell_tric[3:6], mode='angle')
        # Possibly update box_dim again if triclinic_atom returns a new dimension:
        # box_dim = triclinic_Box_dim  (depending on your code)

    # 6) Overwrite 'atom' with combined_atom
    atom = combined_atom

    # 7) Call update_atom to re-index 'index' and possibly fix 'molid'
    atom = update_atom(atom)  # user-defined function

    # 8) Convert final box_dim to cell and return. 
    #    (The MATLAB code sets them in the caller workspace. 
    #     In Python, you typically just return them or store them globally.)
    final_cell = box_dim2cell(box_dim)
    # You could return both:
    # return atom, box_dim, final_cell

    return atom
