import numpy as np

def dist_matrix_atom(atoms, Box_dim):
    """Calculate the full distance matrix between all atoms with periodic boundary conditions for triclinic cells.
    
    Args:
       atoms: list of atom dictionaries, each having 'x', 'y', 'z' coordinates.
       Box_dim: 1x9 list representing the cell vectors in row-major order. If not length 9, a default orthorhombic cell is assumed.
       
    Returns:
       A tuple of four numpy arrays: 
       - A numpy array of shape (N, N) with pairwise distances.
       - Three numpy arrays of shape (N, N) with pairwise x, y, z differences.
    """
    N = len(atoms)
    positions = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])

    # Construct cell matrix
    if Box_dim is not None and len(Box_dim) == 9:
        cell = np.array(Box_dim).reshape((3,3))
    else:
        # Assume orthorhombic if Box_dim not properly provided
        if Box_dim is not None and len(Box_dim) >= 3:
            cell = np.diag(Box_dim[:3])
        else:
            cell = np.eye(3)

    # Compute fractional coordinates for positions
    # Solve cell * frac = pos for each pos, i.e., frac = inv(cell) * pos
    inv_cell = np.linalg.inv(cell)
    pos_frac = positions.dot(inv_cell.T)  # (N, 3)

    # Compute all pairwise differences in fractional coordinates
    diff_frac = pos_frac[:, None, :] - pos_frac[None, :, :]  # shape (N, N, 3)
    # Apply minimum image convention
    diff_frac -= np.round(diff_frac)

    # Convert back to cartesian differences
    diff_cart = diff_frac.dot(cell)  # shape (N, N, 3)

    # Compute the Euclidean distances
    diff_matrix = np.sqrt(np.sum(diff_cart**2, axis=-1))
    dx = diff_cart[:,:,0]
    dy = diff_cart[:,:,1]
    dz = diff_cart[:,:,2]
    return diff_matrix, dx, dy, dz
