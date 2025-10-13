import numpy as np

from atomipy.cell_utils import Cell2Box_dim, Box_dim2Cell

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    has_tqdm = True
except ImportError:
    print("Note: Install tqdm package for progress bars (pip install tqdm)")
    has_tqdm = False

def dist_matrix(atoms, Box):
    """Calculate the distance matrix between atoms following the MATLAB implementation approach.
    
    This function closely mimics the behavior of the MATLAB dist_matrix_MATLAB.m function,
    calculating distances with periodic boundary conditions in a per-atom loop approach.
    
    Args:
        atoms: list of atom dictionaries, each having 'x', 'y', 'z' coordinates.
        Box: a 1x3, 1x6 or 1x9 list representing Cell dimensions (in Angstroms):
            - For orthogonal boxes, a 1x3 list [lx, ly, lz] where Box = Box_dim, and Cell would be [lx, ly, lz, 90, 90, 90]
            - For Cell parameters, a 1x6 list [a, b, c, alpha, beta, gamma] (Cell format)
            - For triclinic boxes, a 1x9 list [lx, ly, lz, 0, 0, xy, 0, xz, yz] (GROMACS Box_dim format)
       
    Returns:
        A tuple of four numpy arrays: 
        - A numpy array of shape (N, N) with pairwise distances.
        - Three numpy arrays of shape (N, N) with pairwise x, y, z differences.
        
    Note:
        This implementation follows the approach in the MATLAB dist_matrix_MATLAB.m function,
        using per-atom iteration.
    """
    
    # Initialize variables
    Cell = None
    Box_dim = None
    
    if Box is None:
        raise ValueError("Box parameter must be provided")
    
    # Determine Box format and convert as needed
    if len(Box) == 9:
        # Triclinic Box in GROMACS format [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        Box_dim = Box
        Cell = Box_dim2Cell(Box_dim)
    elif len(Box) == 6:
        # Cell parameters [a, b, c, alpha, beta, gamma]
        Cell = Box
        Box_dim = Cell2Box_dim(Cell)
    elif len(Box) == 3:
        # Simple orthogonal Box [lx, ly, lz]
        Box_dim = Box
        Cell = list(Box) + [90.0, 90.0, 90.0]
    else:
        raise ValueError("Box must be length 3, 6, or 9")
    
    # Extract Box dimensions
    if len(Box_dim) == 3:
        # Orthogonal Box
        lx, ly, lz = Box_dim
        xy, xz, yz = 0, 0, 0
    elif len(Box_dim) == 9:
        # Triclinic Box in GROMACS format [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        lx, ly, lz = Box_dim[0], Box_dim[1], Box_dim[2]
        xy, xz, yz = Box_dim[5], Box_dim[7], Box_dim[8]
    
    # Extract atomic positions
    n_atoms = len(atoms)
    xyz = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms], dtype=np.float32)

    
    # Initialize output arrays
    distances = np.zeros((n_atoms, n_atoms), dtype=np.float32)
    dx = np.zeros((n_atoms, n_atoms), dtype=np.float32)
    dy = np.zeros((n_atoms, n_atoms), dtype=np.float32)
    dz = np.zeros((n_atoms, n_atoms), dtype=np.float32)

    # Setup progress tracking
    total_distances_processed = 0
    
    # Calculate distance matrix
    if len(Box_dim) == 3:
        # Orthogonal Box approach
        # Setup progress bar
        if has_tqdm:
            atom_iterator = tqdm(range(n_atoms), desc="Finding dists", unit="atom")
        else:
            print("Finding distances...")
            atom_iterator = range(n_atoms)
            last_percent = -1
            
        for i in atom_iterator:
            # Update progress percentage for non-tqdm case
            if not has_tqdm and n_atoms > 100:
                percent = int(100 * i / n_atoms)
                if percent > last_percent and percent % 10 == 0:
                    print(f"  {percent}% complete...")
                    last_percent = percent
            # Calculate distance components
            rx = xyz[i, 0] - xyz[:, 0]
            ry = xyz[i, 1] - xyz[:, 1]
            rz = xyz[i, 2] - xyz[:, 2]
            
            # Apply minimum image convention for orthogonal Box
            rx[rx > lx/2] -= lx
            rx[rx < -lx/2] += lx
            
            ry[ry > ly/2] -= ly
            ry[ry < -ly/2] += ly
            
            rz[rz > lz/2] -= lz
            rz[rz < -lz/2] += lz
            
            # Calculate distances
            r = np.sqrt(rx**2 + ry**2 + rz**2)
            
            # Store results
            distances[:, i] = r
            dx[:, i] = -rx  # Note the negative sign to match MATLAB implementation
            dy[:, i] = -ry
            dz[:, i] = -rz
    else:
        # Triclinic Box approach
        # Setup progress bar
        if has_tqdm:
            atom_iterator = tqdm(range(n_atoms), desc="Finding dists", unit="atom")
        else:
            print("Finding distances...")
            atom_iterator = range(n_atoms)
            last_percent = -1
            
        for i in atom_iterator:
            # Update progress percentage for non-tqdm case
            if not has_tqdm and n_atoms > 100:
                percent = int(100 * i / n_atoms)
                if percent > last_percent and percent % 10 == 0:
                    print(f"  {percent}% complete...")
                    last_percent = percent
            # Calculate initial distance components
            rx = xyz[i, 0] - xyz[:, 0]
            ry = xyz[i, 1] - xyz[:, 1]
            rz = xyz[i, 2] - xyz[:, 2]
            
            # Apply minimum image convention for triclinic Box
            # First handle z-direction
            z_gt_ind = rz > lz/2
            z_lt_ind = rz < -lz/2
            
            rz[z_gt_ind] -= lz
            rz[z_lt_ind] += lz
            
            rx[z_gt_ind] -= xz
            rx[z_lt_ind] += xz
            
            ry[z_gt_ind] -= yz
            ry[z_lt_ind] += yz
            
            # Then handle y-direction
            y_gt_ind = ry > ly/2
            y_lt_ind = ry < -ly/2
            
            ry[y_gt_ind] -= ly
            ry[y_lt_ind] += ly
            
            rx[y_gt_ind] -= xy
            rx[y_lt_ind] += xy
            
            # Finally handle x-direction
            x_gt_ind = rx > lx/2
            x_lt_ind = rx < -lx/2
            
            rx[x_gt_ind] -= lx
            rx[x_lt_ind] += lx
            
            # Calculate distances
            r = np.sqrt(rx**2 + ry**2 + rz**2)
            
            # Store results
            distances[:, i] = r
            dx[:, i] = -rx  # Note the negative sign to match MATLAB implementation
            dy[:, i] = -ry
            dz[:, i] = -rz
    
    # Transpose distances to match MATLAB output format
    distances = distances.T
    dx = dx.T
    dy = dy.T
    dz = dz.T
    
    return distances, dx, dy, dz

def dist_matrix_direct(atoms):
    """Calculate a direct distance matrix between atoms without periodic boundaries.
    
    Args:
        atoms: list of atom dictionaries, each having 'x', 'y', 'z' coordinates.
       
    Returns:
        A tuple of four numpy arrays: 
        - A numpy array of shape (N, N) with pairwise distances.
        - Three numpy arrays of shape (N, N) with pairwise x, y, z differences.
    """
    # Extract atomic positions
    n_atoms = len(atoms)
    positions = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
    
    # Calculate direct distances
    distances = np.zeros((n_atoms, n_atoms))
    dx = np.zeros((n_atoms, n_atoms))
    dy = np.zeros((n_atoms, n_atoms))
    dz = np.zeros((n_atoms, n_atoms))
    
    for i in range(n_atoms):
        rx = positions[i, 0] - positions[:, 0]
        ry = positions[i, 1] - positions[:, 1]
        rz = positions[i, 2] - positions[:, 2]
        
        r = np.sqrt(rx**2 + ry**2 + rz**2)
        
        distances[i, :] = r
        dx[i, :] = rx
        dy[i, :] = ry
        dz[i, :] = rz
    
    return distances, dx, dy, dz

def dist_matrix_hybrid(atoms, Box=None, use_pbc=True, Box_dim=None):
    """Calculate the distance matrix using either direct or PBC approach based on use_pbc flag.
    
    Args:
        atoms: list of atom dictionaries, each having 'x', 'y', 'z' coordinates.
        Box: a 1x3, 1x6 or 1x9 list representing Cell dimensions (in Angstroms):
            - For orthogonal boxes, a 1x3 list [lx, ly, lz]
            - For Cell parameters, a 1x6 list [a, b, c, alpha, beta, gamma] (Cell format)
            - For triclinic boxes, a 1x9 list [lx, ly, lz, 0, 0, xy, 0, xz, yz] (Box_dim format)
        use_pbc: Whether to use periodic boundary conditions (default: True).
        Box_dim: Deprecated. Use 'Box' instead. Maintained for backward compatibility.
       
    Returns:
        A tuple of four numpy arrays: 
        - A numpy array of shape (N, N) with pairwise distances.
        - Three numpy arrays of shape (N, N) with pairwise x, y, z differences.
    """
    # Backward compatibility for Box_dim parameter
    if Box is None and Box_dim is not None:
        Box = Box_dim
        
    if use_pbc:
        if Box is None:
            raise ValueError("Box parameter must be provided when use_pbc=True")
        return dist_matrix(atoms, Box)
    else:
        return dist_matrix_direct(atoms)
