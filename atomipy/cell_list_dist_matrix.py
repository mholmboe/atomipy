import numpy as np
from .cell_utils import Box_dim2Cell

# Optional numba support for JIT compilation
try:
    import numba
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    numba = None

# Decorator that applies numba JIT if available, otherwise does nothing
def optional_jit(func):
    """Apply numba JIT compilation if available, otherwise return the function unchanged."""
    if HAS_NUMBA:
        return numba.jit(nopython=True)(func)
    return func

def cell_list_dist_matrix(atoms, cutoff=2.45, Box_dim=None, Cell=None, rmaxH=1.2, H_type='H'):
    """Calculate a sparse distance matrix using the cell list algorithm for efficiently 
    finding all atom pairs within a cutoff distance. This function closely follows the MATLAB 
    implementation of cell_list_dist_matrix_MATLAB.m.
    
    Args:
        atoms: list of atom dictionaries, each having 'x', 'y', 'z' coordinates and 'type' field.
        cutoff: maximum distance to consider for non-hydrogen bonds (default: 2.45 Å).
        Box_dim: (Optional) Box dimensions in the format [lx, ly, lz, 0, 0, xy, 0, xz, yz] for triclinic cells,
                or [lx, ly, lz] for orthogonal cells. Either Box_dim or Cell must be provided.
        Cell: (Optional) Cell parameters in the format [a, b, c, alpha, beta, gamma]. 
              Either Box_dim or Cell must be provided.
        rmaxH: cutoff distance for bonds involving hydrogen atoms (default: 1.2 Å).
        H_type: atom type string identifying hydrogen atoms (default: 'H').
       
    Returns:
        Tuple containing:
        - dist_matrix: NxN numpy array with distances between atoms
        - bond_list: Nx2 numpy array of atom indices forming bonds
        - dist_list: Nx1 numpy array of distances corresponding to bonds
        - X_dist, Y_dist, Z_dist: NxN numpy arrays with distance vector components
        
    Note:
        This implementation closely follows the MATLAB version, including the same cell grid
        approach and handling of different cutoffs for hydrogen atoms.
    """
    # Parse input & setup
    N = len(atoms)
    positions = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
    
    # Possibly convert Box_dim into [a,b,c,alpha,beta,gamma] form
    if Box_dim is not None:
        if len(Box_dim) == 9:
            # Convert from Box_dim format to Cell format
            if Cell is None:
                Cell = box_dim_to_cell(Box_dim)
        elif len(Box_dim) == 3:
            # Orthogonal box
            Cell = list(Box_dim) + [90.0, 90.0, 90.0]
    
    if Cell is None:
        raise ValueError("Either Box_dim or Cell must be provided")
    
    # Extract cell parameters
    a, b, c = Cell[0], Cell[1], Cell[2]
    
    # For orthogonal, angles = 90
    if len(Cell) == 3:
        alpha_rad = np.radians(90)
        beta_rad = np.radians(90)
        gamma_rad = np.radians(90)
    else:
        alpha_rad = np.radians(Cell[3])
        beta_rad = np.radians(Cell[4])
        gamma_rad = np.radians(Cell[5])
    # Construct triclinic box matrix H
    ax = a
    bx = b * np.cos(gamma_rad)
    by = b * np.sin(gamma_rad)
    cx = c * np.cos(beta_rad)
    cy = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    
    H = np.array([
        [ax, bx, cx],
        [0,  by, cy],
        [0,  0,  cz]
    ])
    Hinv = np.linalg.inv(H)
    
    # Build the cell list
    # Pick a cell size that ensures neighbors can be found
    cellSize = 2 * cutoff
    
    # Bounding box size calculation
    boundingBoxSize = np.array([
        max(abs(np.array([ax, bx, cx]))),
        max(abs(np.array([0, by, cy]))),
        max(abs(np.array([0, 0, cz])))
    ])
    
    # Number of cells along each dimension
    nCells = np.maximum(np.floor(boundingBoxSize / cellSize), 1).astype(int)
    
    # Convert coords -> fractional
    fracCoords = np.zeros((N, 3))
    for i in range(N):
        fracCoords[i] = Hinv @ np.array([positions[i][0], positions[i][1], positions[i][2]])
    
    # Keep coordinates in [0,1) range
    fracCoords = fracCoords - np.floor(fracCoords)
    
    # Compute which cell each atom belongs to
    cellIndex = np.floor(fracCoords * nCells).astype(int)
    
    # Fix boundary cases
    outOfBound = cellIndex[:, 0] >= nCells[0]
    cellIndex[outOfBound, 0] = nCells[0] - 1
    
    outOfBound = cellIndex[:, 1] >= nCells[1]
    cellIndex[outOfBound, 1] = nCells[1] - 1
    
    outOfBound = cellIndex[:, 2] >= nCells[2]
    cellIndex[outOfBound, 2] = nCells[2] - 1
    
    # Convert (ix,iy,iz) -> single linear index
    cellLinIdx = np.ravel_multi_index((cellIndex[:, 0], cellIndex[:, 1], cellIndex[:, 2]), nCells)
    
    # Initialize cell lists
    numCells = np.prod(nCells)
    cellList = [[] for _ in range(numCells)]
    
    # Populate cell lists - vectorized approach with bincount and digitize would be faster,
    # but we need to maintain the list structure for compatibility
    for iAtom in range(N):
        cIdx = cellLinIdx[iAtom]
        cellList[cIdx].append(iAtom)
    
    # Identify neighbor cells (27 total with PBC)
    neighborOffsets = []
    for ix in [-1, 0, 1]:
        for iy in [-1, 0, 1]:
            for iz in [-1, 0, 1]:
                neighborOffsets.append((ix, iy, iz))
    
    # Initialize output data structures with estimated capacity
    # Estimate maximum number of possible bonds based on typical coordination numbers
    # and system density for better memory allocation
    vol = a * b * c if len(Cell) >= 6 else Cell[0] * Cell[1] * Cell[2]
    density = N / vol
    coord_factor = min(12, max(4, int(4.0 * np.pi * (cutoff**3) * density / 3.0)))
    est_bonds = int(N * coord_factor / 2)  # divide by 2 to avoid double counting
    
    # Pre-allocate with estimated capacity (still using list append for compatibility)
    bond_list = []
    dist_list = []
    
    # Initialize NxN output matrices (distances & displacement vectors)
    dist_matrix = np.zeros((N, N), dtype=np.float32)
    X_dist = np.zeros((N, N), dtype=np.float32)
    Y_dist = np.zeros((N, N), dtype=np.float32)
    Z_dist = np.zeros((N, N), dtype=np.float32)
    
    # Loop over cells/neighbors & compute distances
    for cID in range(numCells):
        atomListC = cellList[cID]
        if not atomListC:
            continue
            
        # Get cell indices
        cx, cy, cz = np.unravel_index(cID, nCells)
        
        # Check this cell and all neighboring cells
        for dxIdx, dyIdx, dzIdx in neighborOffsets:
            # Calculate neighbor cell with periodic wrapping
            nx = (cx + dxIdx) % nCells[0]
            ny = (cy + dyIdx) % nCells[1] 
            nz = (cz + dzIdx) % nCells[2]
            
            # Get linear index of neighbor cell
            neighborCellLin = np.ravel_multi_index((nx, ny, nz), nCells)
            atomListN = cellList[neighborCellLin]
            if not atomListN:
                continue
                
            # Avoid double-counting
            if neighborCellLin < cID:
                continue
            # Compute pairwise distances iAtom <-> jAtom
            for i, iAtom in enumerate(atomListC):
                for j, jAtom in enumerate(atomListN):
                    # If in the same cell, only consider jAtom > iAtom to avoid duplicates
                    if neighborCellLin == cID and jAtom <= iAtom:
                        continue
                        
                    # Check hydrogen vs. non-hydrogen
                    isH_i = atoms[iAtom]['type'] == H_type
                    isH_j = atoms[jAtom]['type'] == H_type
                    
                    if isH_i or isH_j:
                        localCutoff = rmaxH  # short cutoff for hydrogen
                    else:
                        localCutoff = cutoff  # default cutoff for non-hydrogen
                    
                    # Calculate distance with minimum image convention
                    # Get difference in fractional coordinates
                    ri = np.array([positions[iAtom][0], positions[iAtom][1], positions[iAtom][2]])
                    rj = np.array([positions[jAtom][0], positions[jAtom][1], positions[jAtom][2]])
                    
                    # Convert to fractional space
                    diffFrac = Hinv @ (rj - ri)
                    
                    # Apply minimum image convention - wrap to [-0.5, 0.5)
                    diffFrac = diffFrac - np.round(diffFrac)
                    
                    # Convert back to Cartesian
                    dVec = H @ diffFrac
                    d = np.linalg.norm(dVec)
                    
                    # Only store if within cutoff
                    if d <= localCutoff:
                        # Store pair & distance
                        bond_list.append([iAtom, jAtom])
                        dist_list.append(d)
                        
                        # Fill NxN distance matrix
                        dist_matrix[iAtom, jAtom] = d
                        dist_matrix[jAtom, iAtom] = d  # symmetric
                        
                        # Fill NxN displacement matrices
                        X_dist[iAtom, jAtom] = dVec[0]
                        X_dist[jAtom, iAtom] = -dVec[0]
                        
                        Y_dist[iAtom, jAtom] = dVec[1]
                        Y_dist[jAtom, iAtom] = -dVec[1]
                        
                        Z_dist[iAtom, jAtom] = dVec[2]
                        Z_dist[jAtom, iAtom] = -dVec[2]
    
    # Clean up very small values that might be numerical errors - single mask for efficiency
    small_vals_mask = np.abs(dist_matrix) <= 1e-7
    dist_matrix[small_vals_mask] = 0
    X_dist[small_vals_mask] = 0
    Y_dist[small_vals_mask] = 0
    Z_dist[small_vals_mask] = 0
    
    # Convert lists to numpy arrays
    if bond_list:
        bond_list = np.array(bond_list)
        dist_list = np.array(dist_list)
    else:
        bond_list = np.zeros((0, 2), dtype=int)
        dist_list = np.zeros(0, dtype=float)
    
    return dist_matrix, X_dist, Y_dist, Z_dist, bond_list, dist_list


def box_dim_to_cell(Box_dim):
    """Convert Box_dim [lx, ly, lz, 0, 0, xy, 0, xz, yz] to Cell [a, b, c, alpha, beta, gamma].
    
    Args:
        Box_dim: Box dimensions in the format [lx, ly, lz, 0, 0, xy, 0, xz, yz]
        
    Returns:
        Cell parameters [a, b, c, alpha, beta, gamma] where angles are in degrees
    """
    # Use the more comprehensive implementation from cell_utils
    return Box_dim2Cell(Box_dim)


@optional_jit
def convert_to_sparse_dict(dist_matrix, X_dist, Y_dist, Z_dist, cutoff):
    """Convert full distance matrices to a sparse dictionary format.
    
    Args:
        dist_matrix: NxN distance matrix
        X_dist, Y_dist, Z_dist: NxN displacement component matrices
        cutoff: Only include distances up to this cutoff
        
    Returns:
        Dictionary mapping (i,j) tuples to (dist, dx, dy, dz) tuples
    """
    N = dist_matrix.shape[0]
    distance_dict = {}
    
    # Vectorized approach - find valid indices in one go
    i_indices, j_indices = np.where((dist_matrix > 0) & (dist_matrix <= cutoff) & (np.triu(np.ones(dist_matrix.shape), k=1) > 0))
    
    # Create dictionary entries for all valid pairs
    for idx in range(len(i_indices)):
        i, j = i_indices[idx], j_indices[idx]
        distance_dict[(i, j)] = (dist_matrix[i, j], X_dist[i, j], Y_dist[i, j], Z_dist[i, j])
    
    return distance_dict


@optional_jit
def get_neighbors(dist_matrix, X_dist, Y_dist, Z_dist, atom_index, r_max=None):
    """Get all neighbors of a specific atom from the distance matrices.
    
    Args:
        dist_matrix: NxN distance matrix
        X_dist, Y_dist, Z_dist: NxN displacement component matrices
        atom_index: Index of the atom to get neighbors for
        r_max: Optional maximum distance for neighbors
               
    Returns:
        List of tuples (neighbor_index, distance, dx, dy, dz)
    """
    # Use vectorized operations to get neighbor indices
    row = dist_matrix[atom_index]
    
    if r_max is None:
        mask = (row > 0) & (np.arange(len(row)) != atom_index)
    else:
        mask = (row > 0) & (row <= r_max) & (np.arange(len(row)) != atom_index)
    
    # Get indices where mask is True
    j_indices = np.where(mask)[0]
    
    # Build neighbor list
    neighbors = []
    for j in j_indices:
        neighbors.append((j, dist_matrix[atom_index, j], 
                         X_dist[atom_index, j], 
                         Y_dist[atom_index, j], 
                         Z_dist[atom_index, j]))
    
    # Sort by distance - use numpy for speed if numba available
    if len(neighbors) > 1:
        # Convert to numpy array for faster sorting
        neighbors_array = np.array(neighbors, dtype=[('idx', int), ('dist', float), 
                                                    ('dx', float), ('dy', float), ('dz', float)])
        # Sort by distance
        neighbors_array.sort(order='dist')
        # Convert back to list of tuples
        neighbors = [tuple(row) for row in neighbors_array]
    
    return neighbors
