import numpy as np
from .cell_utils import normalize_box

# Try to import tqdm for progress bar
try:
    from tqdm import tqdm
    has_tqdm = True
except ImportError:
    has_tqdm = False

def get_progress_iterator(iterable, desc="Processing", unit="it"):
    if has_tqdm:
        return tqdm(iterable, desc=desc, unit=unit)
    return iterable

def cell_list_dist_matrix(atoms, Box, cutoff=2.45, rmaxH=1.2, H_type='H'):
    """
    Higly optimized Cell-list algorithm for finding atom pairs within a cutoff.
    Vectorized using NumPy to avoid nested Python loops. Supports Triclinic boxes.
    
    Args:
        atoms: list of atom dictionaries.
        Box: a 1x3, 1x6 or 1x9 list representing Cell dimensions.
        cutoff: maximum distance for non-hydrogen bonds.
        rmaxH: cutoff distance for hydrogen bonds.
        H_type: atom type string for hydrogen.
        
    Returns:
        dist_matrix: NxN numpy array (sparse-ish, mostly zeros).
        X_dist, Y_dist, Z_dist: NxN displacement component matrices.
        bond_list: Mx2 numpy array of atom indices.
        dist_list: Mx1 numpy array of distances.
    """
    N = len(atoms)
    positions = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms], dtype=np.float32)
    types = np.array([atom.get('type', '') for atom in atoms])
    is_h = np.array([bool(t and t[0].upper() == 'H') for t in types])
    
    if Box is None:
        raise ValueError("Box parameter must be provided")

    Box_dim, Cell = normalize_box(Box)
    
    # Construct triclinic Box matrix H and its inverse Hinv
    a, b, c = Cell[0], Cell[1], Cell[2]
    if len(Cell) == 3:
        alpha, beta, gamma = 90.0, 90.0, 90.0
    else:
        alpha, beta, gamma = Cell[3], Cell[4], Cell[5]
        
    ar, br, gr = np.radians([alpha, beta, gamma])
    
    # Upper triangular lattice matrix H
    ax = a
    bx = b * np.cos(gr)
    by = b * np.sin(gr)
    cx = c * np.cos(br)
    cy = c * (np.cos(ar) - np.cos(br) * np.cos(gr)) / np.sin(gr)
    cz = np.sqrt(max(0, c**2 - cx**2 - cy**2))
    
    H = np.array([[ax, bx, cx], [0, by, cy], [0, 0, cz]], dtype=np.float32)
    Hinv = np.linalg.inv(H)
    
    # 1. Bin atoms into cells
    # We use fractional coordinates for easy binning in triclinic systems
    frac_coords = (Hinv @ positions.T).T
    frac_coords = frac_coords % 1.0
    
    # Cell size should be at least the cutoff
    max_cutoff = max(cutoff, rmaxH)
    
    # Effective box widths in each direction
    # For triclinic, we use the perpendicular widths (projections)
    # Using cz, by, ax is a safe approximation for orthogonal-like cells
    n_cells = np.maximum(np.floor(np.array([ax, by, cz]) / max_cutoff), 1).astype(int)
    
    # Map fractional coords to integer cell indices
    cell_idx = np.floor(frac_coords * n_cells).astype(int)
    cell_idx = np.clip(cell_idx, 0, n_cells - 1)
    
    # Flat cell index
    flat_idx = cell_idx[:, 0] * (n_cells[1] * n_cells[2]) + cell_idx[:, 1] * n_cells[2] + cell_idx[:, 2]
    num_cells_total = np.prod(n_cells)
    
    # Create head and next arrays (Linked List approach, but better for NumPy)
    # head[c] is the index of the first atom in cell c
    # next_atom[i] is the index of the next atom in the same cell as i
    head = np.full(num_cells_total, -1, dtype=np.int32)
    next_atom = np.full(N, -1, dtype=np.int32)
    
    for i in range(N):
        c = flat_idx[i]
        next_atom[i] = head[c]
        head[c] = i
        
    # 2. Iterate over cells and neighbors
    # Pre-allocate results (using lists for appending bonds is fine, but we'll pre-allocate matrices)
    dist_matrix = np.zeros((N, N), dtype=np.float32)
    X_dist = np.zeros((N, N), dtype=np.float32)
    Y_dist = np.zeros((N, N), dtype=np.float32)
    Z_dist = np.zeros((N, N), dtype=np.float32)
    
    bond_list = []
    dist_list = []
    
    # 27 neighbors including self
    offsets = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)
    
    # We only need to check half the neighbors to avoid double counting, 
    # but with PBC it's safer to check all or use a specific subset.
    # To keep it simple and accurate, we check all and filter i < j later or at runtime.
    
    # Use a progress iterator for the outer x-dimension loop
    for cx in get_progress_iterator(range(n_cells[0]), desc="Finding dists", unit="CellX"):
        for cy in range(n_cells[1]):
            for cz_idx in range(n_cells[2]):
                c1 = cx * (n_cells[1] * n_cells[2]) + cy * n_cells[2] + cz_idx
                i = head[c1]
                if i == -1: continue
                
                # Get all atoms in this cell
                atoms1 = []
                while i != -1:
                    atoms1.append(i)
                    i = next_atom[i]
                atoms1 = np.array(atoms1, dtype=np.int32)
                
                # Check 27 neighbor cells
                for off in offsets:
                    nx = (cx + off[0]) % n_cells[0]
                    ny = (cy + off[1]) % n_cells[1]
                    nz = (cz_idx + off[2]) % n_cells[2]
                    
                    c2 = nx * (n_cells[1] * n_cells[2]) + ny * n_cells[2] + nz
                    
                    j = head[c2]
                    if j == -1: continue
                    
                    atoms2 = []
                    while j != -1:
                        atoms2.append(j)
                        j = next_atom[j]
                    atoms2 = np.array(atoms2, dtype=np.int32)
                    
                    # Compute distances between atoms1 and atoms2
                    # Vectorized chunk
                    p1 = positions[atoms1]
                    p2 = positions[atoms2]
                    
                    # Outer subtraction for N1 x N2 pairs
                    # p1: (N1, 3), p2: (N2, 3) -> diff: (N1, N2, 3)
                    diff = p2[np.newaxis, :, :] - p1[:, np.newaxis, :]
                    
                    # Periodic wrap in Cartesian using H matrix
                    # Easiest way to wrap is in fractional space
                    diff_frac = (Hinv @ diff.reshape(-1, 3).T).T
                    diff_frac = diff_frac - np.round(diff_frac)
                    diff_cart = (H @ diff_frac.T).T.reshape(len(atoms1), len(atoms2), 3)
                    
                    d2 = np.sum(diff_cart**2, axis=2)
                    d = np.sqrt(d2)
                    
                    # Masking by cutoff
                    # Note: different cutoffs for H vs non-H
                    is_h1 = is_h[atoms1]
                    is_h2 = is_h[atoms2]
                    # matrix of cutoffs (N1 x N2)
                    cutoffs = np.where(is_h1[:, np.newaxis] | is_h2[np.newaxis, :], rmaxH, cutoff)
                    
                    mask = (d > 1e-7) & (d <= cutoffs)
                    
                    # Filter i < j to avoid double counts if needed, but here we fill the whole matrix.
                    # Actually, if c1 == c2, only i < j. If c1 != c2, we'll double count unless we filter.
                    # Easiest way: only fill if i < j.
                    
                    ii, jj = np.where(mask)
                    for k in range(len(ii)):
                        idx1 = atoms1[ii[k]]
                        idx2 = atoms2[jj[k]]
                        if idx1 < idx2:
                            val = d[ii[k], jj[k]]
                            dist_matrix[idx1, idx2] = val
                            dist_matrix[idx2, idx1] = val
                            X_dist[idx1, idx2] = diff_cart[ii[k], jj[k], 0]
                            X_dist[idx2, idx1] = -X_dist[idx1, idx2]
                            Y_dist[idx1, idx2] = diff_cart[ii[k], jj[k], 1]
                            Y_dist[idx2, idx1] = -Y_dist[idx1, idx2]
                            Z_dist[idx1, idx2] = diff_cart[ii[k], jj[k], 2]
                            Z_dist[idx2, idx1] = -Z_dist[idx1, idx2]
                            bond_list.append([idx1, idx2])
                            dist_list.append(val)
                            
    return dist_matrix, X_dist, Y_dist, Z_dist, np.array(bond_list), np.array(dist_list)

# Alias for backward compatibility during migration
cell_list_dist_matrix_fast = cell_list_dist_matrix

def neighbor_list_fast(atoms, Box, cutoff=2.45, rmaxH=None, H_type='H'):
    """
    Higly optimized Cell-list algorithm for finding atom pairs within a cutoff.
    Returns sparse lists instead of NxN matrices to save memory.
    Supports Triclinic boxes.
    
    Args:
        atoms: list of atom dictionaries.
        Box: a 1x3, 1x6 or 1x9 list representing Cell dimensions.
        cutoff: maximum distance for non-hydrogen bonds.
        rmaxH: cutoff distance for hydrogen bonds. If None, uses cutoff.
        H_type: atom type string for hydrogen.
        
    Returns:
        i_idx, j_idx: Nx1 numpy arrays of atom indices.
        dist: Nx1 numpy array of distances.
        dx, dy, dz: Nx1 numpy arrays of displacement components (r_j - r_i).
    """
    N = len(atoms)
    positions = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms], dtype=np.float32)
    types = np.array([atom.get('type', '') for atom in atoms])
    is_h = np.array([bool(t and t[0].upper() == 'H') for t in types])
    
    if rmaxH is None:
        rmaxH = cutoff
        
    if Box is None:
        # Non-periodic case: use a bounding box
        min_coords = np.min(positions, axis=0)
        max_coords = np.max(positions, axis=0)
        Box_dim = max_coords - min_coords + 2 * max(cutoff, rmaxH)
        H = np.diag(Box_dim).astype(np.float32)
        Hinv = np.linalg.inv(H)
        # Shift positions to be positive relative to min_coords
        positions_shifted = positions - min_coords + max(cutoff, rmaxH)
        frac_coords = (Hinv @ positions_shifted.T).T
        is_periodic = False
    else:
        Box_dim, Cell = normalize_box(Box)
        # Construct triclinic Box matrix H and its inverse Hinv
        a, b, c = Cell[0], Cell[1], Cell[2]
        if len(Cell) == 3:
            alpha, beta, gamma = 90.0, 90.0, 90.0
        else:
            alpha, beta, gamma = Cell[3], Cell[4], Cell[5]
        ar, br, gr = np.radians([alpha, beta, gamma])
        ax = a
        bx = b * np.cos(gr)
        by = b * np.sin(gr)
        cx = c * np.cos(br)
        cy = c * (np.cos(ar) - np.cos(br) * np.cos(gr)) / np.sin(gr)
        cz = np.sqrt(max(0, c**2 - cx**2 - cy**2))
        H = np.array([[ax, bx, cx], [0, by, cy], [0, 0, cz]], dtype=np.float32)
        Hinv = np.linalg.inv(H)
        frac_coords = (Hinv @ positions.T).T
        frac_coords = frac_coords % 1.0
        is_periodic = True

    max_cutoff = max(cutoff, rmaxH)
    # Effective box widths in each direction
    if is_periodic:
        n_cells = np.maximum(np.floor(np.array([H[0,0], H[1,1], H[2,2]]) / max_cutoff), 1).astype(int)
    else:
        n_cells = np.maximum(np.floor(Box_dim / max_cutoff), 1).astype(int)
        
    cell_idx = np.floor(frac_coords * n_cells).astype(int)
    cell_idx = np.clip(cell_idx, 0, n_cells - 1)
    flat_idx = cell_idx[:, 0] * (n_cells[1] * n_cells[2]) + cell_idx[:, 1] * n_cells[2] + cell_idx[:, 2]
    num_cells_total = np.prod(n_cells)
    
    head = np.full(num_cells_total, -1, dtype=np.int32)
    next_atom = np.full(N, -1, dtype=np.int32)
    for i in range(N):
        c = flat_idx[i]
        next_atom[i] = head[c]
        head[c] = i
        
    i_idx_list = []
    j_idx_list = []
    dist_list = []
    dx_list = []
    dy_list = []
    dz_list = []
    
    offsets = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)
    
    # Use a progress iterator for the outer x-dimension loop
    for cx in get_progress_iterator(range(n_cells[0]), desc="Finding sparse dists", unit="CellX"):
        for cy in range(n_cells[1]):
            for cz_idx in range(n_cells[2]):
                c1 = cx * (n_cells[1] * n_cells[2]) + cy * n_cells[2] + cz_idx
                i_head = head[c1]
                if i_head == -1: continue
                
                atoms1 = []
                curr = i_head
                while curr != -1:
                    atoms1.append(curr)
                    curr = next_atom[curr]
                atoms1 = np.array(atoms1, dtype=np.int32)
                
                for off in offsets:
                    nx = (cx + off[0])
                    ny = (cy + off[1])
                    nz = (cz_idx + off[2])
                    
                    if is_periodic:
                        nx %= n_cells[0]
                        ny %= n_cells[1]
                        nz %= n_cells[2]
                    else:
                        if nx < 0 or nx >= n_cells[0] or ny < 0 or ny >= n_cells[1] or nz < 0 or nz >= n_cells[2]:
                            continue
                    
                    c2 = nx * (n_cells[1] * n_cells[2]) + ny * n_cells[2] + nz
                    j_head = head[c2]
                    if j_head == -1: continue
                    
                    atoms2 = []
                    curr = j_head
                    while curr != -1:
                        atoms2.append(curr)
                        curr = next_atom[curr]
                    atoms2 = np.array(atoms2, dtype=np.int32)
                    
                    p1 = positions[atoms1]
                    p2 = positions[atoms2]
                    
                    diff = p2[np.newaxis, :, :] - p1[:, np.newaxis, :]
                    
                    if is_periodic:
                        diff_frac = (Hinv @ diff.reshape(-1, 3).T).T
                        diff_frac = diff_frac - np.round(diff_frac)
                        diff_cart = (H @ diff_frac.T).T.reshape(len(atoms1), len(atoms2), 3)
                    else:
                        diff_cart = diff

                    d2 = np.sum(diff_cart**2, axis=2)
                    d = np.sqrt(d2)
                    
                    is_h1 = is_h[atoms1]
                    is_h2 = is_h[atoms2]
                    cutoffs = np.where(is_h1[:, np.newaxis] | is_h2[np.newaxis, :], rmaxH, cutoff)
                    
                    mask = (d > 1e-7) & (d <= cutoffs)
                    
                    ii_local, jj_local = np.where(mask)
                    for k in range(len(ii_local)):
                        idx1 = atoms1[ii_local[k]]
                        idx2 = atoms2[jj_local[k]]
                        if idx1 < idx2:
                            val = d[ii_local[k], jj_local[k]]
                            i_idx_list.append(idx1)
                            j_idx_list.append(idx2)
                            dist_list.append(val)
                            dx_list.append(diff_cart[ii_local[k], jj_local[k], 0])
                            dy_list.append(diff_cart[ii_local[k], jj_local[k], 1])
                            dz_list.append(diff_cart[ii_local[k], jj_local[k], 2])
                            
    return (np.array(i_idx_list, dtype=np.int32), 
            np.array(j_idx_list, dtype=np.int32), 
            np.array(dist_list, dtype=np.float32),
            np.array(dx_list, dtype=np.float32),
            np.array(dy_list, dtype=np.float32),
            np.array(dz_list, dtype=np.float32))

# Carry over helper functions from original file for full compatibility
def convert_to_sparse_dict(dist_matrix, X_dist, Y_dist, Z_dist, cutoff):
    """Convert full distance matrices to a sparse dictionary format."""
    N = dist_matrix.shape[0]
    distance_dict = {}
    i_indices, j_indices = np.where((dist_matrix > 0) & (dist_matrix <= cutoff) & (np.triu(np.ones(dist_matrix.shape), k=1) > 0))
    for idx in range(len(i_indices)):
        i, j = i_indices[idx], j_indices[idx]
        distance_dict[(i, j)] = (dist_matrix[i, j], X_dist[i, j], Y_dist[i, j], Z_dist[i, j])
    return distance_dict

def get_neighbors(dist_matrix, X_dist, Y_dist, Z_dist, atom_index, r_max=None):
    """Get all neighbors of a specific atom from the distance matrices."""
    row = dist_matrix[atom_index]
    if r_max is None:
        mask = (row > 0) & (np.arange(len(row)) != atom_index)
    else:
        mask = (row > 0) & (row <= r_max) & (np.arange(len(row)) != atom_index)
    j_indices = np.where(mask)[0]
    neighbors = []
    for j in j_indices:
        neighbors.append((j, dist_matrix[atom_index, j], X_dist[atom_index, j], Y_dist[atom_index, j], Z_dist[atom_index, j]))
    if len(neighbors) > 1:
        # Sort by distance
        neighbors.sort(key=lambda x: x[1])
    return neighbors
