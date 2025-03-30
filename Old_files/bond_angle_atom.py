import numpy as np
from dist_matrix_atom import dist_matrix_atom
from radius_atom import radius_atom

def bond_angle_atom(atoms, Box_dim, rmaxH=1.2, rmaxM=2.45):
    """Compute bonds and angles for a given atomic structure.
    
    For each atom, bonds are determined based on a distance threshold:
      - rmaxH (default 1.2 Å) if either atom is hydrogen
      - rmaxM (default 2.45 Å) for bonds between non-hydrogen atoms.
    
    Angles are then computed for each pair of bonds at the central atom using the periodic boundary condition (PBC)
    corrected vectors. The function updates each atom's 'neigh', 'bonds', and 'angles' fields in-place.

    Args:
       atoms: list of atom dictionaries.
       Box_dim: 1x9 list representing triclinic cell dimensions.
       rmaxH: cutoff distance for bonds involving hydrogen (default 1.2 Å).
       rmaxM: cutoff distance for bonds between all other atoms (default 2.45 Å).

    Returns:
       Updated atoms list with 'neigh', 'bonds', and 'angles'.
    """
    # Obtain distance matrix and Cartesian differences accounting for PBC
    dmat, dx, dy, dz = dist_matrix_atom(atoms, Box_dim)
    positions = np.array([[atom['x'], atom['y'], atom['z']] for atom in atoms])
    N = len(atoms)

    for i in range(N):
        atoms[i]['neigh'] = []
        atoms[i]['bonds'] = []
        atoms[i]['angles'] = []
        el_i = atoms[i].get('element','X')
        for j in range(N):
            if i == j:
                continue
            el_j = atoms[j].get('element','X')
            threshold = rmaxH if (el_i == 'H' or el_j == 'H') else rmaxM
            if dmat[i, j] < threshold:
                atoms[i]['neigh'].append(j)
                atoms[i]['bonds'].append((j, dmat[i, j]))
        # Compute angles for every unique pair of neighbors using PBC-adjusted vectors
        for m in range(len(atoms[i]['neigh'])):
            for n in range(m+1, len(atoms[i]['neigh'])):
                j = atoms[i]['neigh'][m]
                k = atoms[i]['neigh'][n]
                vec_ij = np.array([dx[i, j], dy[i, j], dz[i, j]])
                vec_ik = np.array([dx[i, k], dy[i, k], dz[i, k]])
                norm_ij = np.linalg.norm(vec_ij)
                norm_ik = np.linalg.norm(vec_ik)
                if norm_ij == 0 or norm_ik == 0:
                    continue
                cos_angle = np.dot(vec_ij, vec_ik) / (norm_ij * norm_ik)
                # Clip to avoid numerical issues
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle = np.degrees(np.arccos(cos_angle))
                atoms[i]['angles'].append(((j, k), angle))
    return atoms
