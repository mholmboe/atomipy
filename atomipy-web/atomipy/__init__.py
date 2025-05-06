"""
atomipy: The atom Toolbox in Python

This package provides tools for working with molecular structures, particularly
focused on mineral systems with periodic boundary conditions. It supports both
orthogonal and triclinic simulation cells, and provides efficient distance
calculations, bond/angle detection, and structure manipulation functions.
"""

# ===== File I/O functions =====
# Import/export functions
from . import import_conf
import_pdb = import_conf.pdb
import_gro = import_conf.gro
import_xyz = import_conf.xyz
import_auto = import_conf.auto

from . import write_conf
write_pdb = write_conf.pdb
write_gro = write_conf.gro
write_xyz = write_conf.xyz
write_auto = write_conf.auto

# Topology file generation
from . import write_top
write_itp = write_top.itp
write_psf = write_top.psf
write_lmp = write_top.lmp

# ===== Atom property functions =====
from .element import element

from .radius import radius

try:
    from .mass import mass
except ImportError:
    pass

# ===== Structure analysis functions =====
from .dist_matrix import dist_matrix
from .cell_list_dist_matrix import cell_list_dist_matrix
from .bond_angle import bond_angle

# ===== Cell and coordinate transformation functions =====
from . import cell_utils
Box_dim2Cell = cell_utils.Box_dim2Cell
Cell2Box_dim = cell_utils.Cell2Box_dim

# Coordinate transformation module
from . import transform
cartesian_to_fractional = transform.cartesian_to_fractional
fractional_to_cartesian = transform.fractional_to_cartesian
wrap_coordinates = transform.wrap_coordinates
triclinic_to_orthogonal = transform.triclinic_to_orthogonal
orthogonal_to_triclinic = transform.orthogonal_to_triclinic
get_orthogonal_box = transform.get_orthogonal_box
get_cell_vectors = transform.get_cell_vectors
direct_cartesian_to_fractional = transform.direct_cartesian_to_fractional
direct_fractional_to_cartesian = transform.direct_fractional_to_cartesian

from . import replicate
replicate_system = replicate.replicate_system

from . import move
translate = move.translate

from . import add
update = add.update

# ===== Force field functions =====
try:
    from .forcefield import minff, clayff
except ImportError:
    pass

# ===== Charge functions =====
try:
    from .charge import charge_minff, charge_clayff, balance_charges, assign_formal_charges
except ImportError:
    pass

# Version information
__version__ = "0.5.0"

# Expose key functions at the package level
__all__ = [
    'import_pdb', 'import_gro', 'import_xyz', 'import_auto',
    'write_pdb', 'write_gro', 'write_xyz', 'write_auto',
    'write_itp', 'write_psf', 'write_lmp',
    'element', 'radius', 'mass',
    'dist_matrix', 'cell_list_dist_matrix', 'bond_angle',
    'Box_dim2Cell', 'Cell2Box_dim',
    'cartesian_to_fractional', 'fractional_to_cartesian', 'wrap_coordinates',
    'triclinic_to_orthogonal', 'orthogonal_to_triclinic', 'get_orthogonal_box', 'get_cell_vectors',
    'direct_cartesian_to_fractional', 'direct_fractional_to_cartesian',
    'replicate_system', 'translate', 'update',
    'minff', 'clayff',
    'charge_minff', 'charge_clayff', 'balance_charges',
    'assign_formal_charges'
]
