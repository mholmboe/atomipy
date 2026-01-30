"""
atomipy: The atom Toolbox in Python

This package provides tools for working with molecular structures, particularly
focused on mineral systems with periodic boundary conditions. It supports both
orthogonal and triclinic simulation cells, and provides efficient distance
calculations, bond/angle detection, structure manipulation, and X-ray diffraction
simulation functions.

Key features:
- Import/export of common molecular file formats (PDB, GRO, XYZ)
- Topology file generation for molecular dynamics simulations (GROMACS, NAMD, LAMMPS)
- GROMACS n2t file generation for atom name to type conversions (gmx x2top)
  with periodic minimum-image distances and environment deduplication
- Distance calculations with efficient algorithms for periodic systems
- Bond/angle/dihedral/1-4 pair detection
- Coordinate transformations between Cartesian and fractional representations
- Support for both orthogonal and triclinic simulation cells
- Structure manipulation (replication, translation)
- Force field implementations for mineral systems (CLAYFF, MINFF)
- High-performance X-ray diffraction pattern simulation with vectorized calculations
- Professional XRD plotting with Miller indices and multiplicities
- Automatic ionic scattering factor selection and MATLAB compatibility
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
# Itp file import helper
from .import_top import import_itp as import_itp_topology

# ===== Atom property functions =====
from .element import element

from .radius import radius

try:
    from .mass import mass, set_atomic_masses, com
except ImportError:
    pass

# ===== Structure analysis functions =====
from .dist_matrix import dist_matrix
from .cell_list_dist_matrix import cell_list_dist_matrix
from .bond_angle import bond_angle, bond_angle_dihedral

# ===== Cell and coordinate transformation functions =====
from . import cell_utils
normalize_box = cell_utils.normalize_box
Box_dim2Cell = cell_utils.Box_dim2Cell
Cell2Box_dim = cell_utils.Cell2Box_dim

# Coordinate transformation module
from . import transform
cartesian_to_fractional = transform.cartesian_to_fractional
fractional_to_cartesian = transform.fractional_to_cartesian
wrap_coordinates = transform.wrap_coordinates
wrap = transform.wrap
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
rotate = move.rotate
place = move.place
center = move.center

from . import add
update = add.update

# ===== General functions =====
from . import general
scale = general.scale

# ===== Build functions =====
from . import build
substitute = build.substitute
molecule = build.molecule
merge = build.merge
slice = build.slice
ionize = build.ionize
insert = build.insert
add_H_atom = build.add_H_atom
adjust_H_atom = build.adjust_H_atom
adjust_Hw_atom = build.adjust_Hw_atom

# ===== Resname functions =====
from .resname import assign_resname

# ===== Solvent functions =====
from .solvent import find_H2O, solvate

# ===== Force field functions =====
try:
    from .forcefield import minff, clayff, write_n2t, get_structure_stats
except ImportError:
    pass

# ===== Forcefield parameters loading =====
try:
    from .ffparams import load_forcefield, list_blocks as list_ff_blocks, get_ffparams_dir
except ImportError:
    pass

# ===== Charge functions =====
try:
    from .charge import charge_minff, charge_clayff, balance_charges, assign_formal_charges
except ImportError:
    pass

# ===== Bond valence functions =====
try:
    from .bond_valence import compute_bvs, global_instability_index, load_bv_params, load_shannon_radii, bond_valence, summarize_bvs
except ImportError:
    pass
try:
    from .bond_valence import analyze_bvs, conf2bvs
except ImportError:
    pass
try:
    from .size import get_radius, bond_distance
except ImportError:
    pass

# ===== Diffraction functions =====
try:
    from .diffraction import xrd, occupancy_atom, atomic_scattering_factors, calculate_multiplicity, bragg_law
except ImportError:
    pass

# Version information
__version__ = "0.92"

# Expose key functions at the package level
__all__ = [
    'import_pdb', 'import_gro', 'import_xyz', 'import_auto',
    'write_pdb', 'write_gro', 'write_xyz', 'write_auto',
    'write_itp', 'write_psf', 'write_lmp', 'import_itp_topology',
    'element', 'radius', 'mass', 'set_atomic_masses', 'com',
    'dist_matrix', 'cell_list_dist_matrix', 'bond_angle', 'bond_angle_dihedral', 'find_H2O',
    'normalize_box','Box_dim2Cell', 'Cell2Box_dim',
    'cartesian_to_fractional', 'fractional_to_cartesian', 'wrap', 'wrap_coordinates',
    'triclinic_to_orthogonal', 'orthogonal_to_triclinic', 'get_orthogonal_box', 'get_cell_vectors',
    'direct_cartesian_to_fractional', 'direct_fractional_to_cartesian',
    'replicate_system', 'translate', 'rotate', 'place', 'center', 'update', 'scale',
    'substitute', 'molecule', 'merge', 'slice', 'solvate', 'ionize', 'insert',
    'add_H_atom', 'adjust_H_atom', 'adjust_Hw_atom',
    'assign_resname',
    'minff', 'clayff', 'write_n2t',
    'load_forcefield', 'list_ff_blocks',
    'charge_minff', 'charge_clayff', 'balance_charges', 'assign_formal_charges',
    'compute_bvs', 'global_instability_index', 'load_bv_params', 'load_shannon_radii', 'bond_valence', 'summarize_bvs',
    'analyze_bvs', 'conf2bvs',
    'get_radius', 'bond_distance',
    'xrd', 'occupancy_atom', 'atomic_scattering_factors', 'calculate_multiplicity', 'bragg_law'
]
