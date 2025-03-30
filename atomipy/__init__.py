"""
atomipy: A Python toolbox for molecular structure analysis
"""

# File I/O functions
# Import functions from the new combined module
from . import import_module
import_pdb = import_module.pdb
import_gro = import_module.gro
import_auto = import_module.auto
from . import write_module
write_pdb = write_module.pdb
write_gro = write_module.gro
write_auto = write_module.auto
from . import write_itp

# Atom property functions
from .element import element
from .radius import radius
try:
    from .mass import mass
except ImportError:
    pass

# Structure analysis functions
from .dist_matrix import dist_matrix
from .bond_angle import bond_angle

# Cell and coordinate transformation functions
from . import cell_utils
from . import triclinic
from . import ortho
from . import fract
from . import replicate

# Force field functions
try:
    from .minff import minff
except ImportError:
    pass

# Charge functions
try:
    from .charge_minff import charge_minff, balance_charges as charge_balance
except ImportError:
    pass

try:
    from .charge_formal import assign_formal_charges
except ImportError:
    pass


