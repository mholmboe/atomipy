import atomipy as ap
from openmm.app import Simulation
from openmm import LangevinIntegrator, Platform
from openmm.unit import kelvin, picosecond, femtoseconds

struct_path = '/Users/miho0052/Dropbox/Coding/Windsurf/atominpython/Pyrophyllite_GII_0.071.gro'
atoms, Box = ap.import_gro(struct_path)
atoms, Box, _ = ap.replicate_system(atoms, Box, replicate=[2, 2, 1])
atoms = ap.minff(atoms, Box=Box)

ap.write_top(atoms, Box=Box, file_path='system_minff_2x2x1.top')
ap.write_gro(atoms, Box=Box, file_path='system_minff_2x2x1.gro')

inc_dir = ap.get_ffparams_dir()
topology, system, positions = ap.load_minff_into_openmm(
    top_path='system_minff_2x2x1.top',
    gro_path='system_minff_2x2x1.gro',
    defines=['GMINFF_k500', 'OPC3_IOD_LM', 'OPC3'],
    include_dir=inc_dir,
    nonbonded_cutoff_nm=0.4,
    constraints=None,
    rigid_water=True,
)

# Break down energy
for i, force in enumerate(system.getForces()):
    force.setForceGroup(i)

integrator = LangevinIntegrator(298.15*kelvin, 1/picosecond, 1.0*femtoseconds)
simulation = Simulation(topology, system, integrator, Platform.getPlatformByName('Reference'))
simulation.context.setPositions(positions)

print("Energy breakdown (2x2x1):")
for i, force in enumerate(system.getForces()):
    state = simulation.context.getState(getEnergy=True, groups={i})
    print(f"{force.__class__.__name__}: {state.getPotentialEnergy()}")
