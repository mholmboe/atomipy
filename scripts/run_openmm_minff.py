import os
import atomipy as ap

def main():
    print("Loading test structure...")
    struct_path = '/Users/miho0052/Dropbox/Coding/Windsurf/atominpython/Pyrophyllite_GII_0.071.gro'
    if not os.path.exists(struct_path):
        print(f"Test structure {struct_path} not found. Run this from the scripts/ directory.")
        return

    atoms, Box = ap.import_gro(struct_path)
    
    print("Replicating system 6x4x3...")
    atoms, Box, _ = ap.replicate_system(atoms, Box, replicate=[6, 4, 3])
    
    print("Assigning MINFF forcefield...")
    atoms = ap.minff(atoms, Box=Box)

    print("Writing GROMACS topology and coordinates...")
    ap.write_top(atoms, Box=Box, file_path='system_minff.top')
    ap.write_gro(atoms, Box=Box, file_path='system_minff.gro')

    inc_dir = ap.get_ffparams_dir()

    print("Loading into OpenMM...")
    topology, system, positions = ap.load_minff_into_openmm(
        top_path='system_minff.top',
        gro_path='system_minff.gro',
        defines=['GMINFF_k500', 'OPC3_IOD_LM', 'OPC3'],
        include_dir=inc_dir,
        nonbonded_cutoff_nm=0.9, # To ensure it's < half box size
        constraints=None,        # MINFF needs flexible bonds
        rigid_water=True,
    )

    print(f"System has {system.getNumParticles()} particles.")

    print("Setting up simulation...")
    from openmm.app import Simulation, PDBReporter, StateDataReporter
    from openmm import LangevinIntegrator, Platform
    from openmm.unit import kelvin, picosecond, femtoseconds
    
    integrator = LangevinIntegrator(298.15*kelvin, 1/picosecond, 1.0*femtoseconds)
    try:
        platform = Platform.getPlatformByName('CPU')
    except:
        platform = Platform.getPlatformByName('Reference')
        
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    
    state = simulation.context.getState(getEnergy=True)
    print(f"Initial energy: {state.getPotentialEnergy()}")

if __name__ == '__main__':
    main()
