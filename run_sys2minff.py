#!/usr/bin/env python3
"""
run_sys2minff.py - Convert an interface structure to MINFF for GROMACS.

Processes a molecular system containing minerals, ions, and water,
assigns MINFF atom types, and generates GROMACS-compatible output files.

Usage:
    python run_sys2minff.py structure.pdb|.gro|.xyz
"""

import sys
import atomipy as ap

if len(sys.argv) < 2:
    print("Usage: python run_sys2minff.py structure.pdb|.gro|.xyz")
    sys.exit(1)

# Import and prepare structure
atoms, Box_dim = ap.import_auto(sys.argv[1])
atoms = ap.element(atoms)

# Separate water, ions, and minerals
SOL, noSOL = ap.find_H2O(atoms, Box_dim)
noSOL = ap.assign_resname(noSOL)
IONS = [a for a in noSOL if a.get('resname') == 'ION']
MIN = [a for a in noSOL if a.get('resname') == 'MIN']

# Combine into full system for charge assignment
if MIN:
    MIN = ap.molecule(MIN, molid=1, resname='MIN')

System = ap.update(MIN, IONS, SOL)

# Assign MINFF atom types to the entire system
System = ap.minff(System, Box_dim)

# Extract MIN atoms from System for ITP (based on resname)
MIN = [a for a in System if a.get('resname') == 'MIN']

# Write ITP and PSF files for MIN part only
ap.write_itp(MIN, Box=Box_dim, file_path='minff.itp')
ap.write_psf(MIN, Box=Box_dim, file_path='minff.psf')

# Load GMINFF forcefield parameters for LAMMPS Pair Coeffs
ff_params = ap.load_forcefield('GMINFF/gminff_all.json', blocks=['GMINFF_k500', 'OPC3_HFE_LM', 'OPC3'])
ap.write_lmp(System, Box=Box_dim, file_path='minff.data', forcefield=ff_params)

# Write full system GRO and LAMMPS data file
ap.write_gro(System, Box=Box_dim, file_path='preem.gro')
ap.write_pdb(System, Box=Box_dim, file_path='preem.pdb')

# Write structure statistics log
ap.get_structure_stats(System, Box=Box_dim)

print("Written: minff.itp, minff.psf, minff.data, preem.gro, preem.pdb, output.log")
