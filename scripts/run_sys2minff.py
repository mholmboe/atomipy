#!/usr/bin/env python3
"""
run_sys2minff.py - Convert an interface structure to MINFF for GROMACS. 

Processes a molecular system containing minerals, ions, and water,
assigns MINFF atom types, and generates GROMACS-compatible output files.

Usage:
    python scripts/run_sys2minff.py structure.pdb|.gro|.xyz
"""

import sys
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import atomipy as ap

if len(sys.argv) < 2:
    print("Usage: python scripts/run_sys2minff.py structure.pdb|.gro|.xyz")
    sys.exit(1)

# Import and prepare structure
atoms, Box_dim = ap.import_auto(sys.argv[1])
atoms = ap.element(atoms)

# Separate water, ions, and minerals
SOL, noSOL = ap.find_H2O(atoms, Box_dim)
# Standardize resnames
noSOL = ap.assign_resname(noSOL)

# Separate and consolidate mineral atoms (to allow bonding across residues if molid=1)
MIN = [a for a in noSOL if a.get('resname') == 'MIN']
OTHER = [a for a in noSOL if a.get('resname') != 'MIN']

if MIN:
    # Force all mineral atoms into one molecule (molid=1)
    MIN = ap.update(MIN, molid=1)

# Assemble system robustly (preserve all atoms)
# Order: Ions/Other first, then Mineral, then Water
System = ap.update(OTHER, MIN, SOL)

# Assign MINFF atom types to the entire system
System = ap.minff(System, Box_dim)

# Extract MIN atoms from System for ITP (based on resname)
MIN = [a for a in System if a.get('resname') == 'MIN']

# Write ITP files for MIN part only
ap.write_itp(MIN, Box=Box_dim, file_path='minff.itp')

# Write PSF file for the entire system (MIN + ION/Na + SOL/Ow,Hw), with bimodal detection and angle filtering
ap.write_psf(System, Box=Box_dim, file_path='minff.psf', detect_bimodal=True, max_angle=150)

# Load GMINFF forcefield parameters for LAMMPS Pair Coeffs, with and without bimodal detection of atomtype triplet angle terms
ff_params = ap.load_forcefield('GMINFF/gminff_all.json', blocks=['GMINFF_k500', 'OPC3_HFE_LM', 'OPC3'])
ap.write_lmp(System, Box=Box_dim, file_path='minff.data', forcefield=ff_params,detect_bimodal=True)
ap.write_lmp(System, Box=Box_dim, file_path='minff_no150angles.data', forcefield=ff_params, detect_bimodal=True, max_angle=150)

# Write full system GRO and LAMMPS data file
ap.write_gro(System, Box=Box_dim, file_path='preem.gro')
ap.write_pdb(System, Box=Box_dim, file_path='preem.pdb')

# Write structure statistics log
ap.get_structure_stats(System, Box=Box_dim)

print("Written: minff.itp, minff.psf, minff.data, minff_no150angles.data, preem.gro, preem.pdb, output.log")
