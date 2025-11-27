#!/usr/bin/env python3
"""Minimal MINFF workflow for a single structure.

This script loads a PDB structure, assigns MINFF atom types and charges,
logs structure statistics, and writes both a GROMACS ``.itp`` topology file
and a typed PDB.
"""

from atomipy import import_pdb, minff, write_itp, write_pdb

pdb_file = 'Kaolinite_GII_0.0487.pdb'
pdb_file_out = 'minff_Kao.pdb'
itp_file_out = 'minff_Kao.itp'
log_file_out = 'minff_Kao.log'

atoms, cell = import_pdb(file_path=pdb_file)
minff_atoms = minff(atoms, Box=cell, log=True)

write_itp(minff_atoms, Box=cell, file_path=itp_file_out)
write_pdb(minff_atoms, Box=cell, file_path=pdb_file_out)
    
