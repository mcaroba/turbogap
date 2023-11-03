#!/usr/bin/python

from ase.io import read, write
import sys

if len(sys.argv) >= 2:
    atoms_name = sys.argv[1]
else:
    atoms_name = "final_config.xyz"

final_atoms = read( "trajectory_out.xyz", index=-1, format="extxyz" )
write(atoms_name, final_atoms, format = 'extxyz'  )
