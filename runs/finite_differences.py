import numpy as np
from ase.io import read, write

atoms = read("2829_supercell.xyz")

for i in range(10):
    for c3 in range(3):
        atoms.positions[i,c3] = atoms.positions[i,c3] + 0.001
        filename = "2829_%d_p%d.xyz" % (i+1, c3+1)
        print(filename)
        write(filename,atoms)
        atoms.positions[i,c3] = atoms.positions[i,c3] - 0.001
