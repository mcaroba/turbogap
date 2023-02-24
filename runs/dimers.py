import numpy as np
from ase.io import read, write

d_arr = np.linspace(0.1,9,90)

for d in d_arr:
    atoms = read("dimer_1.0.xyz")
    a = 1.+d
    s = "%.1f" % a
    filename = "dimer_" + s + ".xyz"
    atoms.positions[1][0] = atoms.positions[1][0]+d
    write(filename, atoms)
