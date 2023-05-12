import numpy as np
from ase.io import read, write

d_arr = np.linspace(0.1,9.0,90)

for d in d_arr:
    atoms = read("dimer.xyz")
    atoms.positions[1,0] += d
    r = d+1.0
    s = "%.1f" % r
    filename = "dimer_" + s + ".xyz"
    write(filename, atoms)
