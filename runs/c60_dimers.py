import numpy as np
from ase.io import read, write

d_arr = np.linspace(8,15,36)

atoms1 = read("c60.xyz")

for d in d_arr:
    atoms2 = read("c60.xyz")
    atoms2.positions[:,0] += d
    atoms = atoms1+atoms2
    s = "%.1f" % d
    filename = "c60_dimer_" + s + ".xyz"
    write(filename, atoms)
