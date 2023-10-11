import numpy as np
from ase.io import read, write

# Interlayer distance
d_arr = np.linspace(1.6,5.0,35)

A = read("A.xyz")

for d in d_arr:
    B = read("B.xyz")
    B.positions[:,2] += d
    atoms = A+B
    s = "%.1f" % d
    filename = "bilayer_" + s + ".xyz"
    write(filename, atoms)
