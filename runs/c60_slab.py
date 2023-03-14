import numpy as np
from ase.io import read, write

d_arr = np.linspace(0.0,8.0,41)

slab = read("slab_small.xyz")

for d in d_arr:
    c60 = read("c60_slab_cell.xyz")
    c60.positions[:,2] += d
    atoms = slab+c60
    s = "%.1f" % d
    filename = "c60_slab_" + s + ".xyz"
    write(filename, atoms)
