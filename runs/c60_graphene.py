import numpy as np
from ase.io import read, write

d_arr = np.linspace(4.0,8.0,41)

atoms1 = read("graphene_sheet.xyz")

for d in d_arr:
    atoms2 = read("c60_sheet.xyz")
    atoms2.positions[:,2] += d
    atoms = atoms1+atoms2
    s = "%.1f" % d
    filename = "c60_graphene_" + s + ".xyz"
    write(filename, atoms)
