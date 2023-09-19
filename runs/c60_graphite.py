import numpy as np
from ase.io import read, write

atoms1 = read("graphene_2.xyz")

# Distance of c60 center of mass from top layer: 3.3404496666666645 Å

d_arr = np.linspace(0,8,21)

for d in d_arr:
    atoms2 = read("c60_graphite.xyz")
    atoms2.positions[:,2] += d
    atoms = atoms1+atoms2
    s = "%.1f" % d
    filename = "c60_graphite_" + s + ".xyz"
    write(filename, atoms)

