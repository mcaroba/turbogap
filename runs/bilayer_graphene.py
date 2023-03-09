import numpy as np
from ase.io import read, write

d_arr = np.linspace(3.0,4.0,21)

atoms1 = read("graphene_A.xyz")

for d in d_arr:
    atoms2 = read("graphene_B.xyz")
    atoms2.positions[:,2] += d
    atoms = atoms1+atoms2
    s = "%.2f" % d
    filename = "bilayer_graphene_" + s + ".xyz"
    write(filename, atoms)
