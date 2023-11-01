from ase.io import read, write
import numpy as np

A1 = read("A.xyz")

d_arr = np.linspace(1,8,71)

print(d_arr)

for d in d_arr:
    A2 = read("A.xyz")
    s = "%.1f" %  d
    filename = "AA_" + s + ".xyz"
    print(filename)
    A2.positions[:,2] += d
    write(filename,A1+A2)
