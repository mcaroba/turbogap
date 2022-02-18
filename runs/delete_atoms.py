import numpy as np
from ase.io import read, write
from ase.visualize import view

atoms = read("atoms_extended.xyz")

a = atoms.copy()
list = []
full_list = [j for j in range(60)]

for i in range(a.positions.shape[0]):
    if ( a.get_distance(0,i) > 3.2 ):
        list.append(i)
        full_list.remove(i)

print(a.get_distances(0,full_list))

del a[list]

print(list)
print(full_list)

print(a.get_distances(0,[1,2,3,4,5,6,7,8,9,10,11]))

view(a)

write("POSCAR",a)
