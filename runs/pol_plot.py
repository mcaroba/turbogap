import numpy as np
import matplotlib.pyplot as plt

rcut = np.loadtxt("cutoff_scs.dat")
pol = np.loadtxt("pol_scs.dat")
pol2 = np.loadtxt("pol_scs2.dat")

plt.figure(0)
plt.plot(rcut,pol,'b*-')
plt.plot(rcut,pol2,'r*-')
plt.xlabel("SCS cutoff radius [Å]")
plt.ylabel("SCS polarizability of an atom in C60 (Hartree units)")
plt.show()



