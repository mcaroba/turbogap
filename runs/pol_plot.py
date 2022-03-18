import numpy as np
import matplotlib.pyplot as plt

rcut = np.loadtxt("cutoff_scs.dat")
pol = np.loadtxt("pol_scs.dat")
#pol2 = np.loadtxt("pol_scs2.dat")
#pol3 = np.loadtxt("pol_scs_cutoff.dat")

plt.figure(0)
plt.plot(rcut,pol,'b*-')
#plt.plot(rcut,pol2,'r*-')
#plt.plot(rcut[1:],pol3,'g*-')
plt.xlabel("SCS cutoff radius [Å]")
plt.ylabel("SCS polarizability of an atom in C60 (Hartree units)")
plt.show()



