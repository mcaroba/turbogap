import numpy as np
import matplotlib.pyplot as plt

a_VASP = np.loadtxt("alpha_SCS_2829.dat")
a_turbogap = np.loadtxt("alpha_SCS_2829_turbogap.dat")

plt.figure(0)
plt.scatter(a_VASP,a_turbogap,alpha=0.5)
plt.plot([4.5,11],[4.5,11],'r-')
plt.xlim([4.5,11])
plt.ylim([4.5,11])
plt.xlabel("VASP SCS polarizabilities [Bohr$^3$]")
plt.ylabel("Turbogap SCS polarizabilities [Bohr$^3$]")
plt.show()
