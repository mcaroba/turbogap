import numpy as np
import matplotlib.pyplot as plt

pol = np.loadtxt("polynomial_degree.dat")

alpha = 8.1723730292725634

pol = pol[:12,:]

plt.figure(0)
plt.plot([pol[0,0]-0.5,pol[-1,0]+0.5],[alpha,alpha],'r-',label="Exact polarizability")
plt.plot(pol[:,0],pol[:,1],'b*-',label="Polynomial fit")
plt.xlim([pol[0,0]-0.5,pol[-1,0]+0.5])
plt.legend()
plt.xlabel("Degree of polynomial")
plt.ylabel("Polarizibality [Hartree units]")
plt.title("Polarizability of an atom in C60")
plt.show()
