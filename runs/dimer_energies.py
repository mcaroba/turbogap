import numpy as np
import matplotlib.pyplot as plt


r = np.linspace(5.0,8.0,61)
e = np.loadtxt("dimer_energies.dat")

plt.plot(r,e,'b*-')
plt.show()
