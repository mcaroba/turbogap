import numpy as np
import matplotlib.pyplot as plt
from ase.units import Bohr

a_p = np.loadtxt("alpha_p_C60.dat")
a_m = np.loadtxt("alpha_m_C60.dat")
da = np.loadtxt("dalpha_C60.dat")

da = da[:,2]
da_fd = (a_p[:,1]-a_m[:,1])/0.002*Bohr

print(da_fd)
print(da)

plt.figure(0)
plt.plot([-1.5,1.5],[-1.5,1.5],'r-')
plt.scatter(da_fd,da,alpha=0.5)
plt.xlabel("Finite difference polarizability gradients")
plt.ylabel("Analytical polarizability gradients")
plt.title("Polarizability gradients of C60 w.r.t. atom 1 (inc. Hirshfeld)")
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])


RMSE=np.sqrt(((da-da_fd)**2).sum()/da.shape[0])

plt.text(0.4, -1.2, 'RMSE = %f' % RMSE, fontsize = 12, 
         bbox = dict(facecolor = 'white', alpha = 0.5))


plt.show()
