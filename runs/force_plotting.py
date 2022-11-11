import numpy as np
import matplotlib.pyplot as plt

E_p = np.loadtxt("energies_p.dat")
E_m = np.loadtxt("energies_m.dat")
f = np.loadtxt("forces.dat")
f = f[:,2]

f_FD = -(E_p-E_m)/0.002

RMSE = np.sqrt(((f-f_FD)**2).sum()/f.shape)

print(RMSE)

plt.figure(0)
plt.scatter(f_FD,f,alpha=0.5)
plt.plot([-0.55,0.55],[-0.55,0.55],'r-')
plt.xlim([-0.55,0.55])
plt.ylim([-0.55,0.55])
plt.xlabel("Finite difference MBD forces")
plt.ylabel("Analytical MBD forces")
plt.title("Forces with frozen atom approximation")
plt.text(0, -0.4, "RMSE = %f" % RMSE[0], bbox=dict(facecolor='white', alpha=0.5))
plt.show()
