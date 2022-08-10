import numpy as np
import matplotlib.pyplot as plt

a_full = np.loadtxt("alpha_SCS_full.dat")
a_local = np.loadtxt("alpha_SCS_local_inv.dat")
a_appr = np.loadtxt("alpha_SCS_local_pol.dat")

a_full = a_full[:,1]
a_local = a_local[:,1]
a_appr = a_appr[:,1]


RMSE = np.sqrt( ((a_full-a_appr)**2).sum()/a_full.shape  )

print("RMSE", RMSE)

plt.figure(0)
plt.scatter(a_full,a_local,alpha=0.25)
plt.plot([4,12],[4,12],'r-')
plt.xlabel("Exact polarizabilities")
plt.ylabel("Locally inverted polarizabilities")
plt.xlim([4,12])
plt.ylim([4,12])

plt.show()

plt.figure(1)
plt.plot([4,12],[4,12],'r-')
plt.scatter(a_full,a_appr,alpha=0.25)
#plt.plot([4,12],[4,12],'r-')
plt.title("a-C structure SCS polarizabilities",fontsize=12)
plt.xlabel("Exact polarizabilities [a$_0^3$]",fontsize=12)
plt.ylabel("Local approximation [a$_0^3$]",fontsize=12)
plt.xlim([4,12])
plt.ylim([4,12])
plt.text(9, 5, 'RMSE = %.3f' % RMSE, fontsize = 12, 
         bbox = dict(facecolor = 'white', alpha = 0.5))

plt.show()

plt.figure(2)
plt.scatter(a_local,a_appr,alpha=0.25)
plt.plot([4,12],[4,12],'r-')
plt.xlabel("Locally inverted polarizabilities")
plt.ylabel("Locally polynomial approximation")
plt.xlim([4,12])
plt.ylim([4,12])

plt.show()

