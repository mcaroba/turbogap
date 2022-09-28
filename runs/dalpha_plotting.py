import numpy as np
import matplotlib.pyplot as plt
from ase.units import Bohr

a_p = np.loadtxt("alpha_p.dat")
a_m = np.loadtxt("alpha_m.dat")
da = np.loadtxt("dalpha.dat")

print(da.shape[0])

dalpha = np.zeros([1000])

#print(dalpha)
k=0
for i in range(1,1001):
    if (k < da.shape[0]):
        if (int(da[k,0]) == i):
            print(k,i)
            dalpha[i-1] = da[k,2]
            k = k+1


print(dalpha)

dalpha_fd = (a_p[:,1]-a_m[:,1])/0.002*Bohr

print(dalpha_fd)

plt.figure(0)
plt.plot([-0.7,0.9],[-0.7,0.9],'r-')
plt.scatter(dalpha_fd,dalpha,alpha=0.5)
plt.xlabel("Finite difference polarizability gradients")
plt.ylabel("Analytical polarizability gradients")
plt.title("Polarizability gradients of amorphous carbon w.r.t. atom 4")
plt.xlim([-0.7,0.9])
plt.ylim([-0.7,0.9])

RMSE=np.sqrt(((dalpha-dalpha_fd)**2).sum()/da.shape[0])

plt.text(0.2, -0.5, 'RMSE = %f' % RMSE, fontsize = 12, 
         bbox = dict(facecolor = 'white', alpha = 0.5))

plt.show()


print(np.amax(np.abs(dalpha-dalpha_fd)))
print(np.argmax(np.abs(dalpha-dalpha_fd)))

print(da[:,0])

print(dalpha[3])
print(dalpha_fd[3])
