import numpy as np
import matplotlib.pyplot as plt

n_order = np.array([2,3,4,5,6])

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)


ax = plt.subplot(111)

forces = np.zeros([len(n_order),180])

for i in range(len(n_order)):
    #s = "%.1f" % scs[i]
    filename = "c60_forces_" + str(n_order[i]) + ".dat"
    forces[i,:] = np.loadtxt(filename)
    #x = float(i)/len(n_order)

for n in range(180):
    ax.plot(n_order,forces[:,n],color=(0,0,1),linestyle='solid',marker='.')

plt.xlabel("Polynomial order")
plt.ylabel("Forces [eV/Å]")
plt.show()

plt.plot(n_order,np.mean(forces,axis=1))
plt.show()
