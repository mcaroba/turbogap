import numpy as np
import matplotlib.pyplot as plt

d = np.linspace(1,8,36)
d_005 = np.linspace(3,4,21)
ene = np.loadtxt("bilayer_graphene.dat")
ene_ts = np.loadtxt("bilayer_graphene_ts.dat")
ene_005 = np.loadtxt("bilayer_graphene_0.05.dat")
ene_ts_005 = np.loadtxt("bilayer_graphene_ts_0.05.dat")

graphene_energy = -1754.44380807
graphene_energy_ts = -1749.16165927
total_energy = 2*graphene_energy
total_energy_ts = 2*graphene_energy_ts

ene = ene-total_energy
ene_ts = ene_ts - total_energy_ts
ene_005 = ene_005-total_energy
ene_ts_005 = ene_ts_005-total_energy_ts

print(d_005[np.argmin(ene_005)])
print(d_005[np.argmin(ene_ts_005)])

plt.plot(d[8:11],ene[8:11],'b.-',label='MBD')
plt.plot(d[8:11],ene_ts[8:11],'r.-',label='TS')
plt.plot(d_005,ene_005,'b.-')
plt.plot(d_005,ene_ts_005,'r.-')
plt.plot(d[15:],ene[15:],'b.-')
plt.plot(d[15:],ene_ts[15:],'r.-')
plt.axvline(x = d_005[np.argmin(ene_005)], color = 'b', linestyle=':', label = 'MBD minimum')
plt.axvline(x = d_005[np.argmin(ene_ts_005)], color = 'r', linestyle=':', label = 'TS minimum')
plt.axhline(y=0., color='g', linestyle='--')
plt.xlim([d[8],d[-1]])
plt.xlabel("Interlayer distance [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("Interaction energy of graphene sheets in bilayer graphene (5 Å SCS cut-off)")
plt.legend()
plt.show()
