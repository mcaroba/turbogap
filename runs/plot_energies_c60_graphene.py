import numpy as np
import matplotlib.pyplot as plt

d = np.linspace(4,10,31)
ene = np.loadtxt("c60_graphene.dat")
ene_ts = np.loadtxt("c60_graphene_ts.dat")

c60_energy = -525.88659135
graphene_energy = -1754.44380807
c60_energy_ts = -524.53100704
graphene_energy_ts = -1749.16165927
total_energy = c60_energy+graphene_energy
total_energy_ts = c60_energy_ts+graphene_energy_ts

ene = ene-total_energy
ene_ts = ene_ts - total_energy_ts

plt.plot(d[5:],ene[5:],'b.-',label='MBD')
plt.plot(d[5:],ene_ts[5:],'r.-',label='TS')
plt.axhline(y=0., color='g', linestyle='--')
plt.xlim([d[5],d[-1]])
plt.xlabel("Center-sheet distance [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("Interaction energy of a C60 dimer and graphene sheet (5 Å SCS cut-off)")
plt.legend()
plt.show()
