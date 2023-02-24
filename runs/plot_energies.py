import numpy as np
import matplotlib.pyplot as plt

d = np.linspace(8,15,36)
ene = np.loadtxt("c60_energies.dat")
ene_ts = np.loadtxt("c60_energies_ts.dat")

c60_energy = -525.88826598
c60_energy_ts = -524.53100703
total_energy = c60_energy*2
total_energy_ts = c60_energy_ts*2

ene = ene-total_energy
ene_ts = ene_ts - total_energy_ts

plt.plot(d[5:],ene[5:],'b.-',label='MBD')
plt.plot(d[5:],ene_ts[5:],'r.-',label='TS')
plt.axhline(y=0., color='g', linestyle='--')
plt.xlim([d[5],d[-1]])
plt.xlabel("Intermolecular center-center distance [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("Interaction energy of a C60 dimer")
plt.legend()
plt.show()
