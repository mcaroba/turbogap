import numpy as np
import matplotlib.pyplot as plt

d = np.linspace(1,9,81)
ene = np.loadtxt("dimer_vdw_energies.dat")
#ene_ts = np.loadtxt("c60_energies_ts.dat")

#c60_energy = -525.89797570
#c60_energy = -525.79666221
#c60_energy_ts = -524.53100704
#c60_energy_ts = -524.53100703
#total_energy = c60_energy*2
#total_energy_ts = c60_energy_ts*2

#ene = ene-total_energy
#ene_ts = ene_ts - total_energy_ts

plt.plot(d,ene,'b.-')
plt.axhline(y=0., color='g', linestyle='--')
plt.xlim([d[0],d[-1]])
plt.xlabel("Dimer distance [Å]")
plt.ylabel("MBD energy [eV]")
plt.title("MBD energy of a carbon dimer (5 Å SCS cut-off)")
plt.legend()
plt.show()
