import numpy as np
import matplotlib.pyplot as plt

d = np.linspace(8,15,36)
ene_ts = np.loadtxt("c60_energies_ts.dat")
ene_none = np.loadtxt("c60_energies_none.dat")

c60_ref = np.loadtxt("c60_ref.dat")
c60_energy_ts = -524.53100704
c60_energy_none = -521.67745459
total_energy = c60_ref*2
total_energy_ts = c60_energy_ts*2
total_energy_none = c60_energy_none*2

scs = np.linspace(0,17,35)

ene = np.zeros([len(scs),len(d)])

ene_ts = ene_ts-total_energy_ts
ene_none = ene_none-total_energy_none


plt.plot(d[6:],ene_none[6:],'r.-',label="None")
plt.plot(d[6:],ene_ts[6:],'b.-',label="TS")

for i in range(len(scs)):
    s = "%.1f" % scs[i]
    filename = "c60_energies_" + s + ".dat"
    ene[i,:] = np.loadtxt(filename)
    ene[i,:] = ene[i,:]-total_energy[i]
    x = float(i)/len(scs)
    plt.plot(d[6:],ene[i,6:],color=(0,x,1.0-x),linestyle='solid',marker='.',label="SCS cutoff " + str(scs[i]) + " Å")

plt.xlim([d[6],d[-1]])
plt.title("C60 dimer with varying degrees of vdW correction (MBD cutoff = SCS cutoff + 2 Å)")
plt.xlabel("Intermolecular center-center distance [Å]")
plt.ylabel("Interaction energy [eV]")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
plt.show()

t_scs = np.zeros([len(scs),len(d)])

for i in range(len(scs)):
    s = "%.1f" % scs[i]
    filename = "scs_timing_" + s + ".dat"
    t_scs[i,:] = np.loadtxt(filename)
    x = float(i)/len(scs)
    plt.plot(d,t_scs[i,:],color=(0,x,1.0-x),linestyle='solid',marker='.',label="SCS cutoff " + str(scs[i]) + " Å")

plt.xlim([d[0],d[-1]])
plt.ylabel("CPU time [s]")
plt.xlabel("Intermolecular center-center distance [Å]")
plt.title("SCS timings")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
plt.show()

t_com = np.zeros([len(scs),len(d)])

for i in range(len(scs)):
    s = "%.1f" % scs[i]
    filename = "communication_timing_" + s + ".dat"
    t_com[i,:] = np.loadtxt(filename)
    x = float(i)/len(scs)
    plt.plot(d,t_com[i,:],color=(0,x,1.0-x),linestyle='solid',marker='.',label="SCS cutoff " + str(scs[i]) + " Å")

plt.xlim([d[0],d[-1]])
plt.ylabel("CPU time [s]")
plt.xlabel("Intermolecular center-center distance [Å]")
plt.title("Communication timings")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
plt.show()

t_mbd = np.zeros([len(scs),len(d)])

for i in range(len(scs)):
    s = "%.1f" % scs[i]
    filename = "mbd_timing_" + s + ".dat"
    t_mbd[i,:] = np.loadtxt(filename)
    x = float(i)/len(scs)
    plt.plot(d,t_mbd[i,:],color=(0,x,1.0-x),linestyle='solid',marker='.',label="SCS cutoff " + str(scs[i]) + " Å")

plt.xlim([d[0],d[-1]])
plt.ylabel("CPU time [s]")
plt.xlabel("Intermolecular center-center distance [Å]")
plt.title("MBD timings")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
plt.show()

scs_f = np.linspace(0,7,15)
t_scs_f = np.zeros([len(scs_f),len(d)])
t_com_f = np.zeros([len(scs_f),len(d)])
t_mbd_f = np.zeros([len(scs_f),len(d)])
ene_ac = np.zeros(len(scs_f))
t_scs_ac = np.zeros(len(scs_f))
t_com_ac = np.zeros(len(scs_f))
t_mbd_ac = np.zeros(len(scs_f))
t_scs_f_ac = np.zeros(len(scs_f))
t_com_f_ac = np.zeros(len(scs_f))
t_mbd_f_ac = np.zeros(len(scs_f))

for i in range(len(scs_f)):
    s = "%.1f" % scs_f[i]
    filename = "scs_timing_forces_" + s + ".dat"
    t_scs_f[i,:] = np.loadtxt(filename)
    filename = "communication_timing_forces_" + s + ".dat"
    t_com_f[i,:] = np.loadtxt(filename)
    filename = "mbd_timing_forces_" + s + ".dat"
    t_mbd_f[i,:] = np.loadtxt(filename)
    filename = "aC_energies_" + s + ".dat"
    ene_ac[i] = np.loadtxt(filename)
    filename = "aC_scs_timing_" + s + ".dat"
    t_scs_ac[i] = np.loadtxt(filename)
    filename = "aC_communication_timing_" + s + ".dat"
    t_com_ac[i] = np.loadtxt(filename)
    filename = "aC_mbd_timing_" + s + ".dat"
    t_mbd_ac[i] = np.loadtxt(filename)
    filename = "aC_scs_timing_forces_" + s + ".dat"
    t_scs_f_ac[i] = np.loadtxt(filename)
    filename = "aC_communication_timing_forces_" + s + ".dat"
    t_com_f_ac[i] = np.loadtxt(filename)
    filename = "aC_mbd_timing_forces_" + s + ".dat"
    t_mbd_f_ac[i] = np.loadtxt(filename)




plt.plot(scs,np.mean(t_scs,axis=1),'b.-',label="No forces")
plt.plot(scs_f,np.mean(t_scs_f,axis=1),'r.-',label="With forces")
plt.title("C60 dimer, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("Average CPU time for SCS cycle [s]")
plt.xlim([scs[0],scs[-1]])
plt.legend()
plt.show()

plt.plot(scs,np.mean(t_com,axis=1),'b.-',label="No forces")
plt.plot(scs_f,np.mean(t_com_f,axis=1),'r.-',label="With forces")
plt.title("C60 dimer, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("Average CPU time for communication [s]")
plt.xlim([scs[0],scs[-1]])
plt.legend()
plt.show()

plt.plot(scs,np.mean(t_mbd,axis=1),'b.-',label="No forces")
plt.plot(scs_f,np.mean(t_mbd_f,axis=1),'r.-',label="With forces")
plt.title("C60 dimer, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("Average CPU time for MBD calculation [s]")
plt.xlim([scs[0],scs[-1]])
plt.legend()
plt.show()


plt.plot(scs_f,t_scs_ac,'b.-',label="No forces")
plt.plot(scs_f,t_scs_f_ac,'r.-',label="With forces")
plt.title("Amorphous carbon, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("CPU time for SCS cycle [s]")
plt.xlim([scs_f[0],scs_f[-1]])
plt.legend()
plt.show()

plt.plot(scs_f,t_com_ac,'b.-',label="No forces")
plt.plot(scs_f,t_com_f_ac,'r.-',label="With forces")
plt.title("Amorphous carbon, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("CPU time for communication [s]")
plt.xlim([scs_f[0],scs_f[-1]])
plt.legend()
plt.show()

plt.plot(scs_f,t_mbd_ac,'b.-',label="No forces")
plt.plot(scs_f,t_mbd_f_ac,'r.-',label="With forces")
plt.title("Amorphous carbon, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("CPU time for MBD calculation [s]")
plt.xlim([scs_f[0],scs_f[-1]])
plt.legend()
plt.show()

