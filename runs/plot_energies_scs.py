import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

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

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
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

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
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
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

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
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

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

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
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

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

plt.plot(scs,np.mean(t_com,axis=1),'b.-',label="No forces")
plt.plot(scs_f,np.mean(t_com_f,axis=1),'r.-',label="With forces")
plt.title("C60 dimer, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("Average CPU time for communication [s]")
plt.xlim([scs[0],scs[-1]])
plt.legend()
plt.show()

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

plt.plot(scs,np.mean(t_mbd,axis=1),'b.-',label="No forces")
plt.plot(scs_f,np.mean(t_mbd_f,axis=1),'r.-',label="With forces")
plt.title("C60 dimer, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("Average CPU time for MBD calculation [s]")
plt.xlim([scs[0],scs[-1]])
plt.legend()
plt.show()

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

plt.plot(scs_f,ene_ac,'b.-')
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("MBD energy [eV]")
plt.title("Amorphous carbon MBD energy (125 atoms, MBD cutoff = SCS cutoff + 2Å)")
plt.xlim([scs_f[0],scs_f[-1]])
plt.show()

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

plt.plot(scs_f,t_scs_ac,'b.-',label="No forces")
plt.plot(scs_f,t_scs_f_ac,'r.-',label="With forces")
plt.title("Amorphous carbon, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("CPU time for SCS cycle [s]")
plt.xlim([scs_f[0],scs_f[-1]])
plt.legend()
plt.show()

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

plt.plot(scs_f,t_com_ac,'b.-',label="No forces")
plt.plot(scs_f,t_com_f_ac,'r.-',label="With forces")
plt.title("Amorphous carbon, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("CPU time for communication [s]")
plt.xlim([scs_f[0],scs_f[-1]])
plt.legend()
plt.show()

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

plt.plot(scs_f,t_mbd_ac,'b.-',label="No forces")
plt.plot(scs_f,t_mbd_f_ac,'r.-',label="With forces")
plt.title("Amorphous carbon, MBD cutoff = SCS cutoff + 2 Å")
plt.xlabel("SCS cutoff [Å]")
plt.ylabel("CPU time for MBD calculation [s]")
plt.xlim([scs_f[0],scs_f[-1]])
plt.legend()
plt.show()

plt.rcParams.update({'font.size': 22})

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

d = np.linspace(0,8,41)
scs = np.linspace(0,7,8)
ene_c60_slab = np.zeros([len(scs),len(d)])

ax = plt.subplot(111)
c60_energy_slab_ts = -524.53523070
slab_ts = -8112.81007882
c60_slab_ref = np.loadtxt("c60_slab_ref.dat")
slab_ref = np.loadtxt("slab_ref.dat")
ene_c60_slab_ts = np.loadtxt("c60_slab_ts.dat")
ene_c60_slab_ts = ene_c60_slab_ts - c60_energy_slab_ts - slab_ts
ax.plot(d[10:],ene_c60_slab_ts[10:],'k.-',label="TS")

c60_energy_slab_none = -521.67745458
slab_none = -7975.33552343
ene_c60_slab_none = np.loadtxt("c60_slab_none.dat")
ene_c60_slab_none = ene_c60_slab_none - c60_energy_slab_none - slab_none
ax.plot(d[10:],ene_c60_slab_none[10:],'r.-',label="None")

for i in range(len(scs)):
    s = "%.1f" % scs[i]
    filename = "c60_slab_" + s + ".dat"
    ene_c60_slab[i,:] = np.loadtxt(filename)
    ene_c60_slab[i,:] = ene_c60_slab[i,:] - c60_slab_ref[i] - slab_ref[i]
    x = float(i)/len(scs)
    ax.plot(d[10:],ene_c60_slab[i,10:],color=(0,x,1.0-x),linestyle='solid',marker='.',label="SCS " + str(scs[i]) + " Å, MBD " + str(scs[i]+1) + " Å")

plt.xlim([d[10],d[-1]])
plt.xlabel("Displacement from initial configuration [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("C60 molecule interacting with amorphous carbon slab")
plt.legend()

arr_img = plt.imread('image.png')
im = OffsetImage(arr_img, zoom=.25)
ab = AnnotationBbox(im, (0.35, 0.65), xycoords='axes fraction', bboxprops=dict(alpha=0.0))
arr_img2 = plt.imread('image2.png')
im2 = OffsetImage(arr_img2, zoom=.25)
ab2 = AnnotationBbox(im2, (0.5, 0.65), xycoords='axes fraction', bboxprops=dict(alpha=0.0))
ax.add_artist(ab)
ax.add_artist(ab2)

plt.show()


fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)


ax = plt.subplot(111)
ax.plot(d[10:],ene_c60_slab_ts[10:],'k.-',label="TS")
ax.plot(d[10:],ene_c60_slab_none[10:],'r.-',label="None")


c60_slab_ref_loc = np.loadtxt("c60_slab_ref_loc.dat")
slab_ref_loc = np.loadtxt("slab_ref_loc.dat")
ene_c60_slab_loc = np.zeros([len(scs),len(d)])

for i in range(len(scs)):
    s = "%.1f" % scs[i]
    filename = "c60_slab_loc_" + s + ".dat"
    ene_c60_slab_loc[i,:] = np.loadtxt(filename)
    ene_c60_slab_loc[i,:] = ene_c60_slab_loc[i,:] - c60_slab_ref_loc[i] - slab_ref_loc[i]
    x = float(i)/len(scs)
    ax.plot(d[10:],ene_c60_slab_loc[i,10:],color=(0,x,1.0-x),linestyle='solid',marker='.',label="MBD cutoff " + str(scs[i]) + " Å")

ax.plot(d[10:],ene_c60_slab[len(scs)-1,10:],color=(1,0.5,0),linestyle='solid',marker='.',label="MBD 8 Å, SCS 7 Å")

plt.xlim([d[10],d[-1]])
plt.xlabel("Displacement from initial configuration [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("C60 molecule interacting with amorphous carbon slab (SCS cutoff 4 Å)")
plt.legend()


arr_img = plt.imread('image.png')
im = OffsetImage(arr_img, zoom=.2)
ab = AnnotationBbox(im, (0.45, 0.7), xycoords='axes fraction', bboxprops=dict(alpha=0.0))
arr_img2 = plt.imread('image2.png')
im2 = OffsetImage(arr_img2, zoom=.2)
ab2 = AnnotationBbox(im2, (0.55, 0.7), xycoords='axes fraction', bboxprops=dict(alpha=0.0))
#ax.add_artist(ab)
#ax.add_artist(ab2)

plt.show()

fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)


ax = plt.subplot(111)
ax.plot(d[10:],ene_c60_slab_ts[10:],'k.-',label="TS")
ax.plot(d[10:],ene_c60_slab_none[10:],'r.-',label="None")


c60_slab_ref_scs = np.loadtxt("c60_slab_ref_scs.dat")
slab_ref_scs = np.loadtxt("slab_ref_scs.dat")
ene_c60_slab_scs = np.zeros([len(scs),len(d)])

for i in range(len(scs)):
    s = "%.1f" % scs[i]
    filename = "c60_slab_scs_" + s + ".dat"
    ene_c60_slab_scs[i,:] = np.loadtxt(filename)
    ene_c60_slab_scs[i,:] = ene_c60_slab_scs[i,:] - c60_slab_ref_scs[i] - slab_ref_scs[i]
    x = float(i)/len(scs)
    ax.plot(d[10:],ene_c60_slab_scs[i,10:],color=(0,x,1.0-x),linestyle='solid',marker='.',label="SCS cutoff " + str(scs[i]) + " Å")

ax.plot(d[10:],ene_c60_slab[len(scs)-1,10:],color=(1,0.5,0),linestyle='solid',marker='.',label="MBD 8 Å, SCS 7 Å")

plt.xlim([d[10],d[-1]])
plt.xlabel("Displacement from initial configuration [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("C60 molecule interacting with amorphous carbon slab (MBD cutoff 5 Å)")
plt.legend()

plt.show()
