import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

graphene = -36.91333179
graphene_ts = -37.09438465
graphene_mbd = -37.23697460
graphene_mbd2 = -37.24653502
graphene_dft = -37.24994036

ene = np.loadtxt("bilayer.dat")
ene_ts = np.loadtxt("bilayer_ts.dat")
ene_mbd = np.loadtxt("bilayer_mbd.dat")
ene_mbd2 = np.loadtxt("bilayer_mbd2.dat")
ene_dft = np.loadtxt("bilayer_dft.dat")

ene = ene - 2*graphene
ene_ts = ene_ts - 2*graphene_ts
ene_mbd = ene_mbd - 2*graphene_mbd
ene_mbd2 = ene_mbd2 - 2*graphene_mbd2
ene_dft = ene_dft - 2*graphene_dft

d = np.linspace(1.6,5.0,35)

#plt.figure(figsize=(12,6))
plt.plot(d[9:],ene[9:],'r.-',label="None")
plt.plot(d[9:],ene_ts[9:],'b.-',label="TS")
plt.plot(d[9:],ene_mbd[9:],'c.-',label="MBD (10+5)")
plt.plot(d[9:],ene_mbd2[9:],'m.-',label="MBD (10+10)")
plt.plot(d[9:],ene_dft[9:],'k.-',label="DFT+MBD")
plt.xlim([d[9],d[-1]])
plt.legend()
plt.xlabel("Interlayer distance [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("Interaction energy of bilayer AB graphene")
plt.show()

