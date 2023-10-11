import numpy as np
import matplotlib.pyplot as plt

graphene = -36.90674363
graphene_ts = -37.08779649
graphene_mbd = -37.23038936

ene = np.loadtxt("bilayer.dat")
ene_ts = np.loadtxt("bilayer_ts.dat")
ene_mbd = np.loadtxt("bilayer_mbd.dat")

ene = ene - 2*graphene
ene_ts = ene_ts - 2*graphene_ts
ene_mbd = ene_mbd - 2*graphene_mbd

d = np.linspace(1.6,5.0,35)

plt.plot(d[9:],ene[9:],'r.-')
plt.plot(d[9:],ene_ts[9:],'b.-')
plt.plot(d[9:],ene_mbd[9:],'c.-')
plt.show()

