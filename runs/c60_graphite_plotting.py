import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

#N = 444
plt.rcParams['text.usetex'] = True

# Distance of C60 center of mass from top layer
dist = 3.3404496666666676

c60_none = -521.67745458
graphene_none = -3479.16684325
c60_ts = -524.53100703
graphene_ts = -3511.88592008
c60_mbd = -525.69879996
graphene_mbd = -3515.82023359

ene_ts = np.loadtxt("c60_graphite_ts_full.dat")
ene_none = np.loadtxt("c60_graphite_none_full.dat")
ene_mbd = np.loadtxt("c60_graphite_mbd_full.dat")

ene_ts = ene_ts-c60_ts-graphene_ts
ene_none = ene_none-c60_none-graphene_none
ene_mbd = ene_mbd-c60_mbd-graphene_mbd

ene_ts = ene_ts
ene_none = ene_none
ene_mbd = ene_mbd


d1 = np.linspace(2.4,4,17)
d2 = np.linspace(4.4,8,10)
d = np.append(d1,d2)
d = d+dist
d1 = np.linspace(2.4,3.,7)
d1 = np.append(d1,3.05)
d2 = np.linspace(3.1,4,10)
d1 = np.append(d1,d2)
d2 = np.linspace(4.4,8,10)
d_mbd = np.append(d1,d2)
d_mbd = d_mbd+dist

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.plot(d,ene_none,'r.-',label="None")
ax.plot(d_mbd,ene_ts,'b.-',label="TS")
ax.plot(d_mbd,ene_mbd,'g.-',label="MBD")
ax.set_xlim([d_mbd[0],d_mbd[-5]])
ax.set_xlabel("Distance of the C60 CM from the top layer [Å]")
ax.set_ylabel("Interaction energy between C60 and BLG [eV]")
ax.set_title("C60 interacting with bilayer AB graphene")
arr_img = plt.imread('c60_graphite.png')
im = OffsetImage(arr_img, zoom=.13)
ab = AnnotationBbox(im, (0.5, 0.7), xycoords='axes fraction', bboxprops=dict(alpha=0.0))
ax.add_artist(ab)
ax.legend()
plt.show()

