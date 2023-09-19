import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

plt.rcParams['text.usetex'] = True

ene_ts = np.loadtxt("c60_dimer_ts.dat")
ene_none = np.loadtxt("c60_dimer_none.dat")
ene_mbd = np.loadtxt("c60_dimer_mbd.dat")
ene_dft = np.loadtxt("c60_dimer_dft.dat")
ene_dft_edisp = np.loadtxt("c60_dimer_dft_edisp.dat")
c60_ts = -534.21383883
c60_none = -531.36028639
c60_mbd = -535.76915531
c60_dft = -535.74527446
c60_dft_edisp = -4.43104

ene_ts_edisp = np.loadtxt("c60_dimer_ts_edisp.dat")
ene_mbd_edisp = np.loadtxt("c60_dimer_mbd_edisp.dat")

c60_ts_edisp = -2.85355245
c60_mbd_edisp = -4.40887718

ene_ts = ene_ts - 2*c60_ts
ene_none = ene_none -2*c60_none
ene_mbd = ene_mbd -2*c60_mbd
ene_dft = ene_dft-2*c60_dft #+ene_dft_edisp-2*c60_dft_edisp
#ene_dft_edisp = ene_dft_edisp -2*c60_dft_edisp
ene_ts_edisp = ene_ts_edisp - 2*c60_ts_edisp
ene_mbd_edisp = ene_mbd_edisp -2*c60_mbd_edisp
ene_dft_edisp = ene_dft_edisp -2*c60_dft_edisp

d1 = np.linspace(8.,9.6,9)
d1 = np.append(d1,[9.7,9.8,9.9])
d2 = np.linspace(10,12,11)
d = np.append(d1,d2)
d_dft = np.array([9.6,9.7,9.8,9.9,10.0])

plt.figure(figsize=(12,6))
plt.plot(d[6:],ene_none[6:],'r.-',label="GAP+None")
plt.plot(d[6:],ene_ts[6:],'b.-',label="GAP+TS")
plt.plot(d[6:],ene_mbd[6:],'g.-',label="GAP+MBD")
plt.plot(d_dft,ene_dft,'k.-',label="DFT+MBD")
plt.xlim([d[6],d[-1]])
plt.xlabel("Distance between C60 centers [Å]")
plt.ylabel("Interaction energy [eV]")
plt.title("Interaction energy for C60 dimer")
plt.legend()
#plt.plot(d_dft,ene_dft_edisp,'c*-')
plt.show()

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.plot(d[6:],ene_none[6:],'r.-',label="GAP+None")
ax.plot(d[6:],ene_ts[6:],'b.-',label="GAP+TS")
ax.plot(d[6:],ene_mbd[6:],'g.-',label="GAP+MBD")
ax.plot(d_dft,ene_dft,'k.-',label="DFT+MBD")
ax.set_xlim([d[6],d[-1]])
ax.set_xlabel("Distance between C60 centers [Å]")
ax.set_ylabel("Interaction energy [eV]")
ax.set_title("Interaction energy for C60 dimer")
arr_img = plt.imread('c60_dimer.png')
im = OffsetImage(arr_img, zoom=.13)
ab = AnnotationBbox(im, (0.5, 0.7), xycoords='axes fraction', bboxprops=dict(alpha=0.0))
ax.add_artist(ab)
ax.legend()
plt.show()




plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.plot(d[6:],ene_ts_edisp[6:],'b.-',label="TS from GAP")
ax.plot(d[6:],ene_mbd_edisp[6:],'g.-',label="MBD from GAP")
ax.plot(d_dft,ene_dft_edisp,'k.-',label="MBD from DFT")
ax.set_xlim([d[6],d[-1]])
ax.set_xlabel("Distance between C60 centers [Å]")
ax.set_ylabel("Dispersion interaction energy [eV]")
ax.set_title("Dispersion interaction energy for C60 dimer")
arr_img = plt.imread('c60_dimer.png')
im = OffsetImage(arr_img, zoom=.13)
ab = AnnotationBbox(im, (0.5, 0.3), xycoords='axes fraction', bboxprops=dict(alpha=0.0))
ax.add_artist(ab)
ax.legend()
plt.show()
