import numpy as np
from ase.units import Bohr

B_mat_p = np.loadtxt("B_mat_p.dat")
B_mat_m = np.loadtxt("B_mat_m.dat")

dB = (B_mat_p-B_mat_m)/0.002*Bohr

b_der = np.loadtxt("b_der.dat")
indices = np.loadtxt("indices.dat")
alpha = np.loadtxt("alpha_SCS_full.dat")

indices=indices.astype(int)

b_der_fd = np.array([180,3])

#for p in range(60):
#    i = indices(p)-1

alpha_rearr = np.zeros([180,3])

print(alpha)
print(indices-1)
#print(alpha.shape)
for p in range(60):
    #print(p)
    alpha_rearr[3*p:3*p+3,:] = alpha[3*(indices[p]-1):3*(indices[p]-1)+3,:]

print(alpha_rearr)

print(np.amax(np.abs(np.matmul(dB,alpha_rearr)-b_der)))
