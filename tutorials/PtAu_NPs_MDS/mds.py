from ase import cluster
from quippy.descriptors import Descriptor
from quippy.convert import ase_to_quip
import numpy as np
from ase.io import read,write
import cluster_mds as clmds
import kmedoids


# Read in the database
#db = []
#for i in range(1,660+1):
#    atoms = read("trajs/trajs_1.1/%i.xyz" % i, index=-1)
#    db.append(atoms)

db = read("db0.xyz", index=":")

###############################################################################################
# This builds the global kernel for the NP-wise cl-MDS map
n_cl = 22 # number of clusters
cutoff = 5.7; dc = 0.5; sigma = 0.5
zeta = 6
soap = {"Au": 'soap_turbo alpha_max={8 8} l_max=8 rcut_soft=%.4f rcut_hard=%.4f atom_sigma_r={%.4f %.4f} atom_sigma_t={%.4f %.4f} \
               atom_sigma_r_scaling={0. 0.} atom_sigma_t_scaling={0. 0.} radial_enhancement=1 amplitude_scaling={1. 1.} \
               basis="poly3gauss" scaling_mode="polynomial" species_Z={78 79} n_species=2 central_index=2 central_weight={1. 1.} \
               compress_mode=trivial' % (cutoff-dc, cutoff, *(4*[sigma])),
        "Pt": 'soap_turbo alpha_max={8 8} l_max=8 rcut_soft=%.4f rcut_hard=%.4f atom_sigma_r={%.4f %.4f} atom_sigma_t={%.4f %.4f} \
               atom_sigma_r_scaling={0. 0.} atom_sigma_t_scaling={0. 0.} radial_enhancement=1 amplitude_scaling={1. 1.} \
               basis="poly3gauss" scaling_mode="polynomial" species_Z={78 79} n_species=2 central_index=1 central_weight={1. 1.} \
               compress_mode=trivial' % (cutoff-dc, cutoff, *(4*[sigma]))}

d1 = Descriptor(soap["Pt"])
d2 = Descriptor(soap["Au"])

qs = []
for atoms in db:
    at = ase_to_quip(atoms)
    q1 = d1.calc_descriptor(at)
    q2 = d2.calc_descriptor(at)
    if q1.shape == (1,0):
        q = q2
    elif q2.shape == (1,0):
        q = q1
    else:
        q = np.concatenate((q1, q2))
    qs.append(q)

kern = np.zeros([len(qs), len(qs)])
for i in range(0,len(qs)):
    for j in range(i,len(qs)):
        k = 0.
        for i2 in range(0, len(qs[i])):
            for j2 in range(0, len(qs[j])):
                k += np.dot(qs[i][i2], qs[j][j2])**zeta
        kern[i,j] = k
        kern[j,i] = k

kern_n = np.zeros([len(qs), len(qs)])
for i in range(0,len(qs)):
    for j in range(0,len(qs)):
        kern_n[i,j] = kern[i,j] / np.sqrt(kern[i,i] * kern[j,j])

dist = np.zeros([len(qs),len(qs)])
for i in range(0,len(qs)):
    dist[i, i] = 0.
    for j in range(i+1,len(qs)):
        dist[i, j] = np.sqrt(1.-kern_n[i,j])
        dist[j, i] = np.sqrt(1.-kern_n[i,j])

M, C = kmedoids.kMedoids(dist, n_cl, n_iso=n_cl//2, init_Ms="isolated", n_inits=10)
data = clmds.clMDS(dist_matrix = dist)
data.cluster_MDS([n_cl,1], weight_cluster_mds=2., weight_anchor_mds=10., init_medoids=M)
list_sparse = data.sparse_list
XY_sparse = data.sparse_coordinates
C_sparse = data.sparse_cluster_indices
M_sparse = data.sparse_medoids

f = open("mds.dat", "w")
print("# X, Y, cluster, is_medoid, N, E, x in Pt_{x}Au_{1-x}", file=f)
for i in range(0, len(db)):
    db[i].info["clmds_coordinates"] = XY_sparse[i]
    db[i].info["clmds_cluster"] = C_sparse[i]
    db[i].info["clmds_medoid"] = i in M_sparse
    print(*XY_sparse[i], C_sparse[i], i in M_sparse, len(db[i]), db[i].info["energy"], db[i].symbols.count("Pt")/len(db[i]), file=f)

f.close()
write("db1.xyz", db)
###############################################################################################
