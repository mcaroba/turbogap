import sys
from quippy.descriptors import Descriptor
from quippy.convert import ase_to_quip
import numpy as np
from ase.io import read,write
import cluster_mds as clmds
import kmedoids
from ase_tools import surface_list

# Read in the database
db = read("db1.xyz", index=":")


# Create mappings
which_structure = []
which_atom = []
map_back = {}
local_energy = []
k = 0
for i in range(0, len(db)):
    le = db[i].get_array("local_energy")
    for j in range(0, len(db[i])):
        which_structure.append(i)
        which_atom.append(j)
        map_back[i,j] = k
        local_energy.append(le[j])
        k += 1


# Rolling sphere algorithm parameters
r_min = 4.
r_max = 4.5
n_tries = 1000000


# Build a list of surface atoms
surf_list_all = []
is_surface = np.full(len(which_atom), False, dtype=bool)
k = 0
for i in range(0, len(db)):
    surf_list = surface_list(db[i], r_min, r_max, n_tries, cluster=True)
    surf_list_all.append(surf_list)
    for j in range(0,len(db[i])):
        if j in surf_list:
            is_surface[k] = True
        k += 1
    if True:
        sys.stdout.write('\rBuilding list of surface atoms:%6.1f%%' % (float(i)*100./len(db)) )
        sys.stdout.flush()


###############################################################################################
# This builds the kernel for the site-wise cl-MDS map
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

d = {}

d["Pt"] = Descriptor(soap["Pt"])
d["Au"] = Descriptor(soap["Au"])

qs = []
for atoms in db:
    at = ase_to_quip(atoms)
    N = {}
    q = {}
    for i in range(0,len(atoms)):
        symb = atoms.symbols[i]
        if symb not in N:
            N[symb] = 0
            q[symb] = d[symb].calc_descriptor(at)
        this_q = q[symb][N[symb]]
        qs.append(this_q)
        N[symb] += 1

dist = np.zeros([len(qs),len(qs)])
n = 0
for i in range(0,len(qs)):
    arg = 1. - np.dot(qs[:], qs[i])**zeta
    arg2 = np.clip(arg, 0, 1)
    dist[i,:] = np.sqrt(arg2)
#    sys.stdout.write('\rBuilding list of surface atoms:%6.1f%%' % (float(i)*100./len(qs)) )
    sys.stdout.write('\rBuilding distance matrix:%6.1f%%' % (float(i)*100./len(qs)) )
    sys.stdout.flush()

print("    MBytes: ", dist.nbytes/1024/1024)
print("")

# Embedding
M, C = kmedoids.kMedoids(dist, n_cl, n_iso=n_cl//2, init_Ms="isolated", n_inits=10)
data = clmds.clMDS(dist_matrix = dist, sparsify="cur", n_sparse=500)
data.cluster_MDS([n_cl,1], weight_cluster_mds=2., weight_anchor_mds=10., init_medoids=M)
list_sparse = data.sparse_list
XY_sparse = data.sparse_coordinates
C_sparse = data.sparse_cluster_indices
M_sparse = data.sparse_medoids
M = []
for k in M_sparse:
    M.append(list_sparse[k])

M = np.array( M, dtype=int )
indices = np.array( list(range(0,len(dist))), dtype=int )
data.compute_pts_estim_coordinates( indices=indices )
XY = data.estim_coordinates
C = data.estim_cluster_indices
# Transform to usual cluster array
C_usual = []
for i in range(0, n_cl):
    C_usual.append([])

for i in range(0, len(C)):
    C_usual[C[i]].append(i)

C_list = C
C = C_usual
#

k = 0
for atoms in db:
    atoms.set_array("clmds_cluster", C_list[k:k+len(atoms)])
    atoms.set_array("clmds_coordinates", XY[k:k+len(atoms)])
    k += len(atoms)

k = 0
for atoms in db:
    is_medoid = np.full(len(atoms), False)
    is_surface2 = np.full(len(atoms), False)
    for i in range(0, len(atoms)):
        if is_surface[k]:
            is_surface2[i] = True
        if k in M:
            is_medoid[i] = True
        k += 1
    atoms.set_array("clmds_medoid", is_medoid)
    atoms.set_array("surface", is_surface2)


f = open("mds_sites.dat", "w")
print("# X, Y, cluster, is_medoid, is_surface, le, species, N of NP", file=f)
k = 0
for i in range(0, len(db)):
    for j in range(0, len(db[i])):
        print(*XY[k], C_list[k], k in M, is_surface[k], local_energy[k], db[i].symbols[j], len(db[i]), file=f)
        k += 1

f.close()



write("db2.xyz", db)
###############################################################################################
