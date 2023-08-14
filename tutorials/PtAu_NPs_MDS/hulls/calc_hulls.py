from ase.io import read
import numpy as np

db = read("../db0.xyz", index=":")

all = {}
Ns = []
xs = {}

for atoms in db:
    N = len(atoms)
    x = np.round(atoms.symbols.count("Pt")/N, decimals = 8)
    if (N,x) not in all:
        all[N,x] = [atoms.info["energy"]/N]
    else:
        all[N,x].append(atoms.info["energy"]/N)
    if N not in Ns:
        Ns.append(N)
    if N not in xs:
        xs[N] = [x]
    elif x not in xs[N]:
        xs[N].append(x)

for N in Ns:
    x_vals = []
    e_vals = []
    for x in xs[N]:
        for e in all[N,x]:
            if x not in x_vals:
                x_vals.append(x)
                e_vals.append(e)
            else:
                for i in range(0, len(x_vals)):
                    if x == x_vals[i] and e < e_vals[i]:
                        e_vals[i] = e
    if 0. in x_vals and 1. in x_vals:
        i0 = np.where([x == 0. for x in x_vals])[0][0]
        i1 = np.where([x == 1. for x in x_vals])[0][0]
        e0 = e_vals[i0]
        e1 = e_vals[i1]
        idx = np.argsort(x_vals)
        x_vals = np.array(x_vals)[idx]
        e_vals = np.array(e_vals)[idx]
        for i in range(0, len(x_vals)):
            x = x_vals[i]
            e_vals[i] -= e1*x + e0*(1.-x)
        f = open("%i.dat" % N, "w")
        for i in range(0, len(x_vals)):
            if i == 0 or i == len(x_vals) or e_vals[i] <= 0.:
                print(x_vals[i], e_vals[i], file=f)
        print("", file=f)
        print("", file=f)
        for x in xs[N]:
            for e in all[N,x]:
                print(x, e-e1*x-e0*(1.-x), file=f)
        f.close()
