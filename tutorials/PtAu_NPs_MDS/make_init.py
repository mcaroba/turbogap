import numpy as np
from ase.io import write
from ase import Atoms

a = 3.9668
d_min = 2.

append = False
num = 1

nrand = 10
for n in range(25,50+1,5):
    for f_Au in np.arange(0.,1.+1.e-10,0.1):
        for i in range(0, nrand):
            R = (a**3 / 4. * n * 3./4./np.pi)**(1./3.)

            pos = []
            while len(pos) < n:
                this_pos = (2.*np.random.sample([3]) -1) * R
                d = np.dot(this_pos, this_pos)**0.5
                if d > R:
                    continue
                if len(pos) == 0:
                    pos.append(this_pos)
                else:
                    too_close = False
                    for other_pos in pos:
                        dv = this_pos - other_pos
                        d = np.dot(dv, dv)**0.5
                        if d < d_min:
                            too_close = True
                            break
                    if not too_close:
                        pos.append(this_pos)

            n_Au = int(n*f_Au)
            n_Pt = n - n_Au

            atoms = Atoms("Pt%iAu%i" % (n_Pt, n_Au), positions = pos, pbc=True)
            write("init_xyz/all.xyz", atoms, append = append)
            atoms.center(vacuum = 10.)
            write("init_xyz/%i.xyz" % num, atoms)
            num += 1
            append = True
