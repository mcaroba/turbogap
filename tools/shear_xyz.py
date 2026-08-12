#!/usr/bin/env python3
"""Apply an affine shear to an extxyz frame: cell and positions together, so the
fractional coordinates -- and therefore the local geometry up to the strain --
are preserved. Produces a genuinely triclinic starting point for the
variable-cell relaxation to undo."""
import re, sys

src, dst = sys.argv[1], sys.argv[2]
# 4% shear in xy and 3% in xz, applied as x -> S x with S upper triangular
S = [[1.0, 0.04, 0.03],
     [0.0, 1.0, 0.02],
     [0.0, 0.0, 1.0]]

def mul(S, v):
    return [sum(S[i][j]*v[j] for j in range(3)) for i in range(3)]

with open(src) as f:
    lines = f.read().splitlines()

nat = int(lines[0].split()[0])
comment = lines[1]
m = re.search(r'Lattice="([^"]*)"', comment)
lat = [float(x) for x in m.group(1).split()]
# extxyz stores the lattice row-major as a1 a2 a3 b1 b2 b3 c1 c2 c3
vecs = [lat[0:3], lat[3:6], lat[6:9]]
new = [mul(S, v) for v in vecs]
comment = comment[:m.start(1)] + " ".join("%.10f" % x for x in sum(new, [])) + comment[m.end(1):]
# Drop the per-frame properties that no longer describe this structure
for key in ("energy", "pressure", "temperature"):
    comment = re.sub(r'\s%s=[-0-9.eE+]+' % key, '', comment)
comment = re.sub(r'\s(virial|stress)="[^"]*"', '', comment)

out = [str(nat), comment]
for line in lines[2:2+nat]:
    f = line.split()
    pos = mul(S, [float(f[1]), float(f[2]), float(f[3])])
    out.append("%-4s %18.10f %18.10f %18.10f %s" % (f[0], pos[0], pos[1], pos[2], " ".join(f[4:])))

open(dst, "w").write("\n".join(out) + "\n")
print("wrote", dst)
