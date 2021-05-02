# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   TurboGAP
# HND X
# HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
# HND X
# HND X   TurboGAP is published and distributed under the
# HND X      Academic Software License v1.0 (ASL)
# HND X
# HND X   This file, compress_indices.py, is copyright (c) 2019-2021, Miguel A. Caro
# HND X
# HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
# HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
# HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# HND X   ASL for more details.
# HND X
# HND X   You should have received a copy of the ASL along with this program
# HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
# HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
# HND X   http://github.com/gabor1/ASL
# HND X
# HND X   When using this software, please cite the following reference:
# HND X
# HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
# HND X
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import numpy as np

# Change these to actual values for your basis set (nmax must always be a list, even when it only contains one element)
nmax = [6, 6]
lmax = 5

# Choose the mode for compression
mode = "trivial"



if mode == "trivial":
    pivot = np.zeros(len(nmax), dtype=int)
    pivot[0] = 1
    for i in range(len(nmax)-1):
        pivot[i+1] = pivot[i] + nmax[i]
    counter = 0
    k = 1
    indices = []
    for n in range(1, np.sum(nmax)+1):
        for m in range(n, np.sum(nmax)+1):
            for l in range(0, lmax+1):
                if np.any(n == pivot) or np.any(m == pivot):
                    indices.append(k)
                    counter += 1
                k += 1
else:
    print("ERROR: Mode \"%s\" not implemented!!!!!" % mode)
    exit()

print(*nmax, lmax, counter)
for i in indices:
    print(i)
