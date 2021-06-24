# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   TurboGAP
# HND X
# HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
# HND X
# HND X   TurboGAP is published and distributed under the
# HND X      Academic Software License v1.0 (ASL)
# HND X
# HND X   This file, make_gap_files.py, is copyright (c) 2019-2021, Mikhail S. Kuklin
# HND X   and Miguel A. Caro
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

#  !/usr/bin/env python3

import numpy as np
import re
import sys
from pathlib import Path

print("\n		MAKE GAP FILES | v. 2.1, 03.05.2020 | Author(s): Mikhail S. Kuklin, Miguel A. Caro\n")


# transform list to string function
def listToString(s):
    str1 = " "
    return (str1.join(s))


# dictionary of the elements: 1-57, 72-89
elements = {1: ' H', 2: ' He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
            11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
            21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
            31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
            41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
            51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
            72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
            81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac'}


# identify chemical symbol based on the atomic number
def elem(item):
    if item in elements:
        return elements[item]


# search for the keywords from xml file
def rex(keyw):
    k = keyw
    v = re.findall(k, f)
    return v


# search for the keywords from hirshfeld.xml file
def rexh(keyw):
    k = keyw
    v = re.findall(k, f_h)
    return v


# similar to rex but with focus on gp_coordinates part to avoid duplication
def rexcor(keyw):
    k = keyw
    v = re.findall(k, s_gp_coord)
    return v


# similar to rex but with focus on descriptor part to avoid duplication
def rexdesc(keyw):
    k = keyw
    v = re.findall(k, s_desc)
    return v


# open potential xml file as the first argument in the command line
gap_file = sys.argv[1]
file = open(gap_file, 'rt');
f = file.read()
print("Reading " + sys.argv[1] + "...")

if len(sys.argv) > 3:
    hirs = sys.argv[3]
    file_h = open(hirs, 'rt');
    f_h = file_h.read()
    print("Reading " + sys.argv[3] + "...")

# identify particular parts from the file
desc = rex('<descriptor(.*)descriptor>')
s_desc = str(desc)

# extract covariance type
cov_type = {1: 'exp', 4: 'pol'}


# tranfirm covariance type from quip to turbogap format
def ct(item):
    if item in cov_type:
        return cov_type[item]


# count number of different terms
term1 = "distance_2b"
term1_alt = "distance_Nb"
term2 = "angle_3b"
term3 = "soap_turbo"
term4 = "pairpot"
a = 0
a2 = 0
b = 0
c = 0
d = 0

find_t1 = rex('\s*(\S+) distance_2b')
find_t1_alt = rex('\s*(\S+) distance_Nb')
find_t2 = rex('\s*(\S+) angle_3b')
find_t3 = rex('\s*(\S+) soap_turbo')
find_t4 = rex('<potential_pair')
if len(sys.argv) > 3:
    find_t5 = rexh('\s*(\S+) soap_turbo')

file = open(gap_file, 'rt')
for line in file:
    words = line.split()
    for i in words:
        if (i == term1) and (find_t1[a] == '<descriptor>'):
            a = a + 1
        elif (i == term1_alt) and (find_t3[c] == '<descriptor>'):
            a2 = a2 + 1
        elif (i == term2) and (find_t2[b] == '<descriptor>'):
            b = b + 1
        elif (i == term3) and (find_t3[c] == '<descriptor>'):
            c = c + 1
file.close()

if len(sys.argv) > 3:
    file_h = open(gap_file, 'rt')
    for line in file_h:
        words = line.split()
        for i in words:
            if (i == term3) and (find_t5[d] == '<descriptor>'):
                d = d + 1
    file_h.close()

n_2b = a
if n_2b == 0:
    n_2b = a2
print("Number of distance_2b => " + str(n_2b))
n_3b = b
print("Number of angle_3b => " + str(n_3b))
n_mb = c
print("Number of soap_turbo => " + str(n_mb))
n_cp = len(find_t4)
if n_cp > 0:
    print("Core potentials are found")
n_hirs = d
if n_hirs > 0:
    print("Hirshfeld volumes GAP(s) are found => " + str(n_hirs))

# create file with alphas
alpha_s = 'alpha="(.*)" sparse'
as_tmp = re.findall(alpha_s, f)
if as_tmp == []:
    print("Alphas are not found. Check your potential file!")
alphas = np.array(as_tmp)

# create file with sparce cutoff
sparsec_s = 'sparseCutoff="(.*)"'
ss_tmp = re.findall(sparsec_s, f)
if ss_tmp == []:
    print("sparseCutoffs are not found. Check your potential file!")
cutoffs = np.array(ss_tmp)

# extract n_sparse
nspar_tmp = rex('(.*\n{1}).*</gpCoordinates>')

nspar_2b = []
for i in range(n_2b):
    if n_2b == 0:
        exit
    else:
        nspar_2b_tmp1 = re.findall('<sparseX i="(-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)', str(nspar_tmp[i]))
        nspar_2b_tmp2 = list(map(int, nspar_2b_tmp1))
        nspar_2b.append(nspar_2b_tmp2)

nspar_3b = []
for i in range(n_3b):
    if n_3b == 0:
        exit
    else:
        i = n_2b + i
        nspar_3b_tmp1 = re.findall('<sparseX i="(-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)', str(nspar_tmp[i]))
        nspar_3b_tmp2 = list(map(int, nspar_3b_tmp1))
        nspar_3b.append(nspar_3b_tmp2)

nspar_mb = []
for i in range(n_mb):
    if n_mb == 0:
        exit
    else:
        i = n_2b + n_3b + i
        nspar_mb_tmp1 = re.findall('<sparseX i="(-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)', str(nspar_tmp[i]))
        nspar_mb_tmp2 = list(map(int, nspar_mb_tmp1))
        nspar_mb.append(nspar_mb_tmp2)

nsparse1 = nspar_2b + nspar_3b + nspar_mb
nsparse = [val for sublist in nsparse1 for val in sublist]
if nsparse == []:
    print("Number of sparce configurations are not found. Check your potential file!")

k = 0
j = 0
for i in range(0, n_2b):
    n = nsparse[j]
    this_alphas = alphas[k:k + n]
    this_cutoffs = cutoffs[k:k + n]
    t = open("alphas_2b_%i.dat" % (i + 1), "w+")
    for i2 in range(len(this_alphas)):
        print(this_alphas[i2], this_cutoffs[i2], file=t)
    t.close()
    k += n
    j += 1
print("Creating alphas_2b files...")

for i in range(0, n_3b):
    n = nsparse[j]
    this_alphas = alphas[k:k + n]
    this_cutoffs = cutoffs[k:k + n]
    t = open("alphas_3b_%i.dat" % (i + 1), "w+")
    for i2 in range(len(this_alphas)):
        print(this_alphas[i2], this_cutoffs[i2], file=t)
    t.close()
    k += n
    j += 1
print("Creating alphas_3b files...")

for i in range(0, n_mb):
    n = nsparse[j]
    this_alphas = alphas[k:k + n]
    t = open("alphas_mb_%i.dat" % (i + 1), "w+")
    for i2 in range(len(this_alphas)):
        print(this_alphas[i2], file=t)
    t.close()
    k += n
    j += 1
print("Creating alphas_mb files...")

# make alphas for Hirshfeld volumes GAPs
if n_hirs > 0:
    nspar_tmp = rexh('(.*\n{1}).*</gpCoordinates>')
    alpha_s = 'alpha="(.*)" sparse'
    as_tmp = re.findall(alpha_s, f_h)
    if as_tmp == []:
        print("Alphas for Hirshfeld volumes GAP(s) are not found. Check your potential file!")
    alphas = np.array(as_tmp)

    # create file with sparce cutoff
    sparsec_s = 'sparseCutoff="(.*)"'
    ss_tmp = re.findall(sparsec_s, f_h)
    if ss_tmp == []:
        print("sparseCutoffs are not found. Check your potential file!")
    cutoffs = np.array(ss_tmp)

    nspar_hirs = []
    for i in range(n_hirs):
        if n_hirs == 0:
            exit
        else:
            nspar_hirs_tmp1 = re.findall('<sparseX i="(-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)', str(nspar_tmp[i]))
            nspar_hirs_tmp2 = list(map(int, nspar_hirs_tmp1))
            nspar_hirs.append(nspar_hirs_tmp2)

    nsparse1 = nspar_hirs
    nsparse = [val for sublist in nsparse1 for val in sublist]
    if nsparse == []:
        print("Number of sparce configurations are not found. Check your potential file!")

    k = 0
    j = 0
    for i in range(0, n_hirs):
        n = nsparse[j]
        this_alphas = alphas[k:k + n]
        t = open("alphas_hirshfeld_%i.dat" % (i + 1), "w+")
        for i2 in range(len(this_alphas)):
            print(this_alphas[i2], file=t)
        t.close()
        k += n
        j += 1
    print("Creating alphas_hirshfeld files...")

# define arguments to be extracted from xml
gp_coord = rex('<gpCoordinates(.*)">')
s_gp_coord = str(gp_coord)
el1_tmp = rex('Z1=([-+]?\d*\.\d+|\d+)')
el2_tmp = rex('Z2=([-+]?\d*\.\d+|\d+)')
delta_tmp = rex('delta=([-+]?\d*\.\d+|\d+)')
sigma_tmp = rex('theta_uniform=([-+]?\d*\.\d+|\d+)')
ex_kt = rex('<gpCoordinates(.*)">')
kt_tmp = rexcor('covariance_type=\s*(\S+)')
pn = [item.replace('"', '') for item in kt_tmp]
cut_tmp = rex('cutoff=([-+]?\d*\.\d+|\d+)')
potname_tmp = rex('gpCoordinates label="GAP_(.*)" dimensions')
ns_tmp = rex('n_species=(\d)')
compress = rexdesc('compress_file')
# check if hirshfeld.xml is given
if len(sys.argv) > 3:
    potname_h_tmp = rexh('gpCoordinates label="GAP_(.*)" dimensions')
    delta_h_tmp = rexh('delta=([-+]?\d*\.\d+|\d+)')
    lp1 = rexh('local_property0=([-+]?\d*\.\d+|\d+)')

# give name of the output file as the second argument in the command when script is executed
output = sys.argv[2]

p = Path(output)
if p.is_file():
    print(
        "WARNING: you have already output file in the folder with the same name. Remove the file and execute script again")
    exit()

print("Creating " + sys.argv[2] + "...")

for i in range(n_2b):
    if n_2b == 0:
        exit
    else:
        with open(output, 'a') as f3:
            if el1_tmp == []:
                keyw_spec = rexdesc('species_Z=([-+]?\d*\.\d+|\d+)')
                spec_tmp1 = keyw_spec[i]
                spec_tmp2 = int(spec_tmp1)
                tmp = elem(spec_tmp2)
                el1 = str(tmp)
            elif el1_tmp == []:
                keyw_spec = rexdesc('species_Z=([-+]?\d*\.\d+|\d+)')
                spec_tmp1 = keyw_spec[i]
                spec_tmp2 = int(spec_tmp1)
                tmp = elem(spec_tmp2)
                el1 = str(tmp)
            else:
                el1 = elem(int(el1_tmp[i]))
            if el2_tmp == []:
                el2 = el1
            else:
                el2 = elem(int(el2_tmp[i]))
            delta = float(delta_tmp[i])
            sigma = float(sigma_tmp[i])
            kt = ct(int(pn[i]))
            cut = float(cut_tmp[i])
            potname = potname_tmp[i]
            f3.write('gap_beg distance_2b\n')
            f3.write('species1 = %s\n' % el1)
            f3.write('species2 = %s\n' % el2)
            f3.write('delta = %s\n' % delta)
            f3.write('sigma = %s\n' % sigma)
            f3.write('rcut = %s\n' % cut)
            f3.write('desc_sparse = \'gap_files/%s.sparseX.GAP_%s\'\n' % (sys.argv[1], potname))
            u = i + 1
            f3.write('alphas_sparse = \'gap_files/alphas_2b_%u.dat\'\n' % u)
            f3.write('gap_end\n')
            f3.write('\n')
            f3.close()

for i in range(n_3b):
    if n_3b == 0:
        exit
    else:
        z = i + n_2b
        with open(output, 'a') as f3:
            specen_tmp = rex('Z_center=([-+]?\d*\.\d+|\d+)')
            elc = elem(int(specen_tmp[i]))
            el1 = elem(int(el1_tmp[z]))
            el2 = elem(int(el2_tmp[z]))
            delta = float(delta_tmp[z])
            ns = int(ns_tmp[0])
            sigma = float(sigma_tmp[z])
            kt = ct(int(pn[z]))
            if kt == []:
                kt = "exp"
            cut = float(cut_tmp[z])
            potname = potname_tmp[z]
            f3.write('gap_beg angle_3b\n')
            f3.write('species_center = %s\n' % elc)
            f3.write('species1 = %s\n' % el1)
            f3.write('species2 = %s\n' % el2)
            f3.write('delta = %s\n' % delta)
            f3.write('sigma = %s %s %s\n' % (sigma, sigma, sigma))
            f3.write('kernel_type = "%s"\n' % kt)
            f3.write('rcut = %s\n' % cut)
            f3.write('desc_sparse = \'gap_files/%s.sparseX.GAP_%s\'\n' % (sys.argv[1], potname))
            u = i + 1
            f3.write('alphas_sparse = \'gap_files/alphas_3b_%u.dat\'\n' % u)
            f3.write('gap_end\n')
            f3.write('\n')
            f3.close()

for i in range(n_mb):
    if n_mb == 0:
        exit
    else:
        y = i + n_2b + n_3b
        with open(output, 'a') as f3:
            ns = int(ns_tmp[0])
            if ns == []:
                ns = 1
            elif ns == 1:
                keyw_spec = rexdesc('species_Z=([-+]?\d*\.\d+|\d+)')
                spec_tmp1 = keyw_spec[i]
                spec_tmp2 = int(spec_tmp1)
                tmp = elem(spec_tmp2)
                spec = str(tmp)
            else:
                keyw_spec = rexdesc('species_Z={(.*?)}')
                spec_tmp1 = keyw_spec[i].split(' ')
                # spec_tmp2 = spec_tmp1.replace(" ", "")
                # print(spec_tmp2)
                # spec_tmp3 = int(spec_tmp2)
                # print(spec_tmp3)
                spec_tmp4 = [int(x) for x in spec_tmp1]
                spec_tmp5 = []
                for n in spec_tmp4:
                    tmp = elem(n)
                    spec_tmp5.append(tmp)
                spec = listToString(spec_tmp5)
            cents_tmp = rexdesc('central_index=(\d)')
            if cents_tmp == []:
                cents = 1
            else:
                cents = int(cents_tmp[i])

            rcuth_tmp = rexdesc('rcut_hard=([-+]?\d*\.\d+|\d+)')
            rcuth_b = float(rcuth_tmp[i])
            rcuth = (str(rcuth_tmp[i]) + " ") * ns
            rcuts_tmp = rexdesc('rcut_soft=([-+]?\d*\.\d+|\d+)')
            rcuts = float(rcuts_tmp[i])
            buf_tmp = rcuth_b - rcuts
            buf = (str(buf_tmp) + ' ') * ns
            atomsr_tmp = rex('atom_sigma_r={{([-+]?\d*\.\d+|\d+)')
            atomsr = (str(atomsr_tmp[i]) + ' ') * ns
            atomst_tmp = rex('atom_sigma_t={{([-+]?\d*\.\d+|\d+)')
            atomst = (str(atomst_tmp[i]) + ' ') * ns
            atomsrs_tmp = rex('atom_sigma_r_scaling={{([-+]?\d*\.\d+|\d+)')
            atomsrs = (str(atomsrs_tmp[i]) + ' ') * ns
            atomsts_tmp = rex('atom_sigma_t_scaling={{([-+]?\d*\.\d+|\d+)')
            atomsts = (str(atomsts_tmp[i]) + ' ') * ns
            amplsc_tmp = rex('amplitude_scaling={{([-+]?\d*\.\d+|\d+)')
            amplsc = (str(amplsc_tmp[i]) + ' ') * ns
            nmax_tmp = rex('alpha_max={{(\d)')
            nmax = (str(nmax_tmp[i]) + " ") * ns
            lmax_tmp = rexdesc('l_max=(\d)')
            lmax = int(lmax_tmp[i])
            centw_tmp = rexdesc('central_weight={([-+]?\d*\.\d+|\d+)')
            if centw_tmp == []:
                centw_tmp = rexdesc('central_weight=([-+]?\d*\.\d+|\d+)')
                centw = (str(centw_tmp[i]) + ' ') * ns
            else:
                centw = (str(centw_tmp[i]) + ' ') * ns
            scalm_tmp = rexdesc(r'scaling_mode=\s*(\S+)')
            if scalm_tmp == []:
                scalm = "polynomial"
            else:
                scalm = scalm_tmp[i]
            bas_tmp = rexdesc(r'basis=\s*(\S+)')
            if bas_tmp == []:
                bas = "poly"
            else:
                bas = bas_tmp[i]
            radenh_tmp = rexdesc('radial_enhancement=(\d)')
            if radenh_tmp == []:
                radenh_tmp = rexdesc('radial_enhancement={(\d)}')
                radenh = int(radenh_tmp[i])
            elif radenh_tmp == []:
                radenh = 0
            else:
                radenh = int(radenh_tmp[i])
            zeta_tmp = rexdesc('zeta=(\d)')
            zeta = int(zeta_tmp[i])
            delta = float(delta_tmp[y])
            potname = potname_tmp[y]
            nf = (str(4) + ' ') * ns
            f3.write('gap_beg soap_turbo\n')
            f3.write('n_species = %s\n' % ns)
            f3.write('species = %s\n' % spec)
            f3.write('central_species = %s\n' % cents)
            f3.write('rcut = %s\n' % rcuth)
            f3.write('buffer = %s\n' % buf)
            f3.write('atom_sigma_r = %s\n' % atomsr)
            f3.write('atom_sigma_t = %s\n' % atomst)
            f3.write('atom_sigma_r_scaling = %s\n' % atomsrs)
            f3.write('atom_sigma_t_scaling = %s\n' % atomsts)
            f3.write('amplitude_scaling = %s\n' % amplsc)
            f3.write('n_max = %s\n' % nmax)
            f3.write('l_max = %s\n' % lmax)
            f3.write('nf = %s\n' % nf)
            f3.write('central_weight = %s\n' % centw)
            f3.write('scaling mode = %s\n' % scalm)
            f3.write('basis = "%s"\n' % bas)
            f3.write('radial_enhancement = %s\n' % radenh)
            f3.write('zeta = %s\n' % zeta)
            f3.write('delta = %s\n' % delta)
            f3.write('desc_sparse = \'gap_files/%s.sparseX.GAP_%s\'\n' % (sys.argv[1], potname))
            u = i + 1
            f3.write('alphas_sparse = \'gap_files/alphas_mb_%u.dat\'\n' % u)
            if len(compress) > 0:
                f3.write('compress_soap = .true.\n')
                f3.write('file_compress_soap = \'gap_files/compress_%u.dat\'\n' % u)
            # add hirshfeld info is exist
            if n_hirs > 0:
                f3.write('has_vdw = .true.\n')
                potname_h = potname_h_tmp[i]
                f3.write('vdw_qs = \'gap_files/%s.sparseX.GAP_%s\'\n' % (sys.argv[3], potname_h))
                f3.write('vdw_alphas = \'gap_files/alphas_hirshfeld_%u.dat\'\n' % u)
                f3.write('vdw_zeta = %s\n' % zeta)
                delta_h = float(delta_h_tmp[i])
                f3.write('vdw_delta = %s\n' % delta_h)
                if lp1 == []:
                    keyw_lp = rexh('local_property0={(.*?)}')
                    lp_tmp1 = keyw_lp[0]
                    lp_tmp2 = re.sub(r"(?!(?<=\d)\.(?=\d))[^0-9 ]", " ", lp_tmp1)
                    lp_tmp3 = [float(x) for x in lp_tmp2.split()]
                    lp_tmp4 = np.array(lp_tmp3)
                    lp = lp_tmp4[i]
                else:
                    lp = lp1[0]
                f3.write('vdw_v0 = %s\n' % lp)
                f3.write('gap_end\n')
                f3.write('\n')
                f3.close()

            else:
                f3.write('gap_end\n')
                f3.write('\n')
                f3.close()

# add core potentials if exist
if n_cp > 0:
    dat = rex('<point r=(.*)"/>')
    s_dat = str(dat)
    r = re.findall('"(\s*[-+]?\d*\.\d+|\d+)"', s_dat)
    e = re.findall('E="(-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?)', s_dat)
    nup_tmp = rex('<potential_pair num_points="([-+]?\d*\.\d+|\d+)')
    nup = list(map(int, nup_tmp))
    y1 = rex('" y1=\"(.*)" yn')
    yn = rex('" yn="(.*)"')
    s_am = rex('<Glue_params n_types="(\d)"')
    s_nt = rex('<per_type_data (.*?)">')
    s_at = rex('<per_type_data atomic_num="([-+]?\d*\.\d+|\d+)')
    s_tmp1 = rex('<per_pair_data type1="([-+]?\d*\.\d+|\d+)')
    s_tmp2 = rex('" type2="([-+]?\d*\.\d+|\d+)')

    # create dat files for core_pot
    print("Creating " + str(n_cp) + " dat file(s) for core potential(s)...")

    k = 0
    j = 0
    for i in range(0, n_cp):
        nup2 = nup[j]
        r2 = r[k:k + nup2]
        e2 = e[k:k + nup2]
        td = open("core_pot_%i.dat" % (i + 1), 'w+')
        td.write('%s %s %s\n' % (nup[i], y1[i], yn[i]))
        for i2 in range(len(r2)):
            print(r2[i2], e2[i2], file=td)
        td.close()
        k += nup2
        j += 1

    # write to the gap file
    for i in range(n_cp):
        with open(output, 'a') as f3:
            if int(s_am[0]) == 1:
                s = elem(int(s_at[i]))
                f3.write('gap_beg core_pot\n')
                f3.write('species1 = %s\n' % s)
                f3.write('species2 = %s\n' % s)
                f3.write('core_pot_file = \'gap_files/core_pot_%s.dat\'\n' % (i + 1))
                f3.write('gap_end core_pot\n')
                f3.write('\n')
                f3.close()
                exit
            else:
                for u in range(len(s_nt)):
                    ty = re.findall('type="([-+]?\d*\.\d+|\d+)', s_nt[u])
                    tys = listToString(ty)
                    if s_tmp1[i] == tys:
                        s1 = elem(int(s_at[u]))
                for u in range(len(s_nt)):
                    ty = re.findall('type="([-+]?\d*\.\d+|\d+)', s_nt[u])
                    tys = listToString(ty)
                    if s_tmp2[i] == tys:
                        s2 = elem(int(s_at[u]))
                f3.write('gap_beg core_pot\n')
                f3.write('species1 = %s\n' % s1)
                f3.write('species2 = %s\n' % s2)
                f3.write('core_pot_file = \'gap_files/core_pot_%s.dat\'\n' % (i + 1))
                f3.write('gap_end core_pot\n')
                f3.write('\n')
                f3.close()

print("\nDONE! Do not forget to rename compression files to \"compress_1.dat\" format!")
