# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   TurboGAP
# HND X
# HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
# HND X
# HND X   TurboGAP is published and distributed under the
# HND X      Academic Software License v1.0 (ASL)
# HND X
# HND X   This file, make_gap_files.py, is copyright (c) 2019-2021, Mikhail S. Kuklin,
# HND X   Miguel A. Caro and Richard Jana
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

# issues:
    # in theory, each soap could have its own hirshfeld xml, but this is currently not supported, as the necessary information is not in the xml file!
    # What is nf? Where does the value 4 come from?

#  !/usr/bin/env python3

from bs4 import BeautifulSoup
import numpy as np
import os
from pathlib import Path
import shutil
import sys

print("\n		MAKE GAP FILES | v. 3.0, 25.03.2022 | Author(s): Mikhail S. Kuklin, Miguel A. Caro, Richard Jana\n")

# dictionary of the elements: 1-57, 72-89
elements = {1: ' H', 2: ' He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
            11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
            21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
            31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
            41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
            51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
            72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
            81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac'}

# extract covariance type
#cov_type = {'1': 'exp', '4': 'pol'}
cov_type = {'ard_se': 'exp', 'pp': 'pol'}

def write_alphas_file(gpCoordinates, index, hirshfeld=False):
    if hirshfeld == True:
        with open(f"gap_files/alphas_hirshfeld_{index}.dat", "w+") as alpha:
            for sx in gpCoordinates.find_all('sparseX'):
                alpha.write(f"{sx['alpha']}\n")
    else:
        type = gpCoordinates.find('descriptor').get_text().split()[0]
        with open(f"gap_files/alphas_{type}_{index}.dat", "w+") as alpha:
            for sx in gpCoordinates.find_all('sparseX'):
                alpha.write(f"{sx['alpha']} {sx['sparseCutoff']}\n")

    # copy the original xml file to gap_files/
    shutil.copyfile(gpCoordinates['sparseX_filename'], f"gap_files/{gpCoordinates['sparseX_filename']}")

def descriptor_to_dict(descriptor):
    # there can be rogue spaces (e.g. in multi species terms like {8 8}),
    # therefor just splitting at ' ' wouldn't work reliably!
    dict = {'type': descriptor.get_text().split()[0]}
    keys = []
    vals = []
    for string in (descriptor.get_text()+' null').split('='):
        keys.append(string.split()[-1])
        vals.append(' '.join(string.split()[:-1]))
    for key, val in zip(keys[:-1], vals[1:]):
        dict[key] = val

    return dict

def compare_descriptors(descriptor1, descriptor2):
    ignore_fields = ['config_type_n_sparse', 'n_sparse']
    dict1 = descriptor_to_dict(descriptor1)
    dict2 = descriptor_to_dict(descriptor2)
    keys = list(set(dict1.keys()).union(set(dict2.keys())) - set(ignore_fields))

    for key in keys:
        try:
            if dict1[key] != dict2[key]:
                return False
        except:
            return False

    return True

def write_descriptor_to_output(output_file, gpCoordinates, index):
    descriptor = gpCoordinates.find('descriptor')
    desc_dict = descriptor_to_dict(descriptor)

    if desc_dict['type'] == 'distance_2b':
        with open(output_file, 'a') as output:
            output.write("gap_beg distance_2b\n")
            output.write(f"species1 = {elements[int(desc_dict['Z1'])]}\n")
            output.write(f"species2 = {elements[int(desc_dict['Z2'])]}\n")
            output.write(f"delta = {desc_dict['delta']}\n")
            output.write(f"sigma = {desc_dict['theta_uniform']}\n")
            output.write(f"rcut = {desc_dict['cutoff']}\n")
            output.write('desc_sparse = "gap_files/' + gpCoordinates['sparseX_filename'] + '"\n')
            output.write(f'alphas_sparse = "gap_files/alphas_' + f"{desc_dict['type']}_{index}.dat" + '"\n')
            output.write("gap_end\n")
            output.write("\n")

    if desc_dict['type'] == 'angle_3b':
        with open(output_file, 'a') as output:
            output.write("gap_beg angle_3b\n")
            output.write(f"species_center = {elements[int(desc_dict['Z_center'])]}\n")
            output.write(f"species1 = {elements[int(desc_dict['Z1'])]}\n")
            output.write(f"species2 = {elements[int(desc_dict['Z2'])]}\n")
            output.write(f"delta = {desc_dict['delta']}\n")
            output.write(f"sigma = {gpCoordinates.find('theta').get_text().strip()}\n")
            output.write(f"kernel_type = {cov_type[desc_dict['covariance_type']]}\n")
            output.write(f"rcut = {desc_dict['cutoff']}\n")
            output.write('desc_sparse = "gap_files/' + gpCoordinates['sparseX_filename'] + '"\n')
            output.write(f'alphas_sparse = "gap_files/alphas_' + f"{desc_dict['type']}_{index}.dat" + '"\n')
            output.write("gap_end\n")
            output.write("\n")

    if desc_dict['type'] == 'soap_turbo':
        with open(output_file, 'a') as output:
            output.write("gap_beg soap_turbo\n")
            output.write(f"n_species = {desc_dict['n_species']}\n")
            output.write("species =")
            for z in desc_dict['species_Z'].strip('{}').split():
                output.write(f" {elements[int(z)]}")
            output.write("\n")
            output.write(f"central_species = {desc_dict['central_index']}\n")
            output.write("rcut =")
            for i in range(int(desc_dict['n_species'])):
                output.write(f" {desc_dict['rcut_hard']}")
            output.write("\n")
            output.write("buffer =")
            for i in range(int(desc_dict['n_species'])):
                output.write(f" {float(desc_dict['rcut_hard'])-float(desc_dict['rcut_soft'])}")
            output.write("\n")
            output.write(f"atom_sigma_r = {desc_dict['atom_sigma_r'].strip('{}')}\n")
            output.write(f"atom_sigma_t = {desc_dict['atom_sigma_t'].strip('{}')}\n")
            output.write(f"atom_sigma_r_scaling = {desc_dict['atom_sigma_r_scaling'].strip('{}')}\n")
            output.write(f"atom_sigma_t_scaling = {desc_dict['atom_sigma_t_scaling'].strip('{}')}\n")
            output.write(f"amplitude_scaling = {desc_dict['amplitude_scaling'].strip('{}')}\n")
            output.write(f"n_max = {desc_dict['alpha_max'].strip('{}')}\n")
            output.write(f"l_max = {desc_dict['l_max']}\n")
            output.write("nf =")
            for i in range(int(desc_dict['n_species'])):
                output.write(f" {4}")
            output.write("\n")
            output.write(f"central_weight = {desc_dict['central_weight'].strip('{}')}\n")
            output.write(f"scaling mode = {desc_dict['scaling_mode']}\n")
            output.write("basis = " + '"' + desc_dict['basis'] + '"' + "\n")
            output.write(f"radial_enhancement = {desc_dict['radial_enhancement'].strip('{}')}\n")
            output.write(f"zeta = {desc_dict['zeta']}\n")
            output.write(f"delta = {desc_dict['delta']}\n")
            output.write('desc_sparse = "gap_files/' + gpCoordinates['sparseX_filename'] + '"\n')
            output.write(f'alphas_sparse = "gap_files/alphas_' + f"{desc_dict['type']}_{index}.dat" + '"\n')
            if 'compress_file' in desc_dict:
                output.write('compress_soap = .true.\n')
                output.write('file_compress_soap = "gap_files/' + f"{desc_dict['compress_file']}" + '"\n')
                # copy the compression file to gap_files/
                shutil.copyfile(desc_dict['compress_file'], f"gap_files/{desc_dict['compress_file']}")
            if descriptor_counts['hirshfeld'] > 0:
                # find right Hirshfeld descriptor (if any)
                for index_hirsh, gpc in enumerate(hirshfeld_soup.find_all('gpCoordinates')):
                    if compare_descriptors(gpc.find('descriptor'), gpCoordinates.find('descriptor')):
                        hirshfeld_dict = descriptor_to_dict(gpc.find('descriptor'))
                        output.write('has_vdw = .true.\n')
                        output.write('vdw_qs = "gap_files/' + hirshfeld_soup.find('gpCoordinates')['sparseX_filename'] + '"\n')
                        output.write(f'vdw_alphas = "gap_files/alphas_hirshfeld_{index_hirsh+1}.dat"\n')
                        output.write(f"vdw_zeta = {hirshfeld_dict['zeta']}\n")
                        output.write(f"vdw_delta = {hirshfeld_dict['delta']}\n")
                        for word in hirshfeld_soup.find('command_line').get_text().split():
                            try:
                                key, entry = word.split('=')
                                if key == 'local_property0':
                                    local_property0 = entry
                            except:
                                continue
                        output.write(f"vdw_v0 = {local_property0.strip('{}')}\n")
                        break # only include the first matching Hirshfeld
            output.write("gap_end\n")
            output.write("\n")

def write_pairpot_file(per_pair_data, index):
    num_points = per_pair_data.find('potential_pair')['num_points']
    y1 = per_pair_data.find('potential_pair')['y1']
    yn = per_pair_data.find('potential_pair')['yn']
    with open(f"gap_files/core_pot_{index}.dat", "w+") as output:
        output.write(f"{num_points} {y1} {yn}\n")
        for point in per_pair_data.find_all('point'):
            output.write(f"{point['r']} {point['E']}\n")

def write_pairpot_to_output(output_file, per_pair_data, index):
    with open(output_file, 'a') as output:
        output.write("gap_beg core_pot\n")
        output.write(f"species1 = {elements[int(per_type_data[per_pair_data['type1']])]}\n")
        output.write(f"species2 = {elements[int(per_type_data[per_pair_data['type2']])]}\n")
        output.write(f'core_pot_file = "gap_files/core_pot_' + str(index) + '.dat"\n')
        output.write("gap_end\n")
        output.write("\n")


# open potential xml file as the first argument in the command line
gap_file = sys.argv[1]
# give name of the output file as the second argument in the command when script is executed
output_file = f"gap_files/{sys.argv[2]}"
if Path(output_file).is_file():
    print("WARNING: You already have output file in the folder with the same name. Remove the file and execute script again!")
    exit()
if not Path('gap_files').is_dir():
    os.mkdir('gap_files') # make dir for gap files, if it doesn't already exist

# copy the original xml file to gap_files/
shutil.copyfile(gap_file, f"gap_files/{gap_file}")

print(f"Reading {gap_file} ...")
with open(gap_file) as file:
    xml_soup = BeautifulSoup(file,'xml') # 'xml' instead of 'lxml' to retain capital letters

descriptor_counts = {'distance_2b': 0,
                     'distance_Nb': 0,
                     'angle_3b': 0,
                     'soap_turbo': 0,
                     'pairpot': 0,
                     'hirshfeld': 0}

if len(sys.argv) > 3:
    hirshfeld = sys.argv[3]
    print(f"Reading {hirshfeld} ...")
    with open(hirshfeld) as file:
        hirshfeld_soup = BeautifulSoup(file,'xml')

    for gpc in hirshfeld_soup.find_all('gpCoordinates'):
        descriptor_counts['hirshfeld'] += 1
        write_alphas_file(gpc, descriptor_counts['hirshfeld'], hirshfeld=True)

for gpc in xml_soup.find_all('gpCoordinates'):
    type = gpc.find('descriptor').get_text().split()[0]
    descriptor_counts[type] += 1
    write_alphas_file(gpc, descriptor_counts[type])
    write_descriptor_to_output(output_file, gpc, descriptor_counts[type])

if xml_soup.find('pairpot'):
    per_type_data = {} # convert 'type1' etc. into atomic numbers
    for ptd in xml_soup.find_all('per_type_data'):
        per_type_data[ptd['type']] = int(ptd['atomic_num'])
    for pp in xml_soup.find('pairpot').find_all('per_pair_data'):
        descriptor_counts['pairpot'] += 1
        write_pairpot_file(pp, descriptor_counts['pairpot'])
        write_pairpot_to_output(output_file, pp, descriptor_counts['pairpot'])

if descriptor_counts['distance_2b'] == 0:
    descriptor_counts['distance_2b'] = descriptor_counts['distance_Nb']
print(f"Number of distance_2b => {descriptor_counts['distance_2b']}")
print(f"Number of angle_3b => {descriptor_counts['angle_3b']}")
print(f"Number of soap_turbo => {descriptor_counts['soap_turbo']}")
if descriptor_counts['pairpot'] > 0:
    print(f"Core potentials found => {descriptor_counts['pairpot']}")
if descriptor_counts['hirshfeld'] > 0:
    print(f"Hirshfeld volume GAP(s) are found => {descriptor_counts['hirshfeld']}")
