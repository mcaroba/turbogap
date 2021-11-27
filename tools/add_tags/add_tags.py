# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   TurboGAP
# HND X
# HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
# HND X
# HND X   TurboGAP is published and distributed under the
# HND X      Academic Software License v1.0 (ASL)
# HND X
# HND X   This file, add_tags.py, is copyright (c) 2019-2021, Miguel A. Caro
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

from ase.io import read, write
import numpy as np
import sys

#############################################################################################
# User modifiable values:

# Input/output filenames
input_database = "train.xyz"
output_database = "train_tagged.xyz"

# If you want to remove the stress tag and only rely on virial
virial_only = True
# If you want to add some custom processing to your database
custom_code = True

# Energy and virial regularization options. All values in eV/atom! The code will adjust for
# the total number of atoms in each configuration. sigma_e is a dictionary. The code will
# look in your input XYZ file for a config_type tag. You need to provide the list of config
# types for which you want to add custom regularization to the dictionary. If it can't find a
# config type in the dictionary (or no config type at all), it will use the "default" value.
# So "default" *needs* to be an entry in the dictionary. sigma_e and sigma_v do not
# necessarily need to contain all (or the same) config types, only those for which you want
# to add specific regularization other than the default ones.
sigma_e = {"default": 0.001}
sigma_v = {"default": 0.1}

# Force regularization options. The force fraction will determine how many forces are masked
# randomly out of the database.
# must be 0. <= force_fraction <= 1.
force_fraction = {"default": 1.,
                  "unrelaxed_CNTs": 1.}
sigma_f_min = 0.1 # minimum force regularization in eV/Angstrom
f_scale = 0.1 # scaling factor: sigma_f = f_scale * f (if sigma_f > sigma_f_min)
force_parameter_name = "forces" # what is the force array called?

#############################################################################################



#############################################################################################
# Do not modify below this line (only if you have custom code)
#############################################################################################
print("Reading database...")
at = read(input_database, index=":")
print("")
print("... done.")

print("")
print("Adding energy, virial and force-component reg. param. and force mask...")
print("")

print("Adding energy, virial and force-component reg. param. and force mask...")
# The structures are looped in reversed order so that the user can safely remove them using custom code
for i in reversed(range(0, len(at))):
#   Bookkeeping
    at[i].pbc = True
    vol = at[i].get_volume()
    if "config_type" in at[i].info:
        ct = at[i].info["config_type"]
    else:
        ct = "default"
#   Custom code
    if custom_code:
        pass
#   Remove stress and replace by virial (if virial not present and stress is), if requested
    if virial_only:
        if "stress" in at[i].info:
            stress = at[i].info["stress"]
            del at[i].info["stress"]
            del at[i].calc.results["stress"]
            virial = -vol * stress
            if "virial" not in at[i].info:
                at[i].info["virial"] = virial
#   Force part
    if ct in force_fraction:
        ff = force_fraction[ct]
    else:
        ff = force_fraction["default"]
    mask = np.random.choice([False, True], size=len(at[i]), p=[ff, 1.-ff])
    forces = at[i].get_array(force_parameter_name)
    at[i].set_array("force_mask", mask)
    f = np.clip( np.abs(forces) * f_scale, a_min = sigma_f_min, a_max = None )
    at[i].set_array("force_component_sigma", f)
#   Energy and virial part
    n = len(at[i])
    if ct in sigma_e:
        se = np.sqrt(n) * sigma_e[ct]
    else:
        se = np.sqrt(n) * sigma_e["default"]
    at[i].info["energy_sigma"] = se
    if ct in sigma_v:
        sv = np.sqrt(n) * sigma_v[ct]
    else:
        sv = np.sqrt(n) * sigma_v["default"]
    at[i].info["virial_sigma"] = sv
#   Print progress
    sys.stdout.write('\rProgress:%6.1f%%' % ((len(at)-float(i))*100./float(len(at))) )
    sys.stdout.flush()

print("")
print("")

print("Writing to file...")
write(output_database, at)
print("")
print("... done.")
