gap_beg distance_2b
species1 = O
species2 = O
delta = 0.5
sigma = 0.5
rcut = 4.5
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4141"
alphas_sparse = "gap_files/alphas_distance_2b_1.dat"
gap_end

gap_beg distance_2b
species1 = O
species2 = C
delta = 0.5
sigma = 0.5
rcut = 4.5
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4142"
alphas_sparse = "gap_files/alphas_distance_2b_2.dat"
gap_end

gap_beg distance_2b
species1 = C
species2 = C
delta = 0.5
sigma = 0.5
rcut = 4.5
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4143"
alphas_sparse = "gap_files/alphas_distance_2b_3.dat"
gap_end

gap_beg angle_3b
species_center = O
species1 = O
species2 = O
delta = 0.01
sigma = 0.20000000000000000E+001  0.20000000000000000E+001  0.20000000000000000E+001
kernel_type = pol
rcut = 2.0
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4144"
alphas_sparse = "gap_files/alphas_angle_3b_1.dat"
gap_end

gap_beg angle_3b
species_center = O
species1 = O
species2 = C
delta = 0.01
sigma = 0.20000000000000000E+001  0.20000000000000000E+001  0.20000000000000000E+001
kernel_type = pol
rcut = 2.0
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4145"
alphas_sparse = "gap_files/alphas_angle_3b_2.dat"
gap_end

gap_beg angle_3b
species_center = O
species1 = C
species2 = C
delta = 0.01
sigma = 0.20000000000000000E+001  0.20000000000000000E+001  0.20000000000000000E+001
kernel_type = pol
rcut = 2.0
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4146"
alphas_sparse = "gap_files/alphas_angle_3b_3.dat"
gap_end

gap_beg angle_3b
species_center = C
species1 = O
species2 = O
delta = 0.01
sigma = 0.20000000000000000E+001  0.20000000000000000E+001  0.20000000000000000E+001
kernel_type = pol
rcut = 2.0
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4147"
alphas_sparse = "gap_files/alphas_angle_3b_4.dat"
gap_end

gap_beg angle_3b
species_center = C
species1 = O
species2 = C
delta = 0.01
sigma = 0.20000000000000000E+001  0.20000000000000000E+001  0.20000000000000000E+001
kernel_type = pol
rcut = 2.0
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4148"
alphas_sparse = "gap_files/alphas_angle_3b_5.dat"
gap_end

gap_beg angle_3b
species_center = C
species1 = C
species2 = C
delta = 0.01
sigma = 0.20000000000000000E+001  0.20000000000000000E+001  0.20000000000000000E+001
kernel_type = pol
rcut = 2.0
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_4149"
alphas_sparse = "gap_files/alphas_angle_3b_6.dat"
gap_end

gap_beg soap_turbo
n_species = 2
species = C O
central_species = 1
rcut = 4.5 4.5
buffer = 0.5 0.5
atom_sigma_r = 0.2 0.2
atom_sigma_t = 0.2 0.2
atom_sigma_r_scaling = 0.1 0.1
atom_sigma_t_scaling = 0.1 0.1
amplitude_scaling = 2.0 2.0
n_max = 8 8
l_max = 8
nf = 4 4
central_weight = 1.0 1.0
scaling mode = polynomial
basis = "poly3gauss"
radial_enhancement = 1
zeta = 4
delta = 0.1
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_41410"
alphas_sparse = "gap_files/alphas_soap_turbo_1.dat"
compress_soap = .true.
compress_mode = "trivial"
has_local_properties = .true.
n_local_properties = 1
local_property_qs = 'gap_files/core_electron_be.xml.sparseX.GAP_2024_1_26_120_18_36_33_6091'
local_property_alphas = 'gap_files/alphas_core_electron_be_1.dat'
local_property_labels = 'core_electron_be'
local_property_zetas = 2
local_property_deltas = 1.0
local_property_v0s = 290.0
gap_end

gap_beg soap_turbo
n_species = 2
species = C O
central_species = 2
rcut = 4.5 4.5
buffer = 0.5 0.5
atom_sigma_r = 0.2 0.2
atom_sigma_t = 0.2 0.2
atom_sigma_r_scaling = 0.1 0.1
atom_sigma_t_scaling = 0.1 0.1
amplitude_scaling = 2.0 2.0
n_max = 8 8
l_max = 8
nf = 4 4
central_weight = 1.0 1.0
scaling mode = polynomial
basis = "poly3gauss"
radial_enhancement = 1
zeta = 4
delta = 0.1
desc_sparse = "gap_files/CO.xml.sparseX.GAP_2022_5_19_180_7_17_31_41411"
alphas_sparse = "gap_files/alphas_soap_turbo_2.dat"
compress_soap = .true.
compress_mode = "trivial"
has_local_properties = .true.
n_local_properties = 1
local_property_qs = 'gap_files/core_electron_be.xml.sparseX.GAP_2024_1_26_120_18_36_33_6092'
local_property_alphas = 'gap_files/alphas_core_electron_be_2.dat'
local_property_labels = 'core_electron_be'
local_property_zetas = 2
local_property_deltas = 1.0
local_property_v0s = 536.0
gap_end

gap_beg core_pot
species1 = C
species2 = C
core_pot_file = "gap_files/core_pot_1.dat"
gap_end

gap_beg core_pot
species1 = C
species2 = O
core_pot_file = "gap_files/core_pot_2.dat"
gap_end

gap_beg core_pot
species1 = O
species2 = O
core_pot_file = "gap_files/core_pot_3.dat"
gap_end

