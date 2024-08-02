gap_beg distance_2b
species1 = C
species2 = C
delta = 1.
sigma = 0.5
rcut = 4.5
desc_sparse = "gap_files/carbon.xml.sparseX.GAP_2023_10_28_180_1_14_11_7511"
alphas_sparse = "gap_files/alphas_distance_2b_1.dat"
gap_end

gap_beg angle_3b
species_center = C
species1 = C
species2 = C
delta = 0.01
sigma = 0.40000000000000000E+001  0.40000000000000000E+001  0.40000000000000000E+001
kernel_type = pol
rcut = 2.4
desc_sparse = "gap_files/carbon.xml.sparseX.GAP_2023_10_28_180_1_14_11_7512"
alphas_sparse = "gap_files/alphas_angle_3b_1.dat"
gap_end

gap_beg soap_turbo
n_species = 1
species = C
central_species = 1
rcut = 4.5
buffer = 0.5
atom_sigma_r = 0.5
atom_sigma_t = 0.5
atom_sigma_r_scaling = 0.
atom_sigma_t_scaling = 0.
amplitude_scaling = 1.0
n_max = 8
l_max = 8
nf = 4
central_weight = 1.0
scaling mode = polynomial
basis = "poly3gauss"
radial_enhancement = 1
zeta = 6
delta = 0.1
desc_sparse = "gap_files/carbon.xml.sparseX.GAP_2023_10_28_180_1_14_11_7513"
alphas_sparse = "gap_files/alphas_soap_turbo_1.dat"
compress_soap = .true.
compress_mode = "trivial"
has_vdw = .true.
vdw_qs = "gap_files/hirshfeld.xml.sparseX.GAP_2021_3_9_120_9_59_51_5751"
vdw_alphas = "gap_files/alphas_hirshfeld_1.dat"
vdw_zeta = 6
vdw_delta = 1.0
vdw_v0 = 0.9
gap_end

gap_beg core_pot
species1 = C
species2 = C
core_pot_file = "gap_files/core_pot_1.dat"
gap_end

