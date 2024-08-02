gap_beg distance_2b
species1 = P
species2 = P
delta = 1.
sigma = 0.5
rcut = 5.0
desc_sparse = "gap_files/phosphorus.xml.sparseX.GAP_2024_4_8_180_19_35_16_9231"
alphas_sparse = "gap_files/alphas_distance_2b_1.dat"
gap_end

gap_beg soap_turbo
n_species = 1
species = P
central_species = 1
rcut = 5.0
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
scaling_mode = polynomial
basis = "poly3gauss"
radial_enhancement = 1
zeta = 6
delta = 0.1
desc_sparse = "gap_files/phosphorus.xml.sparseX.GAP_2024_4_8_180_19_35_16_9232"
alphas_sparse = "gap_files/alphas_soap_turbo_1.dat"
compress_soap = .true.
compress_mode = "trivial"
has_vdw = .true.
vdw_qs = "gap_files/hirshfeld.xml.sparseX.GAP_2023_12_12_120_17_50_34_2691"
vdw_alphas = "gap_files/alphas_hirshfeld_1.dat"
vdw_zeta = 4
vdw_delta = 1.
vdw_v0 = 0.9
gap_end

gap_beg core_pot
species1 = P
species2 = P
core_pot_file = "gap_files/core_pot_1.dat"
gap_end

