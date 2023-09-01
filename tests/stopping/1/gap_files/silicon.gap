gap_beg distance_2b
species1 = Si
species2 = Si
delta = 1.0
sigma = 1.0
rcut = 5.0
desc_sparse = 'gap_files/silicon.xml.sparseX.GAP_2020_10_15_180_13_48_39_9741'
alphas_sparse = 'gap_files/alphas_2b.dat'
gap_end

gap_beg angle_3b
species_center = Si
species1 = Si
species2 = Si
delta = 0.01
sigma = 4. 4. 4.
rcut = 4.
kernel_type = "pol"
desc_sparse = 'gap_files/silicon.xml.sparseX.GAP_2020_10_15_180_13_48_39_9742'
alphas_sparse = 'gap_files/alphas_3b.dat'
gap_end

gap_beg soap_turbo
n_species = 1
species = Si
central_species = 1
rcut = 5.0
buffer = 1.0
atom_sigma_r = 0.5
atom_sigma_t = 0.5
atom_sigma_r_scaling = 0.
atom_sigma_t_scaling = 0.
amplitude_scaling = 1.
n_max = 8
l_max = 8
nf = 4.
central_weight = 1.
scaling_mode = polynomial
basis = "poly3gauss"
radial_enhancement = 1
zeta = 6
delta = 1.0
desc_sparse = 'gap_files/silicon.xml.sparseX.GAP_2020_10_15_180_13_48_39_9743'
alphas_sparse = 'gap_files/alphas_mb.dat'
compress_soap = .true.
file_compress_soap = 'gap_files/compress.dat'
gap_end

gap_beg core_pot
species1 = Si
species2 = Si
core_pot_file = 'gap_files/corepot.dat'
gap_end core_pot

