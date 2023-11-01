gap_beg distance_2b
species1 = C
species2 = C
delta = 1.
sigma = 0.5
rcut = 4.5
desc_sparse = 'gap_files/carbon.xml.sparseX.GAP_2021_3_12_120_10_53_57_7991'
alphas_sparse = 'gap_files/alphas_2b.dat'
gap_end

gap_beg angle_3b
species_center = C
species1 = C
species2 = C
delta = 0.01
sigma = 4. 4. 4.
rcut = 2.8
kernel_type = "pol"
desc_sparse = 'gap_files/carbon.xml.sparseX.GAP_2021_3_12_120_10_53_57_7992'
alphas_sparse = 'gap_files/alphas_3b.dat'
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
amplitude_scaling = 1.
n_max = 8
l_max = 8
nf = 4.
central_weight = 1.
scaling_mode = polynomial
basis = "poly3gauss"
radial_enhancement = 1
zeta = 6
delta = 0.1
desc_sparse = 'gap_files/carbon.xml.sparseX.GAP_2021_3_12_120_10_53_57_7993'
alphas_sparse = 'gap_files/alphas_mb.dat'
compress_soap = .true.
file_compress_soap = 'gap_files/compress.dat'
has_vdw = .true.
vdw_qs = 'gap_files/hirshfeld.xml.sparseX.GAP_2021_3_9_120_9_59_51_5751'
vdw_alphas = 'gap_files/alphas_hirshfeld.dat'
vdw_zeta = 6
vdw_delta = 1.
vdw_v0 = 0.9
gap_end

! This is the repulsive core potential
gap_beg core_pot
species1 = C
species2 = C
core_pot_file = 'gap_files/corepot.dat'
gap_end core_pot

! This is the dispersion potential
! You can choose between two different versions, roughly equivalent
! to a Tkatchenko-Scheffler approach with the C6 coefficient fixed to
! that of graphene. The two versions differ on the value of sR (the
! TS screening parameter for the London dispersion damping):
! 0.94 for the standard screening
! 0.893 for the screening optimized to C
! You can enable the dispersion correction by uncommenting this block,
! and switch between the two available implementations by commenting
! or uncommenting one of them (if you leave both uncommented, the last
! one will be read.
!gap_beg core_pot
!species1 = C
!species2 = C
!core_pot_file = 'gap_files/disppot_0.893.dat'
!core_pot_file = 'gap_files/disppot_0.940.dat'
!gap_end core_pot

