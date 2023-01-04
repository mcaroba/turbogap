! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, types.f90, is copyright (c) 2019-2021, Miguel A. Caro
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module types

  implicit none

! GAP+descriptor data structure for SOAP
  type soap_turbo
    real*8, allocatable :: nf(:), rcut_hard(:), rcut_soft(:), atom_sigma_r(:), atom_sigma_t(:), &
                           atom_sigma_r_scaling(:), atom_sigma_t_scaling(:), amplitude_scaling(:), &
                           central_weight(:), global_scaling(:), alphas(:), Qs(:,:), cutoff(:), &
                           vdw_Qs(:,:), vdw_alphas(:), vdw_cutoff(:)
    real*8 :: zeta = 2.d0, delta = 1.d0, rcut_max, vdw_zeta, vdw_delta, vdw_V0
    integer, allocatable :: alpha_max(:), compress_soap_indices(:)
    integer :: n_species, central_species = 0, dim, l_max, radial_enhancement = 0, n_max, n_sparse, &
               vdw_n_sparse
    character*1024 :: file_alphas, file_desc, file_compress = "none", file_vdw_alphas, file_vdw_desc
    character*64 :: basis = "poly3", compress_mode = "none"
    character*32 :: scaling_mode = "polynomial"
    character*8, allocatable :: species_types(:)
    logical :: compress_soap = .false., has_vdw = .false.
  end type soap_turbo

! GAP+descriptor data structure for distance_2b
  type distance_2b
    real*8, allocatable :: cutoff(:), alphas(:), Qs(:,:)
    real*8 :: delta = 1.d0, sigma = 1.d0, rcut, buffer = 0.5d0
    integer :: dim = 1, n_sparse
    character*1024 :: file_alphas, file_desc
    character*8 :: species1, species2
  end type distance_2b

! GAP+descriptor data structure for distance_2b
  type angle_3b
    real*8, allocatable :: cutoff(:), alphas(:), Qs(:,:)
    real*8 :: delta = 1.d0, sigma(1:3) = 1.d0, rcut, buffer = 0.5d0
    integer :: dim = 3, n_sparse
    character*1024 :: file_alphas, file_desc
    character*8 :: species_center, species1, species2
    character*3 :: kernel_type = "exp"
  end type angle_3b

! Data structure for core_pot
  type core_pot
    real*8, allocatable :: x(:), V(:), dVdx2(:)
    real*8 :: yp1, ypn
    integer :: n
    character*1024 :: core_pot_file
    character*8 :: species1, species2
  end type core_pot

! These is the type for the input parameters
  type input_parameters
    real*8, allocatable :: masses_types(:), e0(:), vdw_c6_ref(:), vdw_r0_ref(:), vdw_alpha0_ref(:)
    real*8 :: t_beg = 300.d0, t_end = 300.d0, tau_t = 100.d0, md_step = 1.d0, &
              neighbors_buffer = 0.d0, max_GBytes_per_process = 1.d0, e_tol = 1.d-6, &
              vdw_sr = 0.94d0, vdw_d = 20.d0, vdw_rcut = 10.d0, &
              vdw_buffer = 1.d0, vdw_rcut_inner = 0.5d0, vdw_buffer_inner = 0.5d0, &
              tau_p = 1000.d0, p_beg = 1.d0, p_end = 1.d0, gamma_p = 1.d0, &
              box_scaling_factor(3,3) = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3,3]), &
              core_pot_cutoff = 1.d10, core_pot_buffer = 1.d0, tau_dt = 100.d0, target_pos_step, &
              gamma0 = 0.01d0, max_opt_step = 0.1d0, vdw_scs_rcut = 4.5d0, vdw_mbd_rcut = 5.d0, vdw_2b_rcut = 9.d0
    integer :: md_nsteps = 1, write_xyz = 0, write_thermo = 1, which_atom = 0, vdw_mbd_nfreq = 12, vdw_mbd_norder = 4
    character*1024 :: atoms_file
    character*32 :: vdw_type = "none"
    character*8, allocatable :: species_types(:)
    character*2 :: optimize = "vv"
    character*32 :: barostat = "none", thermostat = "none", barostat_sym = "isotropic"
    logical :: do_md = .false., do_prediction = .false., do_forces = .false., do_derivatives = .false., &
               do_derivatives_fd = .false., write_soap = .false., write_derivatives = .false., &
               do_timing = .false., all_atoms = .true., print_progress = .true., scale_box = .false., &
               write_lv = .false., write_forces = .true., write_velocities = .true., write_hirshfeld_v = .true., &
               write_virial = .true., write_pressure = .true., write_stress = .true., &
               write_local_energies = .true., write_property(1:11) = .true., &
               write_array_property(1:8) = .true., write_masses = .false., write_fixes = .true., &
               variable_time_step = .false., vdw_mbd_grad = .false., vdw_hirsh_grad = .false., &
               vdw_polynomial = .false.
  end type input_parameters

end module
