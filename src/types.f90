! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, types.f90, is copyright (c) 2019-2022, Miguel A. Caro and
! HND X   Tigany Zarrouk
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

  type local_property_soap_turbo
     real*8, allocatable :: Qs(:,:), alphas(:), cutoff(:), data(:,:)
     real*8              :: zeta, delta, V0
     character*1024      :: file_data="none", file_alphas, file_desc, label
     integer             :: n_sparse, n_data, dim
     logical             :: has_data = .false., do_derivatives = .false., compute=.true.
  end type local_property_soap_turbo


  type soap_turbo
    real*8, allocatable :: nf(:), rcut_hard(:), rcut_soft(:), atom_sigma_r(:), atom_sigma_t(:), &
                           atom_sigma_r_scaling(:), atom_sigma_t_scaling(:), amplitude_scaling(:), &
                           central_weight(:), global_scaling(:), alphas(:), Qs(:,:), cutoff(:), &
                           vdw_Qs(:,:), vdw_alphas(:), vdw_cutoff(:), compress_P_el(:)
    real*8 :: zeta = 2.d0, delta = 1.d0, rcut_max, vdw_zeta, vdw_delta, vdw_V0
    integer, allocatable :: alpha_max(:), compress_P_i(:), compress_P_j(:)
    integer :: n_species, central_species = 0, dim, l_max, radial_enhancement = 0, n_max, n_sparse, &
               vdw_n_sparse, compress_P_nonzero, n_local_properties=0&
               &, vdw_index=0, core_electron_be_index=0
    character*1024 :: file_alphas, file_desc, file_compress = "none",&
         & file_vdw_alphas, file_vdw_desc
    character*64 :: basis = "poly3", compress_mode = "none"
    character*32 :: scaling_mode = "polynomial"
    character*8, allocatable :: species_types(:)
    logical :: compress_soap = .false., has_vdw = .false.,&
         & has_core_electron_be=.false., has_local_properties = .false.
    type(local_property_soap_turbo), allocatable :: local_property_models(:)
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
     real*8, allocatable :: masses_types(:), e0(:), vdw_c6_ref(:), vdw_r0_ref(:), vdw_alpha0_ref(:), &
          mc_acceptance(:), energy_scales_exp_data(:), radii(:)
    real*8 :: t_beg = 300.d0, t_end = 300.d0, tau_t = 100.d0, md_step = 1.d0, &
              neighbors_buffer = 0.d0, max_GBytes_per_process = 1.d0, e_tol = 1.d-6, &
              vdw_sr = 0.94d0, vdw_d = 20.d0, vdw_rcut = 10.d0, &
              vdw_buffer = 1.d0, vdw_rcut_inner = 0.5d0, vdw_buffer_inner = 0.5d0, &
              tau_p = 1000.d0, p_beg = 1.d0, p_end = 1.d0, gamma_p = 1.d0, &
              box_scaling_factor(3,3) = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3,3]), &
              core_pot_cutoff = 1.d10, core_pot_buffer = 1.d0, tau_dt = 100.d0, target_pos_step, &
              gamma0 = 0.01d0, max_opt_step = 0.1d0, vdw_scs_rcut = 4.d0, f_tol = 0.01d0, p_tol = 0.01d0, &
              max_opt_step_eps = 0.05d0, mc_mu = 0.0d0, t_extra = 0.d0, p_nested = 0.d0, &
              nested_max_strain = 0.d0, nested_max_volume_change = 0.d0, mc_move_max = 1.d0, &
              mc_lnvol_max = 0.01d0, mc_min_dist = 0.2d0, xps_sigma=0.4d0, mc_reverse_lambda = 0.d0

    integer :: md_nsteps = 1, mc_nsteps = 1, write_xyz = 0, write_thermo = 1, which_atom = 0, &
               vdw_mbd_nfreq = 11, n_mc_types = 0, n_nested = 0,&
               & mc_idx = 1, mc_nrelax=0, n_local_properties=0,&
               & xps_n_samples=200, n_moments=0, n_mc_swaps = 0
    integer, allocatable :: mc_swaps_id(:)

    character*1024 :: atoms_file
    character*1024, allocatable :: compute_local_properties(:)
    character*32 :: vdw_type = "none"
    character*32, allocatable ::  mc_types(:)
    character*8, allocatable :: species_types(:), mc_swaps(:)
    character*16 :: optimize = "vv", mc_relax_opt = "gd", mc_hybrid_opt = "vv"
    character*32 :: barostat = "none", thermostat = "none",&
         & barostat_sym = "isotropic", mc_species = "none",&
         & xps_force_type = "similarity", similarity_type = "none"
    logical :: do_md = .false., do_mc = .false., do_prediction =&
         & .false., do_forces = .false., do_derivatives = .false.,&
         & do_derivatives_fd = .false., write_soap = .false.,&
         & write_derivatives = .false., do_timing = .false.,&
         & all_atoms = .true., print_progress = .true., scale_box =&
         & .false., write_lv = .false., write_forces = .true.,&
         & write_velocities = .true., write_hirshfeld_v = .true.,&
         & write_virial = .true., write_pressure = .true.,&
         & write_stress = .true., write_local_energies = .true.,&
         & write_property(1:11) = .true., write_array_property(1:8) =&
         & .true., write_masses = .false., write_fixes = .true.,&
         & variable_time_step = .false., vdw_mbd_grad = .false.,&
         & do_nested_sampling = .false., scale_box_nested = .false.,&
         & mc_write_xyz = .false., mc_relax = .false., mc_opt_spectra&
         &=.false., optimize_exp_data=.false., print_lp_forces&
         &=.false., mc_hamiltonian = .false., accessible_volume =&
         & .false., mc_reverse = .true.
    logical, allocatable :: write_local_properties(:)

  end type input_parameters


! This is a container for atomic images
  type image
    real*8, allocatable :: positions(:,:), positions_prev(:,:), velocities(:,:), masses(:), &
                           forces(:,:), forces_prev(:,:), energies(:), local_properties(:,:)
    real*8 :: a_box(1:3), b_box(1:3), c_box(1:3), energy, e_kin
    integer, allocatable :: species(:), species_supercell(:)
    integer :: n_sites, indices(1:3)
    logical, allocatable :: fix_atom(:,:)
    character*8, allocatable :: xyz_species(:), xyz_species_supercell(:)
  end type image





  contains


!**************************************************************************
! This provides a way to pass all the individual arrays/variables in the main code to an image container
! In time I should make the image data type the default way to store these properties!!!!!!!
  subroutine from_properties_to_image(this_image, positions, velocities, masses, &
                                      forces, a_box, b_box, c_box, energy, energies, e_kin, &
                                      species, species_supercell, n_sites, indices, fix_atom, &
                                      xyz_species, xyz_species_supercell, local_properties)
    implicit none

!   Input variables
    real*8, intent(in) :: positions(:,:), velocities(:,:), masses(:), energies(:), &
         forces(:,:), a_box(1:3), b_box(1:3), c_box(1:3), energy, e_kin
    real*8, allocatable, intent(in) :: local_properties(:,:)
    integer, intent(in) :: species(:), species_supercell(:), n_sites, indices(1:3)
    logical, intent(in) :: fix_atom(:,:)
    character*8, intent(in) :: xyz_species(:), xyz_species_supercell(:)
!   In/out variables
    type(image), intent(inout) :: this_image
!   Internal variables
    integer :: n, n2

    n = size(positions, 2)
    if( allocated( this_image%positions ) )deallocate( this_image%positions )
    allocate( this_image%positions(1:3, 1:n) )
    this_image%positions = positions

    n = size(velocities, 2)
    if( allocated( this_image%velocities ) )deallocate( this_image%velocities )
    allocate( this_image%velocities(1:3, 1:n) )
    this_image%velocities = velocities

    n = size(masses, 1)
    if( allocated( this_image%masses ) )deallocate( this_image%masses )
    allocate( this_image%masses(1:n) )
    this_image%masses = masses

    n = size(energies, 1)
    if( allocated( this_image%energies ) )deallocate( this_image%energies )
    allocate( this_image%energies(1:n) )
    this_image%energies = energies


    n = size(forces, 2)
    if( allocated( this_image%forces ) )deallocate( this_image%forces )
    allocate( this_image%forces(1:3, 1:n) )
    this_image%forces = forces

    this_image%a_box = a_box

    this_image%b_box = b_box

    this_image%c_box = c_box

    this_image%energy = energy

    this_image%e_kin = e_kin

    n = size(species, 1)
    if( allocated( this_image%species ) )deallocate( this_image%species )
    allocate( this_image%species(1:n) )
    this_image%species = species

    n = size(species_supercell, 1)
    if( allocated( this_image%species_supercell ) )deallocate( this_image%species_supercell )
    allocate( this_image%species_supercell(1:n) )
    this_image%species_supercell = species_supercell

    this_image%n_sites = n_sites

    this_image%indices = indices

    n = size(fix_atom, 2)
    if( allocated( this_image%fix_atom ) )deallocate( this_image%fix_atom )
    allocate( this_image%fix_atom(1:3, 1:n) )
    this_image%fix_atom = fix_atom

    n = size(xyz_species, 1)
    if( allocated( this_image%xyz_species ) )deallocate( this_image%xyz_species )
    allocate( this_image%xyz_species(1:n) )
    this_image%xyz_species = xyz_species

    n = size(xyz_species_supercell, 1)
    if( allocated( this_image%xyz_species_supercell ) )deallocate( this_image%xyz_species_supercell )
    allocate( this_image%xyz_species_supercell(1:n) )
    this_image%xyz_species_supercell = xyz_species_supercell

    if(allocated(local_properties))then
       n = size(local_properties, 1)
       n2 = size(local_properties, 2)
       if( allocated( this_image%local_properties ) )deallocate( this_image%local_properties )
       allocate( this_image%local_properties(1:n, 1:n2) )
       this_image%local_properties = local_properties
    end if

  end subroutine
!**************************************************************************





!**************************************************************************
  subroutine from_image_to_properties(this_image, positions, velocities, masses, &
                                      forces, a_box, b_box, c_box, energy, energies, e_kin, &
                                      species, species_supercell, n_sites, indices, fix_atom, &
                                      xyz_species, xyz_species_supercell, local_properties)
    implicit none

!   Input variables
    type(image), intent(in) :: this_image
!   Output variables
    real*8, allocatable, intent(out) :: positions(:,:), velocities(:,:), masses(:), &
         forces(:,:), energies(:)
    real*8, allocatable, intent(out) :: local_properties(:,:)
    real*8, intent(out) :: a_box(1:3), b_box(1:3), c_box(1:3), energy, e_kin
    integer, allocatable, intent(out) :: species(:), species_supercell(:)
    integer, intent(out) :: n_sites, indices(1:3)
    logical, allocatable, intent(out) :: fix_atom(:,:)
    character*8, allocatable, intent(out) :: xyz_species(:), xyz_species_supercell(:)
!   Internal variables
    integer :: n,n2

    n = size(this_image%positions, 2)
    allocate( positions(1:3, 1:n) )
    positions = this_image%positions

    n = size(this_image%velocities, 2)
    allocate( velocities(1:3, 1:n) )
    velocities = this_image%velocities

    n = size(this_image%masses, 1)
    allocate( masses(1:n) )
    masses = this_image%masses

    n = size(this_image%energies, 1)
    allocate( energies(1:n) )
    energies = this_image%energies


    n = size(this_image%forces, 2)
    allocate( forces(1:3, 1:n) )
    forces = this_image%forces

    a_box = this_image%a_box

    b_box = this_image%b_box

    c_box = this_image%c_box

    energy = this_image%energy

    e_kin = this_image%e_kin

    n = size(this_image%species, 1)
    allocate( species(1:n) )
    species = this_image%species

    n = size(this_image%species_supercell, 1)
    allocate( species_supercell(1:n) )
    species_supercell = this_image%species_supercell

    n_sites = this_image%n_sites

    indices = this_image%indices

    n = size(this_image%fix_atom, 2)
    allocate( fix_atom(1:3, 1:n) )
    fix_atom = this_image%fix_atom

    n = size(this_image%xyz_species, 1)
    allocate( xyz_species(1:n) )
    xyz_species = this_image%xyz_species

    n = size(this_image%xyz_species_supercell, 1)
    allocate( xyz_species_supercell(1:n) )
    xyz_species_supercell = this_image%xyz_species_supercell

    if(allocated(this_image%local_properties))then
       n = size(this_image%local_properties, 1)
       n2 = size(this_image%local_properties, 2)
       allocate( local_properties(1:n,1:n2) )
       local_properties = this_image%local_properties
    end if

  end subroutine
!**************************************************************************


end module
