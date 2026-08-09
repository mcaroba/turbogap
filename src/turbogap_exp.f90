! Copyright (c) 2020-2026 by Albert Bartok and Miguel Caro
!
! The XPS spectrum predicted from core-electron binding energies, and the
! energies and forces that fitting it against experiment produces.
!
! Lifted verbatim from turbogap.f90 (lines 1958-2096 at the time of the move).
! This is the GPU branch's copy of the CPU branch's turbogap_exp, with the same
! module name and the same procedure name, so the two can be compared module
! against module.
!
! The one call whose argument list was split across #ifdef _MPIF90 collapses
! here -- inside a procedure the dummy has one name either way, so the caller
! chooses once. The other #ifdef in this block has no #else and guards the
! MPI-only allocation of those same arrays; that one is an ordinary conditional
! over whole statements and survives the move unchanged.
!
! compute_exp_spectra -- the pair-distribution, structure-factor and diffraction
! block -- is now here too, so this module holds the same two procedures as the
! CPU branch's copy and the two can be compared procedure against procedure.
!
! It was held up on the device state (cuBLAS handles, streams, gpu_exp,
! gpu_neigh) appearing to cross the boundary. Most of it does not: the count
! that said 71 was mostly declaration-only crossings of buffers private to the
! block. What genuinely crosses is eight symbols, and they now come from
! gpu_context by USE, which is why the body could move unchanged.
!
! It is 51 arguments against the CPU branch's 37. The difference is the
! batched-device machinery -- the batch decomposition lists, the per-batch
! bounds, and omp_task, which stays an argument because it is meant to be
! OpenMP-private and must not become shared module state.
module turbogap_exp

   use kinds

   use timing
   use types
   use exp_utils
   use exp_interface
   use iso_c_binding
   use F_B_C
   use gpu_context
#ifdef _MPIF90
   use mpi
   use mpi_helper
#endif

   implicit none

   private
   public :: compute_exp_xps
   public :: compute_exp_spectra

contains

!**************************************************************************
   subroutine compute_exp_xps(params, n_sites, n_xyz, xyz, neighbors_list, n_neigh, &
                              local_properties, local_properties_cart_der, soap_turbo_hypers, &
                              a_box, b_box, c_box, indices, i_beg, i_end, j_beg, j_end, rank, &
                              md_istep, mc_istep, valid_xps, xps_idx, core_be_lp_index, &
                              write_condition, overwrite_condition, exp_output, &
                              energies_lp, forces_lp, virial_lp, time)

      implicit none

      type(input_parameters), intent(inout) :: params
      type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)
      integer, intent(in) :: n_sites
      integer, intent(in) :: n_xyz
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: rank
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: xps_idx
      integer, intent(in) :: core_be_lp_index
      logical, intent(in) :: valid_xps
      logical, intent(inout) :: write_condition
      logical, intent(inout) :: overwrite_condition
      character*32, intent(inout) :: exp_output
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in), allocatable :: xyz(:, :)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      real(dp), intent(inout), allocatable :: local_properties(:, :)
      real(dp), intent(inout), allocatable :: local_properties_cart_der(:, :, :)

!   The caller passes the this_-prefixed arrays under MPI and the plain ones
!   otherwise.
      real(dp), allocatable, intent(inout) :: energies_lp(:)
      real(dp), allocatable, intent(inout) :: forces_lp(:, :)
      real(dp), intent(inout) :: virial_lp(1:3, 1:3)
      type(times_t), intent(inout) :: time

!   Was driver state; block-local now.
      real(dp), allocatable :: x_xps(:)
      real(dp), allocatable :: y_xps(:)
      real(dp), allocatable :: y_i_pred_all(:, :)
      real(dp), allocatable :: v_neigh_lp(:)
      real(dp) :: mag
      integer :: i
      integer :: j
      integer :: j2
      integer :: k

      !     Compute core_electron_be energies and forces
      if (any(soap_turbo_hypers(:)%has_core_electron_be) .and. (params%do_prediction) &
          .and. valid_xps) then

         !time%xps(1) = MPI_wtime()
         call time_start(time%xps, "xps")

#ifdef _MPIF90
         allocate (energies_lp(1:n_sites))
         energies_lp = 0.d0
         if (params%do_forces) then
            allocate (forces_lp(1:3, 1:n_sites))
            forces_lp = 0.d0
            virial_lp = 0.d0
         end if
#endif
         allocate (v_neigh_lp(1:j_end - j_beg + 1))
         v_neigh_lp = 0.d0
         k = 0
         do i = i_beg, i_end
            do j = 1, n_neigh(i)
               !           I'm not sure if this is necessary or neighbors_list is already bounded between 1 and n_sites -> CHECK THIS
               j2 = mod(neighbors_list(j_beg + k) - 1, n_sites) + 1
               k = k + 1
               !                 v_neigh_lp(k) = hirshfeld_v(j2)
               v_neigh_lp(k) = local_properties(j2, core_be_lp_index)
            end do
         end do

         call get_xps_spectra(params%exp_data(xps_idx)%data(1, :),&
           & params%exp_data(xps_idx)%data(2, :), params&
           &%xps_sigma, params%exp_data(xps_idx)%n_samples, mag,&
           & params%exp_data(xps_idx)%x, params&
           &%exp_data(xps_idx)%y_pred, y_i_pred_all,&
           & local_properties(1:n_sites, core_be_lp_index),&
           & .true.)

         ! call get_compare_xps_spectra(params%exp_data(xps_idx)%data&
         !      & , local_properties(1:n_sites, core_be_lp_index),&
         !      & params%xps_sigma, params%exp_data(xps_idx) &
         !      &%n_samples, mag, params%exp_data(xps_idx)%similarity&
         !      & , params%exp_data(xps_idx)%x, params &
         !      &%exp_data(xps_idx)%y, params%exp_data(xps_idx) &
         !      &%y_pred, y_i_pred_all, .not. allocated(params &
         !      &%exp_data(xps_idx)%x), params%exp_similarity_type )

         ! print *, params%exp_data(xps_idx)%n_samples, xps_idx
         call get_energy_scale(params%do_md, params%do_mc,&
           & md_istep, params%md_nsteps, mc_istep, params&
           &%mc_nsteps, params &
           &%exp_energy_scales_initial(xps_idx), params &
           &%exp_energy_scales_final(xps_idx), params &
           &%exp_energy_scales(xps_idx))

         call get_exp_pred_spectra_energies_forces(params&
           &%exp_energy_scales(xps_idx),&
           & local_properties(i_beg:i_end, core_be_lp_index),&
           & local_properties_cart_der(1:3, j_beg:j_end,&
           & core_be_lp_index), n_neigh(i_beg:i_end),&
           & neighbors_list(j_beg:j_end), params%xps_sigma,&
           & params%exp_data(xps_idx)%n_samples, mag, params&
           &%exp_data(xps_idx)%x, params%exp_data(xps_idx)%y,&
           & params%exp_data(xps_idx)%y_pred,&
           & y_i_pred_all(i_beg:i_end, 1:params &
           &%exp_data(xps_idx)%n_samples), params%do_forces, &
           & xyz(1:3, j_beg:j_end),&
         & energies_lp(i_beg:i_end), forces_lp, virial_lp, params%exp_similarity_type, rank)

         ! if (rank == 0)then
         !    open(unit=11, file="tg_xps.dat", status="unknown")
         !    do i = 1, params%exp_data(xps_idx)%n_samples
         !       write(11, '(1X,F20.8,1X,F20.8)') params%exp_data(xps_idx)%x(i), params%exp_data(xps_idx)%y_pred(i)
         !    end do
         !    close(11)
         ! end if

         call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

         if (rank == 0 .and. params%write_exp .and. write_condition) then

            call get_overwrite_condition(params%do_mc, params%do_md&
              &, mc_istep, md_istep, params%write_xyz, overwrite_condition)

            call write_exp_datan(params%exp_data(xps_idx)&
              &%x(1:params%exp_data(xps_idx)%n_samples), params&
              &%exp_data(xps_idx)%y_pred(1:params&
              &%exp_data(xps_idx)%n_samples),&
              & overwrite_condition, "xps_prediction.dat",&
              & params%exp_data(xps_idx)%label)

            if (.not. params%exp_data(xps_idx)%wrote_exp) then

               call preprocess_exp_data(params, params%exp_data(xps_idx)%x,&
                 & params%exp_data(xps_idx)%y, params%exp_data(xps_idx)%label,&
                 & n_sites, dot_product(cross_product(a_box,&
                 & b_box), c_box)/(dfloat(indices(1)*indices(2) &
                 &*indices(3))), params%exp_data(xps_idx)%input, exp_output, .true.)

               call write_exp_datan(params%exp_data(xps_idx)&
                 &%x(1:params%exp_data(xps_idx)%n_samples),&
                 & params%exp_data(xps_idx)%y(1:params&
                 &%exp_data(xps_idx)%n_samples),&
                 & overwrite_condition, "xps_exp.dat", params&
                 &%exp_data(xps_idx)%label)
               params%exp_data(xps_idx)%wrote_exp = .true.
            end if

            ! else
            !    call write_exp_data(params%exp_data(xps_idx)%x, params%exp_data(xps_idx)%y_pred, mc_istep == 0, "xps_prediction.dat" )
            !    call write_exp_data(params%exp_data(xps_idx)%x, params%exp_data(xps_idx)%y, mc_istep == 0, "xps_exp.dat" )

         end if

         !deallocate( params%exp_data(xps_idx)%y_pred )
         if (allocated(y_i_pred_all)) deallocate (y_i_pred_all)
         ! sim_exp_pred would be an energy if multiplied by some energy scale \gamma * ( 1 - sim )
         ! sim_exp_pred_der would be the array of forces if multiplied by (- \gamma )
         deallocate (v_neigh_lp)

         !time%xps(2) = MPI_wtime()
         call time_end(time%xps, "xps")
         !           if (rank == 0) print *, rank, " TIME_XPS = ", time%xps(3)

      else if (any(soap_turbo_hypers(:)%has_core_electron_be) .and. params%do_xps) then
         ! Get the linspace of the xps spectrum and then perform the
         ! calculation and write to the prediction file
         !
         if (rank == 0) then
            call get_xps_spectra_standalone(&
                 & params%xps_e_min,&
                 & params%xps_e_max, &
                 & params%xps_sigma, &
                 & params%xps_n_samples,&
                 & x_xps, &
                 & y_xps, &
                 & local_properties(1:n_sites, core_be_lp_index))

            call get_overwrite_condition(params%do_mc, params%do_md&
                 &, mc_istep, md_istep, params%write_xyz, overwrite_condition)

            if (n_xyz > 0) then
               overwrite_condition = (n_xyz == 1)
            end if

            call write_exp_datan(x_xps(1:params%xps_n_samples), &
                 & y_xps(1:params%xps_n_samples), &
                 & overwrite_condition, &
                 &"xps_prediction.dat",&
                 &"core_electron_be xps")
         end if

      end if

   end subroutine compute_exp_xps
!**************************************************************************

!**************************************************************************
! The (partial) pair distribution functions, structure factors and
! diffraction patterns, and the energies and forces that fitting them
! against experiment produces.
!
! Lifted from turbogap.f90 (lines 1979-2680 at the time of the move). The
! only change to the body is that four argument lists split across
! #ifdef _MPIF90 collapsed: the two branches differed only by the this_
! prefix on an energies/forces/virial triple, and in here the dummy has
! one name either way, so the caller chooses once at the single call.
!
! The device context this block shares with the electrostatics and SOAP
! paths -- streams, cuBLAS handles, gpu_exp, gpu_neigh, gpu_batch_storage
! -- comes from gpu_context by USE, which is what lets the body move
! unchanged instead of growing eight more arguments.
!**************************************************************************
   subroutine compute_exp_spectra(params, n_sites, species, rjs, xyz, neighbors_list, &
                                  n_neigh, neighbor_species, indices, a_box, b_box, c_box, i_beg, i_end, j_beg, j_end, &
                                  rank, ntasks, ierr, md_istep, mc_istep, energies_sf, forces_sf, virial_sf, &
                                  energies_xrd, forces_xrd, virial_xrd, energies_nd, forces_nd, virial_nd, time, &
                                  i_beg_list, i_end_list, j_beg_list, j_end_list, &
                                  n_omp, omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end, n_sites_temp, &
                                  n_pairs_temp, write_condition, overwrite_condition, temp_string, &
                                  species_types_actual, v_uc)

      implicit none

!   ---- dummy arguments ----
!   Read-only inputs.
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: rank
      integer, intent(in) :: ntasks
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: species(:)

!   params carries the do_* switches and the exp_data containers, which the
!   exp_interface routines update, so it cannot be intent(in).
      type(input_parameters), intent(inout) :: params
      integer, intent(inout) :: ierr

!   The three contribution families this routine produces. The caller passes
!   the this_-prefixed arrays under MPI and the plain ones otherwise; in here
!   there is one name either way, which is what removes four
!   preprocessor-interrupted argument lists from the driver.
      real(dp), allocatable, intent(inout) :: energies_sf(:)
      real(dp), allocatable, intent(inout) :: forces_sf(:, :)
      real(dp), allocatable, intent(inout) :: energies_xrd(:)
      real(dp), allocatable, intent(inout) :: forces_xrd(:, :)
      real(dp), allocatable, intent(inout) :: energies_nd(:)
      real(dp), allocatable, intent(inout) :: forces_nd(:, :)
      real(dp), intent(inout) :: virial_sf(1:3, 1:3)
      real(dp), intent(inout) :: virial_xrd(1:3, 1:3)
      real(dp), intent(inout) :: virial_nd(1:3, 1:3)

!   Timing buckets, accumulated across snapshots.
      type(times_t), intent(inout) :: time

!   The batch decomposition. Built by the caller for the electrostatics and
!   SOAP paths and consumed here too, which is why it is passed rather than
!   held in gpu_context.
      integer, allocatable, intent(inout) :: i_beg_list(:)
      integer, allocatable, intent(inout) :: i_end_list(:)
      integer, allocatable, intent(inout) :: j_beg_list(:)
      integer, allocatable, intent(inout) :: j_end_list(:)

!   Per-batch bounds and the OpenMP task index. omp_task is intended to be
!   thread-private, so it stays a variable the caller owns rather than moving
!   into gpu_context with the arrays it indexes.
      integer, intent(inout) :: n_omp
      integer, intent(inout) :: omp_task
      integer, intent(inout) :: this_i_beg
      integer, intent(inout) :: this_i_end
      integer, intent(inout) :: this_j_beg
      integer, intent(inout) :: this_j_end
      integer, intent(inout) :: n_sites_temp
      integer, intent(inout) :: n_pairs_temp

!   Output bookkeeping shared with the XPS path.
      logical, intent(inout) :: write_condition
      logical, intent(inout) :: overwrite_condition
      character*1024, intent(inout) :: temp_string
      character*8, allocatable, intent(inout) :: species_types_actual(:)
      real(dp), intent(inout) :: v_uc

!   ---- local variables ----
!   Block-private. All of these were driver declarations that nothing outside
!   this block referenced.
      real(dp), parameter :: pi = acos(-1.0)

!   i, j, k and f were driver scratch. Each is written before it is read here,
!   and the driver rewrites i and j (do-loops) before next reading them; k and f
!   it never reads again. f is only ever used inside this block.
      integer :: i
      integer :: j
      integer :: k
      real(dp) :: f
      integer :: n_dim_idx
      integer :: n_dim_partial
      integer :: n_species_actual
      integer :: q_beg
      integer :: q_end

      real(dp), allocatable, target :: dV(:)
      real(dp), allocatable, target :: n_atoms_of_species(:)
      real(dp), allocatable, target :: prefactor(:)
      real(dp), allocatable, target :: sinc_factor_matrix(:, :)
      real(dp), allocatable, target :: pair_distribution_der(:, :)
      real(dp), allocatable, target :: pair_distribution_partial(:, :)
      real(dp), allocatable, target :: pair_distribution_partial_temp(:, :)
      real(dp), allocatable, target :: pair_distribution_partial_der(:, :, :)
      real(dp), allocatable, target :: pair_distribution_partial_temp_der(:, :, :)
      real(dp), allocatable, target :: structure_factor_partial(:, :)
      real(dp), allocatable, target :: structure_factor_partial_temp(:, :)
      real(dp), allocatable, target :: x_pair_distribution(:)
      real(dp), allocatable, target :: y_pair_distribution(:)
      real(dp), allocatable, target :: y_pair_distribution_temp(:)
      real(dp), allocatable, target :: x_structure_factor(:)
      real(dp), allocatable, target :: x_structure_factor_temp(:)
      real(dp), allocatable, target :: y_structure_factor(:)
      real(dp), allocatable, target :: y_structure_factor_temp(:)
      real(dp), allocatable, target :: x_xrd(:)
      real(dp), allocatable, target :: x_xrd_temp(:)
      real(dp), allocatable, target :: y_xrd(:)
      real(dp), allocatable, target :: y_xrd_temp(:)
      real(dp), allocatable, target :: x_nd(:)
      real(dp), allocatable, target :: x_nd_temp(:)
      real(dp), allocatable, target :: y_nd(:)
      real(dp), allocatable, target :: y_nd_temp(:)

      integer, allocatable :: nk(:)

!   Device-side scratch, private to this block.
      type(c_ptr) :: xpdf_d
      type(c_ptr) :: dV_d
      type(c_ptr) :: prefactor_d
      type(c_ptr) :: sinc_factor_matrix_d
      type(c_ptr), allocatable :: nk_d(:)
      type(c_ptr), allocatable :: k_index_d(:)
      type(c_ptr), allocatable :: j2_index_d(:)
      type(c_ptr), allocatable :: xyz_k_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_der_d(:)
      type(c_ptr), allocatable :: all_scattering_factors_d(:)
      integer(c_size_t) :: st_x_d
      integer(c_size_t) :: st_sinc_factor_matrix_d
      integer(c_size_t) :: st_prefactor_d
      integer(c_size_t), allocatable :: st_nk_d(:)
      integer(c_size_t), allocatable :: st_k_index_d(:)
      integer(c_size_t), allocatable :: st_j2_index_d(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_d(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_der_d(:)

      type(gpu_host_storage_type), allocatable :: gpu_host_exp_storage(:)

      !##############################################################!
      !###---   (Partial) Pair distribution functions and XRD   ---###!
      !##############################################################!

      ! We use these to calculate the (partial) structure factors, which
      ! can be used for X-Ray scattering and (in the future)
      ! neutron scattering.
      !
      ! > We use the formalism which was detailed by
      !   Gutierrez and Johansson, Physical Review B, Volume 65, 104202 (2002)
      !   such that there is consistency between the rdfs and the scattering we calculate.
      !
      ! > Furthermore, calculating the pair distribution function and
      !   the structure factor/XRD becomes much faster, as there is just
      !   a sum over species rather than a double sum over all the
      !   atomic species.
      !
      !   (The ASE implementation of XRD intensity is problematic and
      !   uses the sinc function implemented by DSP
      !   (sin(pi*x)/(pi*x)) which is not what is in the
      !   literature).
      !
      ! ***--- Steps for calculation ---***
      ! 1) We calculate the partial pair distribution functions g_ab
      ! 2) Partial static structure factors S_ab are then calculated from this by Fourier transform
      !    > This calculation can includes a window function ( sin(pi*rij/r_cut) / (pi*rij/r_cut) )
      !      such that termination effects of large sinusoids which come from the cutoff are removed.
      ! 3) If X-Ray Diffraction (xrd) is specified, then the intensity is
      !    calculated from these partial structure factors by the
      !    inclusion of the X-Ray form factors.
      !
      ! ***--- Definitions ---***
      ! There are many definitions of these various functions
      ! in the literature, however, we shall use similar ones
      ! to those in the paper of Gutierrez
      !
      ! Definiton of the structure factors given by Ashcroft and Langreth
      ! N. W. Ashcroft and D. C. Langreth. Phys. Rev., 156(3):685-692 (1967)
      !
      ! PDFs can be smoothed by using a kernel density estimate by a gaussian function when kde_sigma > 0.d0
      !
      ! R(r)    = Radial Distribution Function     === A histogram of atomic distances divided by N, goes as r^2
      ! g(r)    = Pair Distribution Function (PCF) === Scales R(r) by 1/(4 pi r^2) such that it lays flat, converges to 1.
      !         = (N_{r_l < r < r_h} / N) / ( 4 pi r^2 dr  * ( N / V ) )
      !         = n_{r_l < r < r_h} / ( dV * rho_0 )
      !           > n_{r_l < r < r_h} is the average number of particles between r_l and r_h
      !           > rho_0 is the density
      !           > dV is the differential volume between shells
      !
      ! g_ab(r) = Partial Pair Distribution Func.  === Same as above but only for particles a and b
      !         =  (N^{b}_{r_l < r < r_h} / N_b ) / ( 4 pi r^2 dr ) * ( N_a / V )
      !         With kde
      !         =  ( sum_i sum_i/=j exp( - (r - r_ij)^2 / sigma^2 / 2 ) / N_b ) / ( 4 pi r^2 dr ) * ( V / N_a ) &
      !              & * ( (r_max - r_min) / sigma / sqrt(2pi) )
      !
      !           The full pair distribution function is given by the sum (say for the binary system, with a, b atoms)
      !             g(r) = (N_a/N) * g_aa(r) + 2(N_a/N * N_b/N)g_ab + (N_b/N)g_bb
      !
      ! S_ab(q) = delta_ab + 4 pi rho (ca cb)^1/2 int_0^r_cut dr r^2 [ g_ab(r) - 1 ] sin(qr)/(qr) * sin( pi r / R )/ (pi r /R)
      !
      !
      ! XRD(q) = 1/N ( d cross_section/ d Omega )
      !
      ! Total scattering function F^X(q)
      ! F^x(q) = [ XRD(q) - \sum_n c_i f_i(q)^2 ] / [ \sum_n c_i f_i(q) ]^2

      ! First get the number of species in actuality
      n_species_actual = 0
      do i = 1, n_sites
         if (species(i) > n_species_actual) n_species_actual = n_species_actual + 1
      end do

      ! Now find the unique species ids
      if (allocated(species_types_actual)) deallocate (species_types_actual)
      allocate (species_types_actual(1:n_species_actual))

      n_species_actual = 0
      do i = 1, n_sites
         if (species(i) > n_species_actual) then
            n_species_actual = n_species_actual + 1
            species_types_actual(n_species_actual) = params%species_types(species(i))
         end if
      end do

      n_dim_partial = n_species_actual*(n_species_actual + 1)/2

      if (params%gpu_batched .and. (params%do_xrd .or. params%do_nd) &
          .and. params%exp_forces .and. params%do_forces) then

         !           print *, "> Starting batched xrd "
         !           call cpu_time( time%exp_batched(1) )
         call time_start(time%exp_batched, "exp_batched")

!        Size the batch count from the device, unless the input opted out.
!
!        The estimator was written for exactly this and its call site has been
!        commented out ever since -- so gpu_n_batches has always been whatever
!        the input said, with nothing anywhere comparing it to the card. It
!        enumerates the pdf/xrd device buffers one by one; the dominant term is
!        Gk = n_samples x n_pairs x 8 bytes, which is why a large cutoff blows
!        up so much faster than a large cell.
!
!        gpu_batches_for_gb only ever raises the count, and does nothing at all
!        when gpu_mem_fraction is unset, so an input that worked yesterday
!        batches identically today.
         call estimate_max_exp_forces_device_memory_usage(i_end - i_beg + 1, j_end - j_beg + 1, n_dim_partial, &
                                                          params%pair_distribution_n_samples, &
                                                          params%structure_factor_n_samples, gpu_memory_usage, &
                                                          be_verbose=(rank == 0 .and. params%gpu_mem_fraction > 0.d0))

         params%gpu_n_batches = gpu_batches_for_gb(gpu_memory_usage, params%gpu_n_batches, &
                                                   i_end - i_beg + 1, rank, "pdf/xrd batched forces")

         call get_gpu_batches(n_neigh(i_beg:i_end), rjs(j_beg:j_end), params%pair_distribution_rcut, &
                              params%gpu_n_batches, gpu_memory_usage, &
                              params%gpu_max_batch_size, i_beg_list, i_end_list, j_beg_list, j_end_list)

         gpu_memory_usage = 0.d0

         ! Now we want to use open mp to get the partial pdfs for each batch and then collect them
         ! Eventually we will use OMP parallel do to do the loop and submit to different gpu streams

         ! allocate the pdfs
         ! call gpu_check_error()
         call setup_batched_pair_distribution(params%r_range_min, params%r_range_max, &
                                         params%pair_distribution_rcut, params%pair_distribution_n_samples, x_pair_distribution, dV)

         call get_n_atoms_of_species(n_atoms_of_species, n_sites, species, n_species_actual)

         ! call gpu_check_error()
         allocate (gpu_neigh(1:size(i_beg_list)))
         allocate (gpu_exp(1:size(i_beg_list)))
         allocate (gpu_batch_storage(1:size(i_beg_list)))

         st_x_d = params%pair_distribution_n_samples*c_double

         call gpu_malloc_async(xpdf_d, st_x_d, gpu_stream)
         call cpy_htod(c_loc(x_pair_distribution), xpdf_d, st_x_d, gpu_stream)

         call gpu_malloc_async(dV_d, st_x_d, gpu_stream)
         call cpy_htod(c_loc(dV), dV_d, st_x_d, gpu_stream)

         ! > We have n_omp streams created. In general the number
         !   of batches will likely be greater than the number of
         !   omp tasks, but include num_threads(n_omp) to make sure we just use the right number of threads.

         ! ! $OMP SHARED(cublas_handles, gpu_streams, n_neigh, species, neighbor_species, neighbors_list, rjs, xyz, &
         ! ! $OMP xpdf_d, dV_d, dV, n_atoms_of_species, n_species_actual, n_sites, v_uc, gpu_batch_storage, gpu_exp)

         !$OMP PARALLEL DO num_threads(n_omp) DEFAULT(SHARED) &
         !$OMP PRIVATE(i,  omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end, n_sites_temp, n_pairs_temp)

         do i = 1, size(i_beg_list)

            ! > In sequential operation, we just want to use the one stream, and only 1 will be created
            omp_task = gpu_omp_task()

            this_i_beg = i_beg - 1 + i_beg_list(i)
            this_i_end = i_beg - 1 + i_end_list(i)
            this_j_beg = j_beg - 1 + j_beg_list(i)
            this_j_end = j_beg - 1 + j_end_list(i)

            n_sites_temp = this_i_end - this_i_beg + 1
            n_pairs_temp = this_j_end - this_j_beg + 1

            call gpu_malloc_neighbors(gpu_neigh(i), &
                                      n_sites_temp, n_pairs_temp, &
                                      n_neigh(this_i_beg:this_i_end), &
                                      species(this_i_beg:this_i_end), &
                                      neighbor_species(this_j_beg:this_j_end), &
                                      neighbors_list(this_j_beg:this_j_end), &
                                      rjs(this_j_beg:this_j_end), &
                                      xyz(1:3, this_j_beg:this_j_end), &
                                      gpu_streams(omp_task), rank)

            print *, " Total neighbors   = ", j_end
            print *, "   n_batch         = ", i
            print *, "   this neighbors  = ", this_j_end - this_j_beg + 1
            print *, "   Allocating ->   = ", (this_j_end - this_j_beg + 1)*8*5, " bytes"
            print *, "              ->   = ", &
               dfloat((this_j_end - this_j_beg + 1)*8*5)/1024.d0/1024.d0/1024.d0, " Gbytes"

            call total_gpu_memory(dfloat((this_j_end - this_j_beg + 1)*8*5))

            ! call gpu_check_error()

            ! Now as we don't need anything else we can just try and
            ! calculate the electrostatics of this batch directly,
            ! with the forces too

            call calculate_batched_pair_distribution( &
               gpu_exp(i), &
               gpu_batch_storage(i), &
               gpu_neigh(i), &
               x_pair_distribution, &
               dV, &
               n_atoms_of_species, &
               n_species_actual, &
               n_sites, &
               this_i_beg, &
               this_i_end, &
               this_j_beg, &
               this_j_end, &
               params%pair_distribution_n_samples, &
               params%r_range_min, &
               params%r_range_max, &
               params%pair_distribution_rcut, &
               params%pair_distribution_kde_sigma, &
               gpu_streams(omp_task), &
               xpdf_d, &
               dV_d, &
               v_uc, &
               rank)

            ! call gpu_check_error()
            call gpu_free_neighbors(gpu_neigh(i), gpu_streams(omp_task))

            ! call gpu_check_error()
            write (*, '(A,I4,A,I4,A,I4,A,I4,A,I4,A,I4)') "pdf batches finished---Rank ", rank, " ---Thread ", omp_task, &
               " / ", n_omp, " i = ", i, &
               " i_beg = ", this_i_beg, &
               " i_end = ", this_i_end, &
               " j_beg = ", this_j_beg, &
               " j_end = ", this_j_end
            call flush (101)

!           This batch's work is on gpu_streams(omp_task), not on the default
!           stream. Syncing gpu_stream here waited for a stream this iteration
!           never touched -- so it neither ordered anything this loop needs nor,
!           once threads are real, kept a thread from racing ahead of its own
!           batch. The batch results are collected after the loop, and the
!           barrier at END PARALLEL DO is what orders that; this sync only has
!           to make the batch's own device work complete before its host-side
!           bookkeeping below.
            call gpu_stream_sync(gpu_streams(omp_task))

            call gpu_meminfo()
         end do
         !$OMP END PARALLEL DO
         deallocate (gpu_neigh)

         ! Collect the pdf
         !call gpu_
         call collect_batched_pair_distribution(size(i_beg_list), gpu_batch_storage, &
                                                n_dim_partial, params%pair_distribution_n_samples, &
                                                pair_distribution_partial, n_species_actual, n_atoms_of_species, v_uc)

         ! call gpu_check_error()

         ! Write out the partial pair distribution functions
         call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

         if (rank == 0 .and. params%write_pair_distribution .and. write_condition) then
            call get_overwrite_condition(params%do_mc, params%do_md,&
              & mc_istep, md_istep, params%write_xyz,&
              & overwrite_condition)

            n_dim_idx = 1
            outerchk: do j = 1, n_species_actual
               do k = 1, n_species_actual
                  if (j > k) cycle

                  write (temp_string, '(A)')&
                    & 'pair_distribution_'//trim(species_types_actual(j))//'_'//trim(species_types_actual(k))//&
                    & "_prediction.dat"
                  call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
                    & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx),&
                    & overwrite_condition, temp_string, 'pair_distribution')

                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) exit outerchk
               end do
            end do outerchk
         end if

         ! SETTING EXP FORCES TO FALSE
         params%do_forces = .false.
         params%exp_forces = .false.

         if (params%do_structure_factor) then

            !time%sf(1) = MPI_wtime()
            call time_start(time%sf, "structure_factor")

            call calculate_structure_factor(params, x_structure_factor, x_structure_factor_temp,&
              & y_structure_factor, y_structure_factor_temp,&
              & structure_factor_partial, structure_factor_partial_temp,&
              & x_pair_distribution, y_pair_distribution, &
              & pair_distribution_partial, n_species_actual, species_types_actual, n_atoms_of_species,&
              & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
              & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
              & neighbor_species, species, rank, q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
              & nk, nk_d, k_index_d, j2_index_d, xyz_k_d,&
              & pair_distribution_partial_d,&
              & pair_distribution_partial_der_d, st_nk_d,&
              & st_k_index_d, st_j2_index_d,&
              & st_pair_distribution_partial_d,&
              & st_pair_distribution_partial_der_d,&
            & pair_distribution_partial_der, energies_sf, forces_sf&
              &, virial_sf, params%structure_factor_matrix_forces&
              &, cublas_handle, gpu_stream, &
              & gpu_host_exp_storage, params%gpu_low_memory)

            !time%sf(2) = MPI_wtime()
            call time_end(time%sf, "structure_factor")

         end if

         if (params%do_xrd) then

            !time%xrd(1) = MPI_wtime()
            call time_start(time%xrd, "xrd")

            call calculate_xrd(params, x_xrd, x_xrd_temp, y_xrd,&
              & y_xrd_temp, x_structure_factor,&
              & x_structure_factor_temp,&
              & structure_factor_partial,&
              & structure_factor_partial_temp, n_species_actual,&
              & species_types_actual, n_atoms_of_species,&
              & n_sites, a_box, b_box, c_box, indices, md_istep,&
              & mc_istep, i_beg, i_end, j_beg, j_end, ierr, rjs,&
              & xyz, neighbors_list, n_neigh, neighbor_species,&
              & species, rank, q_beg, q_end, ntasks,&
              & sinc_factor_matrix, params%exp_forces, nk, nk_d,&
              & k_index_d, j2_index_d, xyz_k_d,&
              & pair_distribution_partial_d,&
              & pair_distribution_partial_der_d, st_nk_d,&
              & st_k_index_d, st_j2_index_d,&
              & st_pair_distribution_partial_d,&
              & st_pair_distribution_partial_der_d,&
            & pair_distribution_partial_der, energies_xrd,&
              & forces_xrd, virial_xrd, .false., params&
              &%structure_factor_matrix_forces, cublas_handle,&
              & gpu_stream, gpu_host_exp_storage, params&
              &%gpu_low_memory)

            !time%xrd(2) = MPI_wtime()
            call time_end(time%xrd, "xrd")

         end if

         if (params%do_nd) then

            !time%nd(1) = MPI_wtime()
            call time_start(time%nd, "nd")

            call calculate_xrd(params, x_nd, x_nd_temp, y_nd,&
              & y_nd_temp, x_structure_factor,&
              & x_structure_factor_temp,&
              & structure_factor_partial,&
              & structure_factor_partial_temp, n_species_actual,&
              & species_types_actual, n_atoms_of_species,&
              & n_sites, a_box, b_box, c_box, indices, md_istep,&
              & mc_istep, i_beg, i_end, j_beg, j_end, ierr, rjs,&
              & xyz, neighbors_list, n_neigh, neighbor_species,&
              & species, rank, q_beg, q_end, ntasks,&
              & sinc_factor_matrix, params%exp_forces, nk, nk_d,&
              & k_index_d, j2_index_d, xyz_k_d,&
              & pair_distribution_partial_d,&
              & pair_distribution_partial_der_d, st_nk_d,&
              & st_k_index_d, st_j2_index_d,&
              & st_pair_distribution_partial_d,&
              & st_pair_distribution_partial_der_d,&
            & pair_distribution_partial_der, energies_nd, forces_nd&
              &, virial_nd, .true., params &
              &%structure_factor_matrix_forces, cublas_handle,&
              & gpu_stream, gpu_host_exp_storage, params&
              &%gpu_low_memory)

            !time%nd(2) = MPI_wtime()
            call time_end(time%nd, "nd")
         end if

         call time_start(time%xrd, "xrd_batched")

         ! RESETTING EXP FORCES TO FALSE
         params%do_forces = .true.
         params%exp_forces = .true.

         ! Get the scattering factors

         ! call get_all_scattering_factors(all_scattering_factors,&
         !      & n_sites, params%exp_data(params%xrd_idx)%x,&
         !      & species_types_actual, params&
         !      &%structure_factor_n_samples, n_species_actual,&
         !      & x_structure_factor, j,k, params%do_xrd, params&
         !      &%xrd_output, n_atoms_of_species, .false.)

         ! print *, " alloc y_xrd ", allocated(y_xrd), size( y_xrd )
         ! print *, " alloc y_exp ", allocated(params%exp_data(params%xrd_idx)%y), size(params%exp_data(params%xrd_idx)%y)

         allocate (prefactor(1:size(y_xrd)))
         prefactor(1:size(y_xrd)) = (y_xrd(1:size(y_xrd)) - params%exp_data(params%xrd_idx)%y(1:size(y_xrd)))

         st_prefactor_d = size(prefactor, 1)*c_double
         call gpu_malloc_async(prefactor_d, st_prefactor_d, gpu_stream)
         call cpy_htod(c_loc(prefactor), prefactor_d, st_prefactor_d, gpu_stream)

         ! call gpu_check_error()
         ! st_all_scattering_factors_d = size( all_scattering_factors, 1) * c_double
         ! call gpu_malloc_async(all_scattering_factors_d, st_all_scattering_factors_d, gpu_stream)
         ! call cpy_htod(c_loc( all_scattering_factors ), all_scattering_factors_d, st_all_scattering_factors_d, gpu_stream)

         st_sinc_factor_matrix_d = size(sinc_factor_matrix, 1)*size(sinc_factor_matrix, 2)*c_double
         call gpu_malloc_async(sinc_factor_matrix_d, st_sinc_factor_matrix_d, gpu_stream)
         call cpy_htod(c_loc(sinc_factor_matrix), sinc_factor_matrix_d, st_sinc_factor_matrix_d, gpu_stream)

         allocate (all_scattering_factors_d(1:n_dim_partial))

         n_dim_idx = 1
         outerprep: do j = 1, n_species_actual
            do k = 1, n_species_actual
               if (j > k) cycle

               call get_all_scattering_factors(all_scattering_factors_d(n_dim_idx),&
                 & n_sites, x_xrd(1:params%structure_factor_n_samples), species_types_actual, &
                 params%structure_factor_n_samples, n_species_actual,&
                 & x_xrd(1:params%structure_factor_n_samples), j, k, params%do_xrd, params%xrd_output,&
                 & n_atoms_of_species, .false., gpu_stream)

               ! call gpu_check_error()
               n_dim_idx = n_dim_idx + 1
               if (n_dim_idx > n_dim_partial) exit outerprep
            end do
         end do outerprep

         call get_energy_scale(params%do_md, params%do_mc,&
           & md_istep, params%md_nsteps, mc_istep, params &
           &%mc_nsteps, params&
           &%exp_energy_scales_initial(params%xrd_idx), params&
           &%exp_energy_scales_final(params%xrd_idx), params&
           &%exp_energy_scales(params%xrd_idx))

         ! ! $OMP SHARED(cublas_handles, gpu_streams, n_neigh, species, neighbor_species, neighbors_list, rjs, xyz, &
         ! ! $OMP xpdf_d, dV_d, dV, n_atoms_of_species, n_species_actual, n_sites, v_uc, gpu_batch_storage, gpu_exp, x_pair_distribution)

         !$OMP PARALLEL DO num_threads(n_omp) DEFAULT(SHARED) &
         !$OMP PRIVATE(i, omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end) &
         !$OMP PRIVATE(n_sites_temp, n_pairs_temp, n_dim_idx, j, k, f)

         do i = 1, size(i_beg_list)

            omp_task = gpu_omp_task()

            this_i_beg = i_beg - 1 + i_beg_list(i)
            this_i_end = i_beg - 1 + i_end_list(i)
            this_j_beg = j_beg - 1 + j_beg_list(i)
            this_j_end = j_beg - 1 + j_end_list(i)

            n_pairs_temp = this_j_end - this_j_beg + 1

            ! I have the neighbors allocated in in gpu_neigh(i)
            print *, " "
!           I0, not I8, and eight descriptors for eight integers. The format had
!           six, so it wrapped onto a second line, and the pair indices overflowed
!           it into '********' -- a million atoms at a 6 A cutoff gives ~168M
!           pairs, which is nine digits.
            write (*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') "xrd batches---Rank ", rank, " ---Thread ", omp_task, &
               " / ", n_omp, " i ", i, &
               " this_i_beg = ", this_i_beg, &
               " this_i_end = ", this_i_end, &
               " this_j_beg = ", this_j_beg, &
               " this_j_end = ", this_j_end

            n_dim_idx = 1
            outer: do j = 1, n_species_actual
               do k = 1, n_species_actual
                  if (j > k) cycle

                  ! Now calculate the pair_distribution_partial_der_d
                  ! we have the rjs stored

                  call calculate_batched_pair_distribution_der(gpu_exp(i), gpu_batch_storage(i), &
                                                           x_pair_distribution, dV, n_atoms_of_species, n_species_actual, n_sites, &
                                                               this_i_beg, this_i_end, this_j_beg, this_j_end, &
                                                       params%pair_distribution_n_samples, params%r_range_min, params%r_range_max, &
                                                               params%pair_distribution_rcut, params%pair_distribution_kde_sigma, &
                                                               gpu_streams(omp_task), xpdf_d, dV_d, j, k, n_dim_idx, v_uc)

                  ! call gpu_check_error()

                  !                    print * , "starting xrd setup"
                  call setup_gpu_xrd_forces(gpu_exp(i), gpu_batch_storage(i), n_dim_idx, gpu_streams(omp_task))

                  !                    print * , "finished xrd setup"

                  if (j == k) f = 1.d0
                  if (j /= k) f = 2.d0

                  f = 4.d0*pi*f*((n_atoms_of_species(j)*n_atoms_of_species(k))/dfloat(n_sites)/&
                    & dfloat(n_sites))*(dfloat(n_sites)/v_uc)

                  allocate (gpu_batch_storage(i)%host(n_dim_idx)%forces_h(1:3, 1:n_sites))
                  gpu_batch_storage(i)%host(n_dim_idx)%forces_h = 0.d0
                  gpu_batch_storage(i)%host(n_dim_idx)%virial_h = 0.d0

                  ! call gpu_check_error()
                  !                    print * , "starting xrd forces"

                  call get_structure_factor_forces_matrix_original(.true., params%exp_energy_scales(params%xrd_idx), &
                                                                   gpu_batch_storage(i)%host(n_dim_idx)%forces_h, &
                                                                   gpu_batch_storage(i)%host(n_dim_idx)%virial_h, &
                                                            params%pair_distribution_n_samples, params%structure_factor_n_samples, &
                                                                   f, &
                                                                   gpu_exp(i)%nk(n_dim_idx), &
                                                                   gpu_exp(i)%k_index_d(n_dim_idx), &
                                                                   gpu_exp(i)%j2_index_d(n_dim_idx), &
                                                                   gpu_exp(i)%xyz_k_d(n_dim_idx), &
                                                                   gpu_exp(i)%pair_distribution_partial_der_d(n_dim_idx), &
                                                                   sinc_factor_matrix_d, &
                                  all_scattering_factors_d(n_dim_idx), prefactor_d, gpu_streams(omp_task), cublas_handles(omp_task))

                  ! print *, "virial after function exit ", gpu_batch_storage(i) % host( n_dim_idx ) % virial_h
                  call free_gpu_xrd_forces(gpu_exp(i), gpu_batch_storage(i), n_dim_idx, gpu_streams(omp_task))

                  ! call gpu_check_error()
                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) exit outer

               end do
            end do outer
         end do
         !$OMP END PARALLEL DO

         call gpu_free_async(xpdf_d, gpu_stream)
         call gpu_free_async(dV_d, gpu_stream)

         n_dim_idx = 1
         outer_post: do j = 1, n_species_actual
            do k = 1, n_species_actual
               if (j > k) cycle
               call gpu_free_async(all_scattering_factors_d(n_dim_idx), gpu_stream)
               n_dim_idx = n_dim_idx + 1
               if (n_dim_idx > n_dim_partial) exit outer_post
            end do
         end do outer_post
         deallocate (all_scattering_factors_d)

         call gpu_free_async(sinc_factor_matrix_d, gpu_stream)
         call gpu_free_async(prefactor_d, gpu_stream)

         deallocate (prefactor)

         ! call gpu_check_error()
         call gpu_stream_sync(gpu_stream)

         allocate (forces_xrd(1:3, 1:n_sites))
         forces_xrd = 0.d0
         virial_xrd = 0.d0

         call collect_batched_forces(size(i_beg_list), gpu_batch_storage, n_dim_partial, &
                                     forces_xrd, virial_xrd, n_sites)

         call free_host_batches(gpu_batch_storage, params%gpu_n_batches, n_dim_partial)
         call free_exp_batches(gpu_exp, params%gpu_n_batches)

         if (allocated(i_beg_list)) deallocate (i_beg_list, i_end_list, j_beg_list, j_end_list)

         if (allocated(gpu_neigh)) deallocate (gpu_neigh)
         if (allocated(gpu_exp)) deallocate (gpu_exp)
         if (allocated(gpu_batch_storage)) deallocate (gpu_batch_storage)

         call time_end(time%xrd, "xrd_batched")

!        Was a hand-written copy of time_end -- get_time into (2), then the
!        accumulation into (3) -- with the xrd close in between. time_end is
!        that, and only it closes the NVTX range opened at the top of the
!        block, so the hand-written form would leave the range open and nest
!        every later phase inside it.
         call time_end(time%exp_batched, "exp_batched")

         if (rank == 0) print *, " "
         if (rank == 0) print *, "--Rank, ", rank, " --- TIME EXP BATCHED = ", time%exp_batched(3)
         if (rank == 0) print *, " "
      end if

      !################################################################!
      !###---   Compute similarity of experimental predictions   ---###!
      !################################################################!

      if (params%do_exp) then
         do i = 1, params%n_exp
            ! First normalize the spectrum if it matches some type of experimental data

            ! Allocate the prediction data
            ! if (.not. allocated( params%exp_data(i)%y_pred ) )then
            !    allocate( params%exp_data(i)%y_pred( 1:size(params%exp_data(i)%y, 1) ))
            ! end if

            if (trim(params%exp_data(i)%label) == 'pair_distribution') then
               params%exp_data(i)%y_pred = y_pair_distribution
            elseif (trim(params%exp_data(i)%label) == 'structure_factor') then
               params%exp_data(i)%y_pred = y_structure_factor
            elseif (trim(params%exp_data(i)%label) == 'xrd') then
               params%exp_data(i)%y_pred = y_xrd
            elseif (trim(params%exp_data(i)%label) == 'nd') then
               params%exp_data(i)%y_pred = y_nd
            end if

            if (params%exp_data(i)%compute_similarity .and. allocated(params%exp_data(i)%y)) then
               call get_data_similarity(params%exp_data(i)%y, params&
                 &%exp_data(i)%y_pred, params&
                 &%exp_data(i)%similarity, params&
                 &%exp_similarity_type)
            end if

            ! deallocate(params%exp_data(i)%x)
            ! deallocate(params%exp_data(i)%y)
            ! deallocate(params%exp_data(i)%y_pred)

         end do
      end if

      !##############################################!
      !###---   Finalize experimental arrays   ---###!
      !##############################################!

      if (params%do_pair_distribution) then
         call finalize_pair_distribution(params, x_pair_distribution&
           &, y_pair_distribution, y_pair_distribution_temp,&
           & pair_distribution_partial, pair_distribution_partial_temp, params&
           &%exp_forces, pair_distribution_der,&
           & pair_distribution_partial_der,&
           & pair_distribution_partial_temp_der, n_atoms_of_species, rank)

      end if

      ! Now calculate the structure factors
      if (params%do_structure_factor) then

         call finalize_structure_factor(params, x_structure_factor, x_structure_factor_temp,&
           & y_structure_factor, y_structure_factor_temp,&
           & structure_factor_partial, structure_factor_partial_temp,&
           & x_pair_distribution, y_pair_distribution, &
           & pair_distribution_partial, sinc_factor_matrix)

      end if

      if (params%do_xrd) then
         call finalize_xrd(params, x_xrd, x_xrd_temp,&
           & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
           & structure_factor_partial, structure_factor_partial_temp)
      end if

      if (params%do_nd) then
         call finalize_xrd(params, x_nd, x_nd_temp,&
           & y_nd, y_nd_temp, x_structure_factor, x_structure_factor_temp,&
           & structure_factor_partial, structure_factor_partial_temp)
      end if

      deallocate (species_types_actual)

   end subroutine compute_exp_spectra
!**************************************************************************

end module turbogap_exp
