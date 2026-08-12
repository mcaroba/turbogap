! Experimental scattering observables: the pair distribution function, the
! partial structure factors, and the X-ray and neutron diffraction patterns
! derived from them.
!
! Lifted verbatim out of turbogap.f90 (lines 1623-1894 at the time of the
! move). The body is unchanged apart from one thing: four call statements had
! their argument list split across #ifdef _MPIF90, with this_-prefixed arrays
! on one side and un-prefixed ones on the other. Inside a procedure the dummy
! has one name either way, so those conditionals collapse and the caller makes
! the choice once. That is what makes this file parseable by ordinary Fortran
! tooling -- an argument list interrupted by the preprocessor is not.
!
! 28 arrays that were driver state are now locals here: every x_/y_/partial
! array for pdf, sf, xrd and nd, plus sinc_factor_matrix, n_atoms_of_species,
! species_types_actual and the q_beg/q_end range that calculate_structure_factor
! returns.

module turbogap_exp

   use kinds

   use types
   use exp_utils
   use exp_interface
   use timing
#ifdef _MPIF90
   use mpi
   use mpi_helper
#endif

   implicit none

contains

!**************************************************************************
   subroutine compute_exp_spectra(params, n_sites, species, positions, rjs, xyz, neighbors_list, &
                                  n_neigh, neighbor_species, indices, a_box, b_box, c_box, &
                                  i_beg, i_end, j_beg, j_end, rank, ntasks, ierr, md_istep, mc_istep, &
                                  energies_pdf, forces_pdf, virial_pdf, &
                                  energies_sf, forces_sf, virial_sf, &
                                  energies_xrd, forces_xrd, virial_xrd, &
                                  energies_nd, forces_nd, virial_nd, &
                                  time)

      implicit none

!   Passed through to the exp_interface routines, which declare the arrays
!   allocatable, so these must be allocatable too.
      type(input_parameters), intent(inout) :: params
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
      integer, intent(inout) :: ierr
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
!   The Debye route sums over every pair in the cell rather than walking the
!   neighbour list, so it needs the positions themselves. Every rank holds the
!   whole array, so no communication is implied by using it here.
      real(dp), intent(in), allocatable :: positions(:, :)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: species(:)

!   The four contribution families this routine produces. The caller passes the
!   this_-prefixed arrays under MPI and the plain ones otherwise.
      real(dp), allocatable, intent(inout) :: energies_pdf(:)
      real(dp), allocatable, intent(inout) :: forces_pdf(:, :)
      real(dp), allocatable, intent(inout) :: energies_sf(:)
      real(dp), allocatable, intent(inout) :: forces_sf(:, :)
      real(dp), allocatable, intent(inout) :: energies_xrd(:)
      real(dp), allocatable, intent(inout) :: forces_xrd(:, :)
      real(dp), allocatable, intent(inout) :: energies_nd(:)
      real(dp), allocatable, intent(inout) :: forces_nd(:, :)
      real(dp), intent(inout) :: virial_pdf(1:3, 1:3)
      real(dp), intent(inout) :: virial_sf(1:3, 1:3)
      real(dp), intent(inout) :: virial_xrd(1:3, 1:3)
      real(dp), intent(inout) :: virial_nd(1:3, 1:3)
      type(times_t), intent(inout) :: time

!   Was driver state; block-local now.
      real(dp), allocatable :: x_pair_distribution(:)
      real(dp), allocatable :: y_pair_distribution(:)
      real(dp), allocatable :: y_pair_distribution_temp(:)
      real(dp), allocatable :: pair_distribution_partial(:, :)
      real(dp), allocatable :: pair_distribution_partial_temp(:, :)
      real(dp), allocatable :: pair_distribution_der(:, :)
      real(dp), allocatable :: pair_distribution_partial_der(:, :, :)
      real(dp), allocatable :: pair_distribution_partial_temp_der(:, :, :)
      real(dp), allocatable :: n_atoms_of_species(:)
      real(dp), allocatable :: sinc_factor_matrix(:, :)
      real(dp), allocatable :: x_structure_factor(:)
      real(dp), allocatable :: x_structure_factor_temp(:)
      real(dp), allocatable :: y_structure_factor(:)
      real(dp), allocatable :: y_structure_factor_temp(:)
      real(dp), allocatable :: structure_factor_partial(:, :)
      real(dp), allocatable :: structure_factor_partial_temp(:, :)
      real(dp), allocatable :: x_xrd(:)
      real(dp), allocatable :: x_xrd_temp(:)
      real(dp), allocatable :: y_xrd(:)
      real(dp), allocatable :: y_xrd_temp(:)
      real(dp), allocatable :: x_nd(:)
      real(dp), allocatable :: x_nd_temp(:)
      real(dp), allocatable :: y_nd(:)
      real(dp), allocatable :: y_nd_temp(:)
      character*8, allocatable :: species_types_actual(:)
      integer :: i
      integer :: n_species_actual
      integer :: q_beg
      integer :: q_end

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

      if (params%do_pair_distribution) then
         call time_start(time%pdf)

         call calculate_pair_distribution(params, x_pair_distribution&
              &, y_pair_distribution, y_pair_distribution_temp,&
              & pair_distribution_partial, pair_distribution_partial_temp, &
              & n_species_actual, species_types_actual, n_atoms_of_species, n_sites, a_box,&
              & b_box, c_box, indices, md_istep, mc_istep, i_beg, i_end,&
              & j_beg, j_end, ierr, rjs, xyz, neighbors_list,&
              & n_neigh, neighbor_species, species, rank, params%exp_forces, &
              & pair_distribution_der,&
              & pair_distribution_partial_der,&
              & pair_distribution_partial_temp_der, energies_pdf, forces_pdf, virial_pdf)

         call time_end(time%pdf)
         !           if (rank == 0) print *, rank, " TIME_PDF = ", time%pdf(3)

      end if

      ! Now calculate the structure factors
      if (params%do_structure_factor) then
         call time_start(time%sf)
         call calculate_structure_factor(params, x_structure_factor, x_structure_factor_temp,&
              & y_structure_factor, y_structure_factor_temp,&
              & structure_factor_partial, structure_factor_partial_temp,&
              & x_pair_distribution, y_pair_distribution, &
              & pair_distribution_partial, n_species_actual, species_types_actual, n_atoms_of_species,&
              & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
              & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
              & neighbor_species, species, rank, q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
              & pair_distribution_partial_der, energies_sf, forces_sf, virial_sf, &
              & params%structure_factor_matrix_forces)

         call time_end(time%sf)

      end if

      if (params%do_xrd) then
         call time_start(time%xrd)
         if (params%xrd_debye) then
            call calculate_xrd_debye(params, x_xrd, x_xrd_temp, y_xrd, y_xrd_temp,&
                 & n_sites, positions, species, md_istep, mc_istep, i_beg, i_end,&
                 & ierr, rank, params%exp_forces, energies_xrd, forces_xrd,&
                 & virial_xrd, .false.)
         else
            call calculate_xrd(params, x_xrd, x_xrd_temp,&
                 & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
                 & structure_factor_partial, structure_factor_partial_temp,&
                 & n_species_actual, species_types_actual, n_atoms_of_species,&
                 & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
                 & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
                 & neighbor_species, species, rank, q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
                 & pair_distribution_partial_der, energies_xrd,&
                 & forces_xrd, virial_xrd, .false., params&
                 &%structure_factor_matrix_forces)
         end if

         call time_end(time%xrd)

         !           if (rank == 0) print *, rank, " TIME_XRD = ", time%xrd(3)

      end if

      if (params%do_nd) then
         call time_start(time%nd)
         if (params%xrd_debye) then
            call calculate_xrd_debye(params, x_nd, x_nd_temp, y_nd, y_nd_temp,&
                 & n_sites, positions, species, md_istep, mc_istep, i_beg, i_end,&
                 & ierr, rank, params%exp_forces, energies_nd, forces_nd,&
                 & virial_nd, .true.)
         else
            call calculate_xrd(params, x_nd, x_nd_temp,&
                 & y_nd, y_nd_temp, x_structure_factor, x_structure_factor_temp,&
                 & structure_factor_partial, structure_factor_partial_temp,&
                 & n_species_actual, species_types_actual, n_atoms_of_species,&
                 & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
                 & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
                 & neighbor_species, species, rank, q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
                 & pair_distribution_partial_der, energies_nd,&
                 & forces_nd, virial_nd, .true., params&
                 &%structure_factor_matrix_forces)
         end if

         call time_end(time%nd)

         !           if (rank == 0) print *, rank, " TIME_XRD = ", time%xrd(3)

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
         if (params%xrd_debye) then
!           The Debye route never touched the structure-factor arrays, and
!           finalize_xrd would free them out from under a structure_factor
!           calculation the deck asked for in its own right.
            call finalize_xrd_debye(x_xrd, x_xrd_temp, y_xrd, y_xrd_temp)
         else
            call finalize_xrd(params, x_xrd, x_xrd_temp,&
                 & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
                 & structure_factor_partial, structure_factor_partial_temp)
         end if
      end if

      if (params%do_nd) then
         if (params%xrd_debye) then
            call finalize_xrd_debye(x_nd, x_nd_temp, y_nd, y_nd_temp)
         else
            call finalize_xrd(params, x_nd, x_nd_temp,&
                 & y_nd, y_nd_temp, x_structure_factor, x_structure_factor_temp,&
                 & structure_factor_partial, structure_factor_partial_temp)
         end if
      end if

      deallocate (species_types_actual)

   end subroutine compute_exp_spectra
!**************************************************************************

!**************************************************************************
! The XPS spectrum predicted from core-electron binding energies, and the
! energies and forces that fitting it against experiment produces.
!
! Lifted verbatim out of turbogap.f90 (lines 1462-1622 at the time of the
! move). As with compute_exp_spectra, the single call whose argument list was
! split across #ifdef _MPIF90 collapses here and the caller chooses which set
! of arrays to pass. That was the last such statement in the driver.
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
      if (any_has_core_electron_be(soap_turbo_hypers) .and. (params%do_prediction) &
          .and. valid_xps) then
         call time_start(time%xps)

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

         call time_end(time%xps)
         !           if (rank == 0) print *, rank, " TIME_XPS = ", time%xps(3)

      else if (any_has_core_electron_be(soap_turbo_hypers) .and. params%do_xps) then
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

end module turbogap_exp
