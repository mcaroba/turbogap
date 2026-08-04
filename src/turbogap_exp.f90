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
  subroutine compute_exp_spectra( params, n_sites, species, rjs, xyz, neighbors_list, &
       n_neigh, neighbor_species, indices, a_box, b_box, c_box, &
       i_beg, i_end, j_beg, j_end, rank, ntasks, ierr, md_istep, mc_istep, &
       energies_pdf, forces_pdf, virial_pdf, &
       energies_sf,  forces_sf,  virial_sf,  &
       energies_xrd, forces_xrd, virial_xrd, &
       energies_nd,  forces_nd,  virial_nd,  &
       time_pdf, time_sf, time_xrd, time_nd )

    implicit none

!   Passed through to the exp_interface routines, which declare the arrays
!   allocatable, so these must be allocatable too.
    type(input_parameters), intent(inout) :: params
    integer, intent(in) :: n_sites, i_beg, i_end, j_beg, j_end, rank, ntasks
    integer, intent(in) :: indices(1:3), md_istep, mc_istep
    integer, intent(inout) :: ierr
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:), &
         neighbor_species(:), species(:)

!   The four contribution families this routine produces. The caller passes the
!   this_-prefixed arrays under MPI and the plain ones otherwise.
    real*8, allocatable, intent(inout) :: energies_pdf(:), forces_pdf(:,:), &
         energies_sf(:),  forces_sf(:,:),  &
         energies_xrd(:), forces_xrd(:,:), &
         energies_nd(:),  forces_nd(:,:)
    real*8, intent(inout) :: virial_pdf(1:3,1:3), virial_sf(1:3,1:3), &
         virial_xrd(1:3,1:3), virial_nd(1:3,1:3)
    real*8, intent(inout) :: time_pdf(1:3), time_sf(1:3), time_xrd(1:3), time_nd(1:3)

!   Was driver state; block-local now.
    real*8, allocatable :: x_pair_distribution(:), y_pair_distribution(:), &
         y_pair_distribution_temp(:), pair_distribution_partial(:,:), &
         pair_distribution_partial_temp(:,:), pair_distribution_der(:,:), &
         pair_distribution_partial_der(:,:,:), pair_distribution_partial_temp_der(:,:,:), &
         n_atoms_of_species(:), sinc_factor_matrix(:,:), &
         x_structure_factor(:), x_structure_factor_temp(:), &
         y_structure_factor(:), y_structure_factor_temp(:), &
         structure_factor_partial(:,:), structure_factor_partial_temp(:,:), &
         x_xrd(:), x_xrd_temp(:), y_xrd(:), y_xrd_temp(:), &
         x_nd(:), x_nd_temp(:), y_nd(:), y_nd_temp(:)
    character*8, allocatable :: species_types_actual(:)
    integer :: i, n_species_actual, q_beg, q_end

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
    if (allocated( species_types_actual)) deallocate( species_types_actual )
    allocate(species_types_actual(1:n_species_actual))

    n_species_actual = 0
    do i = 1, n_sites
       if (species(i) > n_species_actual)then
          n_species_actual = n_species_actual + 1
          species_types_actual(n_species_actual) = params%species_types( species(i) )
       end if
    end do


    if (params%do_pair_distribution)then
       call get_time(time_pdf(1))

       call calculate_pair_distribution( params, x_pair_distribution&
            &, y_pair_distribution, y_pair_distribution_temp,&
            & pair_distribution_partial, pair_distribution_partial_temp, &
            & n_species_actual, species_types_actual, n_atoms_of_species, n_sites, a_box,&
            & b_box, c_box, indices, md_istep, mc_istep, i_beg, i_end,&
            & j_beg, j_end, ierr , rjs, xyz, neighbors_list,&
            & n_neigh, neighbor_species, species, rank, params%exp_forces, &
            & pair_distribution_der,&
            & pair_distribution_partial_der,&
            & pair_distribution_partial_temp_der, energies_pdf, forces_pdf, virial_pdf)

       call get_time(time_pdf(2))
       time_pdf(3) = time_pdf(3) + time_pdf(2) - time_pdf(1)
       !           if (rank == 0) print *, rank, " TIME_PDF = ", time_pdf(3)



    end if

    ! Now calculate the structure factors
    if (params%do_structure_factor )then
       call get_time(time_sf(1))
       call calculate_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
            & y_structure_factor, y_structure_factor_temp,&
            & structure_factor_partial, structure_factor_partial_temp,&
            & x_pair_distribution, y_pair_distribution, &
            & pair_distribution_partial, n_species_actual, species_types_actual, n_atoms_of_species,&
            & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
            & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
            & neighbor_species, species, rank , q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
            & pair_distribution_partial_der, energies_sf, forces_sf, virial_sf, &
            & params%structure_factor_matrix_forces)


       call get_time(time_sf(2))
       time_sf(3) = time_sf(3) + time_sf(2) - time_sf(1)



    end if

    if ( params%do_xrd )then
       call get_time(time_xrd(1))
       call calculate_xrd( params, x_xrd, x_xrd_temp,&
            & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
            & structure_factor_partial, structure_factor_partial_temp,&
            & n_species_actual, species_types_actual, n_atoms_of_species,&
            & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
            & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
            & neighbor_species, species, rank , q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
            & pair_distribution_partial_der, energies_xrd,&
            & forces_xrd, virial_xrd, .false., params&
            &%structure_factor_matrix_forces)


       call get_time(time_xrd(2))
       time_xrd(3) = time_xrd(3) + time_xrd(2) - time_xrd(1)

       !           if (rank == 0) print *, rank, " TIME_XRD = ", time_xrd(3)

    end if


    if ( params%do_nd )then
       call get_time(time_nd(1))
       call calculate_xrd( params, x_nd, x_nd_temp,&
            & y_nd, y_nd_temp, x_structure_factor, x_structure_factor_temp,&
            & structure_factor_partial, structure_factor_partial_temp,&
            & n_species_actual, species_types_actual, n_atoms_of_species,&
            & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
            & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
            & neighbor_species, species, rank , q_beg, q_end, ntasks, sinc_factor_matrix, params%exp_forces, &
            & pair_distribution_partial_der, energies_nd,&
            & forces_nd, virial_nd, .true., params&
            &%structure_factor_matrix_forces)


       call get_time(time_nd(2))
       time_nd(3) = time_nd(3) + time_nd(2) - time_nd(1)

       !           if (rank == 0) print *, rank, " TIME_XRD = ", time_xrd(3)

    end if


    !################################################################!
    !###---   Compute similarity of experimental predictions   ---###!
    !################################################################!

    if ( params%do_exp )then
       do i = 1, params%n_exp
          ! First normalize the spectrum if it matches some type of experimental data

          ! Allocate the prediction data
          ! if (.not. allocated( params%exp_data(i)%y_pred ) )then
          !    allocate( params%exp_data(i)%y_pred( 1:size(params%exp_data(i)%y, 1) ))
          ! end if

          if ( trim(params%exp_data(i)%label) == 'pair_distribution' )then
             params%exp_data(i)%y_pred = y_pair_distribution
          elseif ( trim(params%exp_data(i)%label) == 'structure_factor' )then
             params%exp_data(i)%y_pred = y_structure_factor
          elseif ( trim(params%exp_data(i)%label) == 'xrd' )then
             params%exp_data(i)%y_pred = y_xrd
          elseif ( trim(params%exp_data(i)%label) == 'nd' )then
             params%exp_data(i)%y_pred = y_nd
          end if

          if ( params%exp_data(i)%compute_similarity .and. allocated(params%exp_data(i)%y) )then
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

    if (params%do_pair_distribution)then
       call finalize_pair_distribution( params, x_pair_distribution&
            &, y_pair_distribution, y_pair_distribution_temp,&
            & pair_distribution_partial, pair_distribution_partial_temp, params&
            &%exp_forces, pair_distribution_der,&
            & pair_distribution_partial_der,&
            & pair_distribution_partial_temp_der, n_atoms_of_species, rank)

    end if

    ! Now calculate the structure factors
    if (params%do_structure_factor )then

       call finalize_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
            & y_structure_factor, y_structure_factor_temp,&
            & structure_factor_partial, structure_factor_partial_temp,&
            & x_pair_distribution, y_pair_distribution, &
            & pair_distribution_partial, sinc_factor_matrix)

    end if

    if ( params%do_xrd )then
       call finalize_xrd( params, x_xrd, x_xrd_temp,&
            & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
            & structure_factor_partial, structure_factor_partial_temp)
    end if

    if ( params%do_nd )then
       call finalize_xrd( params, x_nd, x_nd_temp,&
            & y_nd, y_nd_temp, x_structure_factor, x_structure_factor_temp,&
            & structure_factor_partial, structure_factor_partial_temp)
    end if

    deallocate(species_types_actual)

  end subroutine compute_exp_spectra
!**************************************************************************

end module turbogap_exp
