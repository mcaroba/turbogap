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
! The CPU branch's copy also holds compute_exp_spectra, the pair-distribution,
! structure-factor and diffraction block. That one is not brought across yet:
! here it is 705 lines against 272 there, with 71 boundary variables against 37,
! because the device state -- cuBLAS handles, streams, gpu_exp, gpu_neigh --
! crosses the boundary too. It wants the backend seam, not a straight lift.
module turbogap_exp

  use timing
  use types
  use exp_utils
  use exp_interface
#ifdef _MPIF90
  use mpi
  use mpi_helper
#endif

  implicit none

  private
  public :: compute_exp_xps

contains

!**************************************************************************
  subroutine compute_exp_xps( params, n_sites, xyz, neighbors_list, n_neigh, &
       local_properties, local_properties_cart_der, soap_turbo_hypers, &
       a_box, b_box, c_box, indices, i_beg, i_end, j_beg, j_end, rank, &
       md_istep, mc_istep, valid_xps, xps_idx, core_be_lp_index, &
       write_condition, overwrite_condition, exp_output, &
       energies_lp, forces_lp, virial_lp, time_xps )

    implicit none

    type(input_parameters), intent(inout) :: params
    type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)
    integer, intent(in) :: n_sites, i_beg, i_end, j_beg, j_end, rank
    integer, intent(in) :: indices(1:3), md_istep, mc_istep, xps_idx, core_be_lp_index
    logical, intent(in) :: valid_xps
    logical, intent(inout) :: write_condition, overwrite_condition
    character*32, intent(inout) :: exp_output
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8,  intent(in), allocatable :: xyz(:,:)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:)
    real*8,  intent(inout), allocatable :: local_properties(:,:), local_properties_cart_der(:,:,:)

!   The caller passes the this_-prefixed arrays under MPI and the plain ones
!   otherwise.
    real*8, allocatable, intent(inout) :: energies_lp(:), forces_lp(:,:)
    real*8, intent(inout) :: virial_lp(1:3,1:3), time_xps(1:3)

!   Was driver state; block-local now.
    real*8, allocatable :: y_i_pred_all(:,:), v_neigh_lp(:)
    real*8  :: mag
    integer :: i, j, j2, k

    !     Compute core_electron_be energies and forces
    if( any( soap_turbo_hypers(:)%has_core_electron_be ) .and.( params%do_prediction ) &
      .and. valid_xps )then

      !time_xps(1) = MPI_wtime()
      call get_time( time_xps(1)  )


#ifdef _MPIF90
      allocate( energies_lp(1:n_sites) )
      energies_lp = 0.d0
      if( params%do_forces )then
        allocate( forces_lp(1:3,1:n_sites) )
        forces_lp = 0.d0
        virial_lp = 0.d0
      end if
#endif
      allocate(v_neigh_lp(1:j_end-j_beg+1))
      v_neigh_lp = 0.d0
      k = 0
      do i = i_beg, i_end
        do j = 1, n_neigh(i)
          !           I'm not sure if this is necessary or neighbors_list is already bounded between 1 and n_sites -> CHECK THIS
          j2 = mod(neighbors_list(j_beg + k)-1, n_sites) + 1
          k = k + 1
          !                 v_neigh_lp(k) = hirshfeld_v(j2)
          v_neigh_lp(k) = local_properties(j2, core_be_lp_index)
        end do
      end do

      call get_xps_spectra(params%exp_data(xps_idx)%data(1,:),&
        & params%exp_data(xps_idx)%data(2,:), params&
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
      call get_energy_scale( params%do_md, params%do_mc,&
        & md_istep, params%md_nsteps, mc_istep, params&
        &%mc_nsteps, params &
        &%exp_energy_scales_initial(xps_idx), params &
        &%exp_energy_scales_final(xps_idx), params &
        &%exp_energy_scales(xps_idx) )

      call get_exp_pred_spectra_energies_forces( params&
        &%exp_energy_scales(xps_idx),&
        & local_properties(i_beg:i_end,core_be_lp_index),&
        & local_properties_cart_der(1:3, j_beg:j_end,&
        & core_be_lp_index ), n_neigh(i_beg:i_end),&
        & neighbors_list(j_beg:j_end), params%xps_sigma,&
        & params%exp_data(xps_idx)%n_samples, mag, params&
        &%exp_data(xps_idx)%x, params %exp_data(xps_idx)%y,&
        & params%exp_data(xps_idx) %y_pred,&
        & y_i_pred_all(i_beg:i_end, 1:params &
        &%exp_data(xps_idx)%n_samples), params %do_forces, &
        & xyz(1:3, j_beg:j_end),&
      & energies_lp(i_beg:i_end), forces_lp, virial_lp, params%exp_similarity_type, rank )

      ! if (rank == 0)then
      !    open(unit=11, file="tg_xps.dat", status="unknown")
      !    do i = 1, params%exp_data(xps_idx)%n_samples
      !       write(11, '(1X,F20.8,1X,F20.8)') params%exp_data(xps_idx)%x(i), params%exp_data(xps_idx)%y_pred(i)
      !    end do
      !    close(11)
      ! end if


      call get_write_condition( params%do_mc, params%do_md&
        &, mc_istep, md_istep, params%write_xyz,&
        & write_condition)

      if (rank == 0 .and. params%write_exp .and. write_condition)then

        call get_overwrite_condition( params%do_mc, params%do_md&
          &, mc_istep, md_istep, params%write_xyz, overwrite_condition)

        call write_exp_datan(params%exp_data(xps_idx)&
          &%x(1:params%exp_data(xps_idx)%n_samples), params&
          &%exp_data(xps_idx)%y_pred(1:params&
          &%exp_data(xps_idx)%n_samples),&
          & overwrite_condition, "xps_prediction.dat",&
          & params%exp_data(xps_idx)%label)

        if ( .not.  params%exp_data(xps_idx)%wrote_exp ) then


          call preprocess_exp_data(params, params%exp_data(xps_idx)%x,&
            & params%exp_data(xps_idx)%y, params%exp_data(xps_idx)%label,&
            & n_sites, dot_product( cross_product(a_box,&
            & b_box), c_box ) / (dfloat(indices(1)*indices(2) &
            &*indices(3)) ), params%exp_data(xps_idx)%input, exp_output, .true. )

          call write_exp_datan(params%exp_data(xps_idx)&
            &%x(1:params%exp_data(xps_idx)%n_samples),&
            & params%exp_data(xps_idx)%y(1:params&
            &%exp_data(xps_idx)%n_samples),&
            & overwrite_condition, "xps_exp.dat" , params&
            &%exp_data(xps_idx)%label)
          params%exp_data(xps_idx)%wrote_exp = .true.
        end if

        ! else
        !    call write_exp_data(params%exp_data(xps_idx)%x, params%exp_data(xps_idx)%y_pred, mc_istep == 0, "xps_prediction.dat" )
        !    call write_exp_data(params%exp_data(xps_idx)%x, params%exp_data(xps_idx)%y, mc_istep == 0, "xps_exp.dat" )

      end if


      !deallocate( params%exp_data(xps_idx)%y_pred )
      if (allocated(y_i_pred_all)) deallocate(y_i_pred_all)
      ! sim_exp_pred would be an energy if multiplied by some energy scale \gamma * ( 1 - sim )
      ! sim_exp_pred_der would be the array of forces if multiplied by (- \gamma )
      deallocate(v_neigh_lp)



      !time_xps(2) = MPI_wtime()
      call get_time( time_xps(2)  )

      time_xps(3) = time_xps(3) + time_xps(2) - time_xps(1)
      !           if (rank == 0) print *, rank, " TIME_XPS = ", time_xps(3)

    end if

  end subroutine compute_exp_xps
!**************************************************************************

end module turbogap_exp
