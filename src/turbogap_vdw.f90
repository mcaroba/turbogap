! Copyright (c) 2020-2026 by Albert Bartok and Miguel Caro
!
! van der Waals / many-body dispersion phase of the TurboGAP driver.
!
! This is a verbatim lift of the block that used to sit inline in
! turbogap.f90 between the "Compute vdW energies and forces" comment and its
! matching end if. Nothing about the computation changed; the point of the
! move is that 29 arrays which the main program declared, but which only this
! block ever touched, are now locals here and are gone from the driver's
! ~300-variable shared scope.
!
! The argument list is long. That is not the end state -- it is what remains
! once the block-local variables are removed, and it is what bundling the
! energy/force/virial families into a derived type is meant to collapse. It is
! left explicit here so that this commit is a pure move and can be reviewed
! as one.

module turbogap_vdw

  use timing
  use types
  use vdw
  use misc
#ifdef _MPIF90
  use mpi
  use mpi_helper
#endif

  implicit none

  private
  public :: vdw_state, compute_vdw

! Persistent state of the ts+mbd correction.
!
! These buffers deliberately outlive a single call: the MBD correction is
! recomputed only every mbd_correction_freq steps and reused in between, so
! their contents must survive from one call to the next. As locals of the main
! program that requirement was invisible, and it is exactly what the crash in
! tests/regression/KNOWN_ISSUES.md is -- the allocation and the consumer branch
! on the same modulo() expression, so in predict mode the consumer runs having
! never been allocated. Naming the lifetime is the first step to fixing it.
  type :: vdw_state
     real*8, allocatable :: this_energies_vdw_corr(:)
     real*8, allocatable :: this_forces_vdw_corr(:,:)
     real*8, allocatable :: this_local_virial_vdw_diag_corr(:,:)
     real*8 :: this_virial_vdw_corr(1:3, 1:3) = 0.d0
     real*8 :: virial_vdw_corr(1:3, 1:3) = 0.d0
  end type vdw_state

contains

  subroutine compute_vdw( params, has_vdw, n_sites, n_neigh, neighbors_list, &
       neighbor_species, rjs, xyz, local_properties, local_properties_cart_der, &
       vdw_lp_index, i_beg, i_end, j_beg, j_end, n_atom_pairs_by_rank, &
       site_in_rank, indices, rank, ntasks, md_istep, state, &
       energies_vdw, forces_vdw, virial_vdw, local_virial_vdw_diag, &
       this_energies_vdw, this_forces_vdw, this_virial_vdw, &
       this_local_virial_vdw_diag, energies_vdw_corr, forces_vdw_corr, &
       local_virial_vdw_diag_corr, mbd_ts_scaling, this_mbd_ts_scaling, &
       update_mbd_ts_scaling, time_vdw )

    implicit none

!   Input
    type(input_parameters), intent(in) :: params
    logical, intent(in) :: has_vdw
    integer, intent(in) :: n_sites, vdw_lp_index, i_beg, i_end, j_beg, j_end
    integer, intent(in) :: rank, ntasks, md_istep
    integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
    integer, intent(in) :: n_atom_pairs_by_rank(:), site_in_rank(:), indices(1:3)
    real*8, intent(in) :: rjs(:), xyz(:,:)
    real*8, intent(in) :: local_properties(:,:), local_properties_cart_der(:,:,:)

!   Correction state, persistent across calls
    type(vdw_state), intent(inout) :: state

!   Output
    real*8, allocatable, intent(inout) :: energies_vdw(:), forces_vdw(:,:)
    real*8, allocatable, intent(inout) :: local_virial_vdw_diag(:,:)
    real*8, intent(inout) :: virial_vdw(1:3, 1:3), this_virial_vdw(1:3, 1:3)
    real*8, allocatable, intent(inout) :: this_energies_vdw(:), this_forces_vdw(:,:)
    real*8, allocatable, intent(inout) :: this_local_virial_vdw_diag(:,:)
    real*8, allocatable, intent(inout) :: energies_vdw_corr(:), forces_vdw_corr(:,:)
    real*8, allocatable, intent(inout) :: local_virial_vdw_diag_corr(:,:)
    real*8, allocatable, intent(inout) :: mbd_ts_scaling(:), this_mbd_ts_scaling(:)
    logical, intent(inout) :: update_mbd_ts_scaling
    real*8, intent(inout) :: time_vdw(1:3)

!   Local. Every one of these was previously a variable of the main program
!   that no other part of it referenced.
    real*8, allocatable :: v_neigh_vdw(:)
    real*8, allocatable :: alpha_SCS(:), omega_SCS(:), alpha_SCS_grad(:,:)
    real*8, allocatable :: this_alpha_SCS(:), this_omega_SCS(:)
    real*8, allocatable :: c6_scs(:), r0_scs(:), alpha0_scs(:), S_xyz_inv(:,:)
    real*8, allocatable :: hirshfeld_v_cart_der_send(:,:), hirshfeld_v_cart_der_receive(:,:)
    real*8, allocatable :: this_hirshfeld_v_cart_der_receive(:,:), hirshfeld_v_cart_der_ji(:,:)
    integer, allocatable :: hirshfeld_transfer(:,:), this_hirshfeld_transfer(:)
    integer, allocatable :: i_send(:), j_send(:), k_array(:), k_start(:)
    integer, allocatable :: i_receive(:), j_receive(:), this_i_receive(:), this_j_receive(:)
    integer, allocatable :: hirshfeld_disp(:)
    integer :: i, j, k, i2, j2, k2, n, jx, jy, jz, ierr
    logical :: include_2b
    logical :: is_correction_step
!   Scratch timers for the commented-out stage timings inside the block. The
!   main program's time1/time2 were used for this before; both are written
!   before they are next read there, so nothing depended on the clobbering.
    real*8 :: time1, time2

    if( has_vdw .and. ( params%do_prediction ) &
         .and. ( params%vdw_type == "ts" .or. params%vdw_type == "mbd".or. params%vdw_type == "ts+mbd" ) )then
           call get_time(time_vdw(1))

!          Does this call recompute the ts+mbd correction, or reuse the one a
!          previous call stored? The reuse half of that cycle only means
!          anything during MD, where there was a previous step to store it.
!          predict and mc leave md_istep at -1, and modulo(-1, freq) is
!          non-zero for every freq > 1, so spelling this test as the modulo
!          alone sent them down the reuse path with nothing stored: reading
!          correction buffers the recompute branch had never allocated, and
!          an mbd_ts_scaling nobody had initialised. Outside MD there is no
!          cycle to be part of, so every call recomputes.
!
!          For md_istep >= 0 this is exactly the old expression, so the MD
!          path is unchanged. It is computed once because it used to be
!          spelled out at six sites that all had to agree, and the crash was
!          two of them disagreeing about whether the buffers existed.
           is_correction_step = ( md_istep < 0 ) .or. &
                                ( modulo(md_istep, params%mbd_correction_freq) == 0 )
#ifdef _MPIF90
           allocate( this_energies_vdw(1:n_sites) )
           this_energies_vdw = 0.d0
           if( params%do_forces )then
              allocate( this_forces_vdw(1:3,1:n_sites) )
              allocate( this_local_virial_vdw_diag(1:3, 1:n_sites) )
              this_forces_vdw = 0.d0
              this_virial_vdw = 0.d0
              this_local_virial_vdw_diag = 0.d0
           end if
! TEST COMMENT THIS OUT
!        call mpi_reduce(local_properties(:,vdw_lp_index), this_local_properties(:,vdw_lp_index), n_sites, &
!                        MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!        if( params%do_forces )then
!         I'm not sure if this is necessary at all... CHECK
!          call mpi_reduce(hirshfeld_v_cart_der, this_hirshfeld_v_cart_der, 3*n_atom_pairs, MPI_DOUBLE_PRECISION, MPI_SUM, &
!                          0, MPI_COMM_WORLD, ierr)
!          hirshfeld_v_cart_der = this_hirshfeld_v_cart_der
!        end if
!        local_properties(:,vdw_lp_index) = this_local_properties(:,vdw_lp_index)
!        call mpi_bcast(local_properties(:,vdw_lp_index), n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       HERE WE TRANSFER THE HIRSHFELD VOLUME GRADIENTS IF NEEDED FOR SCS
!       Putting this outside the if condition to avoid segfaults for vdw_hirsh_grad = .false.
        if( allocated(hirshfeld_v_cart_der_ji) )deallocate( hirshfeld_v_cart_der_ji )
          allocate( hirshfeld_v_cart_der_ji(1:3, 1:n_atom_pairs_by_rank(rank+1)) )
          hirshfeld_v_cart_der_ji = 0.d0
!!!!!!!!!!!!!!!!!
        if( params%do_forces .and. params%vdw_hirsh_grad )then
          allocate( hirshfeld_transfer(1:ntasks, 1:ntasks) )
          allocate( this_hirshfeld_transfer(1:ntasks) )
          this_hirshfeld_transfer = 0
          hirshfeld_transfer = 0
          do i = 1, n_atom_pairs_by_rank(rank+1)
!           j is the atom's index in the supercell, use mod() to wrap it back to the unit cell
            j = neighbors_list(i)
            k = site_in_rank( mod(j-1,n_sites)+1 )
!           Transfer only those derivatives that are not zero
            if( .not. all( abs(local_properties_cart_der(1:3,i,vdw_lp_index)) < 1.d-20 ) )then
!             rank = rank sends another derivative to rank = k
              this_hirshfeld_transfer( k+1 ) = this_hirshfeld_transfer( k+1 ) + 1
            end if
          end do
          call mpi_allgather( this_hirshfeld_transfer, ntasks, MPI_INTEGER, hirshfeld_transfer, ntasks, &
                              MPI_INTEGER, MPI_COMM_WORLD, ierr )
!         Now we repeat the operation above but actually populating the array that is to be scattered by this rank.
!         The gradients are stored in hirshfeld_v_cart_der_send in such a way that contiguous blocks of memory
!         are going to be sent to the same rank
!         i_send stores the *central* atom's index; j_send stores the neighbor atom's index
          allocate( hirshfeld_v_cart_der_send(1:3, 1:sum(hirshfeld_transfer(1:ntasks, rank+1))) )
          allocate( i_send(1:sum(hirshfeld_transfer(1:ntasks, rank+1))) )
          allocate( j_send(1:sum(hirshfeld_transfer(1:ntasks, rank+1))) )
          allocate( k_array(1:ntasks) )
          k = 0
!         This points to the part of hirshfeld_v_cart_der_send where gradients are to be put
          k_array = 0
          do i = 2, ntasks
            k_array(i) = k_array(i-1) + this_hirshfeld_transfer(i-1)
          end do
          do i = i_beg, i_end
            do j = 1, n_neigh(i)
              k = k + 1
              if( j == 1 )then
                i2 = neighbors_list(k)
              end if
              j2 = neighbors_list(k)
              k2 = site_in_rank( mod(j2-1, n_sites)+1 )
              if( .not. all( abs(local_properties_cart_der(1:3, k, vdw_lp_index)) < 1.d-20 ) )then
!               Roll the pointer for rank k2 by 1
                k_array(k2+1) = k_array(k2+1) + 1
                hirshfeld_v_cart_der_send(1:3, k_array(k2+1)) = local_properties_cart_der(1:3, k, vdw_lp_index)
!               i2 and j2 are "supercell" indices, although i2 is always also a primitive unit cell index
                i_send(k_array(k2+1)) = i2
                j_send(k_array(k2+1)) = j2
              end if
            end do
          end do
!         Here we allocate the arrays where we're going to put the received data
          allocate( hirshfeld_v_cart_der_receive(1:3, 1:sum(hirshfeld_transfer(rank+1, 1:ntasks))) )
          allocate( i_receive(1:sum(hirshfeld_transfer(rank+1, 1:ntasks))) )
          allocate( j_receive(1:sum(hirshfeld_transfer(rank+1, 1:ntasks))) )
!         Do the communication | THIS IS SLOW, BUT A RELATIVELY STRAIGHTFORWARD WAY TO IMPLEMENT IT
!                                CHECK POSSIBLE WAYS TO SPEED THIS UP (REDUCE COMMUNICATION)
          allocate( hirshfeld_disp(1:ntasks) )
          do i = 1, ntasks
            allocate( this_hirshfeld_v_cart_der_receive(1:3, 1:hirshfeld_transfer(rank+1, i)) )
            allocate( this_i_receive(1:hirshfeld_transfer(rank+1, i)) )
            allocate( this_j_receive(1:hirshfeld_transfer(rank+1, i)) )
            hirshfeld_disp(1) = 0
            do j = 2, ntasks
              hirshfeld_disp(j) = hirshfeld_disp(j-1) + hirshfeld_transfer(j-1, i)
            end do
            call mpi_scatterv(hirshfeld_v_cart_der_send, 3*hirshfeld_transfer(1:ntasks, i), 3*hirshfeld_disp, &
                              MPI_DOUBLE_PRECISION, this_hirshfeld_v_cart_der_receive, &
                              3*hirshfeld_transfer(rank+1, i), MPI_DOUBLE_PRECISION, i-1, &
                              MPI_COMM_WORLD, ierr )
            call mpi_scatterv(i_send, hirshfeld_transfer(1:ntasks, i), hirshfeld_disp, MPI_INTEGER, this_i_receive, &
                              hirshfeld_transfer(rank+1, i), MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr )
            call mpi_scatterv(j_send, hirshfeld_transfer(1:ntasks, i), hirshfeld_disp, MPI_INTEGER, this_j_receive, &
                              hirshfeld_transfer(rank+1, i), MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr )
            hirshfeld_disp(1) = 0
            do j = 2, ntasks
              hirshfeld_disp(j) = hirshfeld_disp(j-1) + hirshfeld_transfer(rank+1, j-1)
            end do
            hirshfeld_v_cart_der_receive(1:3, hirshfeld_disp(i)+1:hirshfeld_disp(i)+hirshfeld_transfer(rank+1, i)) = &
                this_hirshfeld_v_cart_der_receive
            i_receive(hirshfeld_disp(i)+1:hirshfeld_disp(i)+hirshfeld_transfer(rank+1, i)) = this_i_receive
            j_receive(hirshfeld_disp(i)+1:hirshfeld_disp(i)+hirshfeld_transfer(rank+1, i)) = this_j_receive
            deallocate( this_hirshfeld_v_cart_der_receive, this_i_receive, this_j_receive )
          end do
!         Now we do the inverse mapping so that we know where to find grad_i(nu_j). Things to note:
!         *) On this rank, the central atom is i and the neighbor is j
!         *) i_receive and j_receive contain central and neighbor atoms from the point of view of the rank that
!            sent them. We do need to swap i and j
!         *) We need to consider that for this rank j can be > n_sites if a supercell was used to construct the
!            neighbors list. Therefore, local-rank (i,j) can correspond to a different (j',i') on the remote
!            rank. i <= n_sites and j' <= n_sites, but that's only necessarily true for j and i' modulo n_sites.
!            The (i,j) and (j',i') tuples need their two elements to be separated by the same distance to be
!            equivalent (in addition to being equal modulo n_sites)
!
!         First, reduce the indices to those native to the local rank
          k = 0
          do i = 1, ntasks
            do j = hirshfeld_disp(i)+1, hirshfeld_disp(i)+hirshfeld_transfer(rank+1, i)
              k = k + 1
              i2 = i_receive(k)
              j2 = j_receive(k)
              jx = modulo((j2-1)/n_sites, indices(1))
              jy = modulo((j2-1)/(indices(1)*n_sites), indices(2))
              jz = modulo((j2-1)/(indices(1)*indices(2)*n_sites), indices(3))
!             This reduces j2 to the primitive unit cell
              j_receive(k) = j2 - jx*n_sites - jy*indices(1)*n_sites - jz*indices(1)*indices(2)*n_sites
!             Now we need to provide the correct supercell tag for i2 !!!!!!!!!!!!!!!!!!!!! NOT SURE IF THIS IS CORRECT!!!!!
!              i_receive(k) = i2 + jx*n_sites + jy*indices(1)*n_sites + jz*indices(1)*indices(2)*n_sites
              i_receive(k) = i2 + modulo(-jx, indices(1))*n_sites + modulo(-jy, indices(2))*indices(1)*n_sites &
                                + modulo(-jz, indices(3))*indices(1)*indices(2)*n_sites
            end do
          end do
!
!         Now we need to map the position in the big hirshfeld_v_cart_der array to those in the hirshfeld_v_cart_der_receive array
!         NOTE THIS ARRAY INDEX DOES NOT START BY 1
          allocate( k_start(i_beg:i_end) )
          k = 1
          do i = i_beg, i_end
            k_start(i) = k
            k = k + n_neigh(i)
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIX: PUT THE DEALLOCATION AFTER MBD COMPUTATION
          !if( allocated(hirshfeld_v_cart_der_ji) )deallocate( hirshfeld_v_cart_der_ji )
          !allocate( hirshfeld_v_cart_der_ji(1:3, 1:n_atom_pairs_by_rank(rank+1)) )
          !hirshfeld_v_cart_der_ji = 0.d0
          do k2 = 1, size(i_receive)
!           These indices are inverted here
            j = i_receive(k2)
            i = j_receive(k2)
            do k = k_start(i), k_start(i)+n_neigh(i)-1
              j2 = neighbors_list(k)
              if( j == j2 )then
                hirshfeld_v_cart_der_ji(1:3, k) = hirshfeld_v_cart_der_receive(1:3, k2)
              end if
            end do
          end do
!
          deallocate( this_hirshfeld_transfer, hirshfeld_transfer, hirshfeld_v_cart_der_send, i_send, j_send, &
                      k_array, hirshfeld_v_cart_der_receive, i_receive, j_receive, hirshfeld_disp, k_start )
        end if
#endif


if( params%vdw_type == "ts+mbd" .and. is_correction_step )then
!        if( allocated(state%this_energies_vdw_corr) )deallocate( state%this_energies_vdw_corr, this_mbd_ts_scaling )
        if( allocated(state%this_energies_vdw_corr) )deallocate( state%this_energies_vdw_corr )
        allocate( state%this_energies_vdw_corr(1:n_sites) )
        state%this_energies_vdw_corr = 0.d0
        if( params%do_forces )then
          if( allocated(state%this_forces_vdw_corr) )deallocate( state%this_forces_vdw_corr, state%this_local_virial_vdw_diag_corr )
          allocate( state%this_forces_vdw_corr(1:3,1:n_sites) )
          state%this_forces_vdw_corr = 0.d0
          allocate( state%this_local_virial_vdw_diag_corr(1:3, 1:n_sites) )
          state%this_local_virial_vdw_diag_corr = 0.d0
        end if
else if( params%vdw_type == "ts" )then
        if( allocated(this_mbd_ts_scaling) )deallocate( this_mbd_ts_scaling )
        allocate( this_mbd_ts_scaling(1:n_sites) )
        this_mbd_ts_scaling = 1.d0
end if


           if( .not. allocated(v_neigh_vdw) )allocate(v_neigh_vdw(1:j_end-j_beg+1))
           v_neigh_vdw = 0.d0
           k = 0
if( params%vdw_type == "ts" .or. params%vdw_type == "ts+mbd" .or. params%vdw_type == "mbd")then
           do i = i_beg, i_end
              do j = 1, n_neigh(i)
                 !           I'm not sure if this is necessary or neighbors_list is already bounded between 1 and n_sites -> CHECK THIS
                 j2 = mod(neighbors_list(j_beg + k)-1, n_sites) + 1
                 k = k + 1
                 !                 v_neigh_vdw(k) = hirshfeld_v(j2)
                 v_neigh_vdw(k) = local_properties(j2, vdw_lp_index)
              end do
           end do
end if
! TODO: change this back to get_ts_energy_and_forces and implement call for mbd energy
        if( params%vdw_type == "ts" .or. params%vdw_type == "ts+mbd" )then
          call get_ts_energy_and_forces( local_properties(i_beg:i_end, vdw_lp_index), &
                & local_properties_cart_der(1:3, j_beg:j_end, vdw_lp_index), &
                n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                neighbor_species(j_beg:j_end), &
                params%vdw_rcut, params%vdw_buffer, &
                params%vdw_rcut_inner, params%vdw_buffer_inner, &
                rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                params%vdw_sr, params%vdw_d, params%vdw_c6_ref, params%vdw_r0_ref, &
                params%vdw_alpha0_ref, params%do_forces, &
                params%poly_cut_xmin, params%poly_cut_xmax, &
#ifdef _MPIF90
                this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw, &
                this_local_virial_vdw_diag, this_mbd_ts_scaling )
#else
                energies_vdw(i_beg:i_end), forces_vdw, virial_vdw, local_virial_vdw_diag, &
                mbd_ts_scaling )
#endif
          if( .not. (params%vdw_type == "ts+mbd" .and. is_correction_step) )then
            deallocate(v_neigh_vdw)
          end if
! TESTING FOR MERGING
!write(*,*) "local_properties", local_properties(i_beg:i_end,vdw_lp_index)
!write(*,*) "this_energies_vdw", this_energies_vdw
        end if
        if( params%vdw_type == "mbd" .or. &
            (params%vdw_type == "ts+mbd" .and. is_correction_step) )then
          if( params%vdw_type == "ts+mbd" .and. is_correction_step )then
#ifdef _MPIF90
!            state%this_energies_vdw_corr(i_beg:i_end) = this_energies_vdw(i_beg:i_end)
            state%this_energies_vdw_corr = this_energies_vdw
            state%this_forces_vdw_corr = this_forces_vdw
            state%this_virial_vdw_corr = this_virial_vdw
            state%this_local_virial_vdw_diag_corr = this_local_virial_vdw_diag
            this_energies_vdw = 0.d0
            this_forces_vdw = 0.d0
            this_virial_vdw = 0.d0
            this_local_virial_vdw_diag = 0.d0
#else
!            energies_vdw_corr(i_beg:i_end) = energies_vdw(i_beg:i_end)
            energies_vdw_corr = energies_vdw
            forces_vdw_corr = forces_vdw
            state%virial_vdw_corr = virial_vdw
            local_virial_vdw_diag_corr = local_virial_vdw_diag
            energies_vdw = 0.d0
            forces_vdw = 0.d0
            virial_vdw = 0.d0
            local_virial_vdw_diag = 0.d0
#endif
          end if
          allocate( alpha_SCS(1:n_sites) )
          allocate( this_alpha_SCS(1:n_sites) )
          allocate( omega_SCS(1:n_sites) )
          allocate( this_omega_SCS(1:n_sites) )
          alpha_SCS = 0.d0
          omega_SCS = 0.d0
          this_alpha_SCS = 0.d0
          this_omega_SCS = 0.d0
          !allocate( alpha_SCS_grad(j_beg:j_end,1:3) )
          allocate( alpha_SCS_grad(1:n_sites,1:3) )
          allocate( c6_scs(1:j_end-j_beg+1) )
          allocate( r0_scs(1:j_end-j_beg+1) )
          allocate( alpha0_scs(1:j_end-j_beg+1) )
call get_time(time1)
          !write(*,*) "SCS calculation starts here"
          call get_scs_polarizabilities( n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                                         neighbor_species(j_beg:j_end), &
                                         params%vdw_scs_rcut, params%vdw_buffer, &
                                         rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                                         params%vdw_sr_mbd, params%vdw_d_mbd, params%vdw_c6_ref, params%vdw_r0_ref, &
                                         params%vdw_alpha0_ref, &
                                         params%vdw_polynomial, params%vdw_omega_ref, &
                                         alpha_SCS(i_beg:i_end), omega_SCS(i_beg:i_end), &
#ifdef _MPIF90
                                         this_forces_vdw )
#else
                                         forces_vdw )
#endif
call get_time(time2)
!write(*,*) "SCS timing", time2-time1
call get_time(time1)

call mpi_reduce(alpha_SCS, this_alpha_SCS, n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
alpha_SCS = this_alpha_SCS
call mpi_bcast(alpha_SCS, n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call mpi_reduce(omega_SCS, this_omega_SCS, n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
omega_SCS = this_omega_SCS
call mpi_bcast(omega_SCS, n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

call get_time(time2)
!write(*,*) "Communication timing", time2-time1

!write(*,*) "alpha_SCS"
!do i = 1, n_sites
!  write(*,*) i, alpha_SCS(i), omega_SCS(i)
!end do
call get_time(time1)
!write(*,*) "scs timing", time1-time2
if ( params%vdw_2b_rcut > params%vdw_mbd_rcut ) then ! Call 2b version if primary 2b cut-off is larger than primary mbd cut-off
include_2b = .true.
          call get_mbd_energies_and_forces( hirshfeld_v_cart_der_ji(1:3,j_beg:j_end), &
                                         n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                                         neighbor_species(j_beg:j_end), &
                                         params%vdw_scs_rcut, params%vdw_loc_rcut, params%vdw_2b_rcut, &
                                         params%vdw_2b_rcut2, &
                                         params%vdw_buffer, &
                                         rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                                         params%vdw_sr_mbd, params%vdw_d_mbd, params%vdw_c6_ref, params%vdw_r0_ref, &
                                         params%vdw_alpha0_ref, params%vdw_mbd_grad, params%vdw_hirsh_grad, &
                                         params%vdw_polynomial, params%do_nnls, params%vdw_mbd_nfreq, &
                                         2, params%vdw_mbd_cent_appr, &
                                         params%vdw_omega_ref, alpha_SCS, omega_SCS, include_2b, &
                                         params%poly_cut_xmin, params%poly_cut_xmax, &
#ifdef _MPIF90
                                         this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw, &
                                         this_local_virial_vdw_diag )
#else
                                         energies_vdw(i_beg:i_end), forces_vdw, virial_vdw, local_virial_vdw_diag )
#endif
if ( params%vdw_mbd_norder > 2 ) then ! Call 3b+ version if maximum body order > 2, otherwise only 2b with SCS polarizabilities is included
include_2b = .false.
          call get_mbd_energies_and_forces( hirshfeld_v_cart_der_ji(1:3,j_beg:j_end), &
                                         n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                                         neighbor_species(j_beg:j_end), &
                                         params%vdw_scs_rcut, params%vdw_loc_rcut, params%vdw_mbd_rcut, &
                                         params%vdw_mbd_rcut2, &
                                         params%vdw_buffer, &
                                         rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                                         params%vdw_sr_mbd, params%vdw_d_mbd, params%vdw_c6_ref, params%vdw_r0_ref, &
                                         params%vdw_alpha0_ref, params%vdw_mbd_grad, params%vdw_hirsh_grad, &
                                         params%vdw_polynomial, params%do_nnls, params%vdw_mbd_nfreq, &
                                         params%vdw_mbd_norder, params%vdw_mbd_cent_appr, &
                                         params%vdw_omega_ref, alpha_SCS, omega_SCS, include_2b, &
                                         params%poly_cut_xmin, params%poly_cut_xmax, &
#ifdef _MPIF90
                                         this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw, &
                                         this_local_virial_vdw_diag )
#else
                                         energies_vdw(i_beg:i_end), forces_vdw, virial_vdw, local_virial_vdw_diag )
#endif
end if
else ! If 2b cut-off <= mbd cut-off, just execute normally with 2b and 3b+ for same cut-off
include_2b = .true.
          call get_mbd_energies_and_forces( hirshfeld_v_cart_der_ji(1:3,j_beg:j_end), &
                                         n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                                         neighbor_species(j_beg:j_end), &
                                         params%vdw_scs_rcut, params%vdw_loc_rcut, params%vdw_mbd_rcut, &
                                         params%vdw_mbd_rcut2, &
                                         params%vdw_buffer, &
                                         rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                                         params%vdw_sr_mbd, params%vdw_d_mbd, params%vdw_c6_ref, params%vdw_r0_ref, &
                                         params%vdw_alpha0_ref, params%vdw_mbd_grad, params%vdw_hirsh_grad, &
                                         params%vdw_polynomial, params%do_nnls, params%vdw_mbd_nfreq, &
                                         params%vdw_mbd_norder, params%vdw_mbd_cent_appr, &
                                         params%vdw_omega_ref, alpha_SCS, omega_SCS, include_2b, &
                                         params%poly_cut_xmin, params%poly_cut_xmax, &
#ifdef _MPIF90
                                         this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw, &
                                         this_local_virial_vdw_diag )
#else
                                         energies_vdw(i_beg:i_end), forces_vdw, virial_vdw, local_virial_vdw_diag )
#endif
end if
call get_time(time2)
!write(*,*) "MBD timing", time2-time1

!        call get_ts_energy_and_forces( hirshfeld_v(i_beg:i_end), hirshfeld_v_cart_der(1:3, j_beg:j_end), &
!                                       n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
!                                       neighbor_species(j_beg:j_end), &
!                                       params%vdw_rcut, params%vdw_buffer, &
!                                       params%vdw_rcut_inner, params%vdw_buffer_inner, &
!                                       rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
!                                       params%vdw_sr, params%vdw_d, params%vdw_c6_ref, params%vdw_r0_ref, &
!                                       params%vdw_alpha0_ref, c6_scs, r0_scs, alpha0_scs, params%do_forces, &
!#ifdef _MPIF90
!                                       this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw )
!#else
!                                       energies_vdw(i_beg:i_end), forces_vdw, virial_vdw )
!#endif
!write(*,*) "TS_energies", sum(energies)

          !call get_mbd_energies_and_forces( alpha_SCS, alpha_SCS_grad, n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
          !              neighbor_species(j_beg:j_end), params%vdw_rcut, params%vdw_buffer, &
          !              params%vdw_rcut_inner, params%vdw_buffer_inner, rjs(j_beg:j_end), &
          !              xyz(1:3, j_beg:j_end), params%vdw_sr, params%vdw_d, &
          !              params%vdw_c6_ref, params%vdw_r0_ref, params%vdw_alpha0_ref, &
          !              params%vdw_mbd_grad, energies_vdw(i_beg:i_end), forces_vdw, virial_vdw )


          !write(*,*) "vdw forces"
          !do i = 1, n_sites
          !  write(*,*) forces_vdw(1:3,i), this_forces_vdw(1:3,i), energies_vdw(i), this_energies_vdw(i)
          !end do
          deallocate(alpha_SCS, omega_SCS, alpha_SCS_grad, c6_scs, r0_scs, alpha0_scs, &
                     this_alpha_SCS, this_omega_SCS)
          if( .not. (params%vdw_type == "ts+mbd" .and. is_correction_step) )deallocate(v_neigh_vdw) 
        end if
!write(*,*) "This local virial", this_local_virial_vdw_diag
!       Apply the MBD-TS correction if needed/requested
        if( params%vdw_type == "ts+mbd" )then
if( allocated( S_xyz_inv ) )then
deallocate( S_xyz_inv )
end if
allocate( S_xyz_inv(1:3, 1:n_sites) )
S_xyz_inv = 0.d0
k = j_beg - 1
do i = i_beg, i_end
k = k + 1
n = 0
do j = 2, n_neigh(i) ! exclude self term
k = k + 1
if( rjs(k) < params%vdw_rcut )then ! this should be the TS rcut, make sure!
S_xyz_inv(1:3, i) = S_xyz_inv(1:3, i) + screened_one_over_x(xyz(1:3, k))
n = n + 1
end if
end do
if( n >= 1 )then
S_xyz_inv(1:3, i) = S_xyz_inv(1:3, i) / dfloat(n)
end if
end do
!write(*,*) rank, state%this_local_virial_vdw_diag_corr * S_xyz_inv
!         This updates the correction every mbd_correction_freq steps
          if( is_correction_step )then
#ifdef _MPIF90
!            this_mbd_ts_scaling = 1.d0 + (dabs(this_energies_vdw) - dabs(state%this_energies_vdw_corr)) &
!                                  / (dabs(state%this_energies_vdw_corr) + 0.01d0)
!            this_mbd_ts_scaling = 1.d0
!            state%this_energies_vdw_corr(i_beg:i_end) = this_energies_vdw(i_beg:i_end) - state%this_energies_vdw_corr(i_beg:i_end)
!            state%this_energies_vdw_corr = this_energies_vdw - state%this_energies_vdw_corr
            state%this_energies_vdw_corr = this_energies_vdw ! we do this when scaling is used in TS because TS is run twice
!            state%this_forces_vdw_corr = this_forces_vdw - state%this_forces_vdw_corr
!            state%this_virial_vdw_corr = this_virial_vdw - state%this_virial_vdw_corr
            state%this_virial_vdw_corr = this_virial_vdw ! we do this when scaling is used in TS because TS is run twice
!            state%this_local_virial_vdw_diag_corr = this_local_virial_vdw_diag - state%this_local_virial_vdw_diag_corr
            call mpi_reduce(state%this_local_virial_vdw_diag_corr, local_virial_vdw_diag_corr, 3*n_sites, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call mpi_reduce(this_local_virial_vdw_diag, local_virial_vdw_diag, 3*n_sites, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            if( rank == 0 )then
              if( update_mbd_ts_scaling )then
!                this_mbd_ts_scaling = smooth_ratio(sum(local_virial_vdw_diag, 1), sum(local_virial_vdw_diag_corr, 1), 0.1d0, 3.d0)
                this_mbd_ts_scaling = this_mbd_ts_scaling + 0.1d0 * ( &
                                      smooth_ratio(sum(local_virial_vdw_diag, 1), sum(local_virial_vdw_diag_corr, 1), -1.d0, 3.d0) &
                                      - 1.d0 )
                this_mbd_ts_scaling = clip(this_mbd_ts_scaling, 0.1d0, 2.5d0)
              end if
              update_mbd_ts_scaling = .true.
!write(*,*) "mbd_ts_scaling:", this_mbd_ts_scaling
!write(*,*)
!write(*,*) "Tr(virial):", md_istep, sum(local_virial_vdw_diag), sum(local_virial_vdw_diag_corr), &
!sum(this_mbd_ts_scaling)/size(this_mbd_ts_scaling), std(this_mbd_ts_scaling)!, local_virial_vdw_diag/local_virial_vdw_diag_corr
!sum(this_mbd_ts_scaling*sum(local_virial_vdw_diag_corr,1))
            end if
            call mpi_bcast(this_mbd_ts_scaling, n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
!            mbd_ts_scaling = 1.d0 + (dabs(this_energies_vdw) - dabs(state%this_energies_vdw_corr)) &
!                                  / (dabs(state%this_energies_vdw_corr) + 0.01d0)
!            mbd_ts_scaling = 1.d0
!            energies_vdw_corr(i_beg:i_end) = energies_vdw(i_beg:i_end) - energies_vdw_corr(i_beg:i_end)
!            energies_vdw_corr = energies_vdw - energies_vdw_corr
            energies_vdw_corr = energies_vdw ! we do this when scaling is used in TS because TS is run twice
!            forces_vdw_corr = forces_vdw - forces_vdw_corr
!            state%virial_vdw_corr = virial_vdw - state%virial_vdw_corr
            state%virial_vdw_corr = virial_vdw ! we do this when scaling is used in TS because TS is run twice
!            local_virial_vdw_diag_corr = local_virial_vdw_diag - local_virial_vdw_diag_corr
            if( update_mbd_ts_scaling )then
!              mbd_ts_scaling = smooth_ratio(sum(local_virial_vdw_diag, 1), sum(local_virial_vdw_diag_corr, 1), 0.1d0, 3.d0)
              mbd_ts_scaling = mbd_ts_scaling + 0.1d0 * ( &
                               smooth_ratio(sum(local_virial_vdw_diag, 1), sum(local_virial_vdw_diag_corr, 1), -1.d0, 3.d0) &
                               - 1.d0 )
              mbd_ts_scaling = clip(mbd_ts_scaling, 0.1d0, 2.5d0)
            end if
            update_mbd_ts_scaling = .true.
#endif
            call get_ts_energy_and_forces( local_properties(i_beg:i_end, vdw_lp_index), &
                & local_properties_cart_der(1:3, j_beg:j_end, vdw_lp_index), &
                n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                neighbor_species(j_beg:j_end), &
                params%vdw_rcut, params%vdw_buffer, &
                params%vdw_rcut_inner, params%vdw_buffer_inner, &
                rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
                params%vdw_sr, params%vdw_d, params%vdw_c6_ref, params%vdw_r0_ref, &
                params%vdw_alpha0_ref, params%do_forces, &
                params%poly_cut_xmin, params%poly_cut_xmax, &
#ifdef _MPIF90
                this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw, &
                this_local_virial_vdw_diag, this_mbd_ts_scaling )
#else
                energies_vdw(i_beg:i_end), forces_vdw, virial_vdw, local_virial_vdw_diag, &
                mbd_ts_scaling )
#endif
#ifdef _MPIF90
            state%this_energies_vdw_corr = -(this_energies_vdw - state%this_energies_vdw_corr) ! for scaling with TS
            this_energies_vdw = this_energies_vdw + state%this_energies_vdw_corr
            state%this_virial_vdw_corr = -(this_virial_vdw - state%this_virial_vdw_corr) ! for scaling with TS
            this_virial_vdw = this_virial_vdw + state%this_virial_vdw_corr
#else
            energies_vdw_corr = -(energies_vdw - energies_vdw_corr) ! for scaling with TS
            energies_vdw = energies_vdw + energies_vdw_corr
            state%virial_vdw_corr = -(virial_vdw - state%virial_vdw_corr) ! for scaling with TS
            virial_vdw = virial_vdw + state%virial_vdw_corr
#endif
           deallocate(v_neigh_vdw)
          else
!         This applies the correction
#ifdef _MPIF90
!            this_energies_vdw(i_beg:i_end) = this_energies_vdw(i_beg:i_end) + state%this_energies_vdw_corr(i_beg:i_end)
            this_energies_vdw = this_energies_vdw + state%this_energies_vdw_corr
!            this_forces_vdw = this_forces_vdw + state%this_forces_vdw_corr
            this_virial_vdw = this_virial_vdw + state%this_virial_vdw_corr
!            this_local_virial_vdw_diag = this_local_virial_vdw_diag + state%this_local_virial_vdw_diag_corr
!            this_forces_vdw = this_forces_vdw - state%this_local_virial_vdw_diag_corr * S_xyz_inv
#else
!            energies_vdw(i_beg:i_end) = energies_vdw(i_beg:i_end) + energies_vdw_corr(i_beg:i_end)
            energies_vdw = energies_vdw + energies_vdw_corr
!            forces_vdw = forces_vdw + forces_vdw_corr
            virial_vdw = virial_vdw + state%virial_vdw_corr
!            local_virial_trace_vdw = local_virial_trace_vdw + local_virial_vdw_diag_corr
!            forces_vdw = forces_vdw - local_virial_vdw_diag_corr * S_xyz_inv
#endif
          end if
        end if

! REMOVE EVERYTHING UNDER THIS WHEN YOU ARE DONE MERGING MBD

!           call get_ts_energy_and_forces( local_properties(i_beg:i_end, vdw_lp_index), &
!                & local_properties_cart_der(1:3, j_beg:j_end, vdw_lp_index), &
!                n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
!                neighbor_species(j_beg:j_end), &
!                params%vdw_rcut, params%vdw_buffer, &
!                params%vdw_rcut_inner, params%vdw_buffer_inner, &
!                rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), v_neigh_vdw, &
!                params%vdw_sr, params%vdw_d, params%vdw_c6_ref, params%vdw_r0_ref, &
!                params%vdw_alpha0_ref, params%do_forces, &
!                params%poly_cut_xmin, params%poly_cut_xmax, &
!#ifdef _MPIF90
!                this_energies_vdw(i_beg:i_end), this_forces_vdw, this_virial_vdw, &
!                this_local_virial_vdw_diag, this_mbd_ts_scaling )
!#else
!                energies_vdw(i_beg:i_end), forces_vdw, virial_vdw, local_virial_vdw_diag, &
!                mbd_ts_scaling )
!#endif
           call get_time(time_vdw(2))
           time_vdw(3) = time_vdw(2) - time_vdw(1)
!           if( .not. (params%vdw_type == "ts+mbd" .and. modulo(md_istep, params%mbd_correction_freq) == 0) )then
!             deallocate(v_neigh_vdw)
!           end if
    end if

  end subroutine compute_vdw

end module turbogap_vdw
