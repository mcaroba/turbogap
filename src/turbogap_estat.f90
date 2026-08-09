! Electrostatics phase of the main loop.
!
! Lifted verbatim out of turbogap.f90. The GPU branch carries the same module
! and the same procedure name; its version additionally routes the gsf method
! through calculate_batched_electrostatics when params%gpu_batched is set,
! which is the one part of electrostatics that is genuinely device code. The
! driver calls compute_estat on both branches and never learns which it got.
!
! Three of the sixteen argument lists that no Fortran-aware tool can parse
! (docs/refactor_strategy.md 6.4) were in this block: an argument list
! interrupted by #ifdef _MPIF90 with a this_-prefixed triple on one side and
! the plain one on the other. Inside a procedure the dummy has one name either
! way, so the caller chooses once and the conditional disappears from the
! moved code. The remaining #ifdef here has no #else -- it guards an MPI-only
! allocation of whole statements, which is an ordinary conditional and must
! survive.
!
! energies_estat and forces_estat are allocatable dummies because under MPI
! this procedure allocates them (the caller passes its this_ buffers) and the
! driver's reduce block releases them, the same split lifetime compute_vdw
! already has for this_local_virial_vdw_diag.

module turbogap_estat

   use kinds
   use types
   use timing
   use iso_c_binding
   use F_B_C
   use gpu_context
   use neighbors, only: get_gpu_batches
   use exp_interface, only: gpu_malloc_neighbors, gpu_free_neighbors, free_exp_batches
   use electrostatics, only: compute_coulomb_direct, compute_coulomb_dsf, &
                             compute_coulomb_lamichhane, calculate_batched_electrostatics

   implicit none

   private
   public :: compute_estat

contains

   !**************************************************************************
   subroutine compute_estat(params, do_electrostatics, valid_estat_charges, charge_lp_index, &
                            n_sites, n_neigh, neighbors_list, species, neighbor_species, rjs, xyz, &
                            local_properties, local_properties_cart_der, &
                            i_beg, i_end, j_beg, j_end, rank, n_omp, &
                            energies_estat, forces_estat, virial_estat, time)
      implicit none

!     Input variables
      type(input_parameters), intent(in) :: params
      logical, intent(in) :: do_electrostatics
      logical, intent(in) :: valid_estat_charges
      integer, intent(in) :: charge_lp_index
      integer, intent(in) :: n_sites
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: neighbors_list(:)
      integer, intent(in) :: species(:)
      integer, intent(in) :: neighbor_species(:)
      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: i_beg, i_end, j_beg, j_end
      integer, intent(in) :: rank
      integer, intent(in) :: n_omp

!     In/out variables
      real(dp), intent(inout) :: local_properties(:, :)
      real(dp), intent(inout) :: local_properties_cart_der(:, :, :)
      real(dp), allocatable, intent(inout) :: energies_estat(:)
      real(dp), allocatable, intent(inout) :: forces_estat(:, :)
      real(dp), intent(inout) :: virial_estat(1:3, 1:3)
      type(times_t), intent(inout) :: time

!     Internal variables
      real(dp), allocatable :: chg_neigh_estat(:)
      real(dp) :: charge_sum
      integer :: i, j, k, j2

!     The batch decomposition and its scratch. Allocated AND deallocated inside
!     this block, so despite the driver using the same names elsewhere they are
!     block-private and become locals here.
      integer, allocatable, target :: i_beg_list(:), i_end_list(:), j_beg_list(:), j_end_list(:)
      integer :: omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end
      integer :: n_sites_temp, n_pairs_temp
      real(dp), allocatable, target :: charges_temp(:)
      type(c_ptr) :: charges_d
      integer(c_size_t) :: st_charges_d

!     valid_estat_charges is part of the guard, not an afterthought: without it a
!     deck that asks for electrostatics against a GAP with no atomic_charge local
!     property indexes local_properties with an uninitialised charge_lp_index and
!     segfaults. Same shape as the has_vdw/has_local_properties defect.
      if (do_electrostatics .and. (trim(params%estat_method) /= "none") &
          .and. params%do_prediction) then
         if (.not. valid_estat_charges) then
            if (rank == 0) then
               write (*, *) "WARNING: estat_method = "//trim(params%estat_method)// &
                  " but the GAP provides no atomic_charge local property."
               write (*, *) "         Skipping the electrostatics contribution."
            end if
         else
            call time_start(time%estat, "estat")
#ifdef _MPIF90
            allocate (energies_estat(1:n_sites))
            energies_estat = 0.d0
            if (params%do_forces) then
               allocate (forces_estat(1:3, 1:n_sites))
               forces_estat = 0.d0
            end if
#endif

            ! First, make sure that the system is charge neutral if it is periodic
            !
            ! This is not strictly correct, as charges should be equilibrated,
            ! but with the models that we have available, this
            ! is not possible

            charge_sum = 0.0_dp
            do k = 1, n_sites
               charge_sum = charge_sum + local_properties(k, charge_lp_index)
            end do
            charge_sum = charge_sum/real(n_sites, dp)

            do k = 1, n_sites
               local_properties(k, charge_lp_index) = local_properties(k, charge_lp_index) - charge_sum
            end do
            allocate (chg_neigh_estat(1:j_end - j_beg + 1))
            chg_neigh_estat = 0.d0
            k = 0
            do i = i_beg, i_end
               do j = 1, n_neigh(i)
                  !           I'm not sure if this is necessary or neighbors_list is already bounded between 1 and n_sites -> CHECK THIS
                  j2 = mod(neighbors_list(j_beg + k) - 1, n_sites) + 1
                  k = k + 1
                  chg_neigh_estat(k) = local_properties(j2, charge_lp_index)
               end do
            end do
            ! Prepare to call electrostatics subroutine!
            if (trim(params%estat_method) == "direct") then
               call compute_coulomb_direct( &
                  local_properties(i_beg:i_end, charge_lp_index), &
                  local_properties_cart_der(1:3, j_beg:j_end, charge_lp_index), &
                  n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                  params%estat_rcut, params%estat_rcut_inner, params%estat_inner_width, &
                  rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), chg_neigh_estat, &
                  params%do_forces, &
                  energies_estat(i_beg:i_end), forces_estat, virial_estat, params%estat_options)
            else if (trim(params%estat_method) == "dsf") then
               call compute_coulomb_dsf( &
                  local_properties(i_beg:i_end, charge_lp_index), &
                  local_properties_cart_der(1:3, j_beg:j_end, charge_lp_index), &
                  n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                  params%estat_dsf_alpha, params%estat_rcut, &
                  params%estat_rcut_inner, params%estat_inner_width, &
                  rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), chg_neigh_estat, &
                  params%do_forces, &
                  energies_estat(i_beg:i_end), forces_estat, virial_estat, params%estat_options)
            else if (trim(params%estat_method) == "gsf") then
               if (params%gpu_batched) then

                  print *, "> Starting GPU electrostatics, rank = ", rank

                  call get_gpu_batches(n_neigh(i_beg:i_end), rjs(j_beg:j_end), params%estat_rcut, &
                                       params%gpu_n_batches, gpu_memory_usage, &
                                       params%gpu_max_batch_size, i_beg_list, i_end_list, j_beg_list, j_end_list)

                  gpu_memory_usage = 0.d0

                  ! do i = 1, size(i_beg_list)
                  !    write(*,'(A,I8,A,I8,A,I8,A,I8,A,I8,A,I8)') "batches, rank ", rank, &
                  !         " i ", i, " i_beg_i = ", i_beg_list(i), " i_end_i = ", i_end_list(i), &
                  !         " j_beg_i = ", j_beg_list(i), " j_end_i = ", j_end_list(i)
                  ! end do

                  allocate (gpu_neigh(1:size(i_beg_list)))
                  allocate (gpu_exp(1:size(i_beg_list)))
                  allocate (gpu_batch_storage(1:size(i_beg_list)))

                  allocate (charges_temp(1:n_sites))

                  charges_temp = local_properties(:, charge_lp_index)
                  ! Before the loop, allocate the charges
                  st_charges_d = c_double*n_sites
                  call gpu_malloc_async(charges_d, st_charges_d, gpu_stream)
                  call cpy_htod(c_loc(charges_temp), charges_d, st_charges_d, gpu_stream)

!                 These directives are commented out, and must stay that way
!                 until the accumulators below are dealt with. This is not an
!                 oversight to be tidied up by uncommenting them.
!
!                 The batches are independent in their device state --
!                 gpu_neigh(i), gpu_exp(i), gpu_batch_storage(i) are all indexed
!                 by i, exactly as in the pdf loop that IS threaded -- and
!                 energies_estat(this_i_beg:this_i_end) is a disjoint slice per
!                 batch. But forces_estat and virial_estat are passed WHOLE to
!                 every batch and accumulated into. Two threads in this loop
!                 would read-modify-write the same array, and the result would
!                 be wrong by an amount that changes run to run.
!
!                 To thread it, forces and virial have to become per-batch and
!                 be summed after the loop -- which is what the pdf path already
!                 does with collect_batched_pair_distribution, and what makes
!                 that one safe to run in parallel.
!
!                 Note also that no test can currently catch a mistake here:
!                 estat_gsf is the only case that reaches this code, it is an
!                 xfail, and it exhausts device memory before finishing.
!
!                 !$OMP PARALLEL DO num_threads(n_omp) DEFAULT(SHARED) &
!                 !$OMP PRIVATE(i, omp_task, this_i_beg, this_i_end, this_j_beg, this_j_end, n_sites_temp, n_pairs_temp)

                  do i = 1, size(i_beg_list)

                     ! Serial today, so this is stream 1; gpu_omp_task() is what
                     ! it becomes when the loop above is threaded.
                     omp_task = gpu_omp_task()

                     this_i_beg = i_beg - 1 + i_beg_list(i)
                     this_i_end = i_beg - 1 + i_end_list(i)
                     this_j_beg = j_beg - 1 + j_beg_list(i)
                     this_j_end = j_beg - 1 + j_end_list(i)

                     n_sites_temp = this_i_end - this_i_beg + 1
                     n_pairs_temp = this_j_end - this_j_beg + 1

                     print *, " "
                     write (*, '(9(A,I4))') &
                        "estat batches---Rank ", rank, " ---Thread ", omp_task, &
                        " / ", n_omp, " i = ", i, " / ", params%gpu_n_batches, &
                        " i_beg = ", this_i_beg, &
                        " i_end = ", this_i_end, &
                        " j_beg = ", this_j_beg, &
                        " j_end = ", this_j_end

                     call gpu_malloc_neighbors(gpu_neigh(i), &
                                               n_sites_temp, n_pairs_temp, &
                                               n_neigh(this_i_beg:this_i_end), &
                                               species(this_i_beg:this_i_end), &
                                               neighbor_species(this_j_beg:this_j_end), &
                                               neighbors_list(this_j_beg:this_j_end), &
                                               rjs(this_j_beg:this_j_end), &
                                               xyz(1:3, this_j_beg:this_j_end), &
                                               gpu_streams(omp_task), rank)

                     call calculate_batched_electrostatics(gpu_exp(i), gpu_batch_storage(i), &
                                                           gpu_neigh(i), n_sites, &
                                                           this_i_beg, this_i_end, this_j_beg, this_j_end, rank, &
                                                           params%estat_rcut, params%estat_dsf_alpha, &
                                                           charges_temp, charges_d, &
                                                           local_properties_cart_der(1:3, this_j_beg:this_j_end, charge_lp_index), &
                                                           params%do_forces, &
                                                           energies_estat(this_i_beg:this_i_end), forces_estat, virial_estat, &
                                                          params%estat_options, params%estat_rcut_inner, params%estat_inner_width, &
                                                           gpu_streams(omp_task) &
                                                           )

                     call gpu_free_neighbors(gpu_neigh(i), gpu_streams(omp_task))

            write (*, '(A,I4,A,I4,A,I4,A,I4,A,I4,A,I4)') "Electrostatics batches finished---Rank ", rank, " ---Thread ", omp_task, &
                        " / ", n_omp, " i = ", i, &
                        " i_beg = ", this_i_beg, &
                        " i_end = ", this_i_end, &
                        " j_beg = ", this_j_beg, &
                        " j_end = ", this_j_end
                  end do
                  !     !$OMP END PARALLEL DO
                  deallocate (gpu_neigh)

                  call gpu_free_async(charges_d, gpu_stream)
                  deallocate (charges_temp)

                  call free_exp_batches(gpu_exp, params%gpu_n_batches)

                  if (allocated(i_beg_list)) deallocate (i_beg_list, i_end_list, j_beg_list, j_end_list)

                  if (allocated(gpu_neigh)) deallocate (gpu_neigh)
                  if (allocated(gpu_exp)) deallocate (gpu_exp)
                  if (allocated(gpu_batch_storage)) deallocate (gpu_batch_storage)
               else
                  call compute_coulomb_lamichhane( &
                     local_properties(i_beg:i_end, charge_lp_index), &
                     local_properties_cart_der(1:3, j_beg:j_end, charge_lp_index), &
                     n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                     params%estat_dsf_alpha, params%estat_rcut, &
                     rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), chg_neigh_estat, &
                     params%do_forces, &
                     energies_estat(i_beg:i_end), forces_estat, virial_estat, params%estat_options)
               end if

            else ! This really shouldn't happen... but we both know it could
               print("WARNING: Unknown electrostatic method "//params%estat_method)
               write (*, *) "Ignoring..."
            end if
            deallocate (chg_neigh_estat)
            call time_end(time%estat, "estat")
         end if
      end if

   end subroutine compute_estat
   !**************************************************************************

end module turbogap_estat
