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
   use electrostatics, only: compute_coulomb_direct, compute_coulomb_dsf, compute_coulomb_lamichhane

   implicit none

   private
   public :: compute_estat

contains

   !**************************************************************************
   subroutine compute_estat(params, do_electrostatics, valid_estat_charges, charge_lp_index, &
                            n_sites, n_neigh, neighbors_list, rjs, xyz, &
                            local_properties, local_properties_cart_der, &
                            i_beg, i_end, j_beg, j_end, rank, &
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
      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: i_beg, i_end, j_beg, j_end
      integer, intent(in) :: rank

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
            call time_start(time%estat)
#ifdef _MPIF90
            allocate (energies_estat(1:n_sites))
            energies_estat = 0.d0
            if (params%do_forces) then
               allocate (forces_estat(1:3, 1:n_sites))
               forces_estat = 0.d0
            end if
#endif
            allocate (chg_neigh_estat(1:j_end - j_beg + 1))
            chg_neigh_estat = 0.d0

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

            k = 0
            do i = i_beg, i_end
               do j = 1, n_neigh(i)
                  j2 = mod(neighbors_list(j_beg + k) - 1, n_sites) + 1
                  k = k + 1
                  chg_neigh_estat(k) = local_properties(j2, charge_lp_index)
               end do
            end do
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
               call compute_coulomb_lamichhane( &
                  local_properties(i_beg:i_end, charge_lp_index), &
                  local_properties_cart_der(1:3, j_beg:j_end, charge_lp_index), &
                  n_neigh(i_beg:i_end), neighbors_list(j_beg:j_end), &
                  params%estat_dsf_alpha, params%estat_rcut, &
                  rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), chg_neigh_estat, &
                  params%do_forces, &
                  energies_estat(i_beg:i_end), forces_estat, virial_estat, params%estat_options)
            else
               write (*, *) "WARNING: Unknown electrostatic method "//trim(params%estat_method)
               write (*, *) "Ignoring..."
            end if
            deallocate (chg_neigh_estat)
            call time_end(time%estat)
         end if
      end if

   end subroutine compute_estat
   !**************************************************************************

end module turbogap_estat
