!**************************************************************************
!  SOAP-kernel metadynamics: a history-dependent bias on local environments
!  
!  Methodological details: to be continued...
!
!  How it plugs into TurboGAP?
!  --------------------------
!  metad_init          Called once after the input has been read.
!  metad_begin_step    Called every MD step before the descriptors are evaluated.
!                      Zeroes this step's outputs, and every `window` steps
!                      picks a new target atom and empties the bias.
!  metad_add_bias      From inside get_gap_soap, straight after the SOAP vectors
!                      and their derivatives are computed. It resycles descriptors
!                      that GAP computes anyway so no new SOAP evaluations needed.
!                      Energies, forces and virial are communicated into res% after
!                      the compute_soap routine outside of it.
!  metad_deposit       Called after the descriptor loop, on all ranks. Gathers the
!                      deposited vectors on rank 0, grows the sparse set, refits
!                      alphas and broadcasts the result. The new bias is felt
!                      from the next step on. As in standard metadynamics, a
!                      hill is added after the forces of the step that deposits
!                      it.
!  metad_write_output  On rank 0, after the energies and forces are reduced.
!
!  MPI layout
!  ----------
!  Every rank holds the same copy of Qs and alphas, which it needs to evaluate
!  the bias on its own atoms. The visit counts and the Gram matrix are only
!  needed for the update, so they live on rank 0 only.
!
!  The energies, forces and virial in metad_state are partial sums on each rank
!  until the contribution reduction in turbogap.f90 (C_METAD) sums them onto
!  rank 0, the same way as energies_soap and forces_soap.
!**************************************************************************

module enhanced_sampling

   use kinds
   use types, only: input_parameters, soap_turbo
   use grow_array
   use error, only: turbogap_abort
#ifdef _MPIF90
   use mpi
#endif

   implicit none
   private

   public :: metad_state
   public :: metad_init
   public :: metad_begin_step
   public :: metad_add_bias
   public :: metad_deposit
   public :: metad_write_output

   real(dp), parameter :: pi = 3.14159265358979323846d0

!  Output files. Both are written by rank 0 only.
   character(len=*), parameter :: log_file = "metadynamics.dat"
   character(len=*), parameter :: forces_file = "metadynamics_forces.dat"

!  One bias sparse set, attached to one soap_turbo descriptor.
!  metad_state can own multiple sparse sets.
   type metad_set
!     True if this descriptor is biased for the current target atom.
      logical :: active = .false.
!     Length of this descriptor's SOAP vector.
      integer :: n_soap = 0
!     Number of sparse points in use. The arrays below are usually larger,
!     because their capacity grows by doubling (see grow_array.f90).
      integer :: n_sparse = 0
!     The sparse points Qs(1:n_soap, 1:n_sparse) and their kernel weights
!     alphas(1:n_sparse). Every rank holds the same copy.
      real(dp), allocatable :: Qs(:, :)
      real(dp), allocatable :: alphas(:)
!     Rank 0 only. counts(k) is N_k, the number of visits to sparse point k.
!     gram(k, l) = Q_k . Q_l grows by one row and column per new point, so a
!     refit never recomputes old dot products.
      integer, allocatable :: counts(:)
      real(dp), allocatable :: gram(:, :)
!     SOAP vectors this rank collected during the current step that are waiting
!     to be deposited: cand(1:n_soap, 1:n_cand).
      integer :: n_cand = 0
      real(dp), allocatable :: cand(:, :)
   end type metad_set

!  Metadynamics state
   type metad_state
      logical :: on = .false.          ! master switch

!     Input parameters 
      real(dp) :: xi = 2.d0            ! height scale [eV]
      real(dp) :: delta_T = 1.d0       ! well-tempering: larger -> heights saturate later
      real(dp) :: threshold = 0.99d0   ! q_i . q_j above this means "the same environment"
      real(dp) :: reg = 1.d-6          ! diagonal regularisation of the kernel matrix
      integer :: zeta = 10             ! kernel exponent
      integer :: stride = 10           ! deposit every `stride` MD steps
      integer :: window = 400          ! MD steps spent on one target atom

!     The current target atom (global index, 0 until the first one is drawn)
!     and the MD step at which it was drawn.
      integer :: center = 0
      integer :: reset_step = 0
!     True if the environments seen in this step are deposited at its end.
      logical :: deposit_now = .false.

!     One sparse set per soap_turbo descriptor.
      type(metad_set), allocatable :: sets(:)

!     This step's bias contributions. They are partial sums on each rank until
!     the main program's MPI reduction, and totals on rank 0 after it.
      real(dp), allocatable :: energies(:)      ! (1:n_sites)
      real(dp), allocatable :: forces(:, :)     ! (1:3, 1:n_sites)
      real(dp) :: virial(1:3, 1:3) = 0.d0
   end type metad_state

contains

!  Setup metadynamics stuff
   subroutine metad_init(metad, params, soap_turbo_hypers, n_soap_turbo, rank)
      implicit none
      type(metad_state), intent(inout) :: metad
      type(input_parameters), intent(in) :: params
      type(soap_turbo), allocatable, intent(in) :: soap_turbo_hypers(:)
      integer, intent(in) :: n_soap_turbo
      integer, intent(in) :: rank

      integer :: i

      metad%on = params%do_metadynamics
      if (.not. metad%on) return

!     The bias is driven by md_istep and needs SOAP derivatives, so it only
!     makes sense in a plain MD run.
      if (.not. params%do_md .or. params%do_mc .or. params%do_nested_sampling) then
         if (rank == 0) write (*, *) "ERROR: do_metadynamics needs a plain MD run ", &
            "(turbogap md, without MC or nested sampling)"
         call turbogap_abort()
      end if
      if (.not. params%do_forces) then
         if (rank == 0) write (*, *) "ERROR: do_metadynamics needs do_forces = .true."
         call turbogap_abort()
      end if
      if (n_soap_turbo < 1) then
         if (rank == 0) write (*, *) "ERROR: do_metadynamics needs at least one soap_turbo descriptor"
         call turbogap_abort()
      end if
      if (params%metadynamics_stride < 1 .or. params%metadynamics_window < 1 .or. &
          params%metadynamics_zeta < 1 .or. params%metadynamics_dT <= 0.d0) then
         if (rank == 0) write (*, *) "ERROR: metadynamics_stride, metadynamics_window, ", &
            "metadynamics_zeta and metadynamics_dT must be positive"
         call turbogap_abort()
      end if

      metad%xi = params%metadynamics_xi
      metad%delta_T = params%metadynamics_dT
      metad%threshold = params%metadynamics_threshold
      metad%reg = params%metadynamics_reg
      metad%zeta = params%metadynamics_zeta
      metad%stride = params%metadynamics_stride
      metad%window = params%metadynamics_window

      allocate (metad%sets(1:n_soap_turbo))
      do i = 1, n_soap_turbo
         metad%sets(i)%n_soap = soap_turbo_hypers(i)%dim
      end do

   end subroutine metad_init

!  Begin metadynamics (called every md step before the descriptors are calculated)
   subroutine metad_begin_step(metad, md_istep, n_sites, xyz_species, soap_turbo_hypers, rank)
      implicit none
      type(metad_state), intent(inout) :: metad
      integer, intent(in) :: md_istep
      integer, intent(in) :: n_sites
      character*8, intent(in) :: xyz_species(:)
      type(soap_turbo), intent(in) :: soap_turbo_hypers(:)
      integer, intent(in) :: rank

      integer :: i

      if (.not. metad%on) return

!     This step's outputs
      if (allocated(metad%energies)) then
         if (size(metad%energies) /= n_sites) deallocate (metad%energies, metad%forces)
      end if
      if (.not. allocated(metad%energies)) then
         allocate (metad%energies(1:n_sites))
         allocate (metad%forces(1:3, 1:n_sites))
      end if
      metad%energies = 0.d0
      metad%forces = 0.d0
      metad%virial = 0.d0

!     Time for a new target atom?
      if (metad%center == 0 .or. md_istep - metad%reset_step >= metad%window) then
         call new_target(metad, md_istep, n_sites, xyz_species, soap_turbo_hypers, rank)
      end if

!     Deposits are counted from the moment the current target was chosen, so
!     the first step on every target always deposits.
      metad%deposit_now = (modulo(md_istep - metad%reset_step, metad%stride) == 0)

      do i = 1, size(metad%sets)
         metad%sets(i)%n_cand = 0
      end do

   end subroutine metad_begin_step

!  Draw a new target atom on rank 0, share it, and empty every sparse set.
!  The memory is kept for reuse.
   subroutine new_target(metad, md_istep, n_sites, xyz_species, soap_turbo_hypers, rank)
      implicit none
      type(metad_state), intent(inout) :: metad
      integer, intent(in) :: md_istep
      integer, intent(in) :: n_sites
      character*8, intent(in) :: xyz_species(:)
      type(soap_turbo), intent(in) :: soap_turbo_hypers(:)
      integer, intent(in) :: rank

      real(dp) :: u
      integer :: i
#ifdef _MPIF90
      integer :: ierr
#endif

!     Uniformly random over all atoms. random_number is seeded per rank in
!     turbogap_setup, so only rank 0 draws and the others are told.
      if (rank == 0) then
         call random_number(u)
         metad%center = min(1 + int(u*dfloat(n_sites)), n_sites)
      end if
#ifdef _MPIF90
      call mpi_bcast(metad%center, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
      metad%reset_step = md_istep

      do i = 1, size(metad%sets)
         metad%sets(i)%n_sparse = 0
         metad%sets(i)%active = is_biased(soap_turbo_hypers(i), xyz_species(metad%center))
      end do

   end subroutine new_target

   logical function is_biased(hypers, species)
      implicit none
      type(soap_turbo), intent(in) :: hypers
      character*8, intent(in) :: species

      if (hypers%is_dipole_model) then
         is_biased = .false.
      else if (hypers%central_species == 0) then
         is_biased = .true.
      else
         is_biased = (hypers%species_types(hypers%central_species) == species)
      end if
   end function is_biased

!  Main workhorse of metadynamics. Energy, force and virial contributions of the bias 
!  are calculated here and passed into metad%energies, metad%forces etc
   subroutine metad_add_bias(metad, i_set, n_sites, n_neigh, neighbors_list, rjs, xyz, &
                             in_to_out_site, rcut, soap, soap_cart_der)
      implicit none
      type(metad_state), intent(inout) :: metad
      integer, intent(in) :: n_sites               ! central sites in this call
      integer, intent(in) :: i_set
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: neighbors_list(:)
      integer, intent(in) :: in_to_out_site(:)
      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(in) :: rcut                 ! descriptor cutoff = bias window radius
      real(dp), intent(in) :: soap(:, :)           ! (1:n_soap, 1:n_sites)
      real(dp), intent(in) :: soap_cart_der(:, :, :)

!     Per selected site m: its internal site index, its first pair, the pair
!     that points at the target atom, w(r_cm) and dw/dr(r_cm).
      integer, allocatable :: sel(:), sel_pair_beg(:), sel_pair_c(:)
      real(dp), allocatable :: w(:), dw(:)
!     Work arrays for the kernel evaluation.
      real(dp), allocatable :: q(:, :), kern(:, :), kern_der(:, :), grad(:, :), e_site(:)
      real(dp) :: r_min, f(1:3), g(1:3), d(1:3)
      integer :: i, k, l, l_min, m, n_soap, n_sparse, i_glob, j_glob

      if (.not. metad%on) return
      if (.not. metad%sets(i_set)%active) return
      if (n_sites == 0) return

      n_soap = size(soap, 1)

!     -----------------------------------------------------------------------
!     1) Find the sites that lie within rcut of the target atom.
!        This is decided from each site's own neighbour list (the pair list is
!        symmetric), so it needs no communication. If several periodic images
!        of the target are in range, the closest one sets the weight.
!     -----------------------------------------------------------------------
      allocate (sel(1:n_sites), sel_pair_beg(1:n_sites), sel_pair_c(1:n_sites))
      allocate (w(1:n_sites), dw(1:n_sites))
      m = 0
      k = 0                                   ! pairs of sites before site i
      do i = 1, n_sites
         r_min = rcut
         l_min = 0
         do l = k + 1, k + n_neigh(i)
            if (in_to_out_site(neighbors_list(l)) == metad%center .and. rjs(l) < r_min) then
               r_min = rjs(l)
               l_min = l
            end if
         end do
         if (l_min > 0) then
            m = m + 1
            sel(m) = i
            sel_pair_beg(m) = k + 1
            sel_pair_c(m) = l_min
            w(m) = window_fn(r_min, rcut)
            dw(m) = window_der(r_min, rcut)
         end if
         k = k + n_neigh(i)
      end do
      if (m == 0) return

!     -----------------------------------------------------------------------
!     2) Keep these environments for the deposit at the end of the step.
!     -----------------------------------------------------------------------
      if (metad%deposit_now) then
         do i = 1, m
            call push_vector(metad%sets(i_set)%cand, metad%sets(i_set)%n_cand, soap(1:n_soap, sel(i)))
         end do
      end if

!     -----------------------------------------------------------------------
!     3) Evaluate the bias with the current sparse set. The same algebra as
!        get_soap_energy_and_forces in gap.f90, done on the m selected sites
!        only.
!     -----------------------------------------------------------------------
      n_sparse = metad%sets(i_set)%n_sparse
      if (n_sparse == 0) return

      associate (Qs => metad%sets(i_set)%Qs, alphas => metad%sets(i_set)%alphas, zeta => metad%zeta)

         allocate (q(1:n_soap, 1:m))
         allocate (kern(1:m, 1:n_sparse), kern_der(1:m, 1:n_sparse))
         allocate (grad(1:n_soap, 1:m), e_site(1:m))

!        Pack the selected SOAP vectors so that BLAS can work on them.
         do i = 1, m
            q(1:n_soap, i) = soap(1:n_soap, sel(i))
         end do

!        kern(i, k) = q_i . Q_k
         call dgemm("t", "n", m, n_sparse, n_soap, 1.d0, q, n_soap, Qs, size(Qs, 1), 0.d0, kern, m)

!        kern_der(i, k) = alpha_k (q_i . Q_k)^(zeta-1)
!        Then E(q_i) = sum_k kern_der(i,k) kern(i,k)
!        and dE(q_i)/dq_i = zeta sum_k kern_der(i,k) Q_k
         do k = 1, n_sparse
            kern_der(1:m, k) = alphas(k)*kern(1:m, k)**(zeta - 1)
         end do
         e_site(1:m) = sum(kern_der*kern, dim=2)

!        grad(:, i) = -zeta sum_k kern_der(i,k) Q_k = -dE(q_i)/dq_i ...
         call dgemm("n", "t", n_soap, m, n_sparse, -dfloat(zeta), Qs, size(Qs, 1), kern_der, m, &
                    0.d0, grad, n_soap)
!        ... times the window weight, so it is -d[w_i E(q_i)]/dq_i at fixed w_i.
         do i = 1, m
            grad(1:n_soap, i) = w(i)*grad(1:n_soap, i)
         end do

      end associate

!     -----------------------------------------------------------------------
!     4) Accumulate energies, forces and virial.
!     -----------------------------------------------------------------------
      do i = 1, m
         i_glob = in_to_out_site(sel(i))
         metad%energies(i_glob) = metad%energies(i_glob) + w(i)*e_site(i)

!        (a) Through the SOAP vector q_i: every atom j in site i's
!            environment feels  f_j = -w_i dE(q_i)/dq_i . dq_i/dr_j
         do l = sel_pair_beg(i), sel_pair_beg(i) + n_neigh(sel(i)) - 1
            j_glob = in_to_out_site(neighbors_list(l))
            f(1:3) = matmul(soap_cart_der(1:3, 1:n_soap, l), grad(1:n_soap, i))
            metad%forces(1:3, j_glob) = metad%forces(1:3, j_glob) + f(1:3)
            call add_to_virial(metad%virial, f, xyz(1:3, l))
         end do

!        (b) Through the window weight w(r_ci), with d = r_c - r_i:
!            d[w E_i]/dr_c = E_i w'(r) d/r = g,  so  f_c = -g  and  f_i = +g.
!            The target site itself (r = 0) has w' = 0 and is skipped.
         l = sel_pair_c(i)
         if (rjs(l) > 0.d0) then
            d(1:3) = xyz(1:3, l)
            g(1:3) = e_site(i)*dw(i)*d(1:3)/rjs(l)
            metad%forces(1:3, metad%center) = metad%forces(1:3, metad%center) - g(1:3)
            metad%forces(1:3, i_glob) = metad%forces(1:3, i_glob) + g(1:3)
!           The force on the neighbour (the target) times r_c - r_i, the same
!           convention as the GAP virial.
            call add_to_virial(metad%virial, -g, d)
         end if
      end do

   end subroutine metad_add_bias

!  Deposit this step's environments into the bias. 
   subroutine metad_deposit(metad, rank, ntasks)
      implicit none
      type(metad_state), intent(inout) :: metad
      integer, intent(in) :: rank
      integer, intent(in) :: ntasks

      real(dp), allocatable :: cand_all(:, :)
      integer :: i_set, n_cand_all, n_old

      if (.not. metad%on) return
      if (.not. metad%deposit_now) return

      do i_set = 1, size(metad%sets)
!        metad%sets(:)%active is the same on every rank, so all ranks take
!        the same branches and the collectives below match up.
         if (.not. metad%sets(i_set)%active) cycle

         n_old = metad%sets(i_set)%n_sparse
         call gather_candidates(metad%sets(i_set), rank, ntasks, cand_all, n_cand_all)
         if (rank == 0) call grow_sparse_set(metad, i_set, cand_all, n_cand_all)
         call broadcast_sparse_set(metad%sets(i_set), n_old, rank)
      end do

   end subroutine metad_deposit

   subroutine gather_candidates(set, rank, ntasks, cand_all, n_cand_all)
      implicit none
      type(metad_set), intent(inout) :: set
      integer, intent(in) :: rank
      integer, intent(in) :: ntasks
      real(dp), allocatable, intent(out) :: cand_all(:, :)
      integer, intent(out) :: n_cand_all

      integer :: n_soap
#ifdef _MPIF90
      integer, allocatable :: n_by_rank(:), recvcounts(:), displs(:)
      integer :: k, ierr
#endif

      n_soap = set%n_soap

#ifdef _MPIF90
!     Make sure the send buffer exists even on a rank that collected nothing.
      call ensure_capacity_2d(set%cand, n_soap, 1)

      allocate (n_by_rank(1:ntasks), recvcounts(1:ntasks), displs(1:ntasks))
      call mpi_gather(set%n_cand, 1, MPI_INTEGER, n_by_rank, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (rank == 0) then
         n_cand_all = sum(n_by_rank)
         recvcounts = n_soap*n_by_rank
         displs(1) = 0
         do k = 2, ntasks
            displs(k) = displs(k - 1) + recvcounts(k - 1)
         end do
         allocate (cand_all(1:n_soap, 1:max(n_cand_all, 1)))
      else
         n_cand_all = 0
         allocate (cand_all(1:n_soap, 1:1))           ! unused receive buffer
      end if

      call mpi_gatherv(set%cand, n_soap*set%n_cand, MPI_DOUBLE_PRECISION, &
                       cand_all, recvcounts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
      n_cand_all = set%n_cand
      allocate (cand_all(1:n_soap, 1:max(n_cand_all, 1)))
      if (n_cand_all > 0) cand_all(1:n_soap, 1:n_cand_all) = set%cand(1:n_soap, 1:n_cand_all)
#endif

   end subroutine gather_candidates


   subroutine grow_sparse_set(metad, i_set, cand, n_cand)
      implicit none
      type(metad_state), intent(inout) :: metad
      integer, intent(in) :: i_set
      real(dp), intent(in) :: cand(:, :)          ! (1:n_soap, 1:n_cand)
      integer, intent(in) :: n_cand

      real(dp), allocatable :: overlap(:, :), sims(:)
      logical, allocatable :: keep(:)
      integer :: i, j, n, n_soap

      if (n_cand == 0) return

      associate (set => metad%sets(i_set))

         n_soap = set%n_soap

!        a) Drop near-duplicates within this deposit, with the rule of
!           clean_desc in the Python prototype: candidate i goes if any later
!           candidate is more similar to it than the threshold.
         allocate (overlap(1:n_cand, 1:n_cand), keep(1:n_cand))
         call dgemm("t", "n", n_cand, n_cand, n_soap, 1.d0, cand, n_soap, cand, n_soap, 0.d0, overlap, n_cand)
         keep = .true.
         do i = 1, n_cand
            do j = i + 1, n_cand
               if (overlap(i, j) > metad%threshold) then
                  keep(i) = .false.
                  exit
               end if
            end do
         end do

!        b) Each surviving candidate either counts as another visit to the
!           closest sparse point, if that one is similar enough, or becomes a
!           new sparse point. Candidates added earlier in this loop count as
!           sparse points for the later ones.
         allocate (sims(1:set%n_sparse + n_cand))
         do i = 1, n_cand
            if (.not. keep(i)) cycle
            n = set%n_sparse
            if (n > 0) then
!              sims(k) = Q_k . q_i for all sparse points k
               call dgemv("t", n_soap, n, 1.d0, set%Qs, size(set%Qs, 1), cand(1:n_soap, i), 1, &
                          0.d0, sims, 1)
               j = maxloc(sims(1:n), dim=1)
               if (sims(j) >= metad%threshold) then
                  set%counts(j) = set%counts(j) + 1
                  cycle
               end if
            end if
            call append_sparse_point(set, cand(1:n_soap, i), sims(1:n))
         end do

      end associate

!     c) Visit counts have changed, and so have the heights: refit.
      call refit_alphas(metad, i_set)

   end subroutine grow_sparse_set

!  Add point to sparse set
   subroutine append_sparse_point(set, q, sims)
      implicit none
      type(metad_set), intent(inout) :: set
      real(dp), intent(in) :: q(:)
      real(dp), intent(in) :: sims(:)            ! (1:n_sparse)

      integer :: n

      n = set%n_sparse + 1
      call ensure_capacity_2d(set%Qs, set%n_soap, n)
      call ensure_capacity_1d_int(set%counts, n)
      call ensure_capacity_square(set%gram, n)

      set%Qs(1:set%n_soap, n) = q
      set%counts(n) = 1
      set%gram(1:n - 1, n) = sims(1:n - 1)
      set%gram(n, 1:n - 1) = sims(1:n - 1)
      set%gram(n, n) = dot_product(q, q)
      set%n_sparse = n
   end subroutine append_sparse_point

!  Invert the kernel matrix and solve alphas
   subroutine refit_alphas(metad, i_set)
      implicit none
      type(metad_state), intent(inout) :: metad
      integer, intent(in) :: i_set

      real(dp), allocatable :: kmat(:, :), rhs(:)
      integer, allocatable :: ipiv(:)
      integer :: n, k, info

      associate (set => metad%sets(i_set))

         n = set%n_sparse
         allocate (kmat(1:n, 1:n), rhs(1:n))

         call build_kernel_matrix()
         rhs(1:n) = bias_height(set%counts(1:n), metad%xi, metad%delta_T)

!        K + reg I is symmetric positive definite, so Cholesky (dposv) is the
!        cheap way to solve it...
         call dposv("u", n, 1, kmat, n, rhs, n, info)

!        ... unless reg is so small that rounding makes K indefinite. Then
!        fall back to LU, which only needs K to be non-singular.
         if (info /= 0) then
            allocate (ipiv(1:n))
            call build_kernel_matrix()
            rhs(1:n) = bias_height(set%counts(1:n), metad%xi, metad%delta_T)
            call dgesv(n, 1, kmat, n, ipiv, rhs, n, info)
            if (info /= 0) then
               write (*, *) "WARNING: metadynamics kernel matrix is singular; ", &
                  "try a larger metadynamics_reg"
            end if
         end if

         call ensure_capacity_1d(set%alphas, n)
         set%alphas(1:n) = rhs(1:n)

      end associate

   contains

      subroutine build_kernel_matrix()
         kmat(1:n, 1:n) = metad%sets(i_set)%gram(1:n, 1:n)**metad%zeta
         do k = 1, n
            kmat(k, k) = kmat(k, k) + metad%reg
         end do
      end subroutine build_kernel_matrix

   end subroutine refit_alphas

!  Send the updated sparse set from rank 0 to all ranks. Only the newly
!  appended points travel; the old ones are already everywhere.
   subroutine broadcast_sparse_set(set, n_old, rank)
      implicit none
      type(metad_set), intent(inout) :: set
      integer, intent(in) :: n_old             ! n_sparse before this deposit (same on every rank)
      integer, intent(in) :: rank

#ifdef _MPIF90
      integer :: n, ierr

      call mpi_bcast(set%n_sparse, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      n = set%n_sparse
      if (n == 0) return

      if (rank /= 0) then
         call ensure_capacity_2d(set%Qs, set%n_soap, n)
         call ensure_capacity_1d(set%alphas, n)
      end if
      if (n > n_old) then
         call mpi_bcast(set%Qs(1:set%n_soap, n_old + 1:n), set%n_soap*(n - n_old), &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if
!     Every alpha changes on every refit.
      call mpi_bcast(set%alphas(1:n), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

   end subroutine broadcast_sparse_set


!  Rank 0 only, after the energies and forces have been reduced.
!  Write files metadynamics.dat & metadynamics_forces.dat
   subroutine metad_write_output(metad, md_istep, rank, print_forces)
      implicit none
      type(metad_state), intent(in) :: metad
      integer, intent(in) :: md_istep
      integer, intent(in) :: rank
      logical, intent(in) :: print_forces

      integer :: unit, i

      if (.not. metad%on) return
      if (rank /= 0) return

      if (metad%deposit_now) then
         if (md_istep == 0) then
            open (newunit=unit, file=log_file, status="replace")
            write (unit, "(A)", advance="no") "# md_step  target_atom  E_bias[eV]"
            do i = 1, size(metad%sets)
               write (unit, "(A,I0,A)", advance="no") "  n_sparse(", i, ")"
            end do
            write (unit, *)
         else
            open (newunit=unit, file=log_file, status="old", position="append")
         end if
         write (unit, "(I10, 1X, I10, 1X, ES20.10, *(1X, I10))") md_istep, metad%center, &
            sum(metad%energies), (metad%sets(i)%n_sparse, i=1, size(metad%sets))
         close (unit)
      end if

      if (print_forces) then
         open (newunit=unit, file=forces_file, status="replace")
         write (unit, "(A, I0, A, I0)") "# md_step ", md_istep, "  target_atom ", metad%center
         write (unit, "(A)") "# E_bias[eV]  F_x  F_y  F_z [eV/A]"
         do i = 1, size(metad%energies)
            write (unit, "(4(1X, ES20.10))") metad%energies(i), metad%forces(1:3, i)
         end do
         close (unit)
      end if

   end subroutine metad_write_output

!  -----------------
!  Small helpers :=)
!  -----------------


!  Well-tempered height of a sparse point visited n times.
   real(dp) elemental function bias_height(n, xi, delta_T) result(h)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: xi, delta_T
      h = xi*log(1.d0 + dfloat(n)/delta_T)
   end function bias_height

!  Weight of a site at distance r from the target: 1 at r = 0, 0 at r = rcut.
   real(dp) elemental function window_fn(r, rcut) result(w)
      implicit none
      real(dp), intent(in) :: r, rcut
      if (r >= rcut) then
         w = 0.d0
      else
         w = 0.5d0*(1.d0 + cos(pi*r/rcut))
      end if
   end function window_fn

!  dw/dr of window_fn.
   real(dp) elemental function window_der(r, rcut) result(dw)
      implicit none
      real(dp), intent(in) :: r, rcut
      if (r >= rcut) then
         dw = 0.d0
      else
         dw = -0.5d0*pi/rcut*sin(pi*r/rcut)
      end if
   end function window_der

!  virial += sym(f (x) x), the convention of get_soap_energy_and_forces, where
!  f is the force on a neighbour and x its position relative to the site.
   subroutine add_to_virial(virial, f, x)
      implicit none
      real(dp), intent(inout) :: virial(1:3, 1:3)
      real(dp), intent(in) :: f(1:3), x(1:3)
      integer :: k1, k2
      do k2 = 1, 3
         do k1 = 1, 3
            virial(k1, k2) = virial(k1, k2) + 0.5d0*(f(k1)*x(k2) + f(k2)*x(k1))
         end do
      end do
   end subroutine add_to_virial

end module enhanced_sampling

