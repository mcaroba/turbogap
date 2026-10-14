! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, domain_legacy.f90, is copyright (c) 2024-2026, Martti
! HND X   Louhivuori and Tigany Zarrouk
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

!  The legacy spatial decomposition (dd_type = legacy); docs/LEGACY_DD.md.
!
!  The cell is cut into grid(1) x grid(2) x grid(3) slabs of fractional
!  coordinate, at cuts placed so that each slab holds the same number of atoms
!  (or the same volume). Rank 0 bins the atoms and sends each domain its own;
!  the domain roots then trade halos, the atoms within rcut of a face, in six
!  sweeps: down and up along axis 1, then 2, then 3. An atom received along one
!  axis is forwarded along the later ones, which is what fills edges and
!  corners. legacy_reverse runs the sweeps backwards to return what landed on
!  halo atoms, and legacy_gather brings the owned atoms' results to rank 0.
!
!  Between rebuilds the domains keep their atoms: legacy_refresh sends each
!  domain the new positions of the atoms it already holds, as the replicated
!  mode refreshes its neighbour list between rebuilds.
!
!  dd_type rcb replaces the grid by recursive coordinate bisection, weighted
!  by each atom's pair count (rcb_partition): the domains are boxes of any
!  shape and number, so a cluster, a pore or a prime rank count balance too.
!  Each root sends its partners, the domains whose boxes come within rcut of
!  its own, the atoms they need directly (rcb_halos), so forward and reverse
!  are one exchange each, and a migrating atom goes straight to its owner.
!
!  Atoms travel as REC_ATOM doubles: global index, species, position, and the
!  eight characters of the xyz species moved bit for bit.
module domain_legacy

   use kinds, only: dp
   use turbogap_comm, only: comm_t, comm_split, comm_free, comm_cart_create, comm_cart_shift, &
                            comm_cart_coords, comm_cart_rank, comm_sendrecv, comm_scatterv, comm_gatherv, &
                            comm_gather, comm_bcast, comm_sum_all, comm_max_all, comm_alltoall, comm_exchange
   use error, only: turbogap_abort
   use timing, only: time_start, time_end

   implicit none

   private
   public :: legacy_init
   public :: legacy_make_grid
   public :: legacy_auto_grid
   public :: legacy_free
   public :: legacy_refresh
   public :: legacy_forward
   public :: legacy_rebuild_owned
   public :: legacy_migrate
   public :: legacy_check_widths
   public :: legacy_distribute
   public :: legacy_reverse
   public :: legacy_gather
   public :: legacy_frac
   public :: legacy_cell

   integer, parameter :: REC_ATOM = 6
!   An atom that migrates carries, after the REC_ATOM fields: velocity (3),
!   mass, fix_atom (3, as 0 or 1), previous position (3), previous force (3)
!   and the displacement since the last neighbour build (3).
   integer, parameter, public :: REC_DYN = 22
   integer, parameter :: BINS_PER_SLAB = 256
!   RCB: histogram bins per pass of a bisection (two passes, so a cut is
!   placed to 1/RCB_BINS**2 of its box), and how far the heaviest domain may
!   outweigh the average before the cuts are moved at a rebuild.
   integer, parameter :: RCB_BINS = 32
   real(dp), parameter :: RCB_TOLERANCE = 1.05d0
   integer, parameter :: DOWN = 1
   integer, parameter :: UP = 2

!   One halo sweep as it happened: the local sites sent, and where each site
!   received went (0 for a duplicate that was dropped). Replaying it moves any
!   per-site value from owners to their halo copies.
   type :: sweep_t
      integer, allocatable :: sent(:)
      integer, allocatable :: landed(:)
   end type sweep_t

   type, public :: legacy_dd_t
!     Set by legacy_init; legacy_free frees only what it made.
      logical :: active = .false.
      integer :: grid(1:3) = 1
      integer :: n_domains = 1
!     The grid is made once the cell is known when it is chosen automatically.
      logical :: grid_ready = .false.
!     Cut so each slab holds the same number of atoms, or the same volume.
      logical :: balance_count = .true.
!     Fractional coordinate of the cuts along each axis: cuts(0:grid(a), a).
      real(dp), allocatable :: cuts(:, :)
      integer :: n_builds = 0
      integer :: n_refreshes = 0
!     Atoms this root sent to a neighbour in legacy_migrate, over the run.
      integer :: n_migrated = 0
!     The ranks of this domain, and one root per domain ranked by domain.
!     roots has size 0 on every rank that is not a root.
      type(comm_t) :: local
      type(comm_t) :: roots
!     The roots, before the grid is laid over them.
      type(comm_t) :: root_split
      integer :: coords(1:3) = 0
!     roots ranks of the neighbouring domains: (DOWN or UP, axis).
      integer :: nbr(1:2, 1:3) = 0
!     Rank 0 only: the roots rank of the domain at grid coordinates (i, j, k).
      integer, allocatable :: rank_of(:, :, :)
      integer :: n_global = 0
      integer :: n_owned = 0
      integer :: n_local = 0
!     Global index of each local site, owned ones first, then halo.
      integer, allocatable :: ids(:)
!     Roots only: global index -> local site, 0 when absent.
      integer, allocatable :: slot(:)
!     The halo kept from each sweep, as a run of local sites.
      integer :: halo_beg(1:2, 1:3) = 1
      integer :: halo_n(1:2, 1:3) = 0
      type(sweep_t) :: sweeps(1:2, 1:3)
!     Rank 0 only: every domain's local ids, domain after domain, for refresh.
      integer, allocatable :: all_ids(:)
      integer, allocatable :: local_counts(:)
      integer, allocatable :: local_displs(:)
!     dd_type rcb. Every root holds every domain's box in fractional
!     coordinates and the tree that finds an atom's owner: node_child(side, k)
!     is the next node, or minus the domain at a leaf.
      logical :: rcb = .false.
      integer :: n_rebalances = 0
      integer :: n_nodes = 0
      real(dp), allocatable :: box_lo(:, :)
      real(dp), allocatable :: box_hi(:, :)
      integer, allocatable :: node_axis(:)
      real(dp), allocatable :: node_cut(:)
      integer, allocatable :: node_child(:, :)
!     The halo plan, per partner (a roots rank): the owned sites sent to it,
!     send_idx(send_off(p-1)+1:send_off(p)), and the local site each one it
!     sends lands on, land(recv_off(p-1)+1:recv_off(p)), 0 for a duplicate.
      integer, allocatable :: partners(:)
      integer, allocatable :: send_off(:)
      integer, allocatable :: send_idx(:)
      integer, allocatable :: recv_off(:)
      integer, allocatable :: land(:)
   end type legacy_dd_t

contains

!  Group the ranks into domains. block puts consecutive ranks in a domain;
!  cyclic deals them out round-robin. grid = 0 0 0 asks for one domain per rank,
!  laid out by legacy_auto_grid once the cell is known; any other grid is laid
!  out now.
   subroutine legacy_init(ld, world, grid, affinity, balance, rcb)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      integer, intent(in) :: grid(1:3)
      character(len=*), intent(in) :: affinity
      character(len=*), intent(in) :: balance
      logical, intent(in) :: rcb
      integer :: color
      integer :: per_domain

      if (trim(balance) == "count") then
         ld%balance_count = .true.
      else if (trim(balance) == "volume") then
         ld%balance_count = .false.
      else
         if (world%rank == 0) write (*, '(A,A,A)') 'ERROR: dd_balance is "', trim(balance), &
            '"; it must be count or volume'
         call turbogap_abort()
      end if
      ld%grid = grid
      ld%n_domains = product(grid)
      if (all(grid == 0)) ld%n_domains = world%size
      if (any(grid < 1) .and. .not. all(grid == 0) .or. ld%n_domains > world%size) then
         if (world%rank == 0) write (*, '(A,3(1X,I0),A,I0,A)') 'ERROR: dd_grid', grid, ' needs at least ', &
            max(ld%n_domains, 1), ' MPI ranks'
         call turbogap_abort()
      end if
      if (trim(affinity) == "block") then
         if (mod(world%size, ld%n_domains) /= 0) then
            if (world%rank == 0) write (*, '(A,I0,A,I0,A)') 'ERROR: dd_grid_affinity block needs the ', world%size, &
               ' MPI ranks to divide evenly into ', ld%n_domains, ' domains'
            call turbogap_abort()
         end if
         per_domain = world%size/ld%n_domains
         color = world%rank/per_domain
      else if (trim(affinity) == "cyclic") then
         color = mod(world%rank, ld%n_domains)
      else
         if (world%rank == 0) write (*, '(A,A,A)') 'ERROR: dd_grid_affinity is "', trim(affinity), &
            '"; it must be block or cyclic'
         call turbogap_abort()
      end if

      call comm_split(world, color, world%rank, ld%local)
      ld%active = .true.
      if (ld%local%rank == 0) then
         call comm_split(world, 0, world%rank, ld%root_split)
      else
         call comm_split(world, -1, world%rank, ld%root_split)
      end if
      ld%rcb = rcb
      if (rcb) then
         ld%roots = ld%root_split
         ld%grid_ready = .true.
      else if (.not. all(grid == 0)) then
         call legacy_make_grid(ld, world, grid)
      end if
   end subroutine legacy_init

!  Lay the roots out as a periodic grid. Collective over world.
   subroutine legacy_make_grid(ld, world, grid)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      integer, intent(in) :: grid(1:3)
      integer :: axis
      integer :: i
      integer :: j
      integer :: k
      integer :: r

      ld%grid = grid
      if (ld%local%rank == 0) then
         call comm_cart_create(ld%root_split, ld%grid, ld%roots)
         call comm_cart_coords(ld%roots, ld%roots%rank, ld%coords)
         do axis = 1, 3
            call comm_cart_shift(ld%roots, axis, ld%nbr(DOWN, axis), ld%nbr(UP, axis))
         end do
      else
         ld%roots = ld%root_split
      end if

      if (world%rank == 0) then
         allocate (ld%rank_of(0:grid(1) - 1, 0:grid(2) - 1, 0:grid(3) - 1))
         do i = 0, grid(1) - 1
            do j = 0, grid(2) - 1
               do k = 0, grid(3) - 1
                  call comm_cart_rank(ld%roots, [i, j, k], r)
                  ld%rank_of(i, j, k) = r
               end do
            end do
         end do
      end if
      allocate (ld%cuts(0:maxval(grid), 1:3))
      ld%cuts = 1.d0
      do axis = 1, 3
         do i = 0, grid(axis)
            ld%cuts(i, axis) = dble(i)/dble(grid(axis))
         end do
      end do
      ld%grid_ready = .true.
   end subroutine legacy_make_grid

!  The grid of n domains whose halos cost least: the factorisation of n that
!  minimises the local/owned ratio prod(1 + 2 rcut / w_a) over the cut axes,
!  w_a being a domain's perpendicular width. Domains narrower than rcut are not
!  allowed; 0 0 0 if no factorisation fits.
   pure function legacy_auto_grid(n, a_box, b_box, c_box, rcut) result(best)
      integer, intent(in) :: n
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in) :: rcut
      integer :: best(1:3)
      real(dp) :: recip(1:3, 1:3)
      real(dp) :: width(1:3)
      real(dp) :: cost
      real(dp) :: best_cost
      integer :: g(1:3)
      integer :: i
      integer :: j
      integer :: axis
      logical :: fits

      call reciprocal_rows(a_box, b_box, c_box, recip)
      do axis = 1, 3
         width(axis) = 1.d0/norm2(recip(axis, 1:3))
      end do
      best = 0
      best_cost = huge(1.d0)
      do i = 1, n
         if (mod(n, i) /= 0) cycle
         do j = 1, n/i
            if (mod(n/i, j) /= 0) cycle
            g = [i, j, n/(i*j)]
            fits = .true.
            cost = 1.d0
            do axis = 1, 3
               if (g(axis) == 1) cycle
               if (width(axis)/dble(g(axis)) < rcut) fits = .false.
               cost = cost*(1.d0 + 2.d0*rcut*dble(g(axis))/width(axis))
            end do
            if (fits .and. cost < best_cost) then
               best_cost = cost
               best = g
            end if
         end do
      end do
   end function legacy_auto_grid

   subroutine legacy_free(ld)
      type(legacy_dd_t), intent(inout) :: ld

      if (.not. ld%active) return
      if (ld%grid_ready .and. ld%local%rank == 0 .and. .not. ld%rcb) call comm_free(ld%roots)
      call comm_free(ld%root_split)
      call comm_free(ld%local)
      if (allocated(ld%rank_of)) deallocate (ld%rank_of)
      if (allocated(ld%cuts)) deallocate (ld%cuts)
      if (allocated(ld%all_ids)) deallocate (ld%all_ids, ld%local_counts, ld%local_displs)
      if (allocated(ld%ids)) deallocate (ld%ids)
      if (allocated(ld%slot)) deallocate (ld%slot)
      if (allocated(ld%box_lo)) deallocate (ld%box_lo, ld%box_hi, ld%node_axis, ld%node_cut, ld%node_child)
      if (allocated(ld%partners)) deallocate (ld%partners, ld%send_off, ld%send_idx, ld%recv_off, ld%land)
   end subroutine legacy_free

!  Rows of the inverse cell, so recip . r is the fractional coordinate.
   pure subroutine reciprocal_rows(a_box, b_box, c_box, recip)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(out) :: recip(1:3, 1:3)
      real(dp) :: vol

      vol = a_box(1)*(b_box(2)*c_box(3) - b_box(3)*c_box(2)) &
            - a_box(2)*(b_box(1)*c_box(3) - b_box(3)*c_box(1)) &
            + a_box(3)*(b_box(1)*c_box(2) - b_box(2)*c_box(1))
      recip(1, 1:3) = cross(b_box, c_box)/vol
      recip(2, 1:3) = cross(c_box, a_box)/vol
      recip(3, 1:3) = cross(a_box, b_box)/vol
   end subroutine reciprocal_rows

   pure function cross(x, y) result(z)
      real(dp), intent(in) :: x(1:3)
      real(dp), intent(in) :: y(1:3)
      real(dp) :: z(1:3)

      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
   end function cross

!  Fractional coordinate, wrapped into [0, 1).
   pure function legacy_frac(recip, r) result(s)
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: r(1:3)
      real(dp) :: s(1:3)

      s = modulo(matmul(recip, r), 1.d0)
   end function legacy_frac

!  Grid cell of a fractional coordinate, counting from 0: the slab whose cuts
!  bracket it along each axis.
   pure function legacy_cell(s, grid, cuts) result(c)
      real(dp), intent(in) :: s(1:3)
      integer, intent(in) :: grid(1:3)
      real(dp), intent(in) :: cuts(0:, :)
      integer :: c(1:3)
      integer :: axis

      do axis = 1, 3
         c(axis) = count(cuts(1:grid(axis) - 1, axis) <= s(axis))
      end do
   end function legacy_cell

!  Rank 0: cuts that give every slab along an axis the same number of atoms,
!  from a histogram of fractional coordinates, then pushed apart until every
!  slab is at least rcut wide.
   subroutine balance_cuts(ld, positions, n_global, recip, rcut)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: positions(:, :)
      integer, intent(in) :: n_global
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: rcut
      real(dp), allocatable :: hist(:, :)

      if (.not. ld%balance_count) return
      call coordinate_histograms(ld, positions, n_global, recip, hist)
      call cuts_from_histograms(ld, hist, recip, rcut)
   end subroutine balance_cuts

!  How many of the given atoms fall in each bin of fractional coordinate along
!  each axis: BINS_PER_SLAB bins per slab, hist(:, axis).
   subroutine coordinate_histograms(ld, positions, n, recip, hist)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(in) :: positions(:, :)
      integer, intent(in) :: n
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), allocatable, intent(out) :: hist(:, :)
      real(dp) :: s(1:3)
      integer :: axis
      integer :: b
      integer :: i

      allocate (hist(1:BINS_PER_SLAB*maxval(ld%grid), 1:3))
      hist = 0.d0
      do i = 1, n
         s = legacy_frac(recip, positions(1:3, i))
         do axis = 1, 3
            b = min(BINS_PER_SLAB*ld%grid(axis), 1 + int(s(axis)*dble(BINS_PER_SLAB*ld%grid(axis))))
            hist(b, axis) = hist(b, axis) + 1.d0
         end do
      end do
   end subroutine coordinate_histograms

!  Cuts at the atom-count quantiles of the histograms, then pushed apart until
!  every slab is at least rcut wide.
   subroutine cuts_from_histograms(ld, hist, recip, rcut)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: hist(:, :)
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: rcut
      real(dp) :: w_min
      real(dp) :: target
      real(dp) :: below
      real(dp) :: total
      integer :: n_bins
      integer :: axis
      integer :: g
      integer :: b
      integer :: k

      do axis = 1, 3
         g = ld%grid(axis)
         if (g == 1) cycle
         n_bins = BINS_PER_SLAB*g
         total = sum(hist(1:n_bins, axis))
         below = 0.d0
         b = 0
         do k = 1, g - 1
            target = dble(k)*total/dble(g)
            do while (b < n_bins .and. below + hist(min(b + 1, n_bins), axis) <= target)
               b = b + 1
               below = below + hist(b, axis)
            end do
            ld%cuts(k, axis) = dble(b)/dble(n_bins)
            if (b < n_bins .and. hist(b + 1, axis) > 0.d0) then
               ld%cuts(k, axis) = ld%cuts(k, axis) + (target - below)/hist(b + 1, axis)/dble(n_bins)
            end if
         end do
         w_min = rcut*norm2(recip(axis, 1:3))
         do k = 1, g - 1
            ld%cuts(k, axis) = max(ld%cuts(k, axis), ld%cuts(k - 1, axis) + w_min)
         end do
         do k = g - 1, 1, -1
            ld%cuts(k, axis) = min(ld%cuts(k, axis), ld%cuts(k + 1, axis) - w_min)
         end do
      end do
   end subroutine cuts_from_histograms

!  A halo comes from the adjacent domain only if every domain is at least rcut
!  across, measured perpendicular to its faces.
   pure logical function legacy_check_widths(ld, a_box, b_box, c_box, rcut) result(ok)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in) :: rcut
      real(dp) :: recip(1:3, 1:3)
      integer :: axis

      call reciprocal_rows(a_box, b_box, c_box, recip)
      ok = .true.
      if (ld%rcb) return
      do axis = 1, 3
         if (ld%grid(axis) == 1) cycle
         if (1.d0/norm2(recip(axis, 1:3))/dble(ld%grid(axis)) < rcut) ok = .false.
      end do
   end function legacy_check_widths

!  Give every rank its domain's atoms, owned then halo. positions, species and
!  xyz_species are read on rank 0 only. On return every rank of a domain holds
!  the same local arrays and ld%ids, ld%n_owned and ld%n_local. The three
!  buckets time the scatter from rank 0, the halo sweeps and the broadcast
!  inside the domain.
   subroutine legacy_distribute(ld, world, n_global, positions, species, xyz_species, a_box, b_box, c_box, rcut, &
                                pos_l, species_l, xyz_species_l, t_scatter, t_halo, t_share)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      integer, intent(in) :: n_global
      real(dp), intent(in) :: positions(:, :)
      integer, intent(in) :: species(:)
      character*8, intent(in) :: xyz_species(:)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in) :: rcut
      real(dp), allocatable, intent(inout) :: pos_l(:, :)
      integer, allocatable, intent(inout) :: species_l(:)
      character*8, allocatable, intent(inout) :: xyz_species_l(:)
      real(dp), intent(inout) :: t_scatter(3)
      real(dp), intent(inout) :: t_halo(3)
      real(dp), intent(inout) :: t_share(3)
      real(dp), allocatable :: rec(:, :)
      real(dp) :: recip(1:3, 1:3)

      ld%n_global = n_global
      call reciprocal_rows(a_box, b_box, c_box, recip)

      if (ld%local%rank == 0) then
         call time_start(t_scatter, "dd_scatter")
         call scatter_owned(ld, world, n_global, positions, species, xyz_species, recip, rcut, rec)
         call time_end(t_scatter, "dd_scatter")
      end if
      call finish_build(ld, recip, rcut, rec, pos_l, species_l, xyz_species_l, t_halo, t_share)
      call time_start(t_scatter, "dd_scatter")
      call collect_ids(ld, world)
      call time_end(t_scatter, "dd_scatter")
   end subroutine legacy_distribute

!  Domain integration, at a rebuild: the roots already hold their owned atoms,
!  as migrated records (REC_ATOM fields first); halos and the share follow.
   subroutine legacy_rebuild_owned(ld, dyn, a_box, b_box, c_box, rcut, pos_l, species_l, xyz_species_l, &
                                   t_halo, t_share)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: dyn(:, :)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in) :: rcut
      real(dp), allocatable, intent(inout) :: pos_l(:, :)
      integer, allocatable, intent(inout) :: species_l(:)
      character*8, allocatable, intent(inout) :: xyz_species_l(:)
      real(dp), intent(inout) :: t_halo(3)
      real(dp), intent(inout) :: t_share(3)
      real(dp), allocatable :: rec(:, :)
      real(dp) :: recip(1:3, 1:3)

      call reciprocal_rows(a_box, b_box, c_box, recip)
      if (ld%local%rank == 0) then
         allocate (rec(1:REC_ATOM, 1:max(2*ld%n_owned, 16)))
         rec(:, 1:ld%n_owned) = dyn(1:REC_ATOM, 1:ld%n_owned)
         ld%n_local = ld%n_owned
      end if
      call finish_build(ld, recip, rcut, rec, pos_l, species_l, xyz_species_l, t_halo, t_share)
   end subroutine legacy_rebuild_owned

!  From the roots' owned records: the id map, the halos, and every rank of the
!  domain given the same local arrays.
   subroutine finish_build(ld, recip, rcut, rec, pos_l, species_l, xyz_species_l, t_halo, t_share)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: rcut
      real(dp), allocatable, intent(inout) :: rec(:, :)
      real(dp), allocatable, intent(inout) :: pos_l(:, :)
      integer, allocatable, intent(inout) :: species_l(:)
      character*8, allocatable, intent(inout) :: xyz_species_l(:)
      real(dp), intent(inout) :: t_halo(3)
      real(dp), intent(inout) :: t_share(3)
      integer :: sizes(1:2)
      integer :: i

      if (ld%local%rank == 0) then
         call build_slots(ld, rec)
         call time_start(t_halo, "dd_halo")
         if (ld%rcb) then
            call rcb_halos(ld, recip, rcut, rec)
         else
            call exchange_halos(ld, recip, rcut, rec)
         end if
         call time_end(t_halo, "dd_halo")
         sizes = [ld%n_owned, ld%n_local]
      end if

      call time_start(t_share, "dd_share")
      call comm_bcast(ld%local, sizes, 2)
      ld%n_owned = sizes(1)
      ld%n_local = sizes(2)
      if (ld%local%rank /= 0) then
         if (allocated(rec)) deallocate (rec)
         allocate (rec(1:REC_ATOM, 1:max(ld%n_local, 1)))
      end if
      call comm_bcast(ld%local, rec, REC_ATOM*ld%n_local)

      if (allocated(ld%ids)) deallocate (ld%ids)
      allocate (ld%ids(1:ld%n_local))
      if (allocated(pos_l)) deallocate (pos_l)
      if (allocated(species_l)) deallocate (species_l)
      if (allocated(xyz_species_l)) deallocate (xyz_species_l)
      allocate (pos_l(1:3, 1:ld%n_local))
      allocate (species_l(1:ld%n_local))
      allocate (xyz_species_l(1:ld%n_local))
      do i = 1, ld%n_local
         call unpack_atom(rec(1:REC_ATOM, i), ld%ids(i), species_l(i), pos_l(1:3, i), xyz_species_l(i))
      end do
      call time_end(t_share, "dd_share")
      ld%n_builds = ld%n_builds + 1
   end subroutine finish_build

!  Rank 0 learns which atoms every domain holds, owned and halo, in local
!  order, so that legacy_refresh can send them their positions.
   subroutine collect_ids(ld, world)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      real(dp), allocatable :: send(:)
      real(dp), allocatable :: recv(:)
      integer :: d

      if (ld%local%rank /= 0) return
      if (.not. allocated(ld%local_counts)) then
         allocate (ld%local_counts(1:ld%n_domains))
         allocate (ld%local_displs(1:ld%n_domains))
      end if
      call comm_gather(ld%roots, ld%n_local, ld%local_counts)
      allocate (send(1:max(ld%n_local, 1)))
      send(1:ld%n_local) = dble(ld%ids(1:ld%n_local))
      if (world%rank == 0) then
         ld%local_displs(1) = 0
         do d = 2, ld%n_domains
            ld%local_displs(d) = ld%local_displs(d - 1) + ld%local_counts(d - 1)
         end do
         allocate (recv(1:max(sum(ld%local_counts), 1)))
      else
         allocate (recv(1))
      end if
      call comm_gatherv(ld%roots, send, ld%n_local, recv, ld%local_counts, ld%local_displs)
      if (world%rank == 0) then
         if (allocated(ld%all_ids)) deallocate (ld%all_ids)
         allocate (ld%all_ids(1:size(recv)))
         ld%all_ids = nint(recv)
      end if
   end subroutine collect_ids

!  Between rebuilds: every domain keeps its atoms and gets their new
!  positions, in the order it holds them. positions is read on rank 0.
   subroutine legacy_refresh(ld, world, positions, pos_l, t_scatter, t_share)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      real(dp), intent(in) :: positions(:, :)
      real(dp), allocatable, intent(inout) :: pos_l(:, :)
      real(dp), intent(inout) :: t_scatter(3)
      real(dp), intent(inout) :: t_share(3)
      real(dp), allocatable :: send(:)
      real(dp), allocatable :: recv(:)
      integer :: i

      if (ld%local%rank == 0) then
         call time_start(t_scatter, "dd_scatter")
         if (world%rank == 0) then
            allocate (send(1:3*max(size(ld%all_ids), 1)))
            do i = 1, size(ld%all_ids)
               send(3*i - 2:3*i) = positions(1:3, ld%all_ids(i))
            end do
         else
            allocate (send(1:3))
         end if
         allocate (recv(1:3*max(ld%n_local, 1)))
         call comm_scatterv(ld%roots, send, 3*ld%local_counts, 3*ld%local_displs, recv, 3*ld%n_local)
         call time_end(t_scatter, "dd_scatter")
      else
         allocate (recv(1:3*max(ld%n_local, 1)))
      end if
      call time_start(t_share, "dd_share")
      call comm_bcast(ld%local, recv, 3*ld%n_local)
      if (allocated(pos_l)) deallocate (pos_l)
      allocate (pos_l(1:3, 1:ld%n_local))
      pos_l = reshape(recv(1:3*ld%n_local), [3, ld%n_local])
      call time_end(t_share, "dd_share")
      ld%n_refreshes = ld%n_refreshes + 1
   end subroutine legacy_refresh

   pure subroutine pack_atom(id, species, r, xyz_species, rec)
      integer, intent(in) :: id
      integer, intent(in) :: species
      real(dp), intent(in) :: r(1:3)
      character*8, intent(in) :: xyz_species
      real(dp), intent(out) :: rec(1:REC_ATOM)

      rec(1) = dble(id)
      rec(2) = dble(species)
      rec(3:5) = r
      rec(6) = transfer(xyz_species, 0.d0)
   end subroutine pack_atom

   pure subroutine unpack_atom(rec, id, species, r, xyz_species)
      real(dp), intent(in) :: rec(1:REC_ATOM)
      integer, intent(out) :: id
      integer, intent(out) :: species
      real(dp), intent(out) :: r(1:3)
      character*8, intent(out) :: xyz_species

      id = nint(rec(1))
      species = nint(rec(2))
      r = rec(3:5)
      xyz_species = transfer(rec(6), xyz_species)
   end subroutine unpack_atom

!  Rank 0 bins the atoms and each root receives its domain's, in global order.
!  Called on the roots only.
   subroutine scatter_owned(ld, world, n_global, positions, species, xyz_species, recip, rcut, rec)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      integer, intent(in) :: n_global
      real(dp), intent(in) :: positions(:, :)
      integer, intent(in) :: species(:)
      character*8, intent(in) :: xyz_species(:)
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: rcut
      real(dp), allocatable, intent(inout) :: rec(:, :)
      real(dp), allocatable :: send(:, :)
      real(dp), allocatable :: recv(:)
      integer, allocatable :: domain(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: displs(:)
      integer, allocatable :: fill(:)
      integer :: c(1:3)
      integer :: i
      integer :: d

      allocate (counts(1:ld%n_domains))
      allocate (displs(1:ld%n_domains))
      if (ld%rcb) then
         call rcb_from_rank0(ld, world, positions, n_global, recip)
      else
         if (world%rank == 0) call balance_cuts(ld, positions, n_global, recip, rcut)
         call comm_bcast(ld%roots, ld%cuts, size(ld%cuts))
      end if
      if (world%rank == 0) then
         allocate (domain(1:n_global))
         counts = 0
         do i = 1, n_global
            if (ld%rcb) then
               domain(i) = rcb_owner(ld, legacy_frac(recip, positions(1:3, i)))
            else
               c = legacy_cell(legacy_frac(recip, positions(1:3, i)), ld%grid, ld%cuts)
               domain(i) = ld%rank_of(c(1), c(2), c(3)) + 1
            end if
            counts(domain(i)) = counts(domain(i)) + 1
         end do
         displs(1) = 0
         do d = 2, ld%n_domains
            displs(d) = displs(d - 1) + counts(d - 1)
         end do
         allocate (send(1:REC_ATOM, 1:max(n_global, 1)))
         allocate (fill(1:ld%n_domains))
         fill = displs
         do i = 1, n_global
            d = domain(i)
            fill(d) = fill(d) + 1
            call pack_atom(i, species(i), positions(1:3, i), xyz_species(i), send(1:REC_ATOM, fill(d)))
         end do
      else
         allocate (send(1:REC_ATOM, 1))
      end if
      call comm_bcast(ld%roots, counts, ld%n_domains)
      call comm_bcast(ld%roots, displs, ld%n_domains)

      ld%n_owned = counts(ld%roots%rank + 1)
      ld%n_local = ld%n_owned
      allocate (recv(1:REC_ATOM*max(ld%n_owned, 1)))
      call comm_scatterv(ld%roots, reshape(send, [size(send)]), REC_ATOM*counts, REC_ATOM*displs, recv, &
                         REC_ATOM*ld%n_owned)
      if (allocated(rec)) deallocate (rec)
      allocate (rec(1:REC_ATOM, 1:max(2*ld%n_owned, 16)))
      rec(:, 1:ld%n_owned) = reshape(recv(1:REC_ATOM*ld%n_owned), [REC_ATOM, ld%n_owned])
   end subroutine scatter_owned

   subroutine build_slots(ld, rec)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: rec(:, :)
      integer :: i

      if (allocated(ld%slot)) then
         if (size(ld%slot) /= ld%n_global) deallocate (ld%slot)
      end if
      if (.not. allocated(ld%slot)) allocate (ld%slot(1:ld%n_global))
      ld%slot = 0
      do i = 1, ld%n_owned
         ld%slot(nint(rec(1, i))) = i
      end do
   end subroutine build_slots

!  The six sweeps. Only sites already here when an axis starts are sent along
!  it, so nothing bounces back the way it came. A site that arrives twice --
!  both neighbours along an axis are one domain when grid is 2 -- is kept once.
   subroutine exchange_halos(ld, recip, rcut, rec)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: rcut
      real(dp), allocatable, intent(inout) :: rec(:, :)
      real(dp), allocatable :: send(:, :)
      real(dp), allocatable :: recv(:)
      integer, allocatable :: sent(:)
      real(dp) :: width
      real(dp) :: lo
      real(dp) :: hi
      real(dp) :: s(1:3)
      integer :: axis
      integer :: dir
      integer :: n_before
      integer :: n_send
      integer :: n_recv
      integer :: dest
      integer :: source
      integer :: id
      integer :: i
      integer :: k

      ld%halo_beg = ld%n_local + 1
      ld%halo_n = 0
      do axis = 1, 3
         do dir = DOWN, UP
            if (allocated(ld%sweeps(dir, axis)%sent)) deallocate (ld%sweeps(dir, axis)%sent)
            if (allocated(ld%sweeps(dir, axis)%landed)) deallocate (ld%sweeps(dir, axis)%landed)
            allocate (ld%sweeps(dir, axis)%sent(0), ld%sweeps(dir, axis)%landed(0))
         end do
      end do
      allocate (send(1:REC_ATOM, 1:16))
      allocate (sent(1:16))
      do axis = 1, 3
         if (ld%grid(axis) == 1) cycle
         width = 1.d0/norm2(recip(axis, 1:3))
         lo = ld%cuts(ld%coords(axis), axis)
         hi = ld%cuts(ld%coords(axis) + 1, axis)
         n_before = ld%n_local
         do dir = DOWN, UP
            n_send = 0
            do i = 1, n_before
               s = legacy_frac(recip, rec(3:5, i))
               if ((dir == DOWN .and. (s(axis) - lo)*width < rcut) .or. &
                   (dir == UP .and. (hi - s(axis))*width < rcut)) then
                  n_send = n_send + 1
                  call grow(send, n_send)
                  send(1:REC_ATOM, n_send) = rec(1:REC_ATOM, i)
                  if (n_send > size(sent)) sent = [sent, sent]
                  sent(n_send) = i
               end if
            end do
            if (dir == DOWN) then
               dest = ld%nbr(DOWN, axis)
               source = ld%nbr(UP, axis)
            else
               dest = ld%nbr(UP, axis)
               source = ld%nbr(DOWN, axis)
            end if
            call comm_sendrecv(ld%roots, n_send, dest, n_recv, source)
            allocate (recv(1:REC_ATOM*max(n_recv, 1)))
            call comm_sendrecv(ld%roots, reshape(send(:, 1:max(n_send, 1)), [REC_ATOM*max(n_send, 1)]), &
                               REC_ATOM*n_send, dest, recv, REC_ATOM*n_recv, source)
            ld%halo_beg(dir, axis) = ld%n_local + 1
            ld%sweeps(dir, axis)%sent = sent(1:n_send)
            deallocate (ld%sweeps(dir, axis)%landed)
            allocate (ld%sweeps(dir, axis)%landed(1:n_recv))
            ld%sweeps(dir, axis)%landed = 0
            do i = 1, n_recv
               k = (i - 1)*REC_ATOM
               id = nint(recv(k + 1))
               if (ld%slot(id) /= 0) cycle
               ld%n_local = ld%n_local + 1
               call grow(rec, ld%n_local)
               rec(1:REC_ATOM, ld%n_local) = recv(k + 1:k + REC_ATOM)
               ld%slot(id) = ld%n_local
               ld%sweeps(dir, axis)%landed(i) = ld%n_local
            end do
            deallocate (recv)
            ld%halo_n(dir, axis) = ld%n_local + 1 - ld%halo_beg(dir, axis)
         end do
      end do
   end subroutine exchange_halos

!  Domain integration, at a rebuild, on the roots: the cuts are re-balanced
!  from a histogram summed over all domains, then every owned atom whose slab
!  changed moves to the neighbour on the way, one axis at a time, until none
!  is left to move. An atom more than one domain away moves one hop per pass.
!  dyn(:, i) is owned atom i, REC_DYN fields; on return it holds this domain's
!  new owned atoms and ld%n_owned their number. weights(i), its pair count,
!  is what dd_type rcb balances.
   subroutine legacy_migrate(ld, dyn, a_box, b_box, c_box, rcut, weights)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), allocatable, intent(inout) :: dyn(:, :)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in) :: rcut
      real(dp), intent(in) :: weights(:)
      real(dp), allocatable :: hist(:, :)
      real(dp), allocatable :: send(:, :)
      real(dp), allocatable :: keep(:, :)
      real(dp), allocatable :: recv(:)
      real(dp) :: recip(1:3, 1:3)
      real(dp) :: moved
      integer :: c(1:3)
      integer :: n_send(1:2)
      integer :: n_recv
      integer :: n_keep
      integer :: axis
      integer :: dir
      integer :: hop
      integer :: g
      integer :: i
      integer :: k

      call reciprocal_rows(a_box, b_box, c_box, recip)
      if (ld%rcb) then
         call rcb_migrate(ld, dyn, recip, weights)
         return
      end if
      if (ld%balance_count) then
         call coordinate_histograms(ld, dyn(3:5, :), ld%n_owned, recip, hist)
         call comm_sum_all(ld%roots, hist, size(hist))
         call cuts_from_histograms(ld, hist, recip, rcut)
      end if

      do axis = 1, 3
         g = ld%grid(axis)
         if (g == 1) cycle
         do
            allocate (keep(1:REC_DYN, 1:max(ld%n_owned, 1)))
            allocate (send(1:REC_DYN, 1:max(ld%n_owned, 1)))
            n_keep = 0
            moved = 0.d0
            do dir = DOWN, UP
               n_send(dir) = 0
            end do
!           DOWN-bound atoms first, then UP-bound, so each run is contiguous.
            do dir = DOWN, UP
               do i = 1, ld%n_owned
                  c = legacy_cell(legacy_frac(recip, dyn(3:5, i)), ld%grid, ld%cuts)
                  hop = modulo(c(axis) - ld%coords(axis), g)
                  if (hop == 0) then
                     if (dir == DOWN) then
                        n_keep = n_keep + 1
                        keep(:, n_keep) = dyn(:, i)
                     end if
                  else if ((hop > g/2 .and. dir == DOWN) .or. (hop <= g/2 .and. dir == UP)) then
                     n_send(dir) = n_send(dir) + 1
                     send(:, n_send(1) + merge(n_send(2), 0, dir == UP)) = dyn(:, i)
                  end if
               end do
            end do
            moved = dble(n_send(1) + n_send(2))
            ld%n_migrated = ld%n_migrated + n_send(1) + n_send(2)
            call comm_sum_all(ld%roots, moved)
            if (moved == 0.d0) then
               deallocate (keep, send)
               exit
            end if
            do dir = DOWN, UP
               call send_atoms(ld, axis, dir, send, merge(0, n_send(1), dir == DOWN), n_send(dir), recv, n_recv)
               if (n_keep + n_recv > size(keep, 2)) call grow(keep, n_keep + n_recv)
               do k = 1, n_recv
                  keep(:, n_keep + k) = recv((k - 1)*REC_DYN + 1:k*REC_DYN)
               end do
               n_keep = n_keep + n_recv
               deallocate (recv)
            end do
            ld%n_owned = n_keep
            if (allocated(dyn)) deallocate (dyn)
            allocate (dyn(1:REC_DYN, 1:max(n_keep, 1)))
            dyn(:, 1:n_keep) = keep(:, 1:n_keep)
            deallocate (keep, send)
         end do
      end do
   end subroutine legacy_migrate

!  Send n records from send(:, first+1:first+n) to the neighbour in direction
!  dir along axis, and receive what the other neighbour sends this way.
   subroutine send_atoms(ld, axis, dir, send, first, n, recv, n_recv)
      type(legacy_dd_t), intent(in) :: ld
      integer, intent(in) :: axis
      integer, intent(in) :: dir
      real(dp), intent(in) :: send(:, :)
      integer, intent(in) :: first
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: recv(:)
      integer, intent(out) :: n_recv
      real(dp), allocatable :: buf(:)
      integer :: dest
      integer :: source

      if (dir == DOWN) then
         dest = ld%nbr(DOWN, axis)
         source = ld%nbr(UP, axis)
      else
         dest = ld%nbr(UP, axis)
         source = ld%nbr(DOWN, axis)
      end if
      allocate (buf(1:REC_DYN*max(n, 1)))
      if (n > 0) buf(1:REC_DYN*n) = reshape(send(:, first + 1:first + n), [REC_DYN*n])
      call comm_sendrecv(ld%roots, n, dest, n_recv, source)
      allocate (recv(1:REC_DYN*max(n_recv, 1)))
      call comm_sendrecv(ld%roots, buf, REC_DYN*n, dest, recv, REC_DYN*n_recv, source)
   end subroutine send_atoms

!  Make room for column n, doubling.
   subroutine grow(x, n)
      real(dp), allocatable, intent(inout) :: x(:, :)
      integer, intent(in) :: n
      real(dp), allocatable :: tmp(:, :)
      integer :: cap

      if (n <= size(x, 2)) return
      cap = size(x, 2)
      do while (cap < n)
         cap = 2*cap
      end do
      allocate (tmp(1:size(x, 1), 1:cap))
      tmp(:, 1:size(x, 2)) = x
      call move_alloc(tmp, x)
   end subroutine grow

!  Copy per-site values from owners to their halo copies, replaying the
!  sweeps of the last build in order, so values reach edges and corners through
!  the same forwards the atoms did. Called on the roots; values(:, i) belongs to
!  local site i, and only its owned columns need be set on entry.
   subroutine legacy_forward(ld, values)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(inout) :: values(:, :)
      real(dp), allocatable :: send(:)
      real(dp), allocatable :: recv(:)
      integer :: nv
      integer :: axis
      integer :: dir
      integer :: dest
      integer :: source
      integer :: n_send
      integer :: n_recv
      integer :: i
      integer :: j

      if (ld%rcb) then
         call rcb_forward(ld, values)
         return
      end if
      nv = size(values, 1)
      do axis = 1, 3
         if (ld%grid(axis) == 1) cycle
         do dir = DOWN, UP
            if (dir == DOWN) then
               dest = ld%nbr(DOWN, axis)
               source = ld%nbr(UP, axis)
            else
               dest = ld%nbr(UP, axis)
               source = ld%nbr(DOWN, axis)
            end if
            n_send = size(ld%sweeps(dir, axis)%sent)
            n_recv = size(ld%sweeps(dir, axis)%landed)
            allocate (send(1:nv*max(n_send, 1)))
            allocate (recv(1:nv*max(n_recv, 1)))
            do i = 1, n_send
               send((i - 1)*nv + 1:i*nv) = values(1:nv, ld%sweeps(dir, axis)%sent(i))
            end do
            call comm_sendrecv(ld%roots, send, nv*n_send, dest, recv, nv*n_recv, source)
            do i = 1, n_recv
               j = ld%sweeps(dir, axis)%landed(i)
               if (j > 0) values(1:nv, j) = recv((i - 1)*nv + 1:i*nv)
            end do
            deallocate (send, recv)
         end do
      end do
   end subroutine legacy_forward

!  Send what landed on halo sites back to the domain each came from and add
!  it there. Axes and directions run in reverse, so a corner site's share
!  reaches the domain that forwarded it before that domain sends its own back.
!  Called on the roots only; values(:, i) belongs to local site i.
   subroutine legacy_reverse(ld, values)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(inout) :: values(:, :)
      real(dp), allocatable :: send(:, :)
      real(dp), allocatable :: recv(:)
      integer :: nv
      integer :: axis
      integer :: dir
      integer :: dest
      integer :: source
      integer :: n_send
      integer :: n_recv
      integer :: i
      integer :: j
      integer :: k

      if (ld%rcb) then
         call rcb_reverse(ld, values)
         return
      end if
      nv = size(values, 1)
      do axis = 3, 1, -1
         if (ld%grid(axis) == 1) cycle
         do dir = UP, DOWN, -1
!           The halo of this sweep came from source; it goes back there.
            if (dir == DOWN) then
               dest = ld%nbr(UP, axis)
               source = ld%nbr(DOWN, axis)
            else
               dest = ld%nbr(DOWN, axis)
               source = ld%nbr(UP, axis)
            end if
            n_send = ld%halo_n(dir, axis)
            allocate (send(1:nv + 1, 1:max(n_send, 1)))
            do i = 1, n_send
               j = ld%halo_beg(dir, axis) + i - 1
               send(1, i) = dble(ld%ids(j))
               send(2:nv + 1, i) = values(1:nv, j)
            end do
            call comm_sendrecv(ld%roots, n_send, dest, n_recv, source)
            allocate (recv(1:(nv + 1)*max(n_recv, 1)))
            call comm_sendrecv(ld%roots, reshape(send, [size(send)]), (nv + 1)*n_send, dest, recv, &
                               (nv + 1)*n_recv, source)
            do i = 1, n_recv
               k = (i - 1)*(nv + 1)
               j = ld%slot(nint(recv(k + 1)))
               values(1:nv, j) = values(1:nv, j) + recv(k + 2:k + nv + 1)
            end do
            deallocate (send, recv)
         end do
      end do
   end subroutine legacy_reverse

!  Collect the owned sites' values on rank 0, in global order. values is read
!  on the roots, out written on rank 0.
   subroutine legacy_gather(ld, world, values, out)
      type(legacy_dd_t), intent(in) :: ld
      type(comm_t), intent(in) :: world
      real(dp), intent(in) :: values(:, :)
      real(dp), intent(inout) :: out(:, :)
      real(dp), allocatable :: send(:, :)
      real(dp), allocatable :: recv(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: displs(:)
      integer :: nv
      integer :: i
      integer :: k
      integer :: d

      if (ld%local%rank /= 0) return
      nv = size(values, 1)
      allocate (send(1:nv + 1, 1:max(ld%n_owned, 1)))
      do i = 1, ld%n_owned
         send(1, i) = dble(ld%ids(i))
         send(2:nv + 1, i) = values(1:nv, i)
      end do
      allocate (counts(1:ld%n_domains))
      allocate (displs(1:ld%n_domains))
      call comm_gather(ld%roots, (nv + 1)*ld%n_owned, counts)
      if (world%rank == 0) then
         displs(1) = 0
         do d = 2, ld%n_domains
            displs(d) = displs(d - 1) + counts(d - 1)
         end do
         allocate (recv(1:(nv + 1)*max(ld%n_global, 1)))
      else
         allocate (recv(1:nv + 1))
      end if
      call comm_gatherv(ld%roots, reshape(send, [size(send)]), (nv + 1)*ld%n_owned, recv, counts, displs)
      if (world%rank == 0) then
         do i = 1, ld%n_global
            k = (i - 1)*(nv + 1)
            out(1:nv, nint(recv(k + 1))) = recv(k + 2:k + nv + 1)
         end do
      end if
   end subroutine legacy_gather

!  Rank 0 holds every atom: the boxes are drawn over them, weighted by count.
!  Collective over the roots.
   subroutine rcb_from_rank0(ld, world, positions, n_global, recip)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      real(dp), intent(in) :: positions(:, :)
      integer, intent(in) :: n_global
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), allocatable :: s(:, :)
      real(dp), allocatable :: w(:)
      integer :: n
      integer :: i

      n = 0
      if (world%rank == 0) n = n_global
      allocate (s(1:3, 1:max(n, 1)), w(1:max(n, 1)))
      do i = 1, n
         s(1:3, i) = legacy_frac(recip, positions(1:3, i))
      end do
      w = 1.d0
      call rcb_partition(ld, s, w, n, recip, .false.)
   end subroutine rcb_from_rank0

!  Recursive coordinate bisection over the n atoms this root holds, at
!  fractional s with whole-number weights w. A box shared by several domains
!  is cut across its atoms' longest extent so that each side's weight is in
!  proportion to its domains. The weights are summed exactly, so every root
!  draws the same tree. With keep_axes the tree's axes stay and only its cuts
!  move, so the boxes shift rather than reshuffle and few atoms migrate.
!  Collective over the roots.
   subroutine rcb_partition(ld, s, w, n, recip, keep_axes)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: s(:, :)
      real(dp), intent(in) :: w(:)
      integer, intent(in) :: n
      real(dp), intent(in) :: recip(1:3, 1:3)
      logical, intent(in) :: keep_axes
!     A box being cut: lo(3), hi(3); and its domains, first domain, parent
!     node and side of the parent.
      real(dp), allocatable :: box(:, :)
      integer, allocatable :: tag(:, :)
      real(dp), allocatable :: new_box(:, :)
      integer, allocatable :: new_tag(:, :)
      real(dp), allocatable :: ext(:, :)
      real(dp), allocatable :: span(:, :)
      integer, allocatable :: axis_of(:)
      integer, allocatable :: n_dom_of(:)
      integer, allocatable :: old_axis(:)
      real(dp), allocatable :: cut_of(:)
      integer, allocatable :: left_of(:)
      integer, allocatable :: leaf(:)
      real(dp) :: width(1:3)
      real(dp) :: blo(1:3)
      real(dp) :: bhi(1:3)
      logical :: keep
      integer :: n_split
      integer :: n_boxes
      integer :: n_new
      integer :: n_dom
      integer :: n1
      integer :: axis
      integer :: k
      integer :: i
      integer :: d

      do axis = 1, 3
         width(axis) = 1.d0/norm2(recip(axis, 1:3))
      end do
      n_dom = ld%n_domains
      keep = keep_axes .and. ld%n_nodes == n_dom - 1 .and. n_dom > 1
      if (keep) old_axis = ld%node_axis(1:ld%n_nodes)
      if (allocated(ld%box_lo)) deallocate (ld%box_lo, ld%box_hi, ld%node_axis, ld%node_cut, ld%node_child)
      allocate (ld%box_lo(1:3, 1:n_dom), ld%box_hi(1:3, 1:n_dom))
      allocate (ld%node_axis(1:max(n_dom - 1, 1)), ld%node_cut(1:max(n_dom - 1, 1)))
      allocate (ld%node_child(1:2, 1:max(n_dom - 1, 1)))
      ld%n_nodes = 0

      n_boxes = 1
      allocate (box(1:6, 1:1), tag(1:4, 1:1))
      box(1:3, 1) = 0.d0
      box(4:6, 1) = 1.d0
      tag(1:4, 1) = [n_dom, 1, 0, 0]
      allocate (leaf(1:max(n, 1)))
      leaf = 1

      do while (n_boxes > 0)
!        Each box's atoms' extent, as maxima of s and of -s.
         allocate (ext(1:6, 1:n_boxes))
         ext = -huge(1.d0)
         do i = 1, n
            k = leaf(i)
            if (k == 0) cycle
            ext(1:3, k) = max(ext(1:3, k), s(1:3, i))
            ext(4:6, k) = max(ext(4:6, k), -s(1:3, i))
         end do
         call comm_max_all(ld%roots, ext, 6*n_boxes)

         allocate (axis_of(1:n_boxes), span(1:2, 1:n_boxes), cut_of(1:n_boxes), left_of(1:n_boxes))
         axis_of = 0
         n_split = 0
         do k = 1, n_boxes
            if (tag(1, k) == 1) cycle
            n_split = n_split + 1
            if (ext(1, k) >= -ext(4, k)) then
               blo = -ext(4:6, k)
               bhi = ext(1:3, k)
            else
               blo = box(1:3, k)
               bhi = box(4:6, k)
            end if
!           Nodes are made box by box, level by level, so the tree being
!           redrawn numbers them the same way.
            if (keep) then
               axis_of(k) = old_axis(ld%n_nodes + n_split)
            else
               axis_of(k) = maxloc((bhi - blo)*width, 1)
            end if
            span(1:2, k) = [blo(axis_of(k)), bhi(axis_of(k))]
         end do
         n_dom_of = tag(1, 1:n_boxes)
         call rcb_cuts(ld, s, w, n, leaf, axis_of, n_dom_of, span, cut_of)

         n_new = 2*count(axis_of > 0)
         allocate (new_box(1:6, 1:max(n_new, 1)), new_tag(1:4, 1:max(n_new, 1)))
         n_new = 0
         do k = 1, n_boxes
            if (axis_of(k) == 0) then
               d = tag(2, k)
               ld%box_lo(1:3, d) = box(1:3, k)
               ld%box_hi(1:3, d) = box(4:6, k)
               if (tag(3, k) > 0) ld%node_child(tag(4, k), tag(3, k)) = -d
               cycle
            end if
            axis = axis_of(k)
            n1 = tag(1, k)/2
            if (cut_of(k) < 0.d0) cut_of(k) = 0.5d0*(box(axis, k) + box(3 + axis, k))
            cut_of(k) = min(max(cut_of(k), box(axis, k)), box(3 + axis, k))

            ld%n_nodes = ld%n_nodes + 1
            ld%node_axis(ld%n_nodes) = axis
            ld%node_cut(ld%n_nodes) = cut_of(k)
            if (tag(3, k) > 0) ld%node_child(tag(4, k), tag(3, k)) = ld%n_nodes

            left_of(k) = n_new + 1
            new_box(1:6, n_new + 1) = box(1:6, k)
            new_box(3 + axis, n_new + 1) = cut_of(k)
            new_tag(1:4, n_new + 1) = [n1, tag(2, k), ld%n_nodes, 1]
            new_box(1:6, n_new + 2) = box(1:6, k)
            new_box(axis, n_new + 2) = cut_of(k)
            new_tag(1:4, n_new + 2) = [tag(1, k) - n1, tag(2, k) + n1, ld%n_nodes, 2]
            n_new = n_new + 2
         end do

         do i = 1, n
            k = leaf(i)
            if (k == 0) cycle
            if (axis_of(k) == 0) then
               leaf(i) = 0
            else if (s(axis_of(k), i) < cut_of(k)) then
               leaf(i) = left_of(k)
            else
               leaf(i) = left_of(k) + 1
            end if
         end do

         deallocate (ext, axis_of, span, cut_of, left_of, box, tag)
         call move_alloc(new_box, box)
         call move_alloc(new_tag, tag)
         n_boxes = n_new
      end do
   end subroutine rcb_partition

!  Where each box being split (axis_of > 0) is cut along its axis so the weight
!  below is n_dom_of/2 of n_dom_of parts: a histogram over the atoms' span
!  finds the bin, a second over that bin the place within it. -1 for a box
!  with no weight. Collective over the roots.
   subroutine rcb_cuts(ld, s, w, n, leaf, axis_of, n_dom_of, span, cut_of)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(in) :: s(:, :)
      real(dp), intent(in) :: w(:)
      integer, intent(in) :: n
      integer, intent(in) :: leaf(:)
      integer, intent(in) :: axis_of(:)
      integer, intent(in) :: n_dom_of(:)
      real(dp), intent(in) :: span(:, :)
      real(dp), intent(out) :: cut_of(:)
      real(dp), allocatable :: hist(:, :)
      real(dp), allocatable :: fine(:, :)
      real(dp), allocatable :: target(:)
      real(dp), allocatable :: below(:)
      integer, allocatable :: coarse(:)
      real(dp) :: dx
      real(dp) :: f
      integer :: n_boxes
      integer :: k
      integer :: i
      integer :: b

      n_boxes = size(axis_of)
      allocate (hist(1:RCB_BINS, 1:n_boxes), fine(1:RCB_BINS, 1:n_boxes))
      allocate (target(1:n_boxes), below(1:n_boxes), coarse(1:n_boxes))
      hist = 0.d0
      do i = 1, n
         k = leaf(i)
         if (k == 0) cycle
         if (axis_of(k) == 0) cycle
         b = rcb_bin(s(axis_of(k), i), span(1:2, k))
         hist(b, k) = hist(b, k) + w(i)
      end do
      call comm_sum_all(ld%roots, hist, RCB_BINS*n_boxes)

!     The coarse bin holding the cut, and the weight below it.
      coarse = 0
      below = 0.d0
      do k = 1, n_boxes
         if (axis_of(k) == 0) cycle
         target(k) = sum(hist(1:RCB_BINS, k))*dble(n_dom_of(k)/2)/dble(n_dom_of(k))
         do while (coarse(k) < RCB_BINS - 1 .and. below(k) + hist(coarse(k) + 1, k) <= target(k))
            coarse(k) = coarse(k) + 1
            below(k) = below(k) + hist(coarse(k), k)
         end do
      end do

      fine = 0.d0
      do i = 1, n
         k = leaf(i)
         if (k == 0) cycle
         if (axis_of(k) == 0) cycle
         if (rcb_bin(s(axis_of(k), i), span(1:2, k)) /= coarse(k) + 1) cycle
         dx = (span(2, k) - span(1, k))/dble(RCB_BINS)
         b = rcb_bin(s(axis_of(k), i), [span(1, k) + dble(coarse(k))*dx, span(1, k) + dble(coarse(k) + 1)*dx])
         fine(b, k) = fine(b, k) + w(i)
      end do
      call comm_sum_all(ld%roots, fine, RCB_BINS*n_boxes)

      cut_of = -1.d0
      do k = 1, n_boxes
         if (axis_of(k) == 0) cycle
         if (sum(hist(1:RCB_BINS, k)) <= 0.d0 .or. span(2, k) <= span(1, k)) cycle
         b = 0
         do while (b < RCB_BINS .and. below(k) + fine(min(b + 1, RCB_BINS), k) <= target(k))
            b = b + 1
            below(k) = below(k) + fine(b, k)
         end do
         f = 0.d0
         if (b < RCB_BINS) then
            if (fine(b + 1, k) > 0.d0) f = (target(k) - below(k))/fine(b + 1, k)
         end if
         dx = (span(2, k) - span(1, k))/dble(RCB_BINS)
         cut_of(k) = span(1, k) + dble(coarse(k))*dx + (dble(b) + f)/dble(RCB_BINS)*dx
      end do
   end subroutine rcb_cuts

!  The histogram bin of coordinate x over span(1):span(2).
   pure integer function rcb_bin(x, span) result(b)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: span(1:2)

      b = 1
      if (span(2) > span(1)) b = min(RCB_BINS, max(1, 1 + int((x - span(1))/(span(2) - span(1))*dble(RCB_BINS))))
   end function rcb_bin

!  The domain whose box holds fractional coordinate s.
   pure integer function rcb_owner(ld, s) result(d)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(in) :: s(1:3)
      integer :: k
      integer :: c

      d = 1
      if (ld%n_nodes == 0) return
      k = 1
      do
         if (s(ld%node_axis(k)) < ld%node_cut(k)) then
            c = ld%node_child(1, k)
         else
            c = ld%node_child(2, k)
         end if
         if (c < 0) then
            d = -c
            return
         end if
         k = c
      end do
   end function rcb_owner

!  Do the boxes of domains d1 < d2 come within w of each other along every
!  axis of the periodic cell? Always asked with the lower domain first, so
!  both sides of a pair reach the same answer.
   pure logical function rcb_boxes_meet(ld, d1, d2, w) result(yes)
      type(legacy_dd_t), intent(in) :: ld
      integer, intent(in) :: d1
      integer, intent(in) :: d2
      real(dp), intent(in) :: w(1:3)
      real(dp) :: l1
      real(dp) :: l2
      real(dp) :: d
      integer :: a

      yes = .false.
      do a = 1, 3
         l1 = ld%box_hi(a, d1) - ld%box_lo(a, d1)
         l2 = ld%box_hi(a, d2) - ld%box_lo(a, d2)
         if (l1 + l2 + 2.d0*w(a) >= 1.d0) cycle
         d = modulo(ld%box_lo(a, d2) - ld%box_lo(a, d1), 1.d0)
         if (d < l1 + w(a) .or. d + l2 > 1.d0 - w(a)) cycle
         return
      end do
      yes = .true.
   end function rcb_boxes_meet

!  Is fractional s within w of domain d's box, the cell being periodic?
   pure logical function rcb_near(ld, d, s, w) result(yes)
      type(legacy_dd_t), intent(in) :: ld
      integer, intent(in) :: d
      real(dp), intent(in) :: s(1:3)
      real(dp), intent(in) :: w(1:3)
      real(dp) :: l
      integer :: a

      yes = .false.
      do a = 1, 3
         l = ld%box_hi(a, d) - ld%box_lo(a, d) + 2.d0*w(a)
         if (l >= 1.d0) cycle
         if (modulo(s(a) - ld%box_lo(a, d) + w(a), 1.d0) >= l) return
      end do
      yes = .true.
   end function rcb_near

!  The halo, in one exchange: every partner sends the owned atoms within rcut
!  of this box, and they are appended after the owned ones, partner by
!  partner. Called on the roots, with the owned records in rec.
   subroutine rcb_halos(ld, recip, rcut, rec)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: rcut
      real(dp), allocatable, intent(inout) :: rec(:, :)
      real(dp), allocatable :: send(:)
      real(dp), allocatable :: recv(:)
      real(dp), allocatable :: n_out(:)
      real(dp), allocatable :: n_in(:)
      real(dp), allocatable :: s(:, :)
      integer, allocatable :: ones(:)
      integer, allocatable :: steps(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: displs(:)
      integer, allocatable :: r_counts(:)
      integer, allocatable :: r_displs(:)
      real(dp) :: w(1:3)
      integer :: me
      integer :: n_p
      integer :: p
      integer :: d
      integer :: i
      integer :: k
      integer :: id

      me = ld%roots%rank + 1
      do i = 1, 3
         w(i) = rcut*norm2(recip(i, 1:3))
      end do
      if (allocated(ld%partners)) deallocate (ld%partners, ld%send_off, ld%send_idx, ld%recv_off, ld%land)
      allocate (ld%partners(1:ld%n_domains))
      n_p = 0
      do d = 1, ld%n_domains
         if (d == me) cycle
         if (rcb_boxes_meet(ld, min(d, me), max(d, me), w)) then
            n_p = n_p + 1
            ld%partners(n_p) = d - 1
         end if
      end do
      ld%partners = ld%partners(1:n_p)

      allocate (s(1:3, 1:max(ld%n_owned, 1)))
      do i = 1, ld%n_owned
         s(1:3, i) = legacy_frac(recip, rec(3:5, i))
      end do
      allocate (ld%send_off(0:n_p), ld%send_idx(1:16))
      ld%send_off(0) = 0
      k = 0
      do p = 1, n_p
         do i = 1, ld%n_owned
            if (.not. rcb_near(ld, ld%partners(p) + 1, s(1:3, i), w)) cycle
            k = k + 1
            if (k > size(ld%send_idx)) ld%send_idx = [ld%send_idx, ld%send_idx]
            ld%send_idx(k) = i
         end do
         ld%send_off(p) = k
      end do

      allocate (ones(1:max(n_p, 1)), steps(1:max(n_p, 1)), n_out(1:max(n_p, 1)), n_in(1:max(n_p, 1)))
      ones = 1
      do p = 1, n_p
         steps(p) = p - 1
         n_out(p) = dble(ld%send_off(p) - ld%send_off(p - 1))
      end do
      call comm_exchange(ld%roots, ld%partners, n_out, ones, steps, n_in, ones, steps)

      allocate (counts(1:max(n_p, 1)), displs(1:max(n_p, 1)), r_counts(1:max(n_p, 1)), r_displs(1:max(n_p, 1)))
      allocate (ld%recv_off(0:n_p))
      ld%recv_off(0) = 0
      do p = 1, n_p
         counts(p) = REC_ATOM*(ld%send_off(p) - ld%send_off(p - 1))
         displs(p) = REC_ATOM*ld%send_off(p - 1)
         ld%recv_off(p) = ld%recv_off(p - 1) + nint(n_in(p))
         r_counts(p) = REC_ATOM*nint(n_in(p))
         r_displs(p) = REC_ATOM*ld%recv_off(p - 1)
      end do
      allocate (send(1:REC_ATOM*max(ld%send_off(n_p), 1)), recv(1:REC_ATOM*max(ld%recv_off(n_p), 1)))
      do k = 1, ld%send_off(n_p)
         send((k - 1)*REC_ATOM + 1:k*REC_ATOM) = rec(1:REC_ATOM, ld%send_idx(k))
      end do
      call comm_exchange(ld%roots, ld%partners, send, counts, displs, recv, r_counts, r_displs)

      allocate (ld%land(1:max(ld%recv_off(n_p), 1)))
      ld%land = 0
      do k = 1, ld%recv_off(n_p)
         id = nint(recv((k - 1)*REC_ATOM + 1))
         if (ld%slot(id) /= 0) cycle
         ld%n_local = ld%n_local + 1
         call grow(rec, ld%n_local)
         rec(1:REC_ATOM, ld%n_local) = recv((k - 1)*REC_ATOM + 1:k*REC_ATOM)
         ld%slot(id) = ld%n_local
         ld%land(k) = ld%n_local
      end do
   end subroutine rcb_halos

!  The element counts and offsets of nv values per site over the plan's runs:
!  off(p-1)+1:off(p) for partner p.
   subroutine rcb_layout(off, nv, counts, displs)
      integer, intent(in) :: off(0:)
      integer, intent(in) :: nv
      integer, allocatable, intent(out) :: counts(:)
      integer, allocatable, intent(out) :: displs(:)
      integer :: p

      allocate (counts(1:max(ubound(off, 1), 1)), displs(1:max(ubound(off, 1), 1)))
      do p = 1, ubound(off, 1)
         counts(p) = nv*(off(p) - off(p - 1))
         displs(p) = nv*off(p - 1)
      end do
   end subroutine rcb_layout

!  Owners' values to their halo copies. Called on the roots.
   subroutine rcb_forward(ld, values)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(inout) :: values(:, :)
      real(dp), allocatable :: send(:)
      real(dp), allocatable :: recv(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: displs(:)
      integer, allocatable :: r_counts(:)
      integer, allocatable :: r_displs(:)
      integer :: nv
      integer :: n_p
      integer :: k

      nv = size(values, 1)
      n_p = size(ld%partners)
      call rcb_layout(ld%send_off, nv, counts, displs)
      call rcb_layout(ld%recv_off, nv, r_counts, r_displs)
      allocate (send(1:nv*max(ld%send_off(n_p), 1)), recv(1:nv*max(ld%recv_off(n_p), 1)))
      do k = 1, ld%send_off(n_p)
         send((k - 1)*nv + 1:k*nv) = values(1:nv, ld%send_idx(k))
      end do
      call comm_exchange(ld%roots, ld%partners, send, counts, displs, recv, r_counts, r_displs)
      do k = 1, ld%recv_off(n_p)
         if (ld%land(k) > 0) values(1:nv, ld%land(k)) = recv((k - 1)*nv + 1:k*nv)
      end do
   end subroutine rcb_forward

!  What landed on halo copies, back to the owners and added there, partner by
!  partner in a fixed order. Called on the roots.
   subroutine rcb_reverse(ld, values)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(inout) :: values(:, :)
      real(dp), allocatable :: send(:)
      real(dp), allocatable :: recv(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: displs(:)
      integer, allocatable :: r_counts(:)
      integer, allocatable :: r_displs(:)
      integer :: nv
      integer :: n_p
      integer :: k
      integer :: j

      nv = size(values, 1)
      n_p = size(ld%partners)
      call rcb_layout(ld%recv_off, nv, counts, displs)
      call rcb_layout(ld%send_off, nv, r_counts, r_displs)
      allocate (send(1:nv*max(ld%recv_off(n_p), 1)), recv(1:nv*max(ld%send_off(n_p), 1)))
      send = 0.d0
      do k = 1, ld%recv_off(n_p)
         if (ld%land(k) > 0) send((k - 1)*nv + 1:k*nv) = values(1:nv, ld%land(k))
      end do
      call comm_exchange(ld%roots, ld%partners, send, counts, displs, recv, r_counts, r_displs)
      do k = 1, ld%send_off(n_p)
         j = ld%send_idx(k)
         values(1:nv, j) = values(1:nv, j) + recv((k - 1)*nv + 1:k*nv)
      end do
   end subroutine rcb_reverse

!  Domain integration, at a rebuild, on the roots: the boxes are redrawn by
!  weight if the heaviest domain outweighs the average by more than
!  RCB_TOLERANCE, then every owned atom goes straight to the domain whose box
!  holds it. Kept atoms stay first, in order; arrivals follow by sender.
   subroutine rcb_migrate(ld, dyn, recip, weights)
      type(legacy_dd_t), intent(inout) :: ld
      real(dp), allocatable, intent(inout) :: dyn(:, :)
      real(dp), intent(in) :: recip(1:3, 1:3)
      real(dp), intent(in) :: weights(:)
      real(dp), allocatable :: s(:, :)
      real(dp), allocatable :: send(:)
      real(dp), allocatable :: recv(:)
      real(dp), allocatable :: keep(:, :)
      integer, allocatable :: dest(:)
      integer, allocatable :: n_to(:)
      integer, allocatable :: n_from(:)
      integer, allocatable :: partners(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: displs(:)
      integer, allocatable :: r_counts(:)
      integer, allocatable :: r_displs(:)
      integer, allocatable :: fill(:)
      real(dp) :: load(1:1)
      real(dp) :: total
      integer :: me
      integer :: n_p
      integer :: n_keep
      integer :: n_recv
      integer :: p
      integer :: d
      integer :: i

      me = ld%roots%rank + 1
      allocate (s(1:3, 1:max(ld%n_owned, 1)))
      do i = 1, ld%n_owned
         s(1:3, i) = legacy_frac(recip, dyn(3:5, i))
      end do
      load(1) = sum(weights(1:ld%n_owned))
      total = load(1)
      call comm_sum_all(ld%roots, total)
      call comm_max_all(ld%roots, load, 1)
      if (load(1) > RCB_TOLERANCE*total/dble(ld%n_domains)) then
         call rcb_partition(ld, s, weights, ld%n_owned, recip, .true.)
         ld%n_rebalances = ld%n_rebalances + 1
      end if

      allocate (dest(1:max(ld%n_owned, 1)), n_to(1:ld%n_domains), n_from(1:ld%n_domains))
      n_to = 0
      do i = 1, ld%n_owned
         dest(i) = rcb_owner(ld, s(1:3, i))
         if (dest(i) /= me) n_to(dest(i)) = n_to(dest(i)) + 1
      end do
      call comm_alltoall(ld%roots, n_to, n_from)
      ld%n_migrated = ld%n_migrated + sum(n_to)

      n_p = count(n_to > 0 .or. n_from > 0)
      allocate (partners(1:max(n_p, 1)), counts(1:max(n_p, 1)), displs(1:max(n_p, 1)))
      allocate (r_counts(1:max(n_p, 1)), r_displs(1:max(n_p, 1)), fill(1:ld%n_domains))
      p = 0
      fill = 0
      do d = 1, ld%n_domains
         if (n_to(d) == 0 .and. n_from(d) == 0) cycle
         p = p + 1
         partners(p) = d - 1
         counts(p) = REC_DYN*n_to(d)
         r_counts(p) = REC_DYN*n_from(d)
         displs(p) = 0
         r_displs(p) = 0
         if (p > 1) then
            displs(p) = displs(p - 1) + counts(p - 1)
            r_displs(p) = r_displs(p - 1) + r_counts(p - 1)
         end if
         fill(d) = displs(p)
      end do

      allocate (send(1:max(REC_DYN*sum(n_to), 1)), recv(1:max(REC_DYN*sum(n_from), 1)))
      allocate (keep(1:REC_DYN, 1:max(ld%n_owned + sum(n_from), 1)))
      n_keep = 0
      do i = 1, ld%n_owned
         if (dest(i) == me) then
            n_keep = n_keep + 1
            keep(1:REC_DYN, n_keep) = dyn(1:REC_DYN, i)
         else
            send(fill(dest(i)) + 1:fill(dest(i)) + REC_DYN) = dyn(1:REC_DYN, i)
            fill(dest(i)) = fill(dest(i)) + REC_DYN
         end if
      end do
      call comm_exchange(ld%roots, partners(1:n_p), send, counts, displs, recv, r_counts, r_displs)
      n_recv = sum(n_from)
      do i = 1, n_recv
         keep(1:REC_DYN, n_keep + i) = recv((i - 1)*REC_DYN + 1:i*REC_DYN)
      end do
      ld%n_owned = n_keep + n_recv
      call move_alloc(keep, dyn)
   end subroutine rcb_migrate

end module domain_legacy
