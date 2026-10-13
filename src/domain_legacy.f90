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
!  Atoms travel as REC_ATOM doubles: global index, species, position, and the
!  eight characters of the xyz species moved bit for bit.
module domain_legacy

   use kinds, only: dp
   use turbogap_comm, only: comm_t, comm_split, comm_free, comm_cart_create, comm_cart_shift, &
                            comm_cart_coords, comm_cart_rank, comm_sendrecv, comm_scatterv, comm_gatherv, &
                            comm_gather, comm_bcast
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
   public :: legacy_check_widths
   public :: legacy_distribute
   public :: legacy_reverse
   public :: legacy_gather
   public :: legacy_frac
   public :: legacy_cell

   integer, parameter :: REC_ATOM = 6
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
   end type legacy_dd_t

contains

!  Group the ranks into domains. block puts consecutive ranks in a domain;
!  cyclic deals them out round-robin. grid = 0 0 0 asks for one domain per rank,
!  laid out by legacy_auto_grid once the cell is known; any other grid is laid
!  out now.
   subroutine legacy_init(ld, world, grid, affinity, balance)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      integer, intent(in) :: grid(1:3)
      character(len=*), intent(in) :: affinity
      character(len=*), intent(in) :: balance
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
      if (.not. all(grid == 0)) call legacy_make_grid(ld, world, grid)
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
      if (ld%grid_ready .and. ld%local%rank == 0) call comm_free(ld%roots)
      call comm_free(ld%root_split)
      call comm_free(ld%local)
      if (allocated(ld%rank_of)) deallocate (ld%rank_of)
      if (allocated(ld%cuts)) deallocate (ld%cuts)
      if (allocated(ld%all_ids)) deallocate (ld%all_ids, ld%local_counts, ld%local_displs)
      if (allocated(ld%ids)) deallocate (ld%ids)
      if (allocated(ld%slot)) deallocate (ld%slot)
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
      integer, parameter :: BINS_PER_SLAB = 256
      real(dp), allocatable :: s(:, :)
      integer, allocatable :: hist(:)
      real(dp) :: w_min
      real(dp) :: target
      integer :: n_bins
      integer :: axis
      integer :: g
      integer :: b
      integer :: k
      integer :: below
      integer :: i

      if (.not. ld%balance_count) return
      allocate (s(1:3, 1:max(n_global, 1)))
      do i = 1, n_global
         s(1:3, i) = legacy_frac(recip, positions(1:3, i))
      end do
      do axis = 1, 3
         g = ld%grid(axis)
         if (g == 1) cycle
         n_bins = BINS_PER_SLAB*g
         allocate (hist(1:n_bins))
         hist = 0
         do i = 1, n_global
            b = min(n_bins, 1 + int(s(axis, i)*dble(n_bins)))
            hist(b) = hist(b) + 1
         end do
         below = 0
         b = 0
         do k = 1, g - 1
            target = dble(k)*dble(n_global)/dble(g)
            do while (b < n_bins .and. dble(below + hist(min(b + 1, n_bins))) <= target)
               b = b + 1
               below = below + hist(b)
            end do
            ld%cuts(k, axis) = dble(b)/dble(n_bins)
            if (b < n_bins .and. hist(b + 1) > 0) then
               ld%cuts(k, axis) = ld%cuts(k, axis) + (target - dble(below))/dble(hist(b + 1))/dble(n_bins)
            end if
         end do
         deallocate (hist)
         w_min = rcut*norm2(recip(axis, 1:3))
         do k = 1, g - 1
            ld%cuts(k, axis) = max(ld%cuts(k, axis), ld%cuts(k - 1, axis) + w_min)
         end do
         do k = g - 1, 1, -1
            ld%cuts(k, axis) = min(ld%cuts(k, axis), ld%cuts(k + 1, axis) - w_min)
         end do
      end do
   end subroutine balance_cuts

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
      integer :: sizes(1:2)
      integer :: i

      ld%n_global = n_global
      call reciprocal_rows(a_box, b_box, c_box, recip)

      if (ld%local%rank == 0) then
         call time_start(t_scatter, "dd_scatter")
         call scatter_owned(ld, world, n_global, positions, species, xyz_species, recip, rcut, rec)
         call build_slots(ld, rec)
         call time_end(t_scatter, "dd_scatter")
         call time_start(t_halo, "dd_halo")
         call exchange_halos(ld, recip, rcut, rec)
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
      call time_start(t_scatter, "dd_scatter")
      call collect_ids(ld, world)
      call time_end(t_scatter, "dd_scatter")
      ld%n_builds = ld%n_builds + 1
   end subroutine legacy_distribute

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
      if (world%rank == 0) call balance_cuts(ld, positions, n_global, recip, rcut)
      call comm_bcast(ld%roots, ld%cuts, size(ld%cuts))
      if (world%rank == 0) then
         allocate (domain(1:n_global))
         counts = 0
         do i = 1, n_global
            c = legacy_cell(legacy_frac(recip, positions(1:3, i)), ld%grid, ld%cuts)
            domain(i) = ld%rank_of(c(1), c(2), c(3)) + 1
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

end module domain_legacy
