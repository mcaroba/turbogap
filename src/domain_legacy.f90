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
!  coordinate. Rank 0 bins the atoms and sends each domain its own; the domain
!  roots then trade halos, the atoms within rcut of a face, in six sweeps:
!  down and up along axis 1, then 2, then 3. An atom received along one axis is
!  forwarded along the later ones, which is what fills edges and corners.
!  legacy_reverse runs the sweeps backwards to return what landed on halo
!  atoms, and legacy_gather brings the owned atoms' results to rank 0.
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
   public :: legacy_free
   public :: legacy_check_widths
   public :: legacy_distribute
   public :: legacy_reverse
   public :: legacy_gather
   public :: legacy_frac
   public :: legacy_cell

   integer, parameter :: REC_ATOM = 6
   integer, parameter :: DOWN = 1
   integer, parameter :: UP = 2

   type, public :: legacy_dd_t
      integer :: grid(1:3) = 1
      integer :: n_domains = 1
!     The ranks of this domain, and one root per domain ranked by domain.
!     roots has size 0 on every rank that is not a root.
      type(comm_t) :: local
      type(comm_t) :: roots
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
   end type legacy_dd_t

contains

!  Group the ranks into domains and the roots into a periodic grid. block puts
!  consecutive ranks in a domain; cyclic deals them out round-robin.
   subroutine legacy_init(ld, world, grid, affinity)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      integer, intent(in) :: grid(1:3)
      character(len=*), intent(in) :: affinity
      type(comm_t) :: split
      integer :: color
      integer :: per_domain
      integer :: axis
      integer :: i
      integer :: j
      integer :: k
      integer :: r

      ld%grid = grid
      ld%n_domains = product(grid)
      if (any(grid < 1) .or. ld%n_domains > world%size) then
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
      if (ld%local%rank == 0) then
         call comm_split(world, 0, world%rank, split)
      else
         call comm_split(world, -1, world%rank, split)
      end if

      if (ld%local%rank == 0) then
         call comm_cart_create(split, ld%grid, ld%roots)
         call comm_free(split)
         call comm_cart_coords(ld%roots, ld%roots%rank, ld%coords)
         do axis = 1, 3
            call comm_cart_shift(ld%roots, axis, ld%nbr(DOWN, axis), ld%nbr(UP, axis))
         end do
      else
         ld%roots = split
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
   end subroutine legacy_init

   subroutine legacy_free(ld)
      type(legacy_dd_t), intent(inout) :: ld

      call comm_free(ld%local)
      call comm_free(ld%roots)
      if (allocated(ld%rank_of)) deallocate (ld%rank_of)
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

!  Grid cell of a fractional coordinate, counting from 0.
   pure function legacy_cell(s, grid) result(c)
      real(dp), intent(in) :: s(1:3)
      integer, intent(in) :: grid(1:3)
      integer :: c(1:3)

      c = min(grid - 1, int(s*dble(grid)))
   end function legacy_cell

!  A halo comes from the adjacent domain only if every domain is at least rcut
!  across, measured perpendicular to its faces.
   subroutine legacy_check_widths(ld, a_box, b_box, c_box, rcut, rank)
      type(legacy_dd_t), intent(in) :: ld
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(in) :: rcut
      integer, intent(in) :: rank
      real(dp) :: recip(1:3, 1:3)
      real(dp) :: width
      integer :: axis

      call reciprocal_rows(a_box, b_box, c_box, recip)
      do axis = 1, 3
         if (ld%grid(axis) == 1) cycle
         width = 1.d0/norm2(recip(axis, 1:3))/dble(ld%grid(axis))
         if (width < rcut) then
            if (rank == 0) write (*, '(A,I0,A,F0.3,A,F0.3,A)') 'ERROR: dd_grid makes domains ', ld%grid(axis), &
               ' across axis ', width, ' A wide; each must be at least rcut_max = ', rcut, ' A'
            call turbogap_abort()
         end if
      end do
   end subroutine legacy_check_widths

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
         call scatter_owned(ld, world, n_global, positions, species, xyz_species, recip, rec)
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
   end subroutine legacy_distribute

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
   subroutine scatter_owned(ld, world, n_global, positions, species, xyz_species, recip, rec)
      type(legacy_dd_t), intent(inout) :: ld
      type(comm_t), intent(in) :: world
      integer, intent(in) :: n_global
      real(dp), intent(in) :: positions(:, :)
      integer, intent(in) :: species(:)
      character*8, intent(in) :: xyz_species(:)
      real(dp), intent(in) :: recip(1:3, 1:3)
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
      if (world%rank == 0) then
         allocate (domain(1:n_global))
         counts = 0
         do i = 1, n_global
            c = legacy_cell(legacy_frac(recip, positions(1:3, i)), ld%grid)
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
      allocate (send(1:REC_ATOM, 1:16))
      do axis = 1, 3
         if (ld%grid(axis) == 1) cycle
         width = 1.d0/norm2(recip(axis, 1:3))
         lo = dble(ld%coords(axis))/dble(ld%grid(axis))
         hi = dble(ld%coords(axis) + 1)/dble(ld%grid(axis))
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
            do i = 1, n_recv
               k = (i - 1)*REC_ATOM
               id = nint(recv(k + 1))
               if (ld%slot(id) /= 0) cycle
               ld%n_local = ld%n_local + 1
               call grow(rec, ld%n_local)
               rec(1:REC_ATOM, ld%n_local) = recv(k + 1:k + REC_ATOM)
               ld%slot(id) = ld%n_local
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
