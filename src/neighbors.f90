! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, neighbors.f90, is copyright (c) 2019-2022, Miguel A. Caro
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

module neighbors

   use kinds

   use soap_turbo_functions
   use timing, only: get_time
!  int64 for the bin-count product in cell_grid, which overflows the default
!  kind before any real cell needs that many bins.
   use, intrinsic :: iso_fortran_env, only: int64
#ifdef _GPU
   use mpi
#endif

!  Headroom on top of soap_batch_memory_model. The terms there are the
!  allocations that scale with the batch; this covers what does not appear in
!  them, and allocator behaviour. Raising it costs batches; lowering it costs
!  the run.
   real(dp), parameter :: SOAP_BATCH_SAFETY = 1.10d0

!  Floor on the room made for one site's neighbours, on top of a quarter over
!  the mean density. Too small only costs the search a second time.
   integer, parameter :: NEIGHBOR_CAP_SLACK = 16

contains

!
! This subroutine returns the distance between ri and rj under certain
! boundary conditions.
   subroutine get_distance(posi, posj, a, b, c, PBC, dist, d, i_shift)

      implicit none

      real(dp), intent(in) :: posi(1:3)
      real(dp), intent(in) :: posj(1:3)
      real(dp), intent(in) :: a(1:3)
      real(dp), intent(in) :: b(1:3)
      real(dp), intent(in) :: c(1:3)
      integer, intent(out) :: i_shift(1:3)
      logical, intent(in) :: PBC(1:3)
      real(dp), intent(out) :: d
      real(dp), intent(out) :: dist(1:3)
      real(dp) :: d2
      real(dp) :: L(1:3)
      real(dp) :: d_tol = 1.d-6
      real(dp) :: mat(1:3, 1:3)
      real(dp) :: md
      real(dp) :: indices_real(1:3)
      real(dp) :: res
      real(dp) :: res_opt
      real(dp) :: dist_opt(1:3)
      real(dp) :: dist_temp(1:3)
      real(dp), save :: a0(1:3) = 0.d0
      real(dp), save :: b0(1:3) = 0.d0
      real(dp), save :: c0(1:3) = 0.d0
      real(dp), save :: mat_inv(1:3, 1:3) = 0.d0
!     Half the smallest perpendicular width of the cell, squared: below this a
!     pair is on its nearest image already. Zero until a lattice has been seen,
!     which only means the full search runs.
      real(dp), save :: res_near = 0.d0
      real(dp) :: d_near
!     The non-orthorhombic branch below caches the reciprocal cell in these and
!     recomputes it only when the lattice changes. Saved state shared between
!     threads is a race, and the callers of this routine are OpenMP loops, so
!     each thread gets its own cache. It costs one extra inversion per thread
!     per lattice change and nothing at all on the orthorhombic path, which
!     never touches them.
!$OMP THREADPRIVATE(a0, b0, c0, mat_inv, res_near)
      integer :: i
      integer :: j
      integer :: k
      integer :: indices(1:3)
      logical :: lattice_check_a(1:3)
      logical :: lattice_check_b(1:3)
      logical :: lattice_check_c(1:3)

      if (dabs(a(2)) < d_tol .and. dabs(a(3)) < d_tol .and. &
          dabs(b(1)) < d_tol .and. dabs(b(3)) < d_tol .and. &
          dabs(c(1)) < d_tol .and. dabs(c(2)) < d_tol) then
         dist = posj - posi
         i_shift = floor([dist(1)/a(1), dist(2)/b(2), dist(3)/c(3)])
!     Fast solution for orthorhombic cells
         L = (/a(1), b(2), c(3)/)
         d2 = 0.d0
         do i = 1, 3
            if (PBC(i)) then
               dist(i) = modulo(posj(i) - posi(i), L(i))
               if (dist(i) > L(i)/2.d0) then
                  dist(i) = dist(i) - L(i)
               end if
            else
               dist(i) = posj(i) - posi(i)
            end if
            d2 = d2 + dist(i)**2
         end do
         d = dsqrt(d2)
      else if (all(PBC)) then
!     Slow solution for other unit cells
         lattice_check_a = (a /= a0)
         lattice_check_b = (b /= b0)
         lattice_check_c = (c /= c0)
         if (any(lattice_check_a) .or. any(lattice_check_b) .or. any(lattice_check_c)) then
            a0 = a
            b0 = b
            c0 = c
!       We construct our matrix to get the MIC only if the lattice vectors have changed
            mat(1, 1) = dot_product(a, a)
            mat(1, 2) = dot_product(a, b)
            mat(1, 3) = dot_product(a, c)
            mat(2, 1) = mat(1, 2)
            mat(2, 2) = dot_product(b, b)
            mat(2, 3) = dot_product(b, c)
            mat(3, 1) = mat(1, 3)
            mat(3, 2) = mat(2, 3)
            mat(3, 3) = dot_product(c, c)
!       We compute the inverse of this matrix analytically
            md = -mat(1, 3)**2*mat(2, 2) + 2.d0*mat(1, 2)*mat(1, 3)*mat(2, 3) - mat(1, 1)*mat(2, 3)**2 &
                 - mat(1, 2)**2*mat(3, 3) + mat(1, 1)*mat(2, 2)*mat(3, 3)
            mat_inv(1, 1) = mat(2, 2)*mat(3, 3) - mat(2, 3)**2
            mat_inv(1, 2) = mat(1, 3)*mat(2, 3) - mat(1, 2)*mat(3, 3)
            mat_inv(1, 3) = mat(1, 2)*mat(2, 3) - mat(1, 3)*mat(2, 2)
            mat_inv(2, 1) = mat_inv(1, 2)
            mat_inv(2, 2) = mat(1, 1)*mat(3, 3) - mat(1, 3)**2
            mat_inv(2, 3) = mat(1, 2)*mat(1, 3) - mat(1, 1)*mat(2, 3)
            mat_inv(3, 1) = mat_inv(1, 3)
            mat_inv(3, 2) = mat_inv(2, 3)
            mat_inv(3, 3) = mat(1, 1)*mat(2, 2) - mat(1, 2)**2
            mat_inv = mat_inv/md
!       w_k = 1/|g_k| and the diagonal of the inverse Gram matrix is |g_k|^2,
!       so the widths come out of the matrix just built rather than out of
!       three more cross products. The 1.d-10 A keeps the threshold strictly
!       inside the bound it stands for.
            d_near = 0.5d0/dsqrt(maxval((/mat_inv(1, 1), mat_inv(2, 2), mat_inv(3, 3)/))) - 1.d-10
            res_near = d_near*d_near
         end if
         dist = posj - posi
         indices_real = -1.d0*(/dot_product(dist, a), dot_product(dist, b), dot_product(dist, c)/)
         indices_real = matmul(mat_inv, indices_real)
         i_shift = floor(-indices_real)
!     Closest integer solution
         indices(1:3) = nint(indices_real(1:3))
!     Changing any of those integers moves the pair by at least half a cell
!     width -- |dr| >= |ds_k| w_k, and a different integer makes some |ds_k| at
!     least 1/2 -- so a pair already closer than that cannot be beaten by the
!     other 26, and they do not have to be looked at. Every pair in a neighbour
!     list is closer than that on any cell wider than twice the cutoff, which
!     is the case the 27 were costing the most in.
         dist_temp(1:3) = dist(1:3) + dfloat(indices(1))*a(1:3) + dfloat(indices(2))*b(1:3) + dfloat(indices(3))*c(1:3)
         res = dot_product(dist_temp, dist_temp)
         if (res < res_near) then
            dist_opt = dist_temp
         else
!     We bruteforce the integer solution among the 27 points surrounding the real solution
            res_opt = 1.d10
            do i = indices(1) - 1, indices(1) + 1
               do j = indices(2) - 1, indices(2) + 1
                  do k = indices(3) - 1, indices(3) + 1
                     dist_temp(1:3) = dist(1:3) + dfloat(i)*a(1:3) + dfloat(j)*b(1:3) + dfloat(k)*c(1:3)
                     res = dot_product(dist_temp, dist_temp)
                     if (res < res_opt) then
                        res_opt = res
                        dist_opt = dist_temp
                     end if
                  end do
               end do
            end do
         end if
         dist = dist_opt
         d = dsqrt(dot_product(dist, dist))
      else
         write (*, *) "Sorry, non-orthorhombic unit cells only work in combination with full PBC"
         stop
      end if

      return
   end subroutine get_distance

!
! This subroutine returns the number of primitive unit cells required to
! construct a supercell whose unit cell's planes are at least 2*rcut apart.
! This construct allows us to search for all the neighbors within the given
! cutoff. indices(1:3) tell the user how many repetitions are required to
! construct a unit cell with the properties outlined above. (1,1,1) means
! that the primitive unit cell is already enough.
   subroutine number_of_unit_cells_for_given_cutoff(a, b, c, rcut, PBC, indices)

      implicit none

      real(dp), intent(in) :: a(1:3)
      real(dp), intent(in) :: b(1:3)
      real(dp), intent(in) :: c(1:3)
      real(dp), intent(in) :: rcut
      logical, intent(in) :: PBC(1:3)
      integer, intent(out) :: indices(1:3)
      real(dp) :: axb(1:3)
      real(dp) :: bxc(1:3)
      real(dp) :: cxa(1:3)
      real(dp) :: indices_real(1:3)
      real(dp) :: mat(1:3, 1:3)
      real(dp) :: mat_inv(1:3, 1:3)
      real(dp) :: md
      integer :: i

      axb = cross_product(a, b)
      axb = axb/dsqrt(dot_product(axb, axb))
      cxa = cross_product(c, a)
      cxa = cxa/dsqrt(dot_product(cxa, cxa))
      bxc = cross_product(b, c)
      bxc = bxc/dsqrt(dot_product(bxc, bxc))

!   Matrix elements
      mat(1, 1) = dot_product(axb, a)
      mat(1, 2) = dot_product(axb, b)
      mat(1, 3) = dot_product(axb, c)
      mat(2, 1) = dot_product(cxa, a)
      mat(2, 2) = dot_product(cxa, b)
      mat(2, 3) = dot_product(cxa, c)
      mat(3, 1) = dot_product(bxc, a)
      mat(3, 2) = dot_product(bxc, b)
      mat(3, 3) = dot_product(bxc, c)
!   We compute the inverse of this matrix analytically
      md = -mat(1, 3)*mat(2, 2)*mat(3, 1) + mat(1, 2)*mat(2, 3)*mat(3, 1) + mat(1, 3)*mat(2, 1)*mat(3, 2) &
           - mat(1, 1)*mat(2, 3)*mat(3, 2) - mat(1, 2)*mat(2, 1)*mat(3, 3) + mat(1, 1)*mat(2, 2)*mat(3, 3)
      mat_inv(1, 1) = mat(2, 2)*mat(3, 3) - mat(2, 3)*mat(3, 2)
      mat_inv(1, 2) = mat(1, 3)*mat(3, 2) - mat(1, 2)*mat(3, 3)
      mat_inv(1, 3) = mat(1, 2)*mat(2, 3) - mat(1, 3)*mat(2, 2)
      mat_inv(2, 1) = mat(2, 3)*mat(3, 1) - mat(2, 1)*mat(3, 3)
      mat_inv(2, 2) = mat(1, 1)*mat(3, 3) - mat(1, 3)*mat(3, 1)
      mat_inv(2, 3) = mat(1, 3)*mat(2, 1) - mat(1, 1)*mat(2, 3)
      mat_inv(3, 1) = mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1)
      mat_inv(3, 2) = mat(1, 2)*mat(3, 1) - mat(1, 1)*mat(3, 2)
      mat_inv(3, 3) = mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1)
      mat_inv = mat_inv/md

      indices_real = 2.d0*(/rcut, rcut, rcut/)
      indices_real = matmul(mat_inv, indices_real)

      indices_real = dabs(indices_real)
      do i = 1, 3
         indices(i) = int(indices_real(i))
         if (indices_real(i) > dfloat(indices(i))) then
            indices(i) = indices(i) + 1
         end if
      end do

   end subroutine

! This subroutine reads in the XYZ file and builds the lists of neighbors, the spherical
! coordinates, etc.
   subroutine build_neighbors_list(positions, a_box, b_box, c_box, do_timing, &
                                   species_supercell, rcut_max, n_atom_pairs, rjs, &
                                   thetas, phis, xyz, n_neigh, neighbors_list, neighbor_species, &
                                   n_sites, indices, rebuild_neighbors_list, do_list, rank)

      implicit none

      real(dp), intent(in) :: rcut_max
      real(dp), intent(in) :: positions(:, :)
!    integer, intent(in) :: species_multiplicity(:), n_species
      integer, intent(in) :: species_supercell(:)
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: n_sites
      integer, intent(in) :: rank
      logical, intent(in) :: do_timing
      logical, intent(in) :: rebuild_neighbors_list
      logical, intent(in) :: do_list(:)

      integer, intent(out) :: n_atom_pairs
!    logical, allocatable, intent(out) :: mask_species(:,:)

      real(dp), allocatable, intent(inout) :: rjs(:)
      real(dp), allocatable, intent(inout) :: thetas(:)
      real(dp), allocatable, intent(inout) :: phis(:)
      real(dp), allocatable, intent(inout) :: xyz(:, :)
      integer, allocatable, intent(inout) :: neighbor_species(:)
      real(dp), intent(inout) :: a_box(1:3)
      real(dp), intent(inout) :: b_box(1:3)
      real(dp), intent(inout) :: c_box(1:3)
      integer, allocatable, intent(inout) :: neighbors_list(:)
      integer, allocatable, intent(inout) :: n_neigh(:)
!    integer, allocatable, intent(inout) :: species_supercell(:,:)

      real(dp) :: time1
      real(dp) :: time2
      real(dp) :: dist(1:3)
      real(dp) :: d
      real(dp) :: neigh_time
      real(dp) :: time3
      real(dp) :: tol
      real(dp) :: d_tol = 1.d-6
!     Rows that turn a Cartesian position into a fractional one, for the cell
!     list the non-orthorhombic branch bins with, and the fractional coordinate
!     of every position: the binning needs it and so does the cheap rejection
!     in the search, which is per candidate rather than per site.
      real(dp) :: recip(1:3, 1:3)
      real(dp), allocatable :: frac(:, :)
      real(dp) :: w_cell(1:3)
      real(dp) :: ds(1:3)
      real(dp) :: rcut_filter
      real(dp) :: vol
      integer, allocatable :: head(:)
      integer, allocatable :: this_list(:)
!     One column per site, holding that site's neighbours until the offsets
!     exist. The search has to finish before the count is known, and a column
!     of its own is what lets a site run without a counter shared with the
!     others.
      integer, allocatable :: scratch(:, :)
      integer :: cap
!     The largest count any build has needed. A system whose first build did
!     not fit does not search twice on every later one.
      integer, save :: cap_hint = 0
!     Where each site's pairs begin. A prefix sum over n_neigh, so both the
!     copy out of the scratch columns and the per-step geometry loop index
!     rather than count, and neither carries state from one site to the next.
      integer, allocatable :: k2_start(:)
!     Neighbours found for one site.
      integer :: nn
!     First and last site this rank owns.
      integer :: i_lo
      integer :: i_hi
      integer :: i
      integer :: j
      integer :: n_sites_supercell
      integer :: k
      integer :: k2
      integer :: i2
      integer :: j2
      integer :: i3
      integer :: j3
      integer :: k3
      integer :: mx
      integer :: my
      integer :: mz
      integer :: i_shift(1:3)
      logical :: is_box_square
      logical :: is_box_small
      logical, save :: print_cutoff_warning = .true.
      logical, save :: print_shape_warning = .true.

      if (do_timing) then
!        get_time, not cpu_time. cpu_time returns processor time summed over
!        threads, so the moment the loops below became OpenMP it stopped
!        reporting how long they took and started reporting how much work they
!        did -- the build appeared to get slower on more threads while the
!        wall-clock bucket around it fell.
         call get_time(time1)
         time3 = time1
      end if

      n_sites_supercell = size(positions, 2)

!   Very inefficiently build neighbor lists. I should write some routine that performs                              <-- FIX THIS
!   overlapping domain decomposition to make this more efficient
!    if( a_box(2) == 0.d0 .and. a_box(3) == 0.d0 .and. b_box(1) == 0.d0 .and. &
!        b_box(3) == 0.d0 .and. c_box(1) == 0.d0 .and. c_box(2) == 0.d0 )then
!        .and. size(positions,2) == n_sites )then
!   This assumes that if the non-diagonal components of the lattice vectors are approximately zero,
!   it is due to numerical noise, and thus the box is square
      if (dabs(a_box(2)) < d_tol .and. dabs(a_box(3)) < d_tol .and. &
          dabs(b_box(1)) < d_tol .and. dabs(b_box(3)) < d_tol .and. &
          dabs(c_box(1)) < d_tol .and. dabs(c_box(2)) < d_tol) then
         a_box(2) = 0.d0
         a_box(3) = 0.d0
         b_box(1) = 0.d0
         b_box(3) = 0.d0
         c_box(1) = 0.d0
         c_box(2) = 0.d0
         is_box_square = .true.
      else
         is_box_square = .false.
         if (rank == 0 .and. print_shape_warning) then
            print_shape_warning = .false.
            write (*, *) '                                       |'
            write (*, *) 'WARNING: your simulation box is not    |  <-- WARNING'
            write (*, *) 'orthorhombic; this will lead to slower |'
            write (*, *) 'code execution. This warning will be   |'
            write (*, *) 'printed only once. If you have more    |'
            write (*, *) 'than one structure in your XYZ file,   |'
            write (*, *) 'you may also have several instances of |'
            write (*, *) 'non-orthorhombic cells.                |'
            write (*, *) '                                       |'
         end if
      end if
!   This is also inefficient <---------------------------------------- FIX THIS
      is_box_small = .false.
      if (any(indices > 1)) then
         is_box_small = .true.
         if (rank == 0 .and. print_cutoff_warning) then
            print_cutoff_warning = .false.
            write (*, *) '                                       |'
            write (*, *) 'WARNING: your simulation box is smaller|  <-- WARNING'
            write (*, *) 'than a neighbors cutoff sphere; this   |'
            write (*, *) 'will lead to slow code execution due to|'
            write (*, *) 'inefficient neighbor list builds. This |'
            write (*, *) 'warning will be printed only once. If  |'
            write (*, *) 'you have more than one structure in    |'
            write (*, *) 'your XYZ file, you may also have       |'
            write (*, *) 'several instances of non-orthorhombic  |'
            write (*, *) 'cells.                                 |'
            write (*, *) '                                       |'
         end if
      end if
!   The list used to be allocated at a guessed 100 neighbours per atom and grown
!   by 10 whenever that ran out, copying the whole array each time -- GST at a
!   5.5 A cutoff needs 117, so it copied twice on every build. It is now sized
!   exactly, from the counts the one search measures.
      if (rebuild_neighbors_list) then
         allocate (n_neigh(1:n_sites))
         n_neigh = 0
         n_atom_pairs = 0
!        The span of sites this rank owns. The scratch columns cover that and
!        not every site, or on many ranks the buffer would scale with the rank
!        count instead of with the work.
         i_lo = n_sites + 1
         i_hi = 0
         do i = 1, n_sites
            if (do_list(i)) then
               if (i < i_lo) i_lo = i
               i_hi = i
            end if
         end do
      end if
      allocate (k2_start(1:n_sites))
!   We have an efficient algorithm for square boxes and inefficient for non-square boxes (sorry!)
!   Another requirement is that the minimum unit cell length is at least twice the cutoff
!
!   Tolerance in Angstrom for the ratio of positions to lattice vector. We add this because otherwise
!   atoms right at the periodic boundary can become problematic
      tol = 1.d-10
      if (is_box_square .and. (.not. is_box_small) .and. rebuild_neighbors_list) then
         mx = int(a_box(1)/rcut_max)
         my = int(b_box(2)/rcut_max)
         mz = int(c_box(3)/rcut_max)
         allocate (head(1:mx*my*mz))
         head = 0
         allocate (this_list(1:n_sites))
         do i = 1, n_sites
            call get_distance([a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0], positions(1:3, i), &
                              a_box(1:3), b_box(1:3), c_box(1:3), (/.true., .true., .true./), dist, d, i_shift)
!       This is the position within the supercell, we must make sure it really is within the supercell
            dist = dist + [a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0]
            j = 1 + modulo(int(dist(1)/(a_box(1) + tol)*mx), mx) &
                + modulo(int(dist(2)/(b_box(2) + tol)*my), my)*mx &
                + modulo(int(dist(3)/(c_box(3) + tol)*mz), mz)*my*mx
            this_list(i) = head(j)
            head(j) = i
         end do
!        One traversal, each site into its own column. The offsets cannot be
!        known before the search -- that is what the second pass used to be
!        for -- so the columns are sized from the mean density and the search
!        repeated if a site did not fit. That repeat is the old two-pass cost,
!        once, and not again for this system.
         vol = a_box(1)*b_box(2)*c_box(3)
         cap = max(cap_hint, neighbor_capacity(n_sites, vol, rcut_max))
         do
            allocate (scratch(1:cap, i_lo:i_hi))
            !$omp parallel do default(shared) schedule(dynamic, 32) &
            !$omp private(i, j, k, nn, i2, j2, k2, i3, j3, k3, dist, d, i_shift)
            do i = 1, n_sites
               if (.not. do_list(i)) cycle
!              We always count atom i as its own neighbor. This is useful when building the derivatives
               nn = 1
               scratch(1, i) = i
!              Cell coordinates for this atom
               call get_distance([a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0], positions(1:3, i), &
                                 a_box(1:3), b_box(1:3), c_box(1:3), (/.true., .true., .true./), dist, d, i_shift)
               dist = dist + [a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0]
               i2 = 1 + int(dist(1)/(a_box(1) + tol)*mx)
               j2 = 1 + int(dist(2)/(b_box(2) + tol)*my)
               k2 = 1 + int(dist(3)/(c_box(3) + tol)*mz)
!              Look for other atoms in this and neighboring cells
               do k3 = k2 - 1, k2 + 1
                  if (mz == 1 .and. k3 /= 1) cycle
                  if (mz == 2 .and. k2 == 1 .and. k3 == 0) cycle
                  if (mz == 2 .and. k2 == 2 .and. k3 == 3) cycle
                  do j3 = j2 - 1, j2 + 1
                     if (my == 1 .and. j3 /= 1) cycle
                     if (my == 2 .and. j2 == 1 .and. j3 == 0) cycle
                     if (my == 2 .and. j2 == 2 .and. j3 == 3) cycle
                     do i3 = i2 - 1, i2 + 1
                        if (mx == 1 .and. i3 /= 1) cycle
                        if (mx == 2 .and. i2 == 1 .and. i3 == 0) cycle
                        if (mx == 2 .and. i2 == 2 .and. i3 == 3) cycle
                        j = 1 + modulo(i3 - 1, mx) + modulo(j3 - 1, my)*mx + modulo(k3 - 1, mz)*mx*my
                        k = head(j)
                        do while (k /= 0)
                           if (k /= i) then
                              call get_distance(positions(1:3, i), positions(1:3, k), a_box(1:3), b_box(1:3), &
                                                c_box(1:3), (/.true., .true., .true./), dist, d, i_shift)
                              if (d < rcut_max) then
                                 nn = nn + 1
                                 if (nn <= cap) scratch(nn, i) = k
                              end if
                           end if
                           k = this_list(k)
                        end do
                     end do
                  end do
               end do
               n_neigh(i) = nn
            end do
            !$omp end parallel do
            if (maxval(n_neigh) <= cap) exit
!           The count it measured, so the repeat cannot overflow in turn, and
!           the slack so a later build needing one more does not repeat again.
            cap = maxval(n_neigh) + NEIGHBOR_CAP_SLACK
            deallocate (scratch)
         end do
         deallocate (head, this_list)
!   The same cell list for every other cell: any lattice with a non-zero
!   off-diagonal component, and any cell smaller than the cutoff sphere. What
!   this replaces tested every site against every supercell position, twice per
!   rebuild, through a minimum image that brute-forces 27 images per call. The
!   glassy-carbon MAD case is 2912 atoms in a hexagonal cell, so it spent 19.7 s
!   of a 25.4 s GH200 run in here against 3.3 s in the descriptor it feeds.
!
!   Binning is in fractional coordinates, so the wrap is modular arithmetic on
!   bin indices and the cell shape never enters. One bin either way is enough:
!   a pair inside rcut has |ds_k| < rcut/w_k for the perpendicular width w_k,
!   and m_k = floor(w_k/rcut) bins make that less than one bin. Which pairs are
!   accepted does not change -- this only decides which ones are tested.
      else if (rebuild_neighbors_list) then
         call cell_grid(a_box, b_box, c_box, rcut_max, n_sites_supercell, recip, w_cell, vol, mx, my, mz)
!        Slack so the bound below can never reject a pair the cutoff test would
!        have taken. The bound is exact in exact arithmetic; 1.d-10 A is far
!        above the rounding and far below anything physical.
         rcut_filter = rcut_max + 1.d-10
         allocate (head(1:mx*my*mz))
         head = 0
         allocate (this_list(1:n_sites_supercell))
         allocate (frac(1:3, 1:n_sites_supercell))
         do i = 1, n_sites_supercell
            frac(1, i) = modulo(dot_product(recip(1, 1:3), positions(1:3, i)), 1.d0)
            frac(2, i) = modulo(dot_product(recip(2, 1:3), positions(1:3, i)), 1.d0)
            frac(3, i) = modulo(dot_product(recip(3, 1:3), positions(1:3, i)), 1.d0)
            i2 = min(mx, 1 + int(frac(1, i)*dfloat(mx)))
            j2 = min(my, 1 + int(frac(2, i)*dfloat(my)))
            k2 = min(mz, 1 + int(frac(3, i)*dfloat(mz)))
            j = i2 + (j2 - 1)*mx + (k2 - 1)*mx*my
            this_list(i) = head(j)
            head(j) = i
         end do
         cap = max(cap_hint, neighbor_capacity(n_sites, vol, rcut_max))
         do
            allocate (scratch(1:cap, i_lo:i_hi))
            !$omp parallel do default(shared) schedule(dynamic, 32) &
            !$omp private(i, j, k, nn, i2, j2, k2, i3, j3, k3, dist, d, i_shift, ds)
            do i = 1, n_sites
               if (.not. do_list(i)) cycle
!              We always count atom i as its own neighbor. This is useful when building the derivatives
               nn = 1
               scratch(1, i) = i
               i2 = min(mx, 1 + int(frac(1, i)*dfloat(mx)))
               j2 = min(my, 1 + int(frac(2, i)*dfloat(my)))
               k2 = min(mz, 1 + int(frac(3, i)*dfloat(mz)))
               do k3 = k2 - 1, k2 + 1
                  if (mz == 1 .and. k3 /= 1) cycle
                  if (mz == 2 .and. k2 == 1 .and. k3 == 0) cycle
                  if (mz == 2 .and. k2 == 2 .and. k3 == 3) cycle
                  do j3 = j2 - 1, j2 + 1
                     if (my == 1 .and. j3 /= 1) cycle
                     if (my == 2 .and. j2 == 1 .and. j3 == 0) cycle
                     if (my == 2 .and. j2 == 2 .and. j3 == 3) cycle
                     do i3 = i2 - 1, i2 + 1
                        if (mx == 1 .and. i3 /= 1) cycle
                        if (mx == 2 .and. i2 == 1 .and. i3 == 0) cycle
                        if (mx == 2 .and. i2 == 2 .and. i3 == 3) cycle
                        j = 1 + modulo(i3 - 1, mx) + modulo(j3 - 1, my)*mx + modulo(k3 - 1, mz)*mx*my
                        k = head(j)
                        do while (k /= 0)
                           if (k /= i) then
!                             |dr| >= |ds_k| w_k for every k, because a.g1 = w1
!                             while b.g1 = c.g1 = 0, so dr's component along g1
!                             is exactly ds_1 w_1. Rounding ds to the nearest
!                             integer minimises each |ds_k| on its own, which
!                             bounds the minimum over every image from below --
!                             so this rejects only pairs the 27-image search
!                             would also have rejected, at a twentieth of its
!                             cost, and it is the reason the search is affordable
!                             on a cell this shape at all.
                              ds(1:3) = frac(1:3, k) - frac(1:3, i)
                              ds(1:3) = ds(1:3) - dnint(ds(1:3))
                              if (dabs(ds(1))*w_cell(1) < rcut_filter .and. &
                                  dabs(ds(2))*w_cell(2) < rcut_filter .and. &
                                  dabs(ds(3))*w_cell(3) < rcut_filter) then
                                 call get_distance(positions(1:3, i), positions(1:3, k), a_box(1:3), b_box(1:3), &
                                                   c_box(1:3), (/.true., .true., .true./), dist, d, i_shift)
                                 if (d < rcut_max) then
                                    nn = nn + 1
                                    if (nn <= cap) scratch(nn, i) = k
                                 end if
                              end if
                           end if
                           k = this_list(k)
                        end do
                     end do
                  end do
               end do
               n_neigh(i) = nn
!              Ascending index, which is the order the loop over every position
!              produced. Not cosmetic: permuting a site's neighbours
!              reassociates every sum downstream and moves the last digit of
!              every force, so without this the change is not verifiable.
               if (nn > 2 .and. nn <= cap) call sort_ascending(scratch(2:nn, i))
            end do
            !$omp end parallel do
            if (maxval(n_neigh) <= cap) exit
            cap = maxval(n_neigh) + NEIGHBOR_CAP_SLACK
            deallocate (scratch)
         end do
         deallocate (head, this_list, frac)
      end if

!   The offsets, then one copy per site out of its column. Each slice is
!   written exactly once, so there is nothing shared between sites here either.
      if (rebuild_neighbors_list) then
         call size_the_list(n_neigh, n_sites, k2_start, n_atom_pairs, neighbors_list)
         !$omp parallel do default(shared) schedule(static) private(i)
         do i = 1, n_sites
            if (do_list(i)) then
               neighbors_list(k2_start(i):k2_start(i) + n_neigh(i) - 1) = scratch(1:n_neigh(i), i)
            end if
         end do
         !$omp end parallel do
         deallocate (scratch)
         cap_hint = cap
      end if

      if (do_timing) then
         call get_time(time2)
         neigh_time = time2 - time1
         time1 = time2
      end if

!   NOTE on performance: looping over interactions I could have chosen to calculate each interaction only once,
!   i.e., (i,j) = (j,i), because if j is i's neighbor, the reverse is also true. This would reduce the
!   calculation load by half. However, this is a design feature, since the code is easier to parallelize if
!   each atom is accompanied by its full list of neighbors. It also prevents ackward memory access.
!   In any case, I may want to rethink this in the future if I want to further gain extra performance (at the
!   expense of complicating the code, that is).
      n_atom_pairs = 0
      do i = 1, n_sites
         n_atom_pairs = n_atom_pairs + n_neigh(i)
      end do
      if (rebuild_neighbors_list) then
         allocate (rjs(1:n_atom_pairs))
         allocate (xyz(1:3, 1:n_atom_pairs))
         allocate (thetas(1:n_atom_pairs))
         allocate (phis(1:n_atom_pairs))
         allocate (neighbor_species(1:n_atom_pairs))
      end if
!   Offsets for the loop below. Recomputed rather than carried, because on a
!   step that did not rebuild there is nothing to carry them from.
      if (n_sites > 0) then
         k2_start(1) = 1
         do i = 2, n_sites
            k2_start(i) = k2_start(i - 1) + n_neigh(i - 1)
         end do
      end if

!   The geometry of every pair, every step, whether or not the topology above
!   was rebuilt -- which is why this and not the build is the loop that matters
!   once the Verlet skin is doing its job. Each k2 is written exactly once, so
!   with the offsets in hand there is nothing shared between sites.
      !$omp parallel do default(shared) schedule(static) &
      !$omp private(i, j, k, k2, dist, d, i_shift)
      do i = 1, n_sites
         if (do_list(i)) then
            k2 = k2_start(i) - 1
            do k = 1, n_neigh(i)
               k2 = k2 + 1
               j = neighbors_list(k2)
               if (k == 1) then
                  rjs(k2) = 0.d0
                  xyz(1:3, k2) = (/0.d0, 0.d0, 0.d0/)
                  thetas(k2) = 0.d0
                  phis(k2) = 0.d0
                  neighbor_species(k2) = species_supercell(i)
               else
                  call get_distance(positions(1:3, i), positions(1:3, j), a_box(1:3), b_box(1:3), &
                                    c_box(1:3), (/.true., .true., .true./), dist, d, i_shift)
                  rjs(k2) = d
                  xyz(1:3, k2) = dist
!           Avoid numerical artifacts
                  if (dabs(dist(3)) >= d) then
                     if (dist(3) > 0.d0) then
                        thetas(k2) = 0.d0
                     else
                        thetas(k2) = dacos(-1.d0)
                     end if
                  else
                     thetas(k2) = dacos(dist(3)/d)
                  end if
                  phis(k2) = datan2(dist(2), dist(1))
                  neighbor_species(k2) = species_supercell(j)
               end if
            end do
         end if
      end do
      !$omp end parallel do
      deallocate (k2_start)

      if (do_timing) then
         call get_time(time2)
         write (*, *) '                                       |'
         write (*, *) 'Atoms timings (build):                 |'
         write (*, *) '                                       |'
         write (*, '(A, F9.3, A)') '  *) Neighbors build: ', neigh_time, ' seconds |'
         write (*, '(A, F7.3, A)') '  *) Spherical coords.: ', time2 - time1, ' seconds |'
         write (*, '(A, F19.3, A)') '  *) Total: ', time2 - time3, ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if

   end subroutine

!  What one atom pair and one site cost the SOAP descriptor path, in bytes.
!
!  Every term is one allocation in get_soap (src/soap_turbo/src/) or in
!  get_gap_soap (src/gap_interface.f90), with do_derivatives on. Element sizes:
!  complex(dp) 16, real(dp) 8, integer 4, logical 4.
!
!  PER PAIR
!    48*k_max*n_max   cnk_rad_der, cnk_azi_der, cnk_pol_der
!    64*k_max         angular_exp_coeff and its rad/azi/pol derivatives
!    48*n_soap        soap_rad_der, soap_azi_der, soap_pol_der, soap_cart_der
!    16*n_max         radial_exp_coeff, radial_exp_coeff_der
!     4*n_species     mask
!    56               rjs, thetas, phis, xyz, neighbors_list, in_to_out_pairs
!
!  PER SITE
!    16*k_max*n_max   cnk
!     8*n_soap        soap
!    16               sqrt_dot_p, n_neigh, species_multiplicity
!
!  These numbers are this branch's, and differ from the GPU branch's for a
!  reason worth stating: there, get_derivatives is never called, so the
!  coefficient and derivative arrays are dummies and the real allocations are on
!  the device. Here they are live host arrays. The two branches enumerate what
!  they each allocate; the coefficients must not be copied between them.
!
!  Per-site is not a rounding correction. cnk is the same order as the per-pair
!  total divided by the neighbour count, and it decides the answer for a short
!  cutoff, where there are few pairs per site.
!
!  The fractional-coordinate rows, and how many bins fit along each lattice
!  direction. w_k = V/|b x c| and its cyclic partners is the distance between
!  the two cell planes normal to k; for an orthorhombic cell that is a(1), b(2),
!  c(3), so this gives the same counts the orthorhombic branch works out itself.
   subroutine cell_grid(a, b, c, rcut, n_atoms, recip, w, vol, mx, my, mz)

      implicit none

      real(dp), intent(in) :: a(1:3)
      real(dp), intent(in) :: b(1:3)
      real(dp), intent(in) :: c(1:3)
      real(dp), intent(in) :: rcut
      integer, intent(in) :: n_atoms
      real(dp), intent(out) :: recip(1:3, 1:3)
!     Perpendicular widths, one per lattice direction: the distance between the
!     two cell planes normal to it. They set the bin counts below and they are
!     what the cheap rejection in the search is a bound on.
      real(dp), intent(out) :: w(1:3)
      real(dp), intent(out) :: vol
      integer, intent(out) :: mx
      integer, intent(out) :: my
      integer, intent(out) :: mz
      real(dp) :: bxc(1:3)
      real(dp) :: cxa(1:3)
      real(dp) :: axb(1:3)
      integer(int64) :: n_bins

      bxc = cross_product(b, c)
      cxa = cross_product(c, a)
      axb = cross_product(a, b)
      vol = dot_product(a, bxc)
      recip(1, 1:3) = bxc/vol
      recip(2, 1:3) = cxa/vol
      recip(3, 1:3) = axb/vol
!     max(1, ...) rather than an error: a cell thinner than the cutoff is legal
!     and reaches here through is_box_small, and one bin on that axis degrades
!     to testing every position along it, which is what used to happen anyway.
      w(1) = dabs(vol)/dsqrt(dot_product(bxc, bxc))
      w(2) = dabs(vol)/dsqrt(dot_product(cxa, cxa))
      w(3) = dabs(vol)/dsqrt(dot_product(axb, axb))
      mx = max(1, int(w(1)/rcut))
      my = max(1, int(w(2)/rcut))
      mz = max(1, int(w(3)/rcut))
!     A cluster in a large vacuum box asks for one bin per rcut^3 of empty
!     space, and more bins than atoms buys nothing. Halving an axis only makes
!     bins bigger, which the one-bin search is still correct for. int64 because
!     the product of three unclamped counts is what overflows first.
      n_bins = int(mx, int64)*int(my, int64)*int(mz, int64)
      do while (n_bins > 8_int64*int(max(1, n_atoms), int64) .and. max(mx, my, mz) > 1)
         if (mx >= my .and. mx >= mz) then
            mx = max(1, mx/2)
         else if (my >= mz) then
            my = max(1, my/2)
         else
            mz = max(1, mz/2)
         end if
         n_bins = int(mx, int64)*int(my, int64)*int(mz, int64)
      end do

   end subroutine cell_grid

!  Insertion sort, on one site's neighbours: a hundred or so entries at the
!  cutoffs this code runs at, where the setup for anything cleverer costs more
!  than the sort does.
   subroutine sort_ascending(list)

      implicit none

      integer, intent(inout) :: list(:)
      integer :: i
      integer :: j
      integer :: v

      do i = 2, size(list)
         v = list(i)
         j = i - 1
         do while (j >= 1)
            if (list(j) <= v) exit
            list(j + 1) = list(j)
            j = j - 1
         end do
         list(j + 1) = v
      end do

   end subroutine sort_ascending

!  Room for one site's neighbours, before the search has counted any. A site in
!  a uniform density has 4/3 pi rcut^3 n/V of them; a surface, an interface or a
!  cluster in vacuum has more than its box average, which the margin covers
!  when it is small and the repeat covers when it is not.
   integer function neighbor_capacity(n_atoms, vol, rcut) result(cap)

      implicit none

      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: vol
      real(dp), intent(in) :: rcut
      real(dp) :: mean

      mean = 4.d0/3.d0*dacos(-1.d0)*rcut**3*dfloat(max(0, n_atoms))/dmax1(dabs(vol), 1.d-10)
!     At least the one slot the site's own index takes.
      cap = max(1, NEIGHBOR_CAP_SLACK + int(min(1.25d0*mean, 1.d6)))

   end function neighbor_capacity

!  Turn the per-site neighbour counts into the offsets the copy writes at, and
!  size the list to exactly the number of pairs there are. Sites this rank does
!  not own have n_neigh = 0 and take no room.
   subroutine size_the_list(n_neigh, n_sites, k2_start, n_atom_pairs, neighbors_list)

      implicit none

      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: n_sites
      integer, intent(out) :: k2_start(:)
      integer, intent(out) :: n_atom_pairs
      integer, allocatable, intent(inout) :: neighbors_list(:)
      integer :: i

      if (n_sites < 1) then
         n_atom_pairs = 0
      else
         k2_start(1) = 1
         do i = 2, n_sites
            k2_start(i) = k2_start(i - 1) + n_neigh(i - 1)
         end do
         n_atom_pairs = k2_start(n_sites) + n_neigh(n_sites) - 1
      end if
      if (allocated(neighbors_list)) deallocate (neighbors_list)
      allocate (neighbors_list(1:n_atom_pairs))

   end subroutine size_the_list
#ifdef _GPU
!  What one atom pair and one site cost the SOAP descriptor path, in bytes.
!
!  Every term below is one allocation in get_soap (src/soap_turbo/src/) or in
!  get_gap_soap (src/gap_interface.f90), with do_derivatives on. Element sizes:
!  complex(dp) 16, real(dp) 8, integer 4, logical 4.
!
!  PER PAIR
!    48*k_max*n_max   cnk_rad_der_d, cnk_azi_der_d, cnk_pol_der_d
!    64*k_max         angular_exp_coeff_d and its rad/azi/pol derivatives
!    72*n_soap        soap_{rad,azi,pol}_der_d, soap_cart_der_d, and the host
!                     soap_cart_der, which is still allocated at full size
!    40*n_max         radial_exp_coeff_d, radial_exp_coeff_der_d and the three
!                     temporaries (ntemp = maxval(alpha_max) <= n_max, so using
!                     n_max here is the conservative direction)
!     8*n_species     mask_d and its host copy
!    88               rjs/thetas/phis/xyz and the pair index arrays, host+device
!
!  PER SITE
!    16*k_max*n_max   cnk_d
!    16*n_soap        soap_d and its host copy
!    24               sqrt_dot_p_d, n_neigh_d, k2_start_d, i_k2_start_d
!
!  Per-site is not a rounding correction: cnk_d is the same order as the
!  per-pair total divided by the neighbour count, and it is the term that
!  decides the answer for a short cutoff, where there are few pairs per site.
   subroutine soap_batch_memory_model(l_max, n_max, n_soap, n_species, &
                                      bytes_per_pair, bytes_per_site)
      implicit none

      integer, intent(in) :: l_max
      integer, intent(in) :: n_max
      integer, intent(in) :: n_soap
      integer, intent(in) :: n_species
      real(dp), intent(out) :: bytes_per_pair
      real(dp), intent(out) :: bytes_per_site

      integer :: k_max

      k_max = 1 + l_max*(l_max + 1)/2 + l_max

      bytes_per_pair = 48.d0*dfloat(k_max)*dfloat(n_max) &
                       + 64.d0*dfloat(k_max) &
                       + 72.d0*dfloat(n_soap) &
                       + 40.d0*dfloat(n_max) &
                       + 8.d0*dfloat(n_species) &
                       + 88.d0

      bytes_per_site = 16.d0*dfloat(k_max)*dfloat(n_max) &
                       + 16.d0*dfloat(n_soap) &
                       + 24.d0

   end subroutine soap_batch_memory_model

   subroutine get_number_of_atom_pairs_batches(n_batches, n_neigh, rjs, rcut, l_max, n_max, n_soap, n_species, &
                                               max_Gbytes_per_process, &
                                               i_beg_list, i_end_list, j_beg_list, j_end_list)

      implicit none

      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: rcut
      real(dp), intent(in) :: max_Gbytes_per_process
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: l_max
      integer, intent(in) :: n_max
      integer, intent(in) :: n_soap
      integer, intent(in) :: n_species
      integer, intent(in) :: n_batches

      integer, allocatable, intent(out) :: i_beg_list(:)
      integer, allocatable, intent(out) :: i_end_list(:)
      integer, allocatable, intent(out) :: j_beg_list(:)
      integer, allocatable, intent(out) :: j_end_list(:)

      real(dp) :: estimated_memory_in_Gbytes
      real(dp) :: bytes_per_pair
      real(dp) :: bytes_per_site
      real(dp) :: mem_ratio
      real(dp) :: pairs_per_chunk
      integer :: n_sites
      integer :: n_atom_pairs
      integer :: k_max
      integer :: n_chunks
      integer :: i
      integer :: j
      integer :: k
      integer :: k2
      integer :: i_chunk
      integer :: n_atom_pairs_in

#else

      subroutine soap_batch_memory_model(l_max, n_max, n_soap, n_species, &
                                         bytes_per_pair, bytes_per_site)
         implicit none

         integer, intent(in) :: l_max
         integer, intent(in) :: n_max
         integer, intent(in) :: n_soap
         integer, intent(in) :: n_species
         real(dp), intent(out) :: bytes_per_pair
         real(dp), intent(out) :: bytes_per_site

         integer :: k_max

         k_max = 1 + l_max*(l_max + 1)/2 + l_max

         bytes_per_pair = 48.d0*dfloat(k_max)*dfloat(n_max) &
                          + 64.d0*dfloat(k_max) &
                          + 48.d0*dfloat(n_soap) &
                          + 16.d0*dfloat(n_max) &
                          + 4.d0*dfloat(n_species) &
                          + 56.d0

         bytes_per_site = 16.d0*dfloat(k_max)*dfloat(n_max) &
                          + 8.d0*dfloat(n_soap) &
                          + 16.d0

      end subroutine soap_batch_memory_model

      subroutine get_number_of_atom_pairs(n_neigh, rjs, rcut, l_max, n_max, n_soap, n_species, &
                                          max_Gbytes_per_process, &
                                          i_beg_list, i_end_list, j_beg_list, j_end_list)

         implicit none

         real(dp), intent(in) :: rjs(:)
         real(dp), intent(in) :: rcut
         real(dp), intent(in) :: max_Gbytes_per_process
         integer, intent(in) :: n_neigh(:)
         integer, intent(in) :: l_max
         integer, intent(in) :: n_max
         integer, intent(in) :: n_soap
         integer, intent(in) :: n_species

         integer, allocatable, intent(out) :: i_beg_list(:)
         integer, allocatable, intent(out) :: i_end_list(:)
         integer, allocatable, intent(out) :: j_beg_list(:)
         integer, allocatable, intent(out) :: j_end_list(:)

         real(dp) :: estimated_memory_in_Gbytes
         real(dp) :: bytes_per_pair
         real(dp) :: bytes_per_site
         real(dp) :: mem_ratio
         real(dp) :: pairs_per_chunk
         integer :: n_sites
         integer :: n_atom_pairs
         integer :: k_max
         integer :: n_chunks
         integer :: i
         integer :: j
         integer :: k
         integer :: k2
         integer :: i_chunk
         integer :: n_atom_pairs_in

#endif
         n_sites = size(n_neigh)
         n_atom_pairs = size(rjs)

         k = 0
         n_atom_pairs_in = 0
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               k = k + 1
               if (rjs(k) < rcut) then
                  n_atom_pairs_in = n_atom_pairs_in + 1
               end if
            end do
         end do

         k_max = 1 + l_max*(l_max + 1)/2 + l_max
!   What the descriptor path will actually allocate for these pairs and sites,
!   enumerated in soap_batch_memory_model rather than folded into one constant.
         call soap_batch_memory_model(l_max, n_max, n_soap, n_species, bytes_per_pair, bytes_per_site)
         estimated_memory_in_Gbytes = SOAP_BATCH_SAFETY* &
                                      (dfloat(n_atom_pairs_in)*bytes_per_pair &
                                       + dfloat(n_sites)*bytes_per_site)/1024.d0**3
         mem_ratio = estimated_memory_in_Gbytes/max_Gbytes_per_process
#ifdef _GPU
         n_chunks = n_batches
#else
         n_chunks = ceiling(mem_ratio)
#endif
         if (n_chunks > n_sites) then
            n_chunks = n_sites
         end if

         if (n_chunks > 0) &
            pairs_per_chunk = dfloat(n_atom_pairs_in)/dfloat(n_chunks)

         allocate (i_beg_list(1:n_chunks))
         allocate (i_end_list(1:n_chunks))
         allocate (j_beg_list(1:n_chunks))
         allocate (j_end_list(1:n_chunks))

         if (n_chunks == 0) then
            return
         end if

         i_beg_list(1) = 1
         j_beg_list(1) = 1
         i_end_list(n_chunks) = n_sites
         j_end_list(n_chunks) = n_atom_pairs

         if (n_chunks == 1) then
            return
         end if

         k = 0
         k2 = 0
         i_chunk = 1
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               k = k + 1
               if (rjs(k) < rcut) then
                  k2 = k2 + 1
               end if
            end do
            if (k2 >= int(float(i_chunk)*pairs_per_chunk)) then
               i_end_list(i_chunk) = i
               j_end_list(i_chunk) = k
               i_chunk = i_chunk + 1
               i_beg_list(i_chunk) = i + 1
               j_beg_list(i_chunk) = k + 1
               if (i_chunk == n_chunks) then
                  exit
               end if
            end if
         end do
         return

      end subroutine
#ifdef _GPU
      subroutine get_number_of_atom_pairs(n_neigh, rjs, rcut, l_max, n_max, n_soap, n_species, &
                                          max_Gbytes_per_process, &
                                          i_beg_list, i_end_list, j_beg_list, j_end_list)

         implicit none

         real(dp), intent(in) :: rjs(:)
         real(dp), intent(in) :: rcut
         real(dp), intent(in) :: max_Gbytes_per_process
         integer, intent(in) :: n_neigh(:)
         integer, intent(in) :: l_max
         integer, intent(in) :: n_max
         integer, intent(in) :: n_soap
         integer, intent(in) :: n_species

         integer, allocatable, intent(out) :: i_beg_list(:)
         integer, allocatable, intent(out) :: i_end_list(:)
         integer, allocatable, intent(out) :: j_beg_list(:)
         integer, allocatable, intent(out) :: j_end_list(:)

         real(dp) :: estimated_memory_in_Gbytes
         real(dp) :: bytes_per_pair
         real(dp) :: bytes_per_site
         real(dp) :: mem_ratio
         real(dp) :: pairs_per_chunk
         integer :: n_sites
         integer :: n_atom_pairs
         integer :: k_max
         integer :: n_chunks
         integer :: i
         integer :: j
         integer :: k
         integer :: k2
         integer :: i_chunk
         integer :: n_atom_pairs_in

         n_sites = size(n_neigh)
         n_atom_pairs = size(rjs)

         k = 0
         n_atom_pairs_in = 0
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               k = k + 1
               if (rjs(k) < rcut) then
                  n_atom_pairs_in = n_atom_pairs_in + 1
               end if
            end do
         end do

         k_max = 1 + l_max*(l_max + 1)/2 + l_max
!   What the descriptor path will actually allocate for these pairs and sites,
!   enumerated in soap_batch_memory_model rather than folded into one constant.
         call soap_batch_memory_model(l_max, n_max, n_soap, n_species, bytes_per_pair, bytes_per_site)
         estimated_memory_in_Gbytes = SOAP_BATCH_SAFETY* &
                                      (dfloat(n_atom_pairs_in)*bytes_per_pair &
                                       + dfloat(n_sites)*bytes_per_site)/1024.d0**3
         mem_ratio = estimated_memory_in_Gbytes/max_Gbytes_per_process
         n_chunks = ceiling(mem_ratio)
         if (n_chunks > n_sites) then
            n_chunks = n_sites
         end if

         pairs_per_chunk = float(n_atom_pairs_in)/float(n_chunks)

         allocate (i_beg_list(1:n_chunks))
         allocate (i_end_list(1:n_chunks))
         allocate (j_beg_list(1:n_chunks))
         allocate (j_end_list(1:n_chunks))

         if (n_chunks == 0) then
            return
         end if

         i_beg_list(1) = 1
         j_beg_list(1) = 1
         i_end_list(n_chunks) = n_sites
         j_end_list(n_chunks) = n_atom_pairs

         if (n_chunks == 1) then
            return
         end if

         k = 0
         k2 = 0
         i_chunk = 1
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               k = k + 1
               if (rjs(k) < rcut) then
                  k2 = k2 + 1
               end if
            end do
            if (k2 >= int(float(i_chunk)*pairs_per_chunk)) then
               i_end_list(i_chunk) = i
               j_end_list(i_chunk) = k
               i_chunk = i_chunk + 1
               i_beg_list(i_chunk) = i + 1
               j_beg_list(i_chunk) = k + 1
               if (i_chunk == n_chunks) then
                  exit
               end if
            end if
         end do
         return

      end subroutine
#endif

!
! This subroutine returns the fractional coordinates from a list of
! Cartesian positions. This subroutine does NOT carry out unit cell
! wrapping. Wrapped Cartesian coordinates should be provided if wrapped
! fractional coordinates are wanted.
      subroutine get_fractional_coordinates(pos, a, b, c, frac)

         implicit none

         real(dp), intent(in) :: pos(:, :)
         real(dp), intent(in) :: a(1:3)
         real(dp), intent(in) :: b(1:3)
         real(dp), intent(in) :: c(1:3)
         real(dp), intent(out) :: frac(1:3, 1:size(pos, 2))
         real(dp) :: L(1:3)
         real(dp) :: d_tol = 1.d-6
         real(dp) :: mat(1:3, 1:3)
         real(dp) :: md
         real(dp), save :: a0(1:3) = 0.d0
         real(dp), save :: b0(1:3) = 0.d0
         real(dp), save :: c0(1:3) = 0.d0
         real(dp), save :: mat_inv(1:3, 1:3) = 0.d0
         integer :: i
         integer :: atom
         integer :: n_atoms
         logical :: lattice_check_a(1:3)
         logical :: lattice_check_b(1:3)
         logical :: lattice_check_c(1:3)

         n_atoms = size(pos, 2)

         if (dabs(a(2)) < d_tol .and. dabs(a(3)) < d_tol .and. &
             dabs(b(1)) < d_tol .and. dabs(b(3)) < d_tol .and. &
             dabs(c(1)) < d_tol .and. dabs(c(2)) < d_tol) then
!     Fast solution for orthorhombic cells
            L = (/a(1), b(2), c(3)/)
            do atom = 1, n_atoms
               do i = 1, 3
                  frac(i, atom) = pos(i, atom)/L(i)
               end do
            end do
         else
!     Slow solution for other unit cells
            lattice_check_a = (a /= a0)
            lattice_check_b = (b /= b0)
            lattice_check_c = (c /= c0)
            if (any(lattice_check_a) .or. any(lattice_check_b) .or. any(lattice_check_c)) then
               a0 = a
               b0 = b
               c0 = c
!       We construct our matrix to get the MIC only if the lattice vectors have changed
               mat(1:3, 1) = a(1:3)
               mat(1:3, 2) = b(1:3)
               mat(1:3, 3) = c(1:3)
!       We compute the inverse of this matrix analytically
           md = -mat(1,3)*mat(3,1)*mat(2,2) + mat(2,1)*mat(1,3)*mat(3,2) + mat(1,2)*mat(3,1)*mat(2,3) - mat(1,1)*mat(2,3)*mat(3,2) &
                    - mat(1, 2)*mat(2, 1)*mat(3, 3) + mat(1, 1)*mat(2, 2)*mat(3, 3)
               mat_inv(1, 1) = mat(2, 2)*mat(3, 3) - mat(2, 3)*mat(3, 2)
               mat_inv(1, 2) = mat(1, 3)*mat(3, 2) - mat(1, 2)*mat(3, 3)
               mat_inv(1, 3) = mat(1, 2)*mat(2, 3) - mat(1, 3)*mat(2, 2)
               mat_inv(2, 1) = mat(2, 3)*mat(3, 1) - mat(2, 1)*mat(3, 3)
               mat_inv(2, 2) = mat(1, 1)*mat(3, 3) - mat(1, 3)*mat(3, 1)
               mat_inv(2, 3) = mat(1, 3)*mat(2, 1) - mat(1, 1)*mat(2, 3)
               mat_inv(3, 1) = mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1)
               mat_inv(3, 2) = mat(1, 2)*mat(3, 1) - mat(1, 1)*mat(3, 2)
               mat_inv(3, 3) = mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1)
               mat_inv = mat_inv/md
            end if
            do atom = 1, n_atoms
               frac(1:3, atom) = matmul(mat_inv, pos(1:3, atom))
            end do
         end if

         return
      end subroutine get_fractional_coordinates

#ifdef _GPU
      subroutine get_gpu_batches(n_neigh, rjs, rcut, n_chunks, estimated_memory_in_Gbytes, max_Gbytes_per_process, &
                                 i_beg_list, i_end_list, j_beg_list, j_end_list)

         implicit none

         real(dp), intent(in) :: rjs(:)
         real(dp), intent(in) :: rcut
         real(dp), intent(in) :: max_Gbytes_per_process
         real(dp), intent(in) :: estimated_memory_in_Gbytes
         integer, intent(in) :: n_neigh(:)

         integer, allocatable, intent(out) :: i_beg_list(:)
         integer, allocatable, intent(out) :: i_end_list(:)
         integer, allocatable, intent(out) :: j_beg_list(:)
         integer, allocatable, intent(out) :: j_end_list(:)

         real(dp) :: mem_ratio
         real(dp) :: pairs_per_chunk
         integer :: n_sites
         integer :: n_atom_pairs
         integer :: k_max
         integer :: n_chunks
         integer :: i
         integer :: j
         integer :: k
         integer :: k2
         integer :: i_chunk
         integer :: n_atom_pairs_in

         n_sites = size(n_neigh)
         n_atom_pairs = size(rjs)

         k = 0
         n_atom_pairs_in = 0
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               k = k + 1
               if (rjs(k) < rcut) then
                  n_atom_pairs_in = n_atom_pairs_in + 1
               end if
            end do
         end do

         ! !   This is a conservative estimate of the maximum memory that this run will need

         ! !
         ! estimated_memory_in_Gbytes = dfloat(n_atom_pairs_in) * 150.d0 / 1024.d0**3
         mem_ratio = estimated_memory_in_Gbytes/max_Gbytes_per_process
!    n_chunks = ceiling(mem_ratio)
         if (n_chunks > n_sites) then
            n_chunks = n_sites
         end if

         pairs_per_chunk = dfloat(n_atom_pairs_in)/dfloat(n_chunks)

         allocate (i_beg_list(1:n_chunks))
         allocate (i_end_list(1:n_chunks))
         allocate (j_beg_list(1:n_chunks))
         allocate (j_end_list(1:n_chunks))

         if (n_chunks == 0) then
            return
         end if

         i_beg_list(1) = 1
         j_beg_list(1) = 1
         i_end_list(n_chunks) = n_sites
         j_end_list(n_chunks) = n_atom_pairs

         if (n_chunks == 1) then
            return
         end if

         k = 0
         k2 = 0
         i_chunk = 1
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               k = k + 1
               if (rjs(k) < rcut) then
                  k2 = k2 + 1
               end if
            end do
            if (k2 >= int(float(i_chunk)*pairs_per_chunk)) then
               i_end_list(i_chunk) = i
               j_end_list(i_chunk) = k
               i_chunk = i_chunk + 1
               i_beg_list(i_chunk) = i + 1
               j_beg_list(i_chunk) = k + 1
               if (i_chunk == n_chunks) then
                  exit
               end if
            end if
         end do
         return

      end subroutine get_gpu_batches

      subroutine estimate_max_exp_forces_device_memory_usage(n_sites, n_pairs, n_dim_partial, n_samples, n_samples_sf, total, &
                                                             be_verbose)
         implicit none
         integer :: n_sites
         integer :: nk
         integer :: n_samples
         integer :: n_samples_sf
         integer :: n_pairs
         integer :: n_dim_partial
         real(dp), intent(inout) :: total
!     Not initialised here: an initialiser in a declaration is an implicit SAVE,
!     and these are running sums. They are zeroed on entry instead.
         real(dp) :: total_exp
         real(dp) :: total_standard
         real(dp) :: nk_int
         real(dp) :: nk_float
         real(dp) :: Gk
         real(dp) :: dermat
         real(dp) :: sf
         real(dp) :: fi
         real(dp) :: pref
         real(dp) :: xyz
         real(dp) :: forces
         real(dp) :: neigh_list
         real(dp) :: to_gb
!     Off unless asked. This is called once per snapshot on the batched path and
!     prints twenty lines; useful the first time, noise for the rest of an MD
!     run, and multiplied by every MPI rank.
         logical, intent(in), optional :: be_verbose
         logical :: verbose

         verbose = .false.
         if (present(be_verbose)) verbose = be_verbose

         total_exp = 0.d0
         total_standard = 0.d0

         total = 0.d0

         nk = n_pairs

         total = 0.d0

         to_gb = 1/dfloat(1024**3)

         neigh_list = dfloat(n_pairs)*4.d0*to_gb

         forces = dfloat(n_sites)*8.d0*to_gb
         nk_int = dfloat(nk)*4.d0*to_gb
         nk_float = dfloat(nk)*8.d0*to_gb
         sf = dfloat(n_samples*n_samples_sf)*8.d0*to_gb

         pref = dfloat(n_samples_sf)*8.d0*to_gb
         Gk = n_samples*nk_float
         dermat = n_samples_sf*nk_float
         fi = 3*nk_float
         xyz = 3*nk_float

         if (verbose) write (*, '(A)') "> Estimating memory for normal allocations "
         total_standard = total_standard + neigh_list
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device neigh_list_d  ", neigh_list, " Gb"

         total_standard = total_standard + neigh_list
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device neigh_spec_d  ", neigh_list, " Gb"

         total_standard = total_standard + neigh_list*2.d0*3.d0
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device xyz_d         ", neigh_list*2.d0*3.d0, " Gb"

         total_standard = total_standard + neigh_list*2.d0
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device rjs_d         ", neigh_list*2.d0, " Gb"

         total_standard = total_standard + neigh_list*2.d0
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device rjs_d         ", neigh_list*2.d0, " Gb"

         total_standard = total_standard + forces/3.d0
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device species_d     ", forces/3.d0, " Gb"

         if (verbose) write (*, '(A)') "> Estimating memory xrd and pdf calculation"
         total_exp = total_exp + nk_int*n_dim_partial
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device k_index_d     ", nk_int, " Gb"

         total_exp = total_exp + nk_int*n_dim_partial
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device j_index_d     ", nk_int, " Gb"

         total_exp = total_exp + 3.d0*Gk
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device      Gk_d     ", Gk, " Gb"
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device     Gka_d     ", Gk, " Gb"

         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device par_pdf_d     ", Gk, " Gb"

         total_exp = total_exp + dermat
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device  dermat_d     ", dermat, " Gb"

         total_exp = total_exp + fi
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device      fi_d     ", fi, " Gb"

         total_exp = total_exp + xyz*n_dim_partial
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device     xyz_d     ", xyz, " Gb"

         total_exp = total_exp + forces
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device   forces_d    ", forces, " Gb"

         total_exp = total_exp + pref
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device    pref_d     ", pref, " Gb"

         total_exp = total_exp + pref
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device  scat_f_d     ", pref, " Gb"

         total_exp = total_exp + sf
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') ">> Gb/core: device  sinc_f_d     ", sf, " Gb"

         total = total + total_exp + total_standard
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') "--- Total from standard neigh:  ", total_standard, " Gb"
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') "--- Total from exp forces imp:  ", total_exp, " Gb"
         if (verbose) write (*, '(A,1X,F10.6,1X,A)') "--- Total device memory usage:  ", total, " Gb"

      end subroutine estimate_max_exp_forces_device_memory_usage
#endif

      end module neighbors
