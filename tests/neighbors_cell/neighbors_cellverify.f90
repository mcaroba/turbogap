! Does the cell list in build_neighbors_list find the same neighbours the
! all-pairs loop it replaced would have found?
!
! The reference here IS that loop: every site against every position, through
! the module's own get_distance, appending in ascending index. So the two
! differ in one thing only -- which pairs are tested -- and any disagreement
! is a binning error and nothing else. get_distance's own image search is not
! under test and cancels out of the comparison.
!
! Seven cells, chosen for the bin counts they produce rather than for physics:
! m = 1, m = 2 and m >= 3 on each axis all appear, as do a hexagonal cell, two
! degrees of shear, a slab and a cell smaller than its own cutoff.
program neighbors_cellverify

   use kinds, only: dp
   use neighbors, only: build_neighbors_list, get_distance

   implicit none

   integer, parameter :: n_cases = 7
   real(dp) :: cell(1:3, 1:3, 1:n_cases)
   real(dp) :: rcut(1:n_cases)
   integer :: n_atoms(1:n_cases)
   integer :: small(1:n_cases)
   character(len=32) :: label(1:n_cases)
   integer :: c
   integer :: failures

!  Large enough for m >= 3, but forced down the non-orthorhombic branch by
!  indices, which is how a cell smaller than the cutoff sphere gets here.
   cell(1:3, 1, 1) = [20.d0, 0.d0, 0.d0]
   cell(1:3, 2, 1) = [0.d0, 20.d0, 0.d0]
   cell(1:3, 3, 1) = [0.d0, 0.d0, 20.d0]
   rcut(1) = 5.d0
   n_atoms(1) = 400
   small(1) = 2
   label(1) = "orthorhombic, forced"

   cell(1:3, 1, 2) = [20.d0, 0.d0, 0.d0]
   cell(1:3, 2, 2) = [-10.d0, 17.320508075688775d0, 0.d0]
   cell(1:3, 3, 2) = [0.d0, 0.d0, 25.d0]
   rcut(2) = 5.d0
   n_atoms(2) = 500
   small(2) = 1
   label(2) = "hexagonal"

   cell(1:3, 1, 3) = [20.d0, 0.d0, 0.d0]
   cell(1:3, 2, 3) = [3.d0, 20.d0, 0.d0]
   cell(1:3, 3, 3) = [2.d0, 1.5d0, 20.d0]
   rcut(3) = 5.d0
   n_atoms(3) = 500
   small(3) = 1
   label(3) = "sheared"

   cell(1:3, 1, 4) = [18.d0, 0.d0, 0.d0]
   cell(1:3, 2, 4) = [7.d0, 16.d0, 0.d0]
   cell(1:3, 3, 4) = [5.d0, 4.d0, 20.d0]
   rcut(4) = 4.5d0
   n_atoms(4) = 400
   small(4) = 1
   label(4) = "strongly sheared"

!  Thin along c, so mz = 2 and the guard against visiting a bin twice runs.
   cell(1:3, 1, 5) = [30.d0, 0.d0, 0.d0]
   cell(1:3, 2, 5) = [-5.d0, 28.d0, 0.d0]
   cell(1:3, 3, 5) = [0.d0, 0.d0, 9.d0]
   rcut(5) = 4.d0
   n_atoms(5) = 600
   small(5) = 1
   label(5) = "slab, mz = 2"

   cell(1:3, 1, 6) = [9.d0, 0.d0, 0.d0]
   cell(1:3, 2, 6) = [-3.d0, 8.5d0, 0.d0]
   cell(1:3, 3, 6) = [0.d0, 0.d0, 10.d0]
   rcut(6) = 4.d0
   n_atoms(6) = 150
   small(6) = 1
   label(6) = "small, m = 2 throughout"

!  Thinner than its own cutoff on every axis, so every m clamps to 1 and the
!  search degenerates to the all-pairs loop it replaced.
   cell(1:3, 1, 7) = [6.d0, 0.d0, 0.d0]
   cell(1:3, 2, 7) = [-2.d0, 5.6d0, 0.d0]
   cell(1:3, 3, 7) = [0.d0, 0.d0, 7.d0]
   rcut(7) = 5.d0
   n_atoms(7) = 80
   small(7) = 1
   label(7) = "cutoff exceeds the cell"

   failures = 0
   do c = 1, n_cases
      call one_case(cell(1:3, 1:3, c), rcut(c), n_atoms(c), small(c), label(c), failures)
   end do

   write (*, *)
   if (failures == 0) then
      write (*, '(A)') "PASS: the cell list and the all-pairs loop agree on every cell"
   else
      write (*, '(A,I0,A)') "FAIL: ", failures, " case(s) disagree"
      stop 1
   end if

contains

!  Positions from a fixed linear congruential generator rather than
!  random_number, whose sequence is not specified and differs between
!  compilers. A test that cannot be reproduced on another machine is not one.
   subroutine fill_positions(cell_in, n, positions)

      real(dp), intent(in) :: cell_in(1:3, 1:3)
      integer, intent(in) :: n
      real(dp), intent(out) :: positions(1:3, 1:n)
      integer(kind=8) :: state
      integer :: i
      integer :: k
      real(dp) :: s(1:3)

      state = 20260916_8
      do i = 1, n
         do k = 1, 3
            state = modulo(6364136223846793005_8*state + 1442695040888963407_8, 9223372036854775807_8)
            s(k) = dble(modulo(state/65536_8, 1000000_8))/1000000.d0
         end do
         positions(1:3, i) = s(1)*cell_in(1:3, 1) + s(2)*cell_in(1:3, 2) + s(3)*cell_in(1:3, 3)
      end do

   end subroutine fill_positions

!  The all-pairs loop, written out: this is what build_neighbors_list did for
!  these cells before it had a cell list.
   subroutine reference_list(positions, a, b, c_vec, rcut_in, n, n_neigh_ref, list_ref, n_pairs_ref)

      real(dp), intent(in) :: positions(:, :)
      real(dp), intent(in) :: a(1:3)
      real(dp), intent(in) :: b(1:3)
      real(dp), intent(in) :: c_vec(1:3)
      real(dp), intent(in) :: rcut_in
      integer, intent(in) :: n
      integer, intent(out) :: n_neigh_ref(1:n)
      integer, intent(out) :: list_ref(:)
      integer, intent(out) :: n_pairs_ref
      integer :: i
      integer :: j
      integer :: nn
      integer :: at
      integer :: i_shift(1:3)
      real(dp) :: dist(1:3)
      real(dp) :: d

      at = 0
      do i = 1, n
         nn = 1
         at = at + 1
         list_ref(at) = i
         do j = 1, size(positions, 2)
            if (j == i) cycle
            call get_distance(positions(1:3, i), positions(1:3, j), a, b, c_vec, &
                              [.true., .true., .true.], dist, d, i_shift)
            if (d < rcut_in) then
               nn = nn + 1
               at = at + 1
               list_ref(at) = j
            end if
         end do
         n_neigh_ref(i) = nn
      end do
      n_pairs_ref = at

   end subroutine reference_list

   subroutine one_case(cell_in, rcut_in, n, small_in, label_in, failures_io)

      real(dp), intent(in) :: cell_in(1:3, 1:3)
      real(dp), intent(in) :: rcut_in
      integer, intent(in) :: n
      integer, intent(in) :: small_in
      character(len=*), intent(in) :: label_in
      integer, intent(inout) :: failures_io
      real(dp), allocatable :: positions(:, :)
      real(dp), allocatable :: rjs(:)
      real(dp), allocatable :: thetas(:)
      real(dp), allocatable :: phis(:)
      real(dp), allocatable :: xyz(:, :)
      real(dp), allocatable :: rjs_ref(:)
      integer, allocatable :: neighbor_species(:)
      integer, allocatable :: neighbors_list(:)
      integer, allocatable :: n_neigh(:)
      integer, allocatable :: n_neigh_ref(:)
      integer, allocatable :: list_ref(:)
      integer, allocatable :: species_supercell(:)
      logical, allocatable :: do_list(:)
      real(dp) :: a(1:3)
      real(dp) :: b(1:3)
      real(dp) :: c_vec(1:3)
      real(dp) :: dist(1:3)
      real(dp) :: d
      integer :: indices(1:3)
      integer :: n_atom_pairs
      integer :: n_pairs_ref
      integer :: i
      integer :: k
      integer :: k2
      integer :: i_shift(1:3)
      integer :: bad_count
      integer :: bad_geom

      allocate (positions(1:3, 1:n))
      call fill_positions(cell_in, n, positions)
      allocate (species_supercell(1:n))
      species_supercell = 1
      allocate (do_list(1:n))
      do_list = .true.
      allocate (n_neigh_ref(1:n))
      allocate (list_ref(1:n*n))
      a = cell_in(1:3, 1)
      b = cell_in(1:3, 2)
      c_vec = cell_in(1:3, 3)
      call reference_list(positions, a, b, c_vec, rcut_in, n, n_neigh_ref, list_ref, n_pairs_ref)

      a = cell_in(1:3, 1)
      b = cell_in(1:3, 2)
      c_vec = cell_in(1:3, 3)
      indices = [small_in, 1, 1]
      call build_neighbors_list(positions, a, b, c_vec, .false., species_supercell, rcut_in, &
                                n_atom_pairs, rjs, thetas, phis, xyz, n_neigh, neighbors_list, &
                                neighbor_species, n, indices, .true., do_list, 0)

      bad_count = 0
      do i = 1, n
         if (n_neigh(i) /= n_neigh_ref(i)) bad_count = bad_count + 1
      end do
      if (bad_count == 0) then
         do k = 1, n_pairs_ref
            if (neighbors_list(k) /= list_ref(k)) bad_count = bad_count + 1
         end do
      end if

!     The distances too, not only the topology: a list in the right order built
!     from the wrong pair would still index correctly.
      allocate (rjs_ref(1:n_atom_pairs))
      k2 = 0
      bad_geom = 0
      do i = 1, n
         do k = 1, n_neigh(i)
            k2 = k2 + 1
            if (k == 1) then
               rjs_ref(k2) = 0.d0
            else
               call get_distance(positions(1:3, i), positions(1:3, neighbors_list(k2)), a, b, c_vec, &
                                 [.true., .true., .true.], dist, d, i_shift)
               rjs_ref(k2) = d
            end if
            if (rjs(k2) /= rjs_ref(k2)) bad_geom = bad_geom + 1
         end do
      end do

      write (*, '(A,A24,A,I6,A,I8,A)') "  ", label_in, ": ", n, " atoms, ", n_atom_pairs, " pairs"
      if (bad_count == 0 .and. bad_geom == 0 .and. n_atom_pairs == n_pairs_ref) then
         write (*, '(A)') "      identical to the all-pairs list, distances included"
      else
         write (*, '(A,I0,A,I0,A,I0,A,I0)') "      MISMATCH: entries ", bad_count, &
            ", distances ", bad_geom, ", pairs ", n_atom_pairs, " against ", n_pairs_ref
         failures_io = failures_io + 1
      end if

      deallocate (positions, species_supercell, do_list, n_neigh_ref, list_ref, rjs_ref)
      deallocate (rjs, thetas, phis, xyz, neighbor_species, neighbors_list, n_neigh)

   end subroutine one_case

end program neighbors_cellverify
