! Standalone checks on src/topology.f90.
!
! Every assertion compares against a structure whose answer is known because
! this program built it: put four ethanols in a box and four is the answer. No
! potential and no GAP file is involved -- topology is geometry and graph
! theory, so it can be checked without any of the machinery that normally
! surrounds it, and a failure here names the topology code and nothing else.
!
! The system is the one the feature is for: ethanol intercalated into an
! oxygenated and hydrogenated carbon sheet. That is the case where a naive
! match goes wrong, because the sheet's own -OH and -H groups are made of the
! same elements as the guest.
program topologyverify

   use kinds
   use topology

   implicit none

   real(dp), parameter :: SCALE = 1.2d0
   integer :: failures

   failures = 0

   call check_reference_fingerprint(failures)
   call check_isolated_molecules(failures)
   call check_intercalated_sheet(failures)
   call check_bonded_guest_is_not_matched(failures)
   call check_across_periodic_boundary(failures)
   call check_bond_cutoff_boundary(failures)

   write (*, *)
   if (failures == 0) then
      write (*, '(A)') "topology: all checks passed"
   else
      write (*, '(A,I0,A)') "topology: ", failures, " check(s) FAILED"
      stop 1
   end if

contains

!  Ethanol, CH3-CH2-OH, at roughly its real geometry: C-C 1.53, C-O 1.43,
!  O-H 0.97, C-H 1.09 Angstrom.
   subroutine ethanol(positions, species)

      implicit none

      real(dp), allocatable, intent(out) :: positions(:, :)
      character*8, allocatable, intent(out) :: species(:)

      allocate (positions(1:3, 1:9), species(1:9))

      species = [character*8 :: "C", "C", "O", "H", "H", "H", "H", "H", "H"]

      positions(1:3, 1) = [0.000d0, 0.000d0, 0.000d0]     ! methyl C
      positions(1:3, 2) = [1.530d0, 0.000d0, 0.000d0]     ! methylene C
      positions(1:3, 3) = [2.245d0, 1.238d0, 0.000d0]     ! O
      positions(1:3, 4) = [-0.360d0, 1.030d0, 0.000d0]    ! 3 H on the methyl C
      positions(1:3, 5) = [-0.360d0, -0.510d0, 0.890d0]
      positions(1:3, 6) = [-0.360d0, -0.510d0, -0.890d0]
      positions(1:3, 7) = [1.890d0, -0.510d0, 0.890d0]    ! 2 H on the methylene C
      positions(1:3, 8) = [1.890d0, -0.510d0, -0.890d0]
      positions(1:3, 9) = [3.145d0, 1.598d0, 0.000d0]     ! hydroxyl H

   end subroutine ethanol

   subroutine reference(fp)

      implicit none

      type(topology_fingerprint), intent(out) :: fp
      real(dp), allocatable :: positions(:, :)
      character*8, allocatable :: species(:)
      character(len=256) :: message
      logical :: ok

      call ethanol(positions, species)
      call topology_reference(positions, species, 9, SCALE, fp, ok, message)
      if (.not. ok) then
         write (*, '(A)') "  could not build the ethanol reference: "//trim(message)
         stop 1
      end if

   end subroutine reference

   subroutine report(name, got, want, failures)

      implicit none

      character(len=*), intent(in) :: name
      integer, intent(in) :: got, want
      integer, intent(inout) :: failures

      if (got == want) then
         write (*, '(A,A,A,I0)') "  ok    ", name, ": ", got
      else
         write (*, '(A,A,A,I0,A,I0)') "  FAIL  ", name, ": got ", got, ", expected ", want
         failures = failures + 1
      end if

   end subroutine report

!  The fingerprint of ethanol, worked out by hand: 9 atoms and 8 bonds, both
!  carbons four-coordinated, the oxygen two, every hydrogen one.
   subroutine check_reference_fingerprint(failures)

      implicit none

      integer, intent(inout) :: failures
      type(topology_fingerprint) :: fp
      integer :: n_ch, i

      write (*, '(A)') "ethanol reference"
      call reference(fp)
      call report("atoms", fp%n_atoms, 9, failures)
      call report("bonds", fp%n_bonds, 8, failures)

!     Six hydrogens with one bond each: key 1*100 + 1.
      call report("one-coordinate H", count(fp%atom_key == 101), 6, failures)
!     Two carbons with four bonds each: 6*100 + 4.
      call report("four-coordinate C", count(fp%atom_key == 604), 2, failures)
!     One oxygen with two bonds: 8*100 + 2.
      call report("two-coordinate O", count(fp%atom_key == 802), 1, failures)

!     Five C-H, one C-C, one C-O, one O-H.
      n_ch = 0
      do i = 1, fp%n_bonds
         if (fp%bond_key(i) == 106) n_ch = n_ch + 1
      end do
      call report("C-H bonds", n_ch, 5, failures)
      call report("C-C bonds", count(fp%bond_key == 606), 1, failures)
      call report("C-O bonds", count(fp%bond_key == 608), 1, failures)
      call report("O-H bonds", count(fp%bond_key == 108), 1, failures)

      call topology_free(fp)

   end subroutine check_reference_fingerprint

!  Four ethanols, far enough apart not to bond to each other. Four is the
!  answer because this put four there.
   subroutine check_isolated_molecules(failures)

      implicit none

      integer, intent(inout) :: failures
      type(topology_fingerprint) :: fp
      real(dp), allocatable :: mol(:, :), positions(:, :)
      character*8, allocatable :: mol_species(:), species(:)
      integer, allocatable :: mol_id(:)
      real(dp) :: a(1:3), b(1:3), c(1:3)
      integer :: i, k, n_found

      write (*, '(A)') "four isolated ethanols"
      call reference(fp)
      call ethanol(mol, mol_species)

      allocate (positions(1:3, 1:36), species(1:36), mol_id(1:36))
      mol_id = 0
      do k = 0, 3
         do i = 1, 9
            positions(1:3, k*9 + i) = mol(1:3, i) + [dfloat(k)*8.d0, 0.d0, 0.d0]
            species(k*9 + i) = mol_species(i)
         end do
      end do

      a = [40.d0, 0.d0, 0.d0]
      b = [0.d0, 20.d0, 0.d0]
      c = [0.d0, 0.d0, 20.d0]

      call topology_find(positions, species, 36, a, b, c, SCALE, fp, 1, mol_id, n_found)
      call report("molecules found", n_found, 4, failures)
      call report("atoms labelled", count(mol_id > 0), 36, failures)
      call report("distinct labels", maxval(mol_id), 4, failures)

      call topology_free(fp)

   end subroutine check_isolated_molecules

!  The case the feature is for. A carbon sheet decorated with -OH and -H, and
!  two ethanols sitting between the sheets without bonding to them.
!
!  The sheet carries the same elements as the guest, in -OH and -H groups whose
!  local bonding looks like parts of ethanol, so anything matching on
!  composition alone would find molecules that are not there. The sheet is one
!  big component, so nothing in it has ethanol's 9 atoms.
   subroutine check_intercalated_sheet(failures)

      implicit none

      integer, intent(inout) :: failures
      type(topology_fingerprint) :: fp
      real(dp), allocatable :: mol(:, :), positions(:, :)
      character*8, allocatable :: mol_species(:), species(:)
      integer, allocatable :: mol_id(:)
      real(dp) :: a(1:3), b(1:3), c(1:3)
      integer :: n, i, j, k, n_found, n_sheet

      write (*, '(A)') "ethanol intercalated in an oxygenated, hydrogenated carbon sheet"
      call reference(fp)
      call ethanol(mol, mol_species)

!     A 6x6 square lattice of carbons at 1.42 A, with an -OH on one and an -H
!     on another. 36 C + 2 O + 2 H(hydroxyl) + 1 H(on C) = 41 atoms.
      n_sheet = 36 + 2 + 2 + 1
      n = n_sheet + 18
      allocate (positions(1:3, 1:n), species(1:n), mol_id(1:n))
      mol_id = 0

      k = 0
      do i = 0, 5
         do j = 0, 5
            k = k + 1
            positions(1:3, k) = [dfloat(i)*1.42d0, dfloat(j)*1.42d0, 0.d0]
            species(k) = "C"
         end do
      end do
!     Two hydroxyls, above carbons 1 and 8.
      k = k + 1; positions(1:3, k) = positions(1:3, 1) + [0.d0, 0.d0, 1.43d0]; species(k) = "O"
      k = k + 1; positions(1:3, k) = positions(1:3, k - 1) + [0.d0, 0.90d0, 0.36d0]; species(k) = "H"
      k = k + 1; positions(1:3, k) = positions(1:3, 8) + [0.d0, 0.d0, 1.43d0]; species(k) = "O"
      k = k + 1; positions(1:3, k) = positions(1:3, k - 1) + [0.d0, 0.90d0, 0.36d0]; species(k) = "H"
!     One bare hydrogen on a carbon.
      k = k + 1; positions(1:3, k) = positions(1:3, 20) + [0.d0, 0.d0, 1.09d0]; species(k) = "H"

!     Two ethanols well above the sheet, and apart from each other.
      do j = 0, 1
         do i = 1, 9
            k = k + 1
            positions(1:3, k) = mol(1:3, i) + [dfloat(j)*9.d0, 1.d0, 6.d0]
            species(k) = mol_species(i)
         end do
      end do

      a = [30.d0, 0.d0, 0.d0]
      b = [0.d0, 30.d0, 0.d0]
      c = [0.d0, 0.d0, 30.d0]

      call topology_find(positions, species, n, a, b, c, SCALE, fp, 1, mol_id, n_found)
      call report("ethanols found", n_found, 2, failures)
!     Nothing in the sheet is labelled: every atom of it is in one component of
!     41 atoms, which is not 9.
      call report("sheet atoms labelled", count(mol_id(1:n_sheet) > 0), 0, failures)
      call report("guest atoms labelled", count(mol_id(n_sheet + 1:n) > 0), 18, failures)

      call topology_free(fp)

   end subroutine check_intercalated_sheet

!  An ethanol chemisorbed on the sheet is not a free molecule, and must not be
!  removable. It ends up in the sheet's component, so its atom count is wrong
!  and it cannot match.
   subroutine check_bonded_guest_is_not_matched(failures)

      implicit none

      integer, intent(inout) :: failures
      type(topology_fingerprint) :: fp
      real(dp), allocatable :: mol(:, :), positions(:, :)
      character*8, allocatable :: mol_species(:), species(:)
      integer, allocatable :: mol_id(:)
      real(dp) :: a(1:3), b(1:3), c(1:3)
      integer :: i, n_found

      write (*, '(A)') "an ethanol bonded to the surface"
      call reference(fp)
      call ethanol(mol, mol_species)

!     One carbon, with an ethanol whose methyl carbon is a bond length away
!     from it: 1 + 9 = 10 atoms in one component.
      allocate (positions(1:3, 1:10), species(1:10), mol_id(1:10))
      mol_id = 0
      positions(1:3, 1) = [0.d0, 0.d0, -1.53d0]
      species(1) = "C"
      do i = 1, 9
         positions(1:3, i + 1) = mol(1:3, i)
         species(i + 1) = mol_species(i)
      end do

      a = [20.d0, 0.d0, 0.d0]
      b = [0.d0, 20.d0, 0.d0]
      c = [0.d0, 0.d0, 20.d0]

      call topology_find(positions, species, 10, a, b, c, SCALE, fp, 1, mol_id, n_found)
      call report("molecules found", n_found, 0, failures)

      call topology_free(fp)

   end subroutine check_bonded_guest_is_not_matched

!  A molecule straddling the cell boundary is still one molecule. Without the
!  minimum image it would be found as two fragments and matched as neither.
   subroutine check_across_periodic_boundary(failures)

      implicit none

      integer, intent(inout) :: failures
      type(topology_fingerprint) :: fp
      real(dp), allocatable :: mol(:, :), positions(:, :)
      character*8, allocatable :: mol_species(:), species(:)
      integer, allocatable :: mol_id(:)
      real(dp) :: a(1:3), b(1:3), c(1:3)
      integer :: i, n_found

      write (*, '(A)') "a molecule across the periodic boundary"
      call reference(fp)
      call ethanol(mol, mol_species)

      allocate (positions(1:3, 1:9), species(1:9), mol_id(1:9))
      mol_id = 0
!     Shifted so the molecule sits astride x = 0 of a 20 A cell.
      do i = 1, 9
         positions(1:3, i) = mol(1:3, i) + [19.d0, 5.d0, 5.d0]
         species(i) = mol_species(i)
      end do

      a = [20.d0, 0.d0, 0.d0]
      b = [0.d0, 20.d0, 0.d0]
      c = [0.d0, 0.d0, 20.d0]

      call topology_find(positions, species, 9, a, b, c, SCALE, fp, 1, mol_id, n_found)
      call report("molecules found", n_found, 1, failures)

      call topology_free(fp)

   end subroutine check_across_periodic_boundary

!  The criterion is a threshold, so it is worth knowing that it is where it is
!  claimed to be. Two carbons at 0.99 and at 1.01 of scale*(r_i + r_j).
   subroutine check_bond_cutoff_boundary(failures)

      implicit none

      integer, intent(inout) :: failures
      type(topology_fingerprint) :: fp
      real(dp) :: positions(1:3, 1:2), cut
      character*8 :: species(1:2)
      character(len=256) :: message
      logical :: ok

      write (*, '(A)') "the bond cutoff is where it says it is"

!     Covalent radius of carbon is 0.76 A, so the cutoff is 1.2*1.52 = 1.824.
      cut = SCALE*(0.76d0 + 0.76d0)
      species = [character*8 :: "C", "C"]
      positions(1:3, 1) = [0.d0, 0.d0, 0.d0]

!     Just inside: one connected pair, so the reference builds.
      positions(1:3, 2) = [0.99d0*cut, 0.d0, 0.d0]
      call topology_reference(positions, species, 2, SCALE, fp, ok, message)
      if (ok) then
         write (*, '(A)') "  ok    just inside the cutoff: bonded"
      else
         write (*, '(A)') "  FAIL  just inside the cutoff: not bonded"
         failures = failures + 1
      end if
      call topology_free(fp)

!     Just outside: two fragments, which topology_reference refuses.
      positions(1:3, 2) = [1.01d0*cut, 0.d0, 0.d0]
      call topology_reference(positions, species, 2, SCALE, fp, ok, message)
      if (.not. ok) then
         write (*, '(A)') "  ok    just outside the cutoff: not bonded"
      else
         write (*, '(A)') "  FAIL  just outside the cutoff: bonded"
         failures = failures + 1
      end if
      call topology_free(fp)

   end subroutine check_bond_cutoff_boundary

end program topologyverify
