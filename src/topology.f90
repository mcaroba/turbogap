! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, topology.f90, is copyright (c) 2019-2026, Miguel A. Caro
! HND X   and Tigany Zarrouk
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

! Which atoms form a molecule, decided from the bonding rather than from a
! record of what was inserted.
!
! Grand-canonical Monte Carlo labels the molecules it inserts, and can only
! remove those: an atom of a molecule that was in the starting structure carries
! no label, so a solvated system can never come to equilibrium with its
! reservoir -- the run can only take back what it put in. This finds molecules
! by their connectivity instead, so a pre-existing one is as removable as an
! inserted one.
!
! Two atoms are BONDED when their separation is below
!
!     scale * (r_cov(Z_i) + r_cov(Z_j))
!
! with covalent radii from elements.f90 and scale about 1.2. Not van der Waals
! radii: those are roughly twice as large, they describe contact rather than
! bonding, and with them every atom in a condensed phase is bonded to every
! neighbour -- the structure becomes one component and nothing matches.
!
! A molecule is then a connected component of that graph, and it MATCHES the
! reference when it has the same fingerprint: the same number of atoms, the same
! multiset of (element, coordination), and the same multiset of bonded element
! pairs. That is not a graph isomorphism, and two graphs can share it without
! being isomorphic -- but they have to be the same size, the same composition,
! the same degree sequence AND the same bond composition to do so, which no
! molecule anyone runs GCMC on manages by accident.
!
! What it deliberately does not match is a molecule COVALENTLY BOUND to the rest
! of the structure. That lands in a bigger component and is not a free molecule;
! removing it would break bonds the potential is holding.
module topology

   use kinds
   use elements, only: symbol_to_z, covalent_radius

   implicit none

   private
   public :: topology_fingerprint, topology_reference, topology_find, topology_free
   public :: topology_describe

!  The size-ordered invariants of one molecular graph. Sorted, so that two
!  fingerprints compare with a single array equality regardless of the order the
!  atoms happened to be listed in.
   type :: topology_fingerprint
      integer :: n_atoms = 0
      integer :: n_bonds = 0
!     Z*100 + coordination, ascending.
      integer, allocatable :: atom_key(:)
!     min(Z_i,Z_j)*100 + max(Z_i,Z_j), ascending.
      integer, allocatable :: bond_key(:)
      logical :: valid = .false.
   end type topology_fingerprint

contains

!  The fingerprint of an isolated molecule, from its own coordinates.
!
!  No periodic images: a molecule file is one molecule, and if two of its atoms
!  are far enough apart to need one then the file is wrong in a way this would
!  hide. O(n^2), which for the handful of atoms in a molecule is the cheapest
!  thing that can be written.
   subroutine topology_reference(positions, xyz_species, n_atoms, scale, fp, ok, message)

      implicit none

      real(dp), intent(in) :: positions(:, :)
      character*8, intent(in) :: xyz_species(:)
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: scale
      type(topology_fingerprint), intent(out) :: fp
      logical, intent(out) :: ok
      character(len=*), intent(out) :: message

      integer, allocatable :: z(:), degree(:), parent(:)
      integer, allocatable :: bond_i(:), bond_j(:)
      real(dp) :: d, cut, r_i, r_j
      integer :: i, j, n_bonds

      ok = .false.
      message = ""

      if (n_atoms < 1) then
         message = "the molecule has no atoms"
         return
      end if

      allocate (z(1:n_atoms), degree(1:n_atoms), parent(1:n_atoms))
      allocate (bond_i(1:n_atoms*(n_atoms - 1)/2 + 1))
      allocate (bond_j(1:n_atoms*(n_atoms - 1)/2 + 1))

      do i = 1, n_atoms
         z(i) = symbol_to_z(xyz_species(i))
         if (covalent_radius(z(i)) < 0.d0) then
            message = "no covalent radius for element "//trim(xyz_species(i))
            return
         end if
      end do

      degree = 0
      n_bonds = 0
      do i = 1, n_atoms
         parent(i) = i
      end do
      do i = 1, n_atoms - 1
         r_i = covalent_radius(z(i))
         do j = i + 1, n_atoms
            r_j = covalent_radius(z(j))
            cut = scale*(r_i + r_j)
            d = dsqrt(sum((positions(1:3, i) - positions(1:3, j))**2))
            if (d < cut) then
               n_bonds = n_bonds + 1
               bond_i(n_bonds) = i
               bond_j(n_bonds) = j
               degree(i) = degree(i) + 1
               degree(j) = degree(j) + 1
               call union(parent, i, j)
            end if
         end do
      end do

!     A reference that falls into two pieces would match any two pieces
!     anywhere, and the removal would take atoms that are not a molecule.
      do i = 2, n_atoms
         if (find(parent, i) /= find(parent, 1)) then
            message = "the molecule is not connected at this bond scale"
            return
         end if
      end do

      call build_fingerprint(z, degree, bond_i, bond_j, n_atoms, n_bonds, fp)
      ok = .true.

   end subroutine topology_reference

!  Label every molecule in the structure that matches the reference.
!
!  mol_id(i) is set to a distinct positive number per matched molecule, counting
!  up from first_id, and left untouched for atoms that are not in one, so this
!  composes with whatever labelling the caller already has. n_found returns how
!  many were matched.
!
!  O(N) in the atom count: the candidate pairs come from a cell list on the
!  fractional coordinates, which is why this can run inside a Monte Carlo loop
!  rather than only at setup.
   subroutine topology_find(positions, xyz_species, n_sites, a_box, b_box, c_box, &
                            scale, fp, first_id, mol_id, n_found)

      implicit none

      real(dp), intent(in) :: positions(:, :)
      character*8, intent(in) :: xyz_species(:)
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
      real(dp), intent(in) :: scale
      type(topology_fingerprint), intent(in) :: fp
      integer, intent(in) :: first_id
      integer, intent(inout) :: mol_id(:)
      integer, intent(out) :: n_found

      integer, allocatable :: z(:), degree(:), parent(:), head(:), next(:)
      integer, allocatable :: bond_i(:), bond_j(:), cell_of(:, :)
      integer, allocatable :: member(:), comp_size(:), comp_root(:)
      real(dp), allocatable :: frac(:, :)
      real(dp) :: r_max, cut, d, lat(1:3, 1:3)
      integer :: n_cells(1:3), i, j, k, n_bonds, root, n_comp, id
      integer :: ia, ib, ic, ja, jb, jc, da, db, dc, cell, other

      n_found = 0
      if (.not. fp%valid .or. n_sites < fp%n_atoms) return

      allocate (z(1:n_sites), degree(1:n_sites), parent(1:n_sites))
      do i = 1, n_sites
         z(i) = symbol_to_z(xyz_species(i))
         parent(i) = i
      end do
      degree = 0

!     The widest bond any pair present can make. Sized from the elements that
!     are actually here, not from the table's maximum, so a cell list for an
!     organic system is built at 1.8 A rather than at 6 A.
      r_max = 0.d0
      do i = 1, n_sites
         if (covalent_radius(z(i)) > r_max) r_max = covalent_radius(z(i))
      end do
      if (r_max <= 0.d0) return
      cut = 2.d0*scale*r_max

      lat(1:3, 1) = a_box
      lat(1:3, 2) = b_box
      lat(1:3, 3) = c_box
      call bin_atoms(positions, n_sites, lat, cut, frac, n_cells, head, next, cell_of)

      allocate (bond_i(1:16*n_sites), bond_j(1:16*n_sites))
      n_bonds = 0

!     Each atom against its own cell and the 26 around it. Every bond is found
!     twice that way, so only i < j is kept.
      do i = 1, n_sites
         ia = cell_of(1, i); ib = cell_of(2, i); ic = cell_of(3, i)
         do da = -1, 1
            do db = -1, 1
               do dc = -1, 1
                  ja = modulo(ia + da, n_cells(1))
                  jb = modulo(ib + db, n_cells(2))
                  jc = modulo(ic + dc, n_cells(3))
                  cell = 1 + ja + n_cells(1)*(jb + n_cells(2)*jc)
                  other = head(cell)
                  do while (other > 0)
                     if (other > i) then
                        d = minimum_image_distance(positions(1:3, i), positions(1:3, other), lat)
                        if (d < scale*(covalent_radius(z(i)) + covalent_radius(z(other)))) then
                           n_bonds = n_bonds + 1
                           if (n_bonds > size(bond_i)) call grow(bond_i, bond_j)
                           bond_i(n_bonds) = i
                           bond_j(n_bonds) = other
                           degree(i) = degree(i) + 1
                           degree(other) = degree(other) + 1
                           call union(parent, i, other)
                        end if
                     end if
                     other = next(other)
                  end do
               end do
            end do
         end do
      end do

!     Components, by root. Only those of exactly the reference's size can match,
!     and testing that first is one integer comparison against building a
!     fingerprint for every component in the cell.
      allocate (comp_size(1:n_sites), comp_root(1:n_sites), member(1:fp%n_atoms))
      comp_size = 0
      do i = 1, n_sites
         root = find(parent, i)
         comp_size(root) = comp_size(root) + 1
      end do

      n_comp = 0
      do i = 1, n_sites
         if (comp_size(i) == fp%n_atoms) then
            n_comp = n_comp + 1
            comp_root(n_comp) = i
         end if
      end do

      id = first_id
      do k = 1, n_comp
         j = 0
         do i = 1, n_sites
            if (find(parent, i) == comp_root(k)) then
               j = j + 1
               member(j) = i
            end if
         end do
         if (component_matches(z, degree, bond_i, bond_j, n_bonds, member, fp)) then
            n_found = n_found + 1
            do i = 1, fp%n_atoms
               mol_id(member(i)) = id
            end do
            id = id + 1
         end if
      end do

   end subroutine topology_find

!  One line saying what the reference is, for the run's echo of its own setup.
   subroutine topology_describe(fp, text)

      implicit none

      type(topology_fingerprint), intent(in) :: fp
      character(len=*), intent(out) :: text

      write (text, '(A,I0,A,I0,A)') "topology: ", fp%n_atoms, " atoms, ", fp%n_bonds, " bonds"

   end subroutine topology_describe

   subroutine topology_free(fp)

      implicit none

      type(topology_fingerprint), intent(inout) :: fp

      if (allocated(fp%atom_key)) deallocate (fp%atom_key)
      if (allocated(fp%bond_key)) deallocate (fp%bond_key)
      fp%valid = .false.
      fp%n_atoms = 0
      fp%n_bonds = 0

   end subroutine topology_free

! ---------------------------------------------------------------- internals

   function component_matches(z, degree, bond_i, bond_j, n_bonds, member, fp) result(same)

      implicit none

      integer, intent(in) :: z(:), degree(:), bond_i(:), bond_j(:), n_bonds, member(:)
      type(topology_fingerprint), intent(in) :: fp
      logical :: same
      type(topology_fingerprint) :: candidate
      integer, allocatable :: local_z(:), local_deg(:), local_bi(:), local_bj(:), where_in(:)
      integer :: n, i, k, count_bonds

      same = .false.
      n = size(member)
      allocate (local_z(1:n), local_deg(1:n), where_in(1:maxval(member)))
      where_in = 0
      do i = 1, n
         where_in(member(i)) = i
         local_z(i) = z(member(i))
         local_deg(i) = degree(member(i))
      end do

      allocate (local_bi(1:n_bonds + 1), local_bj(1:n_bonds + 1))
      count_bonds = 0
      do k = 1, n_bonds
         if (bond_i(k) <= size(where_in) .and. bond_j(k) <= size(where_in)) then
            if (where_in(bond_i(k)) > 0 .and. where_in(bond_j(k)) > 0) then
               count_bonds = count_bonds + 1
               local_bi(count_bonds) = where_in(bond_i(k))
               local_bj(count_bonds) = where_in(bond_j(k))
            end if
         end if
      end do

      call build_fingerprint(local_z, local_deg, local_bi, local_bj, n, count_bonds, candidate)
      if (candidate%n_atoms == fp%n_atoms .and. candidate%n_bonds == fp%n_bonds) then
         if (all(candidate%atom_key == fp%atom_key) .and. all(candidate%bond_key == fp%bond_key)) same = .true.
      end if
      call topology_free(candidate)

   end function component_matches

   subroutine build_fingerprint(z, degree, bond_i, bond_j, n_atoms, n_bonds, fp)

      implicit none

      integer, intent(in) :: z(:), degree(:), bond_i(:), bond_j(:), n_atoms, n_bonds
      type(topology_fingerprint), intent(out) :: fp
      integer :: i

      fp%n_atoms = n_atoms
      fp%n_bonds = n_bonds
      allocate (fp%atom_key(1:max(n_atoms, 1)))
      allocate (fp%bond_key(1:max(n_bonds, 1)))
      fp%atom_key = 0
      fp%bond_key = 0

      do i = 1, n_atoms
         fp%atom_key(i) = 100*z(i) + degree(i)
      end do
      do i = 1, n_bonds
         fp%bond_key(i) = 100*min(z(bond_i(i)), z(bond_j(i))) + max(z(bond_i(i)), z(bond_j(i)))
      end do

      call sort_int(fp%atom_key, n_atoms)
      call sort_int(fp%bond_key, n_bonds)
      fp%valid = .true.

   end subroutine build_fingerprint

!  Insertion sort. n is a molecule's atom count, so it is small and this is
!  faster than anything with a better exponent.
   subroutine sort_int(a, n)

      implicit none

      integer, intent(inout) :: a(:)
      integer, intent(in) :: n
      integer :: i, j, key

      do i = 2, n
         key = a(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= key) exit
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = key
      end do

   end subroutine sort_int

   recursive function find(parent, i) result(root)

      implicit none

      integer, intent(inout) :: parent(:)
      integer, intent(in) :: i
      integer :: root

      if (parent(i) /= i) then
         parent(i) = find(parent, parent(i))
      end if
      root = parent(i)

   end function find

   subroutine union(parent, i, j)

      implicit none

      integer, intent(inout) :: parent(:)
      integer, intent(in) :: i, j
      integer :: ri, rj

      ri = find(parent, i)
      rj = find(parent, j)
      if (ri /= rj) parent(rj) = ri

   end subroutine union

   subroutine grow(bond_i, bond_j)

      implicit none

      integer, allocatable, intent(inout) :: bond_i(:), bond_j(:)
      integer, allocatable :: tmp(:)
      integer :: n

      n = size(bond_i)
      allocate (tmp(1:2*n))
      tmp(1:n) = bond_i
      call move_alloc(tmp, bond_i)
      allocate (tmp(1:2*n))
      tmp(1:n) = bond_j
      call move_alloc(tmp, bond_j)

   end subroutine grow

!  Minimum image, by brute force over the 27 neighbouring images.
!
!  Correct for any cell including a strongly triclinic one, where the usual
!  round-the-fractional-coordinate shortcut picks the wrong image. 27 distance
!  evaluations per candidate pair, and the cell list has already reduced the
!  candidates to a handful per atom.
   function minimum_image_distance(r1, r2, lat) result(d)

      implicit none

      real(dp), intent(in) :: r1(1:3), r2(1:3), lat(1:3, 1:3)
      real(dp) :: d, d2, best, v(1:3), shift(1:3)
      integer :: i, j, k

      best = huge(1.d0)
      do i = -1, 1
         do j = -1, 1
            do k = -1, 1
               shift = dfloat(i)*lat(1:3, 1) + dfloat(j)*lat(1:3, 2) + dfloat(k)*lat(1:3, 3)
               v = r1 - r2 + shift
               d2 = sum(v*v)
               if (d2 < best) best = d2
            end do
         end do
      end do
      d = dsqrt(best)

   end function minimum_image_distance

!  Sort the atoms into a linked-cell grid on the fractional coordinates, so the
!  binning is right for a triclinic cell without any of them being wrapped.
   subroutine bin_atoms(positions, n_sites, lat, cut, frac, n_cells, head, next, cell_of)

      implicit none

      real(dp), intent(in) :: positions(:, :)
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: lat(1:3, 1:3), cut
      real(dp), allocatable, intent(out) :: frac(:, :)
      integer, intent(out) :: n_cells(1:3)
      integer, allocatable, intent(out) :: head(:), next(:), cell_of(:, :)
      real(dp) :: inv(1:3, 1:3), width(1:3), cross(1:3), volume
      integer :: i, c, total, axis

      call invert_3x3_local(lat, inv, volume)

!     The perpendicular width of the cell along each axis: the spacing of the
!     lattice planes, which is what a cell has to be at least `cut` across.
      call cross_product_local(lat(1:3, 2), lat(1:3, 3), cross)
      width(1) = dabs(volume)/dsqrt(sum(cross*cross))
      call cross_product_local(lat(1:3, 3), lat(1:3, 1), cross)
      width(2) = dabs(volume)/dsqrt(sum(cross*cross))
      call cross_product_local(lat(1:3, 1), lat(1:3, 2), cross)
      width(3) = dabs(volume)/dsqrt(sum(cross*cross))

      do axis = 1, 3
         n_cells(axis) = max(1, int(width(axis)/cut))
!        Fewer than three cells on an axis means the 27-image sweep already
!        covers it, and a grid that small only costs bookkeeping.
         if (n_cells(axis) < 3) n_cells(axis) = 1
      end do

      total = n_cells(1)*n_cells(2)*n_cells(3)
      allocate (frac(1:3, 1:n_sites), head(1:total), next(1:n_sites), cell_of(1:3, 1:n_sites))
      head = 0
      next = 0

      do i = 1, n_sites
         frac(1:3, i) = matmul(inv, positions(1:3, i))
         frac(1:3, i) = frac(1:3, i) - dfloat(floor(frac(1:3, i)))
         do axis = 1, 3
            cell_of(axis, i) = min(n_cells(axis) - 1, int(frac(axis, i)*dfloat(n_cells(axis))))
            if (cell_of(axis, i) < 0) cell_of(axis, i) = 0
         end do
         c = 1 + cell_of(1, i) + n_cells(1)*(cell_of(2, i) + n_cells(2)*cell_of(3, i))
         next(i) = head(c)
         head(c) = i
      end do

   end subroutine bin_atoms

   subroutine cross_product_local(u, v, w)

      implicit none

      real(dp), intent(in) :: u(1:3), v(1:3)
      real(dp), intent(out) :: w(1:3)

      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)

   end subroutine cross_product_local

   subroutine invert_3x3_local(m, inv, det)

      implicit none

      real(dp), intent(in) :: m(1:3, 1:3)
      real(dp), intent(out) :: inv(1:3, 1:3), det

      det = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) &
            - m(1, 2)*(m(2, 1)*m(3, 3) - m(2, 3)*m(3, 1)) &
            + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))

      inv(1, 1) = (m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2))/det
      inv(1, 2) = (m(1, 3)*m(3, 2) - m(1, 2)*m(3, 3))/det
      inv(1, 3) = (m(1, 2)*m(2, 3) - m(1, 3)*m(2, 2))/det
      inv(2, 1) = (m(2, 3)*m(3, 1) - m(2, 1)*m(3, 3))/det
      inv(2, 2) = (m(1, 1)*m(3, 3) - m(1, 3)*m(3, 1))/det
      inv(2, 3) = (m(1, 3)*m(2, 1) - m(1, 1)*m(2, 3))/det
      inv(3, 1) = (m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))/det
      inv(3, 2) = (m(1, 2)*m(3, 1) - m(1, 1)*m(3, 2))/det
      inv(3, 3) = (m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1))/det

   end subroutine invert_3x3_local

end module topology
