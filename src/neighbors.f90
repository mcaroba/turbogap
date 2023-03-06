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


  use soap_turbo_functions


  contains


!**************************************************************************
!
! This subroutine returns the distance between ri and rj under certain
! boundary conditions.
!
  subroutine get_distance(posi, posj, a, b, c, PBC, dist, d, i_shift)

    implicit none

    real*8, intent(in) :: posi(1:3), posj(1:3), a(1:3), b(1:3), c(1:3)
    integer, intent(out) :: i_shift(1:3)
    logical, intent(in) :: PBC(1:3)
    real*8, intent(out) :: d
    real*8, intent(out) :: dist(1:3)
    real*8 :: d2, L(1:3), d_tol = 1.d-6
    real*8 :: mat(1:3,1:3), md, indices_real(1:3), res, res_opt, dist_opt(1:3), dist_temp(1:3)
    real*8, save :: a0(1:3) = 0.d0, b0(1:3) = 0.d0, c0(1:3) = 0.d0, mat_inv(1:3,1:3) = 0.d0
    integer :: i, j, k, indices(1:3)
    logical :: lattice_check_a(1:3), lattice_check_b(1:3), lattice_check_c(1:3)

    if( dabs(a(2)) < d_tol .and. dabs(a(3)) < d_tol .and. &
        dabs(b(1)) < d_tol .and. dabs(b(3)) < d_tol .and. &
        dabs(c(1)) < d_tol .and. dabs(c(2)) < d_tol )then
      dist = posj - posi
      i_shift = floor([dist(1)/a(1), dist(2)/b(2), dist(3)/c(3)])
!     Fast solution for orthorhombic cells
      L = (/ a(1), b(2), c(3) /)
      d2 = 0.d0
      do i = 1, 3
        if( PBC(i) )then
          dist(i) = modulo(posj(i) - posi(i), L(i))
          if( dist(i) > L(i)/2.d0 )then
            dist(i) = dist(i) - L(i)
          end if
        else
          dist(i) = posj(i) - posi(i)
        end if
        d2 = d2 + dist(i)**2
      end do
      d = dsqrt(d2)
    else if( all( PBC ) )then
!     Slow solution for other unit cells
      lattice_check_a = ( a /= a0 )
      lattice_check_b = ( b /= b0 )
      lattice_check_c = ( c /= c0 )
      if( any(lattice_check_a) .or. any(lattice_check_b) .or. any(lattice_check_c) )then
        a0 = a
        b0 = b
        c0 = c
!       We construct our matrix to get the MIC only if the lattice vectors have changed
        mat(1,1) = dot_product(a, a)
        mat(1,2) = dot_product(a, b)
        mat(1,3) = dot_product(a, c)
        mat(2,1) = mat(1,2)
        mat(2,2) = dot_product(b, b)
        mat(2,3) = dot_product(b, c)
        mat(3,1) = mat(1,3)
        mat(3,2) = mat(2,3)
        mat(3,3) = dot_product(c, c)
!       We compute the inverse of this matrix analytically
        md = -mat(1,3)**2*mat(2,2) + 2.d0*mat(1,2)*mat(1,3)*mat(2,3) - mat(1,1)*mat(2,3)**2 &
             - mat(1,2)**2*mat(3,3) + mat(1,1)*mat(2,2)*mat(3,3)
        mat_inv(1,1) = mat(2,2)*mat(3,3) - mat(2,3)**2
        mat_inv(1,2) = mat(1,3)*mat(2,3) - mat(1,2)*mat(3,3)
        mat_inv(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
        mat_inv(2,1) = mat_inv(1,2)
        mat_inv(2,2) = mat(1,1)*mat(3,3) - mat(1,3)**2
        mat_inv(2,3) = mat(1,2)*mat(1,3) - mat(1,1)*mat(2,3)
        mat_inv(3,1) = mat_inv(1,3)
        mat_inv(3,2) = mat_inv(2,3)
        mat_inv(3,3) = mat(1,1)*mat(2,2) - mat(1,2)**2
        mat_inv = mat_inv / md
      end if
      dist = posj - posi
      indices_real = -1.d0 * (/ dot_product(dist, a), dot_product(dist, b), dot_product(dist,c) /)
      indices_real = matmul(mat_inv, indices_real)
      i_shift = floor(-indices_real)
!     Closest integer solution
      indices(1:3) = nint(indices_real(1:3))
!     We bruteforce the integer solution among the 27 points surrounding the real solution
      res_opt = 1.d10
      do i = indices(1)-1, indices(1)+1
        do j = indices(2)-1, indices(2)+1
          do k = indices(3)-1, indices(3)+1
            dist_temp(1:3) = dist(1:3) + dfloat(i)*a(1:3) + dfloat(j)*b(1:3) + dfloat(k)*c(1:3)
            res = dot_product(dist_temp, dist_temp)
            if( res < res_opt )then
              res_opt = res
              dist_opt = dist_temp
            end if
          end do
        end do
      end do
      dist = dist_opt
      d = dsqrt( dot_product(dist,dist) )
    else
      write(*,*) "Sorry, non-orthorhombic unit cells only work in combination with full PBC"
      stop
    end if

  return
  end subroutine get_distance
!**************************************************************************







!**************************************************************************
!
! This subroutine returns the number of primitive unit cells required to
! construct a supercell whose unit cell's planes are at least 2*rcut apart.
! This construct allows us to search for all the neighbors within the given
! cutoff. indices(1:3) tell the user how many repetitions are required to
! construct a unit cell with the properties outlined above. (1,1,1) means
! that the primitive unit cell is already enough.
!
  subroutine number_of_unit_cells_for_given_cutoff(a, b, c, rcut, PBC, indices)

    implicit none

    real*8, intent(in) :: a(1:3), b(1:3), c(1:3), rcut
    logical, intent(in) :: PBC(1:3)
    integer, intent(out) :: indices(1:3)
    real*8 :: axb(1:3), bxc(1:3), cxa(1:3), indices_real(1:3)
    real*8 :: mat(1:3,1:3), mat_inv(1:3,1:3), md
    integer :: i

    axb = cross_product(a, b)
    axb = axb / dsqrt( dot_product(axb, axb) )
    cxa = cross_product(c, a)
    cxa = cxa / dsqrt( dot_product(cxa, cxa) )
    bxc = cross_product(b, c)
    bxc = bxc / dsqrt( dot_product(bxc, bxc) )

!   Matrix elements
    mat(1,1) = dot_product(axb, a)
    mat(1,2) = dot_product(axb, b)
    mat(1,3) = dot_product(axb, c)
    mat(2,1) = dot_product(cxa, a)
    mat(2,2) = dot_product(cxa, b)
    mat(2,3) = dot_product(cxa, c)
    mat(3,1) = dot_product(bxc, a)
    mat(3,2) = dot_product(bxc, b)
    mat(3,3) = dot_product(bxc, c)
!   We compute the inverse of this matrix analytically
    md = -mat(1,3)*mat(2,2)*mat(3,1) + mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(2,1)*mat(3,2) &
         -mat(1,1)*mat(2,3)*mat(3,2) - mat(1,2)*mat(2,1)*mat(3,3) + mat(1,1)*mat(2,2)*mat(3,3)
    mat_inv(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
    mat_inv(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
    mat_inv(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
    mat_inv(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
    mat_inv(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
    mat_inv(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)
    mat_inv(3,1) = mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
    mat_inv(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
    mat_inv(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
    mat_inv = mat_inv / md

    indices_real = 2.d0 * (/ rcut, rcut, rcut /)
    indices_real = matmul(mat_inv, indices_real)

    indices_real = dabs(indices_real)
    do i = 1, 3
      indices(i) = int( indices_real(i) )
      if( indices_real(i) > dfloat(indices(i)) )then
        indices(i) = indices(i) + 1
      end if
    end do


  end subroutine
!**************************************************************************






!**************************************************************************
! This subroutine reads in the XYZ file and builds the lists of neighbors, the spherical
! coordinates, etc.
!
  subroutine build_neighbors_list(positions, a_box, b_box, c_box, do_timing, &
                                  species_supercell, rcut_max, n_atom_pairs, rjs, &
                                  thetas, phis, xyz, n_neigh, neighbors_list, neighbor_species, &
                                  n_sites, indices, rebuild_neighbors_list, do_list, rank )

    implicit none

!   Input variables
    real*8, intent(in) :: rcut_max, positions(:,:)
!    integer, intent(in) :: species_multiplicity(:), n_species
    integer, intent(in) :: species_supercell(:), indices(1:3), n_sites, rank
    logical, intent(in) :: do_timing, rebuild_neighbors_list, do_list(:)

!   Output variables
    integer, intent(out) :: n_atom_pairs
!    logical, allocatable, intent(out) :: mask_species(:,:)

!   In and out variables
    real*8, allocatable, intent(inout) :: rjs(:), thetas(:), phis(:), xyz(:,:)
    integer, allocatable, intent(inout) :: neighbor_species(:)
    real*8, intent(inout) :: a_box(1:3), b_box(1:3), c_box(1:3)
    integer, allocatable, intent(inout) :: neighbors_list(:), n_neigh(:)
!    integer, allocatable, intent(inout) :: species_supercell(:,:)

!   Internal variables
    real*8 :: time1, time2, dist(1:3), d, neigh_time, time3, tol, d_tol = 1.d-6
    integer, allocatable :: neighbors_list_temp(:), head(:), this_list(:)
    integer :: n_neigh_max, i, j, n_sites_supercell, &
               k, k2, i2, j2, i3, j3, k3, mx, my, mz, i_shift(1:3)
    logical :: is_box_square, is_box_small
    logical, save :: print_cutoff_warning = .true., print_shape_warning = .true.


    if( do_timing )then
      call cpu_time(time1)
      time3 = time1
    end if

    n_sites_supercell = size(positions,2)

!   Very inefficiently build neighbor lists. I should write some routine that performs                              <-- FIX THIS
!   overlapping domain decomposition to make this more efficient
!    if( a_box(2) == 0.d0 .and. a_box(3) == 0.d0 .and. b_box(1) == 0.d0 .and. &
!        b_box(3) == 0.d0 .and. c_box(1) == 0.d0 .and. c_box(2) == 0.d0 )then
!        .and. size(positions,2) == n_sites )then
!   This assumes that if the non-diagonal components of the lattice vectors are approximately zero,
!   it is due to numerical noise, and thus the box is square
    if( dabs(a_box(2)) < d_tol .and. dabs(a_box(3)) < d_tol .and. &
        dabs(b_box(1)) < d_tol .and. dabs(b_box(3)) < d_tol .and. &
        dabs(c_box(1)) < d_tol .and. dabs(c_box(2)) < d_tol )then
      a_box(2) = 0.d0
      a_box(3) = 0.d0
      b_box(1) = 0.d0
      b_box(3) = 0.d0
      c_box(1) = 0.d0
      c_box(2) = 0.d0
      is_box_square = .true.
    else
      is_box_square = .false.
      if( rank == 0 .and. print_shape_warning )then
        print_shape_warning = .false.
        write(*,*)'                                       |'
        write(*,*)'WARNING: your simulation box is not    |  <-- WARNING'
        write(*,*)'orthorhombic; this will lead to slower |'
        write(*,*)'code execution. This warning will be   |'
        write(*,*)'printed only once. If you have more    |'
        write(*,*)'than one structure in your XYZ file,   |'
        write(*,*)'you may also have several instances of |'
        write(*,*)'non-orthorhombic cells.                |'
        write(*,*)'                                       |'
      end if
    end if
!
!   This is also inefficient <---------------------------------------- FIX THIS
    is_box_small = .false.
    if( any(indices > 1) )then
      is_box_small = .true.
      if( rank == 0 .and. print_cutoff_warning )then
        print_cutoff_warning = .false.
        write(*,*)'                                       |'
        write(*,*)'WARNING: your simulation box is smaller|  <-- WARNING'
        write(*,*)'than a neighbors cutoff sphere; this   |'
        write(*,*)'will lead to slow code execution due to|'
        write(*,*)'inefficient neighbor list builds. This |'
        write(*,*)'warning will be printed only once. If  |'
        write(*,*)'you have more than one structure in    |'
        write(*,*)'your XYZ file, you may also have       |'
        write(*,*)'several instances of non-orthorhombic  |'
        write(*,*)'cells.                                 |'
        write(*,*)'                                       |'
      end if
    end if
!
!   This is an initial guess for the maximum number of neighbors that we expect within a cutoff
!   It is used for the purpose of memory allocation. At the moment we are making this big, but
!   we should find a way to increase this value automatically whenever a bigger array is needed <---- fix
!   A way to fix it is, whenever an atom has more neighbors than n_neigh_max, to have store all
!   the neighbors lists into a temporary array, increase the size of the lists array, and then
!   store the lists back to it
    if( rebuild_neighbors_list )then
      n_neigh_max = 100
      allocate( neighbors_list(1:n_neigh_max*n_sites) )
      neighbors_list = 0
      allocate( n_neigh(1:n_sites) )
      n_neigh = 0
      n_atom_pairs = 0
    end if
!   We have an efficient algorithm for square boxes and inefficient for non-square boxes (sorry!)
!   Another requirement is that the minimum unit cell length is at least twice the cutoff
!
!   Tolerance in Angstrom for the ratio of positions to lattice vector. We add this because otherwise
!   atoms right at the periodic boundary can become problematic
    tol = 1.d-10
    if( is_box_square .and. (.not. is_box_small) .and. rebuild_neighbors_list )then
      mx = int( a_box(1) / rcut_max )
      my = int( b_box(2) / rcut_max )
      mz = int( c_box(3) / rcut_max )
      allocate( head(1:mx*my*mz) )
      head = 0
      allocate( this_list(1:n_sites) )
      do i = 1, n_sites
        call get_distance( [a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0], positions(1:3, i), &
                           a_box(1:3), b_box(1:3), c_box(1:3), (/ .true., .true., .true. /), dist, d, i_shift)
!       This is the position within the supercell, we must make sure it really is within the supercell
        dist = dist + [a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0]
        j = 1 + modulo(int( dist(1) / (a_box(1)+tol) * mx ), mx) &
              + modulo(int( dist(2) / (b_box(2)+tol) * my ), my) * mx &
              + modulo(int( dist(3) / (c_box(3)+tol) * mz ), mz) * my * mx
        this_list(i) = head(j)
        head(j) = i
      end do
      do i = 1, n_sites
        if( do_list(i) )then
!         We always count atom i as its own neighbor. This is useful when building the derivatives
          n_neigh(i) = n_neigh(i) + 1
          n_atom_pairs = n_atom_pairs + 1
          if( n_atom_pairs > n_neigh_max*n_sites )then
            allocate( neighbors_list_temp( 1:(n_neigh_max+10)*n_sites ) )
            neighbors_list_temp(1:n_atom_pairs-1) = neighbors_list(1:n_atom_pairs-1)
            deallocate( neighbors_list )
            allocate( neighbors_list( 1:(n_neigh_max+10)*n_sites ) )
            neighbors_list = neighbors_list_temp
            deallocate( neighbors_list_temp )
            n_neigh_max = n_neigh_max+10
          end if
          neighbors_list(n_atom_pairs) = i
!         Cell coordinates for this atom
          call get_distance( [a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0], positions(1:3, i), &
                             a_box(1:3), b_box(1:3), c_box(1:3), (/ .true., .true., .true. /), dist, d, i_shift)
          dist = dist + [a_box(1)/2.d0, b_box(2)/2.d0, c_box(3)/2.d0]
          i2 = 1 + int( dist(1) / (a_box(1)+tol) * mx )
          j2 = 1 + int( dist(2) / (b_box(2)+tol) * my )
          k2 = 1 + int( dist(3) / (c_box(3)+tol) * mz )
!         Look for other atoms in this and neighboring cells
          do k3 = k2-1, k2+1
            if( mz == 1 .and. k3 /= 1 )cycle
            if( mz == 2 .and. k2 == 1 .and. k3 == 0 )cycle
            if( mz == 2 .and. k2 == 2 .and. k3 == 3 )cycle
            do j3 = j2-1, j2+1
              if( my == 1 .and. j3 /= 1 )cycle
              if( my == 2 .and. j2 == 1 .and. j3 == 0 )cycle
              if( my == 2 .and. j2 == 2 .and. j3 == 3 )cycle
              do i3 = i2-1, i2+1
                if( mx == 1 .and. i3 /= 1 )cycle
                if( mx == 2 .and. i2 == 1 .and. i3 == 0 )cycle
                if( mx == 2 .and. i2 == 2 .and. i3 == 3 )cycle
                j = 1 + modulo(i3-1,mx) + modulo(j3-1,my)*mx + modulo(k3-1,mz)*mx*my
                k = head(j)
                do while( k /= 0 )
                  if( k /= i )then
                    call get_distance(positions(1:3, i), positions(1:3, k), a_box(1:3), b_box(1:3), &
                                      c_box(1:3), (/ .true., .true., .true. /), dist, d, i_shift)
                    if( d < rcut_max )then
                      n_neigh(i) = n_neigh(i) + 1
                      n_atom_pairs = n_atom_pairs + 1
                      if( n_atom_pairs > n_neigh_max*n_sites )then
                        allocate( neighbors_list_temp( 1:(n_neigh_max+10)*n_sites ) )
                        neighbors_list_temp(1:n_atom_pairs-1) = neighbors_list(1:n_atom_pairs-1)
                        deallocate( neighbors_list )
                        allocate( neighbors_list( 1:(n_neigh_max+10)*n_sites ) )
                        neighbors_list = neighbors_list_temp
                        deallocate( neighbors_list_temp )
                        n_neigh_max = n_neigh_max+10
                      end if
                      neighbors_list(n_atom_pairs) = k
                    end if
                  end if
                  k = this_list(k)
                end do
              end do
            end do
          end do
        end if
      end do
      deallocate(head, this_list)
!   Very inefficient algorithm for non-square boxes
    else if( rebuild_neighbors_list )then
      do i = 1, n_sites
        if( do_list(i) )then
!         We always count atom i as its own neighbor. This is useful when building the derivatives
          n_neigh(i) = n_neigh(i) + 1
          n_atom_pairs = n_atom_pairs + 1
          if( n_atom_pairs > n_neigh_max*n_sites )then
              allocate( neighbors_list_temp( 1:(n_neigh_max+10)*n_sites ) )
              neighbors_list_temp(1:n_atom_pairs-1) = neighbors_list(1:n_atom_pairs-1)
              deallocate( neighbors_list )
              allocate( neighbors_list( 1:(n_neigh_max+10)*n_sites ) )
              neighbors_list = neighbors_list_temp
              deallocate( neighbors_list_temp )
              n_neigh_max = n_neigh_max+10
          end if
          neighbors_list(n_atom_pairs) = i
          do j = 1, n_sites_supercell
            if( j /= i )then
              call get_distance(positions(1:3, i), positions(1:3, j), a_box(1:3), b_box(1:3), &
                                c_box(1:3), (/ .true., .true., .true. /), dist, d, i_shift)
              if( d < rcut_max )then
                n_neigh(i) = n_neigh(i) + 1
                n_atom_pairs = n_atom_pairs + 1
                if( n_atom_pairs > n_neigh_max*n_sites )then
                  allocate( neighbors_list_temp( 1:(n_neigh_max+10)*n_sites ) )
                  neighbors_list_temp(1:n_atom_pairs-1) = neighbors_list(1:n_atom_pairs-1)
                  deallocate( neighbors_list )
                  allocate( neighbors_list( 1:(n_neigh_max+10)*n_sites ) )
                  neighbors_list = neighbors_list_temp
                  deallocate( neighbors_list_temp )
                  n_neigh_max = n_neigh_max+10
                end if
!                j2 = mod(j-1, n_sites) + 1 
!                neighbors_list(n_atom_pairs) = j2
                neighbors_list(n_atom_pairs) = j
              end if
            end if
          end do
        end if
      end do
    end if


    if( rebuild_neighbors_list )then
      allocate( neighbors_list_temp(1:n_atom_pairs) )
      neighbors_list_temp = neighbors_list(1:n_atom_pairs)
      deallocate( neighbors_list )
      allocate( neighbors_list(1:n_atom_pairs) )
      neighbors_list = neighbors_list_temp
      deallocate( neighbors_list_temp )
    end if


    if( do_timing )then
      call cpu_time(time2)
      neigh_time = time2 - time1
      time1 = time2
    end if


!   NOTE on performance: looping over interactions I could have chosen to calculate each interaction only once,
!   i.e., (i,j) = (j,i), because if j is i's neighbor, the reverse is also true. This would reduce the
!   calculation load by half. However, this is a design feature, since the code is easier to parallelize if
!   each atom is accompanied by its full list of neighbors. It also prevents ackward memory access.
!   In any case, I may want to rethink this in the future if I want to further gain extra performance (at the
!   expense of complicating the code, that is).
!
    n_atom_pairs = 0
    do i = 1, n_sites
      n_atom_pairs = n_atom_pairs + n_neigh(i)
    end do
!    allocate( mask_species(1:n_atom_pairs, 1:n_species) )
!    mask_species = .false.
    if( rebuild_neighbors_list )then
      allocate( rjs(1:n_atom_pairs) )
      allocate( xyz(1:3, 1:n_atom_pairs) )
      allocate( thetas(1:n_atom_pairs) )
      allocate( phis(1:n_atom_pairs) )
      allocate( neighbor_species(1:n_atom_pairs) )
    end if
    k2 = 0
    do i = 1, n_sites
      if( do_list(i) )then
        do k = 1, n_neigh(i)
          k2 = k2 + 1
          j = neighbors_list(k2)
          if( k == 1 )then
            rjs(k2) = 0.d0
            xyz(1:3, k2) = (/ 0.d0, 0.d0, 0.d0 /)
            thetas(k2) = 0.d0
            phis(k2) = 0.d0
            neighbor_species(k2) = species_supercell(i)
!            do i2 = 1, species_multiplicity(i)
!              mask_species(k2, species_supercell(i2, j)) = .true.
!            end do
          else
            call get_distance(positions(1:3,i), positions(1:3,j), a_box(1:3), b_box(1:3), &
                              c_box(1:3), (/ .true., .true., .true. /), dist, d, i_shift)
            rjs(k2) = d
            xyz(1:3, k2) = dist
!           Avoid numerical artifacts
            if( dabs(dist(3)) >= d )then
              if( dist(3) > 0.d0 )then
                thetas(k2) = 0.d0
              else
                thetas(k2) = dacos(-1.d0)
              end if
            else
              thetas(k2) = dacos( dist(3) / d )
            end if
            phis(k2) = datan2( dist(2), dist(1) )
            neighbor_species(k2) = species_supercell(j)
!            do i2 = 1, species_multiplicity(i)
!              mask_species(k2, species_supercell(i2, j)) = .true.
!            end do
          end if
        end do
      else
        k2 = k2 + n_neigh(i)
      end if
    end do

    if( do_timing )then
      call cpu_time(time2)
      write(*,*)'                                       |'
      write(*,*)'Atoms timings (build):                 |'
      write(*,*)'                                       |'
      write(*,'(A, F9.3, A)') '  *) Neighbors build: ', neigh_time, ' seconds |'
      write(*,'(A, F7.3, A)') '  *) Spherical coords.: ', time2-time1, ' seconds |'
      write(*,'(A, F19.3, A)') '  *) Total: ', time2-time3, ' seconds |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
    end if

  end subroutine
!**************************************************************************




!**************************************************************************
  subroutine get_number_of_atom_pairs( n_neigh, rjs, rcut, l_max, n_max, max_Gbytes_per_process, &
                                       i_beg_list, i_end_list, j_beg_list, j_end_list )

    implicit none

!   Input variables
    real*8, intent(in) :: rjs(:), rcut, max_Gbytes_per_process
    integer, intent(in) :: n_neigh(:), l_max, n_max

!   Output variables
    integer, allocatable, intent(out) :: i_beg_list(:), i_end_list(:), j_beg_list(:), j_end_list(:)

!   Internal variables
    real*8 :: estimated_memory_in_Gbytes, mem_ratio, pairs_per_chunk
    integer :: n_sites, n_atom_pairs, k_max, n_chunks, i, j, k, k2, i_chunk, n_atom_pairs_in


    n_sites = size(n_neigh)
    n_atom_pairs = size(rjs)

    k = 0
    n_atom_pairs_in = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
        if( rjs(k) < rcut )then
          n_atom_pairs_in = n_atom_pairs_in + 1
        end if
      end do
    end do

    k_max = 1 + l_max*(l_max+1)/2 + l_max
!   This is a conservative estimate of the maximum memory that this run will need
    estimated_memory_in_Gbytes = dfloat(n_atom_pairs_in) * dfloat(n_max) * dfloat(k_max) * 150.d0 / 1024.d0**3
    mem_ratio = estimated_memory_in_Gbytes/max_Gbytes_per_process
    n_chunks = ceiling(mem_ratio)
    if( n_chunks > n_sites )then
      n_chunks = n_sites
    end if

    pairs_per_chunk = float(n_atom_pairs_in) / float(n_chunks)

    allocate( i_beg_list(1:n_chunks) )
    allocate( i_end_list(1:n_chunks) )
    allocate( j_beg_list(1:n_chunks) )
    allocate( j_end_list(1:n_chunks) )

    i_beg_list(1) = 1
    j_beg_list(1) = 1
    i_end_list(n_chunks) = n_sites
    j_end_list(n_chunks) = n_atom_pairs

    if( n_chunks == 1 )then
      return
    end if

    k = 0
    k2 = 0
    i_chunk = 1
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
        if( rjs(k) < rcut )then
          k2 = k2 + 1
        end if
      end do
      if( k2 >= int( float(i_chunk)*pairs_per_chunk ) )then
        i_end_list(i_chunk) = i
        j_end_list(i_chunk) = k
        i_chunk = i_chunk + 1
        i_beg_list(i_chunk) = i+1
        j_beg_list(i_chunk) = k+1
        if( i_chunk == n_chunks )then
          exit
        end if
      end if
    end do
    return

  end subroutine
!**************************************************************************





!**************************************************************************
!
! This subroutine returns the fractional coordinates from a list of
! Cartesian positions. This subroutine does NOT carry out unit cell
! wrapping. Wrapped Cartesian coordinates should be provided if wrapped
! fractional coordinates are wanted.
!
  subroutine get_fractional_coordinates(pos, a, b, c, frac)

    implicit none

!   Input variables
    real*8, intent(in) :: pos(:,:), a(1:3), b(1:3), c(1:3)
!   Output variables
    real*8, intent(out) :: frac(1:3,1:size(pos,2))
!   Internal variables
    real*8 :: L(1:3), d_tol = 1.d-6
    real*8 :: mat(1:3,1:3), md
    real*8, save :: a0(1:3) = 0.d0, b0(1:3) = 0.d0, c0(1:3) = 0.d0, mat_inv(1:3,1:3) = 0.d0
    integer :: i, atom, n_atoms
    logical :: lattice_check_a(1:3), lattice_check_b(1:3), lattice_check_c(1:3)

    n_atoms = size(pos, 2)

    if( dabs(a(2)) < d_tol .and. dabs(a(3)) < d_tol .and. &
        dabs(b(1)) < d_tol .and. dabs(b(3)) < d_tol .and. &
        dabs(c(1)) < d_tol .and. dabs(c(2)) < d_tol )then
!     Fast solution for orthorhombic cells
      L = (/ a(1), b(2), c(3) /)
      do atom = 1, n_atoms
        do i = 1, 3
          frac(i, atom) = pos(i, atom) / L(i)
        end do
      end do
    else
!     Slow solution for other unit cells
      lattice_check_a = ( a /= a0 )
      lattice_check_b = ( b /= b0 )
      lattice_check_c = ( c /= c0 )
      if( any(lattice_check_a) .or. any(lattice_check_b) .or. any(lattice_check_c) )then
        a0 = a
        b0 = b
        c0 = c
!       We construct our matrix to get the MIC only if the lattice vectors have changed
        mat(1:3,1) = a(1:3)
        mat(1:3,2) = b(1:3)
        mat(1:3,3) = c(1:3)
!       We compute the inverse of this matrix analytically
        md = -mat(1,3)*mat(3,1)*mat(2,2) + mat(2,1)*mat(1,3)*mat(3,2) + mat(1,2)*mat(3,1)*mat(2,3) - mat(1,1)*mat(2,3)*mat(3,2) &
             - mat(1,2)*mat(2,1)*mat(3,3) + mat(1,1)*mat(2,2)*mat(3,3)
        mat_inv(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
        mat_inv(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
        mat_inv(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
        mat_inv(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
        mat_inv(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
        mat_inv(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)
        mat_inv(3,1) = mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
        mat_inv(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
        mat_inv(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
        mat_inv = mat_inv / md
      end if
      do atom = 1, n_atoms
        frac(1:3, atom)  = matmul(mat_inv, pos(1:3, atom))
      end do
    end if

    return
  end subroutine get_fractional_coordinates
!**************************************************************************


end module neighbors
