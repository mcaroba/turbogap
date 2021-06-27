! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, gap.f90, is copyright (c) 2019-2021, Miguel A. Caro
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

module gap

  use splines

  contains

  subroutine get_soap_energy_and_forces(soap, soap_der, alphas, delta, zeta0, e0, Qs, &
                                        n_neigh, neighbors_list, xyz, do_forces, do_timing, &
                                        energies, forces, virial)
!   **********************************************
!   soap(1:n_soap, 1:n_sites)

    implicit none

    real*8, intent(in) :: soap(:,:), soap_der(:,:,:), alphas(:), delta, Qs(:,:), e0, zeta0, xyz(:,:)
    real*8, intent(out) :: energies(:), forces(:,:), virial(1:3,1:3)
    integer, intent(in) :: n_neigh(:), neighbors_list(:)
    logical, intent(in) :: do_forces, do_timing
    real*8, allocatable :: kernels(:,:), kernels_der(:,:), Qss(:,:), Qs_copy(:,:), this_Qss(:), &
                           kernels_copy(:,:)
    real*8 :: time1, time2, time3, energies_time, forces_time, zeta, this_force(1:3)
    integer :: n_sites, n_sparse, n_soap, i, j, k, l, j2, zeta_int, n_sites0, k1, k2
    logical :: is_zeta_int = .false.
!    integer, allocatable :: neighbors_beg(:), neighbors_end(:)


    if( dabs(zeta0-dfloat(int(zeta0))) < 1.d-5 )then
      is_zeta_int = .true.
      zeta_int = int(zeta0)
      zeta = dfloat(zeta_int)
    else
      zeta = zeta0
    end if  

!   Energies
    if( do_timing )then
      call cpu_time(time1)
      time3 = time1
    end if
    n_sparse = size(alphas)
    n_soap = size(soap, 1)
    n_sites = size(soap, 2)
    n_sites0 = size(forces, 2)

    allocate( kernels(1:n_sites, 1:n_sparse) )
    kernels = 0.d0
    allocate( kernels_copy(1:n_sites, 1:n_sparse) )
    call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
              kernels, n_sites)
!   We copy the kernels because it makes the matmul() operation (WHICH SHOULD BY THE WAY BE WRITTEN
!   USING LAPACK ROUTINES) a lot faster
    if( is_zeta_int )then
      kernels_copy = kernels**zeta_int
      energies = matmul(kernels_copy, alphas)
    else
      kernels_copy = kernels**zeta
      energies = matmul(kernels_copy, alphas)
    end if
    energies = delta**2 * energies + e0
    if( do_timing )then
      call cpu_time(time2)
      energies_time = time2 - time1
    end if



!   Forces
    if( do_forces )then
      if( do_timing )then
        call cpu_time(time1)
      end if
      allocate( kernels_der(1:n_sites, 1:n_sparse) )
      allocate( Qss(1:n_sites, 1:n_soap) )
      Qss = 0.d0
      allocate( Qs_copy(1:n_soap, 1:n_sparse) )
      allocate(this_Qss(1:n_soap))

      Qs_copy = Qs

      if( is_zeta_int )then
        kernels_der = kernels**(zeta_int-1)
      else
        kernels_der = kernels**(zeta-1.d0)
      end if
      if( n_sites < n_soap )then
        do i = 1, n_sites
          kernels_der(i,:) = kernels_der(i,:)*alphas(:)
        end do
      else
        do i = 1, n_soap
          Qs_copy(i,:) = Qs(i,:)*alphas(:)
        end do
      end if

      call dgemm("n", "t", n_sites, n_soap, n_sparse, -zeta*delta**2, kernels_der, n_sites, &
                 Qs_copy, n_soap, 0.d0, Qss, n_sites)

! EXPERIMENTAL CODE
!      call cpu_time(time1)
!      allocate( neighbors_beg(1:n_sites) )
!      allocate( neighbors_end(1:n_sites) )
!      l = 0
!      do i = 1, n_sites
!        neighbors_beg(i) = l + 1
!        do j = 1, n_neigh(i)
!          l = l + 1
!        end do
!        neighbors_end(i) = l
!      end do
! END EXPERIMENTAL CODE

      virial = 0.d0
      forces = 0.d0
!!$OMP parallel do private(i,j,l,j2,this_Qss)
      l = 0
      do i = 1, n_sites
        this_Qss = Qss(i,1:n_soap)
        do j = 1, n_neigh(i)
          l = l + 1
!         do l = neighbors_beg(i), neighbors_end(i)
          j2 = mod(neighbors_list(l)-1, n_sites0) + 1
          do k = 1, 3
            this_force(k) = dot_product(this_Qss, soap_der(k,:,l))
            forces(k, j2) = forces(k, j2) + this_force(k)
          end do
!         This is a many body potential, so there's no factor of 1/2 here
!          virial = virial + dot_product( this_force(1:3), xyz(1:3,l) )
          do k1 = 1, 3
            do k2 =1, 3
              virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,l) + this_force(k2)*xyz(k1,l))
            end do
          end do
        end do
      end do
!!$OMP end parallel do

! EXPERIMENTAL CODE
!      deallocate( neighbors_beg, neighbors_end )
!      call cpu_time(time2)
!      write(*,*) time2-time1
! END EXPERIMENTAL CODE

      if( do_timing )then
        call cpu_time(time2)
        forces_time = time2 - time1
      end if
    end if



!   Wrap it up
    deallocate(kernels, kernels_copy)
    if( do_forces )then
      deallocate(kernels_der, Qs_copy, Qss, this_Qss)
    end if

    if( do_timing )then
      call cpu_time(time2)
      write(*,*)'                                       |'
      write(*,*)'Prediction timings (SOAP):             |'
      write(*,*)'                                       |'
      write(*,'(A, F7.3, A)') '  *) Energy prediction: ', energies_time, ' seconds |'
      if( do_forces )then
        write(*,'(A, F7.3, A)') '  *) Forces prediction: ', forces_time, ' seconds |'
      end if
      write(*,'(A, F8.3, A)') '  *) Total prediction: ', time2-time3, ' seconds |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
    end if

  end subroutine







  subroutine get_2b_energy_and_forces(rjs, xyz, alphas, cutoff, rcut, buffer, delta, sigma, e0, Qs, &
                                      n_neigh, do_forces, do_timing, species, neighbor_species, &
                                      species1, species2, species_types, energies, forces, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: rjs(:), xyz(:,:), alphas(:), cutoff(:), delta, sigma, e0, Qs(:), rcut, buffer
    integer, intent(in) :: n_neigh(:), species(:), neighbor_species(:)
    character*8, intent(in) :: species_types(:), species1, species2
    logical, intent(in) :: do_forces, do_timing

!   Output variables
    real*8, intent(out) :: energies(:), forces(:,:), virial(1:3,1:3)

!   Internal variables
    real*8 :: time1, time2, fcut, pi, dfcut, this_force(1:3)
    integer :: n_sparse, i, j, k, n_sites, n_atom_pairs, s, sp1, sp2, n_sites0, k1, k2

    if( do_timing )then
      call cpu_time(time1)
    end if

    pi = dacos(-1.d0)

    n_sparse = size(alphas)
    n_sites = size(n_neigh)
    n_atom_pairs = size(rjs)
    n_sites0 = size(forces, 2)

!   Map species to index
    do i = 1, size(species_types)
      if( species1 == species_types(i) )then
        sp1 = i
        exit
      end if
    end do
    do i = 1, size(species_types)
      if( species2 == species_types(i) )then
        sp2 = i
        exit
      end if
    end do

!   Energy calculation
    energies = e0
    k = 0
    do i = 1, n_sites
      if( species(i) /= sp1 .and. species(i) /= sp2 )then
        k = k + n_neigh(i)
        cycle
      end if
      k = k + 1
      do j = 2, n_neigh(i)
        k = k + 1
        if( (species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
            (species(i) == sp2 .and. neighbor_species(k) == sp1) )then
          continue
        else
          cycle
        end if
        if( rjs(k) < rcut )then
          if( rjs(k) < rcut - buffer )then
            fcut = 1.d0
          else
            fcut = ( dcos( pi*(rjs(k) - rcut + buffer) / buffer ) + 1.d0 ) / 2.d0
          end if
          do s = 1, n_sparse
            energies(i) = energies(i) + delta**2 * alphas(s) * cutoff(s) * fcut * &
                          dexp(-0.5d0 * (rjs(k) - Qs(s))**2 / sigma**2)
          end do
        end if
      end do
    end do

!   Force calculation
    if( do_forces )then
      forces = 0.d0
      virial = 0.d0
      k = 0
      do i = 1, n_sites
        if( species(i) /= sp1 .and. species(i) /= sp2 )then
          k = k + n_neigh(i)
          cycle
        end if
        k = k + 1
        do j = 2, n_neigh(i)
          k = k + 1
          if( (species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
              (species(i) == sp2 .and. neighbor_species(k) == sp1) )then
            continue
          else
            cycle
          end if
          if( rjs(k) < rcut )then
            if( rjs(k) < rcut - buffer )then
              fcut = 1.d0
              dfcut = 0.d0
            else
              fcut = ( dcos( pi*(rjs(k) - rcut + buffer) / buffer ) + 1.d0 ) / 2.d0
              dfcut = pi / 2.d0 / buffer * dsin( pi*(rjs(k) - rcut + buffer) / buffer )
            end if
            do s = 1, n_sparse
              this_force(1:3) = - 2.d0 * delta**2 * alphas(s) * cutoff(s) * &
                              dexp(-0.5d0 * (rjs(k) - Qs(s))**2 / sigma**2) * &
                              xyz(1:3, k) / rjs(k) * ( (rjs(k) - Qs(s)) / sigma**2 * fcut + dfcut )
              forces(1:3,i) = forces(1:3,i) + this_force(1:3)
!              virial = virial - dot_product( this_force(1:3), xyz(1:3,k) )
              do k1 = 1, 3
                do k2 =1, 3
                  virial(k1, k2) = virial(k1, k2) - 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
                end do
              end do
            end do
          end if
        end do
      end do
!     Half contribution to the virial for pair potentials
      virial = 0.5d0 * virial
    end if

    if( do_timing )then
      call cpu_time(time2)
      write(*,*)'                                       |'
      write(*,*)'Prediction timings (2b):               |'
      write(*,*)'                                       |'
      write(*,'(A, F8.3, A)') '  *) Total prediction: ', time2-time1, ' seconds |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
    end if

  end subroutine








  subroutine get_core_pot_energy_and_forces(rjs, xyz, x, V, yp1, ypn, dVdx2, n_neigh, do_forces, do_timing, species, &
                                            neighbor_species, species1, species2, species_types, energies, forces, &
                                            virial )

    implicit none

!   Input variables
    real*8, intent(in) :: rjs(:), xyz(:,:), x(:), V(:), yp1, ypn, dVdx2(:)
    integer, intent(in) :: n_neigh(:), species(:), neighbor_species(:)
    character*8, intent(in) :: species_types(:), species1, species2
    logical, intent(in) :: do_forces, do_timing

!   Output variables
    real*8, intent(out) :: energies(:), forces(:,:), virial(1:3,1:3)

!   Internal variables
!   There are two ways of doing the core_pot interpolation; most efficient probably depends on
!   whether a small subset or big subset of the total number of atom pairs has a core potential
!   term associated to it. The current implementation is fast going over pairs, but slow computing
!   the splines. This will give the best performance for systems with many atom types
    real*8 :: V_int(1:1), dV_int(1:1), this_force(1:3)
!    real*8, allocatable :: V_int(:), dV_int(:)
    real*8 :: time1, time2, rcut
    integer :: n_sparse, i, j, k, n_sites, n_atom_pairs, s, sp1, sp2, n_sites0, k1, k2

    if( do_timing )then
      call cpu_time(time1)
    end if

    n_sparse = size(x)
    n_sites = size(n_neigh)
    n_atom_pairs = size(rjs)
    n_sites0 = size(forces, 2)

    rcut = maxval(x)

!   Whether we do arrays or not depends on use case. Maybe we can improve the code in the future
!   to make some choices leading to faster execution (i.e., figure out when it pays off to use
!   arrays, or use masks and index assignments to avoid unecessary computations with arrays)
!    allocate( V_int(1:n_atom_pairs) )
!    if( do_forces )then
!      allocate( dV_int(1:n_atom_pairs) )
!    end if

!   Map species to index
    do i = 1, size(species_types)
      if( species1 == species_types(i) )then
        sp1 = i
        exit
      end if
    end do
    do i = 1, size(species_types)
      if( species2 == species_types(i) )then
        sp2 = i
        exit
      end if
    end do

!   Energy calculation
    energies = 0.d0
    k = 0
    do i = 1, n_sites
      if( species(i) /= sp1 .and. species(i) /= sp2 )then
        k = k + n_neigh(i)
        cycle
      end if
      k = k + 1
      do j = 2, n_neigh(i)
        k = k + 1
        if( (species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
            (species(i) == sp2 .and. neighbor_species(k) == sp1) )then
          continue
        else
          cycle
        end if
        if( rjs(k) < rcut )then
          V_int = spline(x, V, dVdx2, yp1, ypn, rjs(k:k), rcut)
          energies(i) = energies(i) + 0.5d0 * V_int(1)
        end if
      end do
    end do

!   Force calculation
    if( do_forces )then
      forces = 0.d0
      virial = 0.d0
      k = 0
      do i = 1, n_sites
        if( species(i) /= sp1 .and. species(i) /= sp2 )then
          k = k + n_neigh(i)
          cycle
        end if
        k = k + 1
        do j = 2, n_neigh(i)
          k = k + 1
          if( (species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
              (species(i) == sp2 .and. neighbor_species(k) == sp1) )then
            continue
          else
            cycle
          end if
          if( rjs(k) < rcut )then
            dV_int = spline_der(x, V, dVdx2, yp1, ypn, rjs(k:k), rcut)
            this_force(1:3) = dV_int(1) * xyz(1:3, k) / rjs(k)
            forces(1:3,i) = forces(1:3,i) + this_force(1:3)
!            virial = virial - dot_product( this_force(1:3), xyz(1:3, k) )
            do k1 = 1, 3
              do k2 =1, 3
                virial(k1, k2) = virial(k1, k2) - 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
              end do
            end do
          end if
        end do
      end do
!     Half contribution to the virial for pair potentials
      virial = 0.5d0 * virial
    end if

    if( do_timing )then
      call cpu_time(time2)
      write(*,*)'                                       |'
      write(*,*)'Prediction timings (core_pot):         |'
      write(*,*)'                                       |'
      write(*,'(A, F8.3, A)') '  *) Total prediction: ', time2-time1, ' seconds |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
    end if

  end subroutine









!**************************************************************************
  subroutine get_3b_energy_and_forces( rjs, xyz, alphas, cutoff, rcut, buffer, delta, sigma, e0, Qs, &
                                       n_neigh, neighbors_list, do_forces, do_timing, kernel_type, &
                                       species, neighbor_species, species_center, species1, species2, &
                                       species_types, energies, forces, virial )

    implicit none

!   Input variables
    real*8, intent(in) :: rjs(:), xyz(:,:), alphas(:), cutoff(:), rcut, buffer, delta, sigma(:), e0, Qs(:,:)
    integer, intent(in) :: n_neigh(:), neighbors_list(:), species(:), neighbor_species(:)
    logical, intent(in) :: do_timing, do_forces
    character*3, intent(in) :: kernel_type
    character*8, intent(in) :: species_center, species1, species2, species_types(:)

!   Output variables 
    real*8, intent(out) :: energies(:), forces(:,:), virial(1:3,1:3)

!   Internal variables
    real*8 :: time1, time2, fcut, pi, r12, r13, r23, xyz12(1:3), xyz13(1:3), &
              xyz23(1:3), q(1:3), fcut12, fcut13, dfcut12(1:3), dfcut13(1:3), &
              force1(1:3), force2(1:3), force3(1:3), dfcut(1:3), &
              xyz12_red(1:3), xyz13_red(1:3), xyz23_red(1:3), this_force(1:3)
    real*8, allocatable :: r(:), drdq(:,:), kernel(:), drdx1(:,:), drdx2(:,:), drdx3(:,:), pref(:), &
                           kernel_der(:)
    integer :: n_sparse, i, j, k, k2, n_sites, n_atom_pairs, s, j2, i3, j3, k3, l, sp0, sp1, sp2, n_sites0, k1, k4

    if( do_timing )then
      call cpu_time(time1)
    end if

    pi = dacos(-1.d0)

    n_sparse = size(alphas)
    n_sites = size(n_neigh)
    n_atom_pairs = size(rjs)
    n_sites0 = size(forces, 2)

!   Map species to index
    do i = 1, size(species_types)
      if( species_center == species_types(i) )then
        sp0 = i
        exit
      end if
    end do
    do i = 1, size(species_types)
      if( species1 == species_types(i) )then
        sp1 = i
        exit
      end if
    end do
    do i = 1, size(species_types)
      if( species2 == species_types(i) )then
        sp2 = i
        exit
      end if
    end do

    allocate( r(1:n_sparse) )
    allocate( kernel(1:n_sparse) )
    allocate( pref(1:n_sparse) )
    if( do_forces )then
      allocate( drdq(1:n_sparse, 1:3) )
      allocate( drdx1(1:n_sparse, 1:3) )
      allocate( drdx2(1:n_sparse, 1:3) )
      allocate( drdx3(1:n_sparse, 1:3) )
      if( kernel_type == "pol" )then
        allocate( kernel_der(1:n_sparse) )
      end if
    end if

!   NOTE: these loops can be made between 2 and 3 times faster by summing over triplets instead of
!   over sites, however it may be a better strategy with parallelization in mind to sum over sites
!   (or maybe not...)

    pref(1:n_sparse) = 2.d0 * delta**2 * alphas(1:n_sparse) * cutoff(1:n_sparse)

!   Energy and force calculation
    energies = e0
    if( do_forces )then
      forces = 0.d0
      virial = 0.d0
    end if
    k = 0
    do i = 1, n_sites
      if( species(i) /= sp0 )then
        k = k + n_neigh(i)
        cycle
      end if
      k = k + 1
!      i3 = neighbors_list(k)
      i3 = mod(neighbors_list(k)-1, n_sites0) + 1
      do j = 2, n_neigh(i)
        k = k + 1
        r12 = rjs(k)
        if( (neighbor_species(k) /= sp1 .and. neighbor_species(k) /= sp2) .or. r12 > rcut )then
          cycle
        end if
!        j3 = neighbors_list(k)
        j3 = mod(neighbors_list(k)-1, n_sites0) + 1
        xyz12 = xyz(1:3, k)
        xyz12_red = xyz12 / r12
        k2 = k
        do j2 = j+1, n_neigh(i)
          k2 = k2 + 1
          r13 = rjs(k2)
          if( ( (neighbor_species(k) == sp1 .and. neighbor_species(k2) == sp2) .or. &
              (neighbor_species(k) == sp2 .and. neighbor_species(k2) == sp1) ) .and. r13 < rcut )then
            continue
          else
            cycle
          end if
!          k3 = neighbors_list(k2)
          k3 = mod(neighbors_list(k2)-1, n_sites0) + 1
          xyz13 = xyz(1:3, k2)
          xyz13_red = xyz13 / r13
          xyz23 = xyz13 - xyz12 
          r23 = dsqrt(sum(xyz23**2))
          xyz23_red = xyz23 / r23
!         It would be nice that if r23 < rcut, we only evaluate the expression if k3 > i3 and j3 > i3,
!         however that gets messy because the cutoff functions are not the same for each of the 3
!         evaluations of the triplet, i.e., atom 1 sees two cutoff functions, atom 2 sees two *different*
!         cutoff functions and 3 sees another 2 differetn cutoff functions. I actually tried and it
!         reduced the calculation time by about 40% for the large random C system with a 3 Angstrom cutoff
          if( r12 < rcut .and. r13 < rcut )then
            if( r12 < rcut - buffer )then
              fcut12 = 1.d0
              dfcut12(1:3) = 0.d0
            else
              fcut12 = ( dcos( pi*(r12 - rcut + buffer) / buffer ) + 1.d0 ) / 2.d0
              dfcut12(1:3) = dsin( pi*(r12 - rcut + buffer) / buffer ) / 2.d0 * pi / buffer * xyz12_red(1:3)
            end if
            if( r13 < rcut - buffer )then
              fcut13 = 1.d0
              dfcut13(1:3) = 0.d0
            else
              fcut13 = ( dcos( pi*(r13 - rcut + buffer) / buffer ) + 1.d0 ) / 2.d0
              dfcut13(1:3) = dsin( pi*(r13 - rcut + buffer) / buffer ) / 2.d0 * pi / buffer * xyz13_red(1:3)
            end if
            fcut = fcut12 * fcut13
            dfcut(1:3) = fcut12 * dfcut13(1:3) + fcut13 * dfcut12(1:3)
!           This builds the actual descriptor
            q = [ r12+r13, (r12-r13)**2, r23 ]
!           This gets the Euclidean distance between q and all the qs
            r(1:n_sparse) = (q(1) - Qs(1:n_sparse, 1))**2 / sigma(1)**2
            r(1:n_sparse) = r(1:n_sparse) + (q(2) - Qs(1:n_sparse, 2))**2 / sigma(2)**2
            r(1:n_sparse) = r(1:n_sparse) + (q(3) - Qs(1:n_sparse, 3))**2 / sigma(3)**2
            r(1:n_sparse) = dsqrt(r(1:n_sparse))
!           This gets the derivatives of the distances wrt the descriptors (the 1/r factor is not included)
            if( do_forces )then
              drdq(1:n_sparse, 1) = (q(1) - Qs(1:n_sparse, 1)) / sigma(1)**2
              drdq(1:n_sparse, 2) = (q(2) - Qs(1:n_sparse, 2)) / sigma(2)**2
              drdq(1:n_sparse, 3) = (q(3) - Qs(1:n_sparse, 3)) / sigma(3)**2
            end if
!           Evaluate the kernels
            if( kernel_type == "exp" )then
!             This kernel already contains the prefactor
              kernel(1:n_sparse) = pref(1:n_sparse) * dexp(-0.5d0 * r(1:n_sparse)**2)
              energies(i) = energies(i) + fcut * sum( kernel(1:n_sparse) )
              if( do_forces )then
                force1 = 0.d0
                force2 = 0.d0
                force3 = 0.d0
                do l = 1, 3
!                 For atom 1
                  drdx1(1:n_sparse, l) = drdq(1:n_sparse, 1) * ( xyz12_red(l) + xyz13_red(l) ) &
                                         + drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * ( xyz12_red(l) - xyz13_red(l) )
!                 For atom 2
                  drdx2(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz12_red(l) &
                                         - drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * xyz12_red(l) &
                                         + drdq(1:n_sparse, 3) * xyz23_red(l)
!                 For atom 3
                  drdx3(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz13_red(l) &
                                         - drdq(1:n_sparse, 2) * 2.d0 * (r13-r12) * xyz13_red(l) &
                                         - drdq(1:n_sparse, 3) * xyz23_red(l)
                  force1(l) = - sum( kernel(1:n_sparse) * (dfcut(l) + fcut*drdx1(1:n_sparse, l)) )
                  force2(l) = - sum( kernel(1:n_sparse) * (-fcut13 * dfcut12(l) + fcut*drdx2(1:n_sparse, l)) )
                  force3(l) = - sum( kernel(1:n_sparse) * (-fcut12 * dfcut13(l) + fcut*drdx3(1:n_sparse, l)) )
                end do
                forces(1:3, i3) = forces(1:3, i3) + force1(1:3)
                forces(1:3, j3) = forces(1:3, j3) + force2(1:3)
                forces(1:3, k3) = forces(1:3, k3) + force3(1:3)
!               force1 acting on i3 does not contribute to the virial
!                virial = virial + dot_product( force2(1:3), xyz12(1:3) )
!                virial = virial + dot_product( force3(1:3), xyz13(1:3) )
                do k1 = 1, 3
                  do k4 =1, 3
                    virial(k1, k4) = virial(k1, k4) + 0.5d0 * (force2(k1)*xyz12(k4) + force2(k4)*xyz12(k1))
                    virial(k1, k4) = virial(k1, k4) + 0.5d0 * (force3(k1)*xyz13(k4) + force3(k4)*xyz13(k1))
                  end do
                end do
              end if
            else if( kernel_type == "pol" )then
!             This kernel already contains the prefactor
              kernel(1:n_sparse) = pref(1:n_sparse) * cov_pp(r(1:n_sparse), 3, 1)
              energies(i) = energies(i) + fcut * sum( kernel(1:n_sparse) )
              if( do_forces )then
                force1 = 0.d0
                force2 = 0.d0
                force3 = 0.d0
!               This derivative contains the 1/r term
                kernel_der(1:n_sparse) = pref(1:n_sparse) * cov_pp_der(r(1:n_sparse), 3, 1)
                do l = 1, 3
!                 For atom 1
                  drdx1(1:n_sparse, l) = drdq(1:n_sparse, 1) * ( xyz12_red(l) + xyz13_red(l) ) &
                                         + drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * ( xyz12_red(l) - xyz13_red(l) )
!                 For atom 2
                  drdx2(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz12_red(l) &
                                         - drdq(1:n_sparse, 2) * 2.d0 * (r12-r13) * xyz12_red(l) &
                                         + drdq(1:n_sparse, 3) * xyz23_red(l)
!                 For atom 3
                  drdx3(1:n_sparse, l) = - drdq(1:n_sparse, 1) * xyz13_red(l) &
                                         - drdq(1:n_sparse, 2) * 2.d0 * (r13-r12) * xyz13_red(l) &
                                         - drdq(1:n_sparse, 3) * xyz23_red(l)
                  force1(l) = - sum( kernel(1:n_sparse) * dfcut(l) &
                                     - kernel_der(1:n_sparse) * fcut*drdx1(1:n_sparse, l) )
                  force2(l) = - sum( - kernel(1:n_sparse) * fcut13 * dfcut12(l) &
                                     - kernel_der(1:n_sparse) * fcut*drdx2(1:n_sparse, l) )
                  force3(l) = - sum( - kernel(1:n_sparse) * fcut12 * dfcut13(l) &
                                     - kernel_der(1:n_sparse) * fcut*drdx3(1:n_sparse, l) )
                end do
                forces(1:3, i3) = forces(1:3, i3) + force1(1:3)
                forces(1:3, j3) = forces(1:3, j3) + force2(1:3)
                forces(1:3, k3) = forces(1:3, k3) + force3(1:3)
!               force1 acting on i3 does not contribute to the virial
!                virial = virial + dot_product( force2(1:3), xyz12(1:3) )
!                virial = virial + dot_product( force3(1:3), xyz13(1:3) )
                do k1 = 1, 3
                  do k4 =1, 3
                    virial(k1, k4) = virial(k1, k4) + 0.5d0 * (force2(k1)*xyz12(k4) + force2(k4)*xyz12(k1))
                    virial(k1, k4) = virial(k1, k4) + 0.5d0 * (force3(k1)*xyz13(k4) + force3(k4)*xyz13(k1))
                 end do
               end do
              end if
            end if
          end if
        end do
      end do
    end do


    deallocate( kernel, pref, r )
    if( do_forces )then
      deallocate( drdq, drdx1, drdx2, drdx3 )
      if( kernel_type == "pol" )then
        deallocate( kernel_der )
      end if
    end if

    if( do_timing )then
      call cpu_time(time2)
      write(*,*)'                                       |'
      write(*,*)'Prediction timings (3b):               |'
      write(*,*)'                                       |'
      write(*,'(A, F8.3, A)') '  *) Total prediction: ', time2-time1, ' seconds |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
    end if

  end subroutine
!**************************************************************************







!**************************************************************************
!
! This is two functions
!
  function cov_pp(r, d, q) result(cov)

    implicit none

    real*8, intent(in) :: r(:)
    integer, intent(in) :: d, q
    real*8, dimension(1:size(r)) :: cov

    real*8 :: j
    integer :: j_int, i

    j_int = d/2 + q + 1
    j = dfloat(j_int)

    if( d == 3 .and. q == 1 )then
      do i = 1, size(r)
        if( r(i) >= 1.d0 )then
          cov(i) = 0.d0
        else if( r(i) <= 0.d0 )then
          cov(i) = 1.d0
        else
          cov(i) = (1.d0 - r(i))**(j_int+1) * ((j+1.d0)*r(i) + 1.d0)
        end if
      end do
    else
      write(*,*) "ERROR: This combination of input parameters is not currently supported for cov_pp!"
      stop
    end if

  end function
  function cov_pp_der(r, d, q) result(cov_der)

    implicit none

    real*8, intent(in) :: r(:)
    integer, intent(in) :: d, q
    real*8, dimension(1:size(r)) :: cov_der

    real*8 :: j
    integer :: j_int, i

    j_int = d/2 + q + 1
    j = dfloat(j_int)

    if( d == 3 .and. q == 1 )then
      do i = 1, size(r)
        if( r(i) >= 1.d0 .or. r(i) <= 0.d0)then
          cov_der(i) = 0.d0
        else
!         This expression contains the 1/r factor
          cov_der(i) = ( -(j+1.d0) * (1.d0 - r(i))**j_int * ((j+1.d0)*r(i) + 1.d0) + &
                         (1.d0 - r(i))**(j_int+1) * (j + 1.d0) ) / r(i)
        end if
      end do
    else
      write(*,*) "ERROR: This combination of input parameters is not currently supported for cov_pp!"
      stop
    end if

  end function
!**************************************************************************





end module gap
