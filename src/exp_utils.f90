! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, exp_utils.f90, is copyright (c) 2019-2023,
! HND X   Miguel A. Caro and Tigany Zarrouk
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

module exp_utils
  use types

contains





  subroutine get_pair_distribution_forces(  n_sites0, energy_scale, y_exp, forces0, virial,  &
       & neighbors_list, n_neigh, neighbor_species, rjs, xyz, r_min, r_max, n_samples,&
       & pair_distribution, r_cut,&
       &  species_1, species_2,  pair_distribution_der, partial_rdf, kde_sigma, c_factor )
    implicit none
    real*8,  intent(in) :: rjs(:), xyz(:,:), y_exp(:), energy_scale,  kde_sigma, c_factor
    integer, intent(in) :: n_sites0
    real*8 :: r_min, r_max, r_cut
    integer, intent(in) :: neighbors_list(:), n_neigh(:), neighbor_species(:)
    integer, intent(in) :: n_samples, species_1, species_2
    integer :: n_sites, n_pairs, count, count_species_1
    integer :: i, j, k, ki, k1, k2,  i2, j2, l, ii, jj, kk, i3, j3, i4,  species_i, species_j
    real*8,  intent(in) :: pair_distribution(1:n_samples)
    real*8,  intent(in) :: pair_distribution_der(:,:)
    real*8,  intent(inout) :: forces0(:,:), virial(1:3,1:3)
    real*8, allocatable ::  prefactor(:)
    real*8 :: r, n_pc, this_force(1:3), f
    real*8, parameter :: pi = acos(-1.0)
    logical, intent(in) :: partial_rdf
    logical :: species_in_list, counted_1=.false.
    ! First allocate the pair correlation function array

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)

    allocate( prefactor( 1:n_samples ) )

    prefactor =  ( pair_distribution - y_exp )

    ! Not reinitalizing the forces as these will be added to for each partial
    !    forces0 = 0.d0
    k = 0
    do i = 1, n_sites
       ! k is the index which keeps a number of the atom pairs
       ! i2 is the index of a particular atom
       i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
       species_i = neighbor_species(k+1)
       do j = 1, n_neigh(i)
          ! Loop through the neighbors of atom i
          ! j2 is the index of a neighboring atom to i
          k = k + 1
          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1

          species_j = neighbor_species(k)

          if (partial_rdf)then
             if (species_i /= species_1) cycle
             if (species_j /= species_2) cycle
          end if

          r = rjs(k) ! atom pair distance

          if (r < 1e-3 .or. r > r_cut) cycle
          if (r < r_min) cycle
          if (r > r_max + kde_sigma*6.d0) cycle


          if ( .not. all( xyz( 1:3, k ) == 0.d0 ) )then
             ! Actual derivative of the pair distribution function given here, no normalisation needed I think

             this_force(1) = dot_product( - 2.d0 * xyz( 1, k ) *&
                  & pair_distribution_der(1:n_samples,  k),&
                  & prefactor(1:n_samples))

             this_force(2) = dot_product( - 2.d0 * xyz( 2, k ) *&
                  & pair_distribution_der(1:n_samples,  k),&
                  & prefactor(1:n_samples))

             this_force(3) = dot_product( - 2.d0 * xyz( 3, k ) *&
                  & pair_distribution_der(1:n_samples,  k),&
                  & prefactor(1:n_samples))

             if (species_1 == species_2) f = 1.d0
             if (species_1 /= species_2) f = 2.d0

             this_force(1:3) =  - f * energy_scale * this_force(1:3)

             forces0(1:3, j2) = forces0(1:3, j2) + this_force(1:3)

             !             Sign is plus because this force is acting on j2. Factor of one is because this is
             !             derived from a local energy
             !              virial = virial + dot_product(this_force(1:3), xyz(1:3,k))
             do k1 = 1, 3
                do k2 =1, 3
                   virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
                end do
             end do
          end if

       end do
    end do


    deallocate(prefactor)

  end subroutine get_pair_distribution_forces


  subroutine get_pair_distribution(  n_sites0, &
       & neighbors_list, n_neigh, neighbor_species, rjs, xyz, r_min, r_max, n_samples, x,&
       & pair_distribution, r_cut, return_histogram,&
       & partial_rdf, species_1, species_2, kde_sigma, rho, do_derivatives, pair_distribution_der, n_dim_idx, j_beg, j_end )
    implicit none
    real*8,  intent(in) :: rjs(:), kde_sigma, rho, xyz(:,:)
    real*8 :: r_min, r_max, r_cut
    integer, intent(in) :: neighbors_list(:), n_neigh(:), neighbor_species(:)
    integer, intent(in) :: n_sites0, n_samples, species_1, species_2, n_dim_idx, j_beg, j_end
    integer :: n_sites, n_pairs, count, count_species_1
    integer :: i, j, k, ki,  i2, j2, l, ii, jj, kk, i3, j3, i4, ind_bin_l, ind_bin_h, species_i, species_j
    real*8,  intent(inout) :: pair_distribution(1:n_samples), x(1:n_samples)
    real*8,  intent(inout), allocatable :: pair_distribution_der(:,:,:)
    real*8, allocatable :: bin_edges(:), dV(:), kde(:)
    real*8 :: r, n_pc, ri_vec(1:3), rj_vec(1:3)
    real*8, parameter :: pi = acos(-1.0)
    logical, intent(in) :: partial_rdf, do_derivatives
    logical :: return_histogram, species_in_list, counted_1=.false.
    ! First allocate the pair correlation function array

    allocate(bin_edges(1:n_samples+1))
    if ( .not. return_histogram ) allocate(dV(1:n_samples))

    if (kde_sigma > 0.d0)then
       allocate(kde(1:n_samples))
       kde = 0.d0
    end if

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)


    ! ! For debuggind the derivatives
    ! do_derivatives = .true.
    ! if ( do_derivatives )then
    !    allocate( pair_distribution_der(1:n_samples, 1:3, 1:n_pairs) )
    ! end if



    x = 0.d0
    bin_edges = 0.d0
    pair_distribution = 0.d0

    ! check that r_min is less than r_max
    if (r_min > r_max)then
       print *, "!!! Given r_min is less than r_max! swapping these values!"
       r = r_max
       r_max = r_min
       r_min = r
    end if


    do i = 1, n_samples + 1
       bin_edges(i) = r_min  +  ( real( (i-1) ) /  real(n_samples) ) * (r_max - r_min)
    end do

    do i = 1, n_samples
       x(i) = (bin_edges(i) + bin_edges(i+1)) / 2.d0

       if ( .not. return_histogram )then
          dV(i) = 4.d0 * pi * (bin_edges(i)**2 * (bin_edges(i+1) - bin_edges(i)) )
       end if

    end do

    count = 0
    count_species_1 = 0
    k = 0
    do i = 1, n_sites
       ! k is the index which keeps a number of the atom pairs
       ! i2 is the index of a particular atom
       i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
       species_i = neighbor_species(k+1)
       do j = 1, n_neigh(i)
          ! Loop through the neighbors of atom i
          ! j2 is the index of a neighboring atom to i
          k = k + 1
          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1

          species_j = neighbor_species(k)

          if (partial_rdf)then
             if (species_i /= species_1) then
                cycle
             elseif(j == 1)then
                count_species_1 = count_species_1 + 1
             end if

          end if

          if (partial_rdf)then
             if (species_j /= species_2) cycle
          end if

          r = rjs(k) ! atom pair distance

          if ( r > r_cut )then
             ! if i2 == j2 or r out of range
             cycle
          end if

          ! edge cases where r is not in range
          if (r < r_min)then
             !            print *, "pair_distribution_function: Given r is less than r_min! Continuing loop "
             cycle
          end if

          if (r > r_max + kde_sigma*6.d0 )then
             !            print *, "pair_distribution_function: Given r is more than r_max! Continuing loop "
             cycle
          end if

          if (kde_sigma > 0.d0)then
             ! do a kernel density estimate
             count = count + 1
             kde = 0.d0
             do l = 1,n_samples
                kde(l) = kde(l) +  exp( -( (x(l) - r) / kde_sigma )**2 / 2.d0 )
             end do
             if (k /= 1) pair_distribution(1:n_samples) = pair_distribution(1:n_samples) + &
                  & kde(1:n_samples)

             if (do_derivatives)then
                ! reuse kde to get the derivatives
                do l = 1,n_samples
                   kde(l) = kde(l)  *  ( (x(l) - r) / kde_sigma**2 ) / r / dV(l)
                end do

                ! the first term is ( delta_jk - delta_ik )*( r_j^alpha - r_i^alpha )
                ! which can be double counted as per the vdw paper as
                ! 2 * ( - delta_ik )*( r_j^alpha - r_i^alpha )
                ! therefore we only have a contributon at i == k

                ! Calculating the derivatives this way takes up too much memory if the cutoff is very large.
                ! There must be some other way
                ! do i4 = 1,3
                !    pair_distribution_der(1:n_samples, i4, k) = - 2.d0 * xyz( i4, k ) * kde
                ! end do

                ! Calculating the derivatives this way simply takes up too much memory
                ! There must be some other way
                ! do i4 = 1,3
                !    pair_distribution_der(1:n_samples, i4, k) = - 2.d0 * xyz( i4, k ) * kde
                ! end do

                pair_distribution_der(1:n_samples, n_dim_idx,  k) = &
                     & pair_distribution_der(1:n_samples, n_dim_idx, &
                     & k) + kde(1:n_samples)


             end if
          else

             call binary_search_real( r, bin_edges, 1, n_samples+1, ind_bin_l, ind_bin_h )
             pair_distribution(ind_bin_l) = pair_distribution(ind_bin_l) +  1.d0

          end if

       end do
    end do


    ! The factor of 1 / ( n_samples * kde_sigma * sqrt(2 pi) )
    !   is for the kernel density estimate.
    ! Not sure exactly why multiplying by (r_max - r_min) works
    !   pdf = kde * dr / norm
    if (kde_sigma>0.d0)then
       pair_distribution(1:n_samples) =&
            & pair_distribution(1:n_samples) * ( ( r_max - r_min) / dfloat(n_samples) ) &
            & / ( sqrt( 2.d0 * pi) * kde_sigma)

       if ( do_derivatives )then
          pair_distribution_der(1:n_samples,n_dim_idx, 1:n_pairs) =&
               & pair_distribution_der(1:n_samples, n_dim_idx, 1:n_pairs) * ( ( r_max - r_min) / dfloat(n_samples) ) &
               & / ( sqrt( 2.d0 * pi) * kde_sigma)
       end if

    end if

    !( dfloat(count_species_1) )

    ! Remember that this is done for each process, so then the final
    ! result should be a reduction summation of all of these
    ! pair_distribution arrays multiplied by the 1/density factor ( pair_distribution *  (V / n_sites) )

    if ( .not. return_histogram )then
       ! Instead of giving the Radial distribution function (which goes as r^2), we give the pair distribution function
       pair_distribution =  pair_distribution / dV
       deallocate(dV)
    end if

    deallocate(bin_edges)
    if (kde_sigma>0.d0) deallocate(kde)

  end subroutine get_pair_distribution


  subroutine binary_search_real( xs, x, lin, hin, l, h )
    ! give delimiting indices (lout, hout) of array of x (which is a set of delimited bins) where xs (x search) is in range
    real*8, allocatable, intent(in) :: x(:)
    real*8, intent(in) :: xs
    real*8 :: xi
    integer, intent(in) :: lin, hin
    integer, intent(out) :: l, h
    integer :: m
    logical :: found = .false.

    l = lin
    h = hin

    ! Check that lower index, l is less than upper index, h
    if (l > h)then
       print *, "binary_search_real: lower bound of search is greater than higher bound, swapping"
       m = h
       h = l
       l = h
    end if

    m = int( real(l + h) / 2.d0 )

    found = .false.
    do while( .not. found  )
       xi = x(m)

       if ( xs < xi )then
          ! its in the lower segment
          h = m
          m = int( real(h + l) / 2.d0 )
       end if

       if ( xs > xi )then
          ! its in the upper segment
          l = m
          m = int( real(h + l) / 2.d0 )
       end if


       !         print *, " h - l = ", h - l
       if ( h - l == 1 ) then
          ! terminate the search
          found = .true.
       end if
    end do

  end subroutine binary_search_real


  subroutine get_structure_factor_explicit( n_sites0, &
       & neighbors_list, n_neigh, neighbor_species, rjs, n_samples, q_list,&
       & structure_factor, r_min, r_max, r_cut,&
       & partial, species_1, species_2, window )
    implicit none
    real*8,  intent(in) :: rjs(:)
    real*8 :: r_min, r_max, r_cut
    integer, intent(in) :: neighbors_list(:), n_neigh(:), neighbor_species(:)
    integer, intent(in) :: n_sites0, n_samples, species_1, species_2
    integer :: n_sites, n_pairs, count, count_species_1
    integer :: i, j, k, i2, j2, l, ind_bin_l, ind_bin_h, species_i, species_j
    real*8,  intent(inout) :: structure_factor(1:n_samples), q_list(1:n_samples)
    real*8 :: r, n_pc, w
    real*8, parameter :: pi = acos(-1.0)
    logical, intent(in) :: partial, window
    logical ::  species_in_list, counted_1=.false.
    ! First allocate the pair correlation function array

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)

    structure_factor = 0.d0

    count = 0
    count_species_1 = 0
    k = 0
    w = 1.d0
    do i = 1, n_sites
       ! k is the index which keeps a number of the atom pairs
       ! i2 is the index of a particular atom
       i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
       species_i = neighbor_species(k+1)

       do j = 1, n_neigh(i)
          ! Loop through the neighbors of atom i
          ! j2 is the index of a neighboring atom to i
          k = k + 1
          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1

          species_j = neighbor_species(k)

          if (partial)then
             if (species_i /= species_1) then
                cycle
             elseif(j == 1)then
                count_species_1 = count_species_1 + 1
             end if

          end if

          if (partial)then
             if (species_j /= species_2) cycle
          end if

          r = rjs(k) ! atom pair distance

          if ( r > r_cut)then
             cycle
          end if

          ! edge cases where r is not in range
          if (r < r_min)then
             !            print *, "pair_distribution_function: Given r is less than r_min! Continuing loop "
             cycle
          end if

          if (r > r_max)then
             !            print *, "pair_distribution_function: Given r is more than r_max! Continuing loop "
             cycle
          end if

          do l = 1, n_samples
             if (window) w = sinc( pi * r / r_cut )
             structure_factor(l) = structure_factor(l) +  sinc(&
                  & 2.d0 * pi * q_list(l) * r  ) * w /&
                  & dfloat(n_sites0)
          end do


       end do
    end do


  end subroutine get_structure_factor_explicit



  subroutine get_structure_factor_from_pdf( q_beg, q_end, structure_factor&
       &, pair_distribution, q_list , rs, r_cut,&
       & n_samples_pc, n_samples_sf, n_species, n_atoms_of_species,&
       & n_sites, rho, window)
    implicit none
    real*8 , intent(in), allocatable :: n_atoms_of_species(:)
    real*8, intent(out) :: structure_factor(:)
    real*8 , intent(in) ::  pair_distribution(:), q_list(:), rs(:), r_cut
    integer, intent(in) :: n_samples_pc, n_samples_sf, n_species, n_sites, q_beg, q_end
    real*8 :: r, q, dr, ca, cb, cabh, w, rho
    integer :: i, j, k, l, idx,  n
    real*8, parameter :: pi = acos(-1.0)
    logical, intent(in) :: window

    ! S_ab(q) = delta_ab + 4 pi rho (ca cb)^0.5
    !                  * int_0^r_cut dr r^2 [ g_ab(r) - 1 ] sin(qr)/(qr) * sin( pi r / R )/ (pi r /R)
    ! window corresponds to sin(pi r / R) / (pi r / R)
    ! rs has size n_samples_pc

    structure_factor = 0.d0
    ! n = q_end - q_beg + 1 !size(q_list)

    dr = rs(2) - rs(1)

    w = 1.d0

    ! Iterate through each q
    do k = q_beg, q_end
       q = q_list(k) * 2.d0 * pi

       do l = 1, n_samples_pc
          ! do the integral
          if (window) w = sinc( pi * rs(l) / r_cut )

          structure_factor(k) = structure_factor(k) + &
               & dr * rs(l)**2 &
               & * ( pair_distribution(l) - 1.d0 ) &
               & * sinc( q * rs(l) ) * w
       end do

       structure_factor(k) = 4.d0 * pi * rho * structure_factor(k)

    end do

  end subroutine get_structure_factor_from_pdf

  subroutine get_sinc_factor_matrix( q_beg, q_end, q_list, rs, r_cut,&
       & n_samples_pc, n_samples_sf, window, sinc_factor_matrix )
    implicit none
    real*8, allocatable, intent(in) :: q_list(:), rs(:)
    real*8, intent(in) :: r_cut
    integer, intent(in) :: n_samples_pc, n_samples_sf, q_beg, q_end
    real*8, allocatable, intent(out) :: sinc_factor_matrix(:,:)
    logical, intent(in) :: window
    real*8 :: dr, q, w=1.d0
    integer :: i, j
    real*8, parameter :: pi = acos(-1.0)
    ! This matrix can be used to get the derivatives. We can parallelize the generation of it on each process
    ! The only thing we have to worry about is memory if n_samples_pc/sf are too large
    allocate( sinc_factor_matrix( 1:n_samples_sf, 1:n_samples_pc ) )
    sinc_factor_matrix = 0.d0

    dr = rs(2) - rs(1)

    do j = 1, n_samples_pc
       if (window) w = sinc( pi * rs(j) / r_cut )

       do i = q_beg, q_end

          q = q_list(i) * 2.d0 * pi
          sinc_factor_matrix( i, j ) = dr * rs(j)**2 * sinc( q * rs(j) ) * w

       end do
    end do
  end subroutine get_sinc_factor_matrix


    



  ! subroutine get_partial_structure_factor_derivatives( &
  !      & pair_distribution_partial_der , q_list , rs, r_cut,&
  !      & n_samples_pc, n_samples_sf, n_species, n_dim_partial,  n_atoms_of_species,&
  !      & n_sites, rho,  window, j_beg, j_end, n_dim_partial,  structure_factor_partial_der)
  !   implicit none
  !   real*8 , intent(in), allocatable :: n_atoms_of_species(:)
  !   real*8, intent(in) :: pair_distribution_partial_der(:,:,:), q_list(:), rs(:), r_cut
  !   integer, intent(in) :: n_samples_pc, n_samples_sf, n_species, n_sites, q_beg, q_end, n_dim_partial
  !   real*8, intent(out) :: structure_factor_partial_der(1:n_samples_sf, 1:n_dim_partial, 1:3,  j_beg:j_end)
  !   real*8,  intent(in) :: rjs(:)
  !   real*8 :: r_min, r_max, r_cut
  !   integer, intent(in) :: neighbors_list(:), n_neigh(:), neighbor_species(:)
  !   integer, intent(in) :: n_sites0, n_samples, species_1, species_2
  !   real*8 :: r, q, dr, ca, cb, cabh, w
  !   integer :: i, j, k, l, i2, j2, i3, j3,  idx,  n, n_dim_idx
  !   integer, allocatable :: indexes_ndim(:,:)
  !   real*8, parameter :: pi = acos(-1.0)
  !   real*8, intent(in) ::  rho
  !   logical, intent(in) :: window

  !   ! Here we actually don't parallelize over q, we just go over all q
  !   ! values, and use the in-built parallelism of the splitting of the rdf gradients.
  !   structure_factor_partial_der = 0.d0

  !   allocate( indexes_ndim( 1:n_species, 1:n_species ) )
  !   indexes_ndim = 0

  !   n_dim_idx = 1
  !   outer: do i = 1, n_species
  !      do j = 1, n_species

  !         if (i > j) cycle

  !         indexes_ndim(i,j) =  n_dim_idx
  !         indexes_ndim(j,i) =  n_dim_idx

  !         n_dim_idx = n_dim_idx + 1

  !         if (n_dim_idx > n_dim_partial) exit outer

  !      end do
  !   end do outer

  !   dr = rs(2) - rs(1)
  !   w = 1.d0

  !   ! Now we integrate over r for every q for each k value in the

  !   k = 0
  !   do i = 1, n_sites
  !      ! k is the index which keeps a number of the atom pairs
  !      ! i2 is the index of a particular atom
  !      i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
  !      species_i = neighbor_species(k+1)
  !      ca = n_atoms_of_species(i) / dfloat( n_sites )
  !      do j = 1, n_neigh(i)
  !         ! Loop through the neighbors of atom i
  !         ! j2 is the index of a neighboring atom to i
  !         k = k + 1
  !         j2 = modulo(neighbors_list(k)-1, n_sites0) + 1

  !         species_j = neighbor_species(k)

  !         cb = n_atoms_of_species(j) / dfloat( n_sites )

  !         cabh = ( ca * cb )**(0.5)

  !         n_dim_idx = indexes_ndim( species_i, species_j )


  !         r = rjs(k) ! atom pair distance

  !         if ( r < 1e-3 .or. r > r_cut)then
  !            cycle
  !         end if

  !         ! edge cases where r is not in range
  !         if (r < r_min)then
  !            !            print *, "pair_distribution_function: Given r is less than r_min! Continuing loop "
  !            cycle
  !         end if

  !         if (r > r_max)then
  !            !            print *, "pair_distribution_function: Given r is more than r_max! Continuing loop "
  !            cycle
  !         end if


  !         do n = 1, n_samples_sf
  !            q = q_list(n) * 2.d0 * pi

  !            do l = 1, n_samples_pc
  !               ! do the integral
  !               if (window) w = sinc( pi * rs(l) / r_cut )

  !               structure_factor_partial_der(n, n_dim_idx, 1:3, k) = &
  !                    & structure_factor_partial_der(n, n_dim_idx, 1:3,  k) + &
  !                    & dr * rs(l)**2 &
  !                    & * ( - 2.d0 * xyz(1:3, l) *  pair_distribution_partial_der(l, n_dim_idx, k)  ) &
  !                    & * sinc( q * rs(l) ) * w
  !            end do

  !            structure_factor_partial(k,n_dim_idx) = 4.d0 * pi * cabh * rho * structure_factor_partial(k,n_dim_idx)

  !            if (i == j)then
  !               structure_factor_partial(k,n_dim_idx)  = structure_factor_partial(k,n_dim_idx)  + 1.d0
  !            end if

  !         end do
  !      end if

  !   end do
  ! end subroutine get_partial_structure_factor_derivatives


  subroutine get_partial_structure_factor( q_beg, q_end, &
       & pair_distribution_partial, q_list , rs, r_cut,&
       & n_samples_pc, n_samples_sf, n_species, n_dim_partial,  n_atoms_of_species,&
       & n_sites, rho,  window,  structure_factor_partial)
    implicit none
    real*8 , intent(in), allocatable :: n_atoms_of_species(:)
    real*8, intent(in) :: pair_distribution_partial(:,:), q_list(:), rs(:), r_cut
    integer, intent(in) :: n_samples_pc, n_samples_sf, n_species, n_sites, q_beg, q_end, n_dim_partial
    real*8, intent(out) :: structure_factor_partial(1:n_samples_sf,1:n_dim_partial)
    real*8 :: r, q, dr, ca, cb, cabh, w
    integer :: i, j, k, l, idx,  n, n_dim_idx
    real*8, parameter :: pi = acos(-1.0)
    real*8, intent(in) ::  rho
    logical, intent(in) :: window

    ! S_ab(q) = delta_ab + 4 pi rho (ca cb)^0.5
    !                  * int_0^r_cut dr r^2 [ g_ab(r) - 1 ] sin(qr)/(qr) * sin( pi r / R )/ (pi r /R)
    ! window corresponds to sin(pi r / R) / (pi r / R)
    ! rs has size n_samples_pc
    ! if (allocated(structure_factor_partial))  deallocate(structure_factor_partial)
    ! allocate( structure_factor_partial(1:params&
    !      &%structure_factor_n_samples, 1 : n_species , 1 :&
    !      & n_species) )
    structure_factor_partial = 0.d0

    ! n = q_end - q_beg + 1 !size(q_list)

    dr = rs(2) - rs(1)

    w = 1.d0

    n_dim_idx = 1
    outers: do i = 1, n_species
       ca = n_atoms_of_species(i) / dfloat( n_sites )
       do j = 1, n_species

          if (i > j) cycle

          cb = n_atoms_of_species(j) / dfloat( n_sites )
          cabh = ( ca * cb )**(0.5)

          ! Iterate through each q
          do k = q_beg, q_end
             q = q_list(k) * 2.d0 * pi


             do l = 1, n_samples_pc
                ! do the integral
                if (window) w = sinc( pi * rs(l) / r_cut )

                structure_factor_partial(k,n_dim_idx) = structure_factor_partial(k,n_dim_idx) + &
                     & dr * rs(l)**2 &
                     & * ( pair_distribution_partial(l, n_dim_idx) - 1.d0 ) &
                     & * sinc( q * rs(l) ) * w
             end do

             structure_factor_partial(k,n_dim_idx) = 4.d0 * pi * cabh * rho * structure_factor_partial(k,n_dim_idx)

             if (i == j)then
                structure_factor_partial(k,n_dim_idx)  = structure_factor_partial(k,n_dim_idx)  + 1.d0
             end if

          end do

          n_dim_idx = n_dim_idx + 1

          if ( n_dim_idx > n_dim_partial )then
             exit outers
          end if

       end do
    end do outers

  end subroutine get_partial_structure_factor



  subroutine get_xrd_explicit( n_sites0, species_types, n_species,  &
       & neighbors_list, n_neigh, neighbor_species, rjs, n_samples, q_list,&
       & y, r_min, r_max, r_cut,&
       & partial, species_1, species_2, window)
    implicit none
    real*8,  intent(in) :: rjs(:)
    real*8 :: r_min, r_max, r_cut
    character*8, allocatable :: species_types(:)
    integer, intent(in) :: neighbors_list(:), n_neigh(:), neighbor_species(:)
    integer, intent(in) :: n_sites0, n_samples, species_1, species_2, n_species
    integer :: n_sites, n_pairs, count, count_species_1
    integer :: i, j, k, i2, j2, l, ind_bin_l, ind_bin_h, species_i, species_j
    real*8,  intent(out) :: y(1:n_samples)
    real*8,  intent(in) :: q_list(1:n_samples)
    real*8 :: r, n_pc, w, wfaci, wfacj, mag
    real*8, allocatable :: sf_parameters(:,:)
    real*8, parameter :: pi = acos(-1.0)
    logical, intent(in) :: partial, window
    logical ::  species_in_list, counted_1=.false.
    ! First allocate the pair correlation function array

    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)

    allocate( sf_parameters(1:9,1:n_species) )


    sf_parameters = 0.d0
    do i = 1, n_species
       call get_scattering_factor_params(species_types(i), sf_parameters(1:9,i))
    end do

    y = 0.d0

    count = 0
    count_species_1 = 0
    k = 0
    w = 1.d0
    do i = 1, n_sites
       ! k is the index which keeps a number of the atom pairs
       ! i2 is the index of a particular atom
       i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
       species_i = neighbor_species(k+1)


       do j = 1, n_neigh(i)
          ! Loop through the neighbors of atom i
          ! j2 is the index of a neighboring atom to i
          k = k + 1
          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1

          species_j = neighbor_species(k)

          if (partial)then
             if (species_i /= species_1) then
                cycle
             elseif(j == 1)then
                count_species_1 = count_species_1 + 1
             end if

          end if

          if (partial)then
             if (species_j /= species_2) cycle
          end if

          r = rjs(k) ! atom pair distance

          if ( r > r_cut)then
             cycle
          end if

          ! edge cases where r is not in range
          if (r < r_min)then
             !            print *, "pair_distribution_function: Given r is less than r_min! Continuing loop "
             cycle
          end if

          if (r > r_max)then
             !            print *, "pair_distribution_function: Given r is more than r_max! Continuing loop "
             cycle
          end if


          if (window) w = sinc( pi * r / r_cut )

          do l = 1, n_samples
             call get_scattering_factor(wfaci, sf_parameters(1:9,species_i), q_list(l)/2.d0)
             call get_scattering_factor(wfacj, sf_parameters(1:9,species_j), q_list(l)/2.d0)

             y(l) = y(l) +  wfaci * wfacj * sinc(2.d0 * pi * q_list(l) * r) * w  
          end do


       end do
    end do
    y = y / dfloat(n_sites0)

    deallocate(sf_parameters)

  end subroutine get_xrd_explicit


  subroutine get_xrd_from_partial_structure_factors(q_beg, q_end, &
       & structure_factor_partial, n_species, species_types,&
       & species, wavelength, damping, alpha, method, use_iwasa, output, &
       & x, y, n_atoms_of_species, y_sub )
    implicit none
    real*8, intent(in) :: damping, wavelength, alpha, structure_factor_partial(:,:)
    integer, intent(in) :: species(:)
    logical, intent(in) :: use_iwasa
    integer, intent(in) :: n_species, q_beg, q_end
    character*8, allocatable :: species_types(:)
    character*32, intent(in) :: method, output
    real*8 :: prefactor, p, c, c2, rij, diff(1:3), mag, sth, wfaci, wfacj, wfac_n, ntot, f, delta
    integer :: i, j, l, n, n_dim_idx, n_dim_partial
    real*8 , intent(in), allocatable :: n_atoms_of_species(:)
    real*8, allocatable :: sf_parameters(:,:)
    real*8, intent(in) :: x(:)
    real*8, intent(out) :: y(:), y_sub(:)
    real*8, parameter :: pi = acos(-1.0)

    n_dim_partial = n_species * ( n_species + 1 ) / 2

    y = 0.d0
    y_sub = 0.d0
    n = q_end - q_beg + 1 !size(x)

    allocate( sf_parameters(1:9,1:n_species) )

    sf_parameters = 0.d0
    do i = 1, n_species
       call get_scattering_factor_params(species_types(i), sf_parameters(1:9,i))
    end do

    ntot = sum(n_atoms_of_species)

    do l = q_beg, q_end

       n_dim_idx = 1
       outer: do i = 1, n_species
          call get_scattering_factor(wfaci, sf_parameters(1:9,i), x(l)/2.d0)

          do j = 1, n_species
             call get_scattering_factor(wfacj, sf_parameters(1:9,j), x(l)/2.d0)
             !               call get_scattering_factor(species_types(j), x(l)/2.d0, wfacj)

             if ( i > j ) cycle

             if ( i == j )  f= 1.d0
             if ( i /= j )  f= 2.d0

             if ( i == j )  delta= 1.d0
             if ( i /= j )  delta= 0.d0

             y(l) = y(l) + f * ( wfaci * wfacj )  * (( n_atoms_of_species(i) / ntot ) * ( n_atoms_of_species(j) / ntot ))**0.5 &
                  & * ( structure_factor_partial(l, n_dim_idx) - delta )

             n_dim_idx = n_dim_idx + 1
             if ( n_dim_idx > n_dim_partial ) exit outer

          end do
       end do outer

       ! if ( trim(output) == "xrd" )then
       !    !          default !
          do i = 1, n_species
             call get_scattering_factor(wfaci, sf_parameters(1:9,i), x(l)/2.d0)
             y_sub(l) = ( n_atoms_of_species(i) / ntot ) * wfaci * wfaci
             y(l) = y(l) + y_sub(l)
          end do
       ! elseif ( trim(output) == "q*i(q)" .or. trim(output) == "q*F(q)")then
       !    ! we now handlw this case in preprocess experimental data
       !    ! output q * i(q) === q * F_x(q)
       !    y(l) = x(l) * y(l)
       ! elseif( trim(output) == "F(q)" .or. trim(output) == "i(q)")
       !    ! do nothing,
       !    ! Output the total scattering functon, i(q) === F_x(q)

       ! end if


    end do

    deallocate(sf_parameters)

  end subroutine get_xrd_from_partial_structure_factors


  !#############################################################!
  !###---   Experimental Interpolation and Similarities   ---###!
  !#############################################################!

  subroutine get_all_similarities( n_exp, exp_data, energy_scales, s_tot )
    implicit none
    integer, intent(in) :: n_exp
    type(exp_data_container), allocatable, intent(in) :: exp_data(:)
    real*8, allocatable :: energy_scales(:)
    integer :: i
    real*8, intent(out) :: s_tot

    s_tot = 0.d0
    do i = 1, n_exp
       if (exp_data(i)%compute_similarity)then

          write(*,'(A,1X,A,1X,F12.6,1X,F12.6)') " Exp similarity&
               & (label, escale, similarity) ", trim(exp_data(i)&
               &%label), energy_scales(i) ,  exp_data(i)%similarity

          s_tot = s_tot +  energy_scales(i) * exp_data(i)%similarity
       end if
    end do

  end subroutine get_all_similarities

  subroutine check_species_in_list( species_i, allowed_species, species_in_list )
    implicit none
    integer, allocatable, intent(in) :: allowed_species(:)
    integer, intent(in) :: species_i
    logical, intent(out) :: species_in_list
    integer :: l
    species_in_list = .false.
    do l = 1, size(allowed_species)
       if (allowed_species(l) == species_i) species_in_list = .true.
    end do
  end subroutine check_species_in_list



  !**************************************************************************
  subroutine calculate_exp_interpolation(x, y, n_samples, data)
    implicit none
    real*8, allocatable, intent(in) ::  data(:,:)
    integer, intent(in) :: n_samples
    real*8, allocatable, intent(inout) :: x(:), y(:)
    real*8 :: dx

    if (.not. allocated(x) )then
       allocate(x(1:n_samples))
       allocate(y(1:n_samples))

       x = 0.d0
       y = 0.d0

       call interpolate_data( x, y, data(1,:), data(2,:), n_samples, dx )
    else
       dx = x(2) - x(1)
    end if


  end subroutine calculate_exp_interpolation



  subroutine get_data_similarity( y, y_pred, sim_exp_pred, exp_similarity_type)
    implicit none
    real*8, allocatable, intent(in) ::  y(:), y_pred(:)
    character*32, intent(in) :: exp_similarity_type
    real*8, intent(out) :: sim_exp_pred

    if (exp_similarity_type == "squared_diff")then
       sim_exp_pred = - 0.5 * dot_product(y - y_pred, y - y_pred)
    else
       sim_exp_pred =  dot_product(y, y_pred)
    end if
  end subroutine get_data_similarity


  !########################################!
  !###---   XPS spectrum utilities   ---###!
  !########################################!


  subroutine get_compare_xps_spectra(data, core_electron_be, &
       & sigma, n_samples, mag, sim_exp_pred, &
       & x_i_exp, y_i_exp,  y_i_pred, y_i_pred_all, &
       & get_exp, exp_similarity_type )
    implicit none
    real*8, allocatable, intent(in) :: data(:,:)
    real*8, intent(in) :: sigma, core_electron_be(:)
    real*8, allocatable, intent(inout) :: x_i_exp(:), y_i_exp(:), &
         & y_i_pred(:), y_i_pred_all(:,:)
    real*8, intent(out) :: sim_exp_pred, mag
    integer, intent(in) :: n_samples
    logical, intent(in) :: get_exp
    character*32, intent(in) :: exp_similarity_type

    if (get_exp)then
       ! Get interpolated exp spectra
       call get_xps_spectra(data(1,:), data(2,:), sigma, n_samples, mag,&
            & x_i_exp,  y_i_exp, y_i_pred_all, core_electron_be,&
            &  .false.)
    end if

    ! Get interpolated predicted spectra
    call get_xps_spectra(data(1,:), data(2,:), sigma, n_samples, mag,&
         & x_i_exp, y_i_pred, y_i_pred_all, core_electron_be,&
         & .true.)

    call get_data_similarity( y_i_exp, y_i_pred, sim_exp_pred, exp_similarity_type)

  end subroutine get_compare_xps_spectra


  subroutine broaden_spectrum(x, x0, y, y_all, idx, sigma)
    implicit none
    real*8, allocatable, intent(in) :: x(:)
    real*8,  intent(inout) :: y(:), y_all(:,:)
    real*8, intent(in) :: x0, sigma
    real*8 :: norm_fac
    integer :: i,idx
    norm_fac = 1.d0 / ( sqrt(2.d0 * 3.14159265359) * sigma )
    do i = 1, size(x)
       y_all(idx, i) = exp( -( x(i) - x0 )**2 / (2.d0 * sigma**2) )
       y(i) = y(i) + y_all(idx, i) !  exp( -( x(i) - x0 )**2 / (2.d0 * sigma**2) )
    end do

  end subroutine broaden_spectrum


  subroutine get_xps_spectra(xi, yi, sigma, n_samples, mag, x, y, y_all,&
       & core_electron_be, broaden)
    ! This just gets the broadened spectra
    ! xi are the predicted core electron binding energies and x is
    ! the one
    implicit none
    real*8, intent(in) :: xi(:), yi(:), core_electron_be(:)
    integer, intent(in) :: n_samples
    real*8, allocatable, intent(out) :: x(:), y(:), y_all(:,:)
    real*8, intent(in) :: sigma
    integer :: i
    real*8 ::  x_min, x_max, x_range, t, dx
    real*8, intent(out) :: mag
    logical, intent(in) :: broaden

    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(y_all)) deallocate(y_all)

    allocate(x(1:n_samples))
    allocate(y(1:n_samples))
    allocate(y_all(1:size(core_electron_be), 1:n_samples))

    x = 0.d0
    y = 0.d0
    y_all = 0.d0

    x_min = xi(1)
    x_max = xi(size(xi))
    x_range = x_max - x_min
    dx = x_range / real(n_samples)

    do i = 1, n_samples
       t = ((i-1) / real(n_samples-1))
       x(i) = (1.d0 - t) * x_min  +  t * x_max !range
    end do

    ! Interpolate
    call lerp(x, y, xi, yi)

    if (broaden) then
       y=0.d0

       do i = 1, size(core_electron_be)
          call broaden_spectrum(x, core_electron_be(i), y, y_all, i, sigma)
       end do

    end if

    mag = sqrt(dot_product(y, y) * dx)  !sum(y * dx) !
    y = y / mag
    y_all = y_all / mag

  end subroutine get_xps_spectra



  subroutine broaden_spectrum_derivative(x, x0, x0_der, y, sigma, y_exp, y_tot, mag, dx)
    implicit none
    real*8, intent(in) :: x(:)
    real*8, intent(inout) :: y(:, :)
    real*8, intent(in) ::  y_exp(:)
    real*8, intent(in) :: x0, x0_der(:), sigma, mag, dx
    real*8 :: f
    real*8, intent(out) :: y_tot(1:3)
    integer :: i
    real*8 :: norm_fac
    norm_fac = 1.d0 / ( sqrt(2.d0 * 3.14159265359) * sigma )

    y = 0.d0

    do i = 1, size(x)
       f = ( ( x(i) - x0 ) / (sigma**2) ) * exp( -( x(i) - x0 )**2 / (2*sigma**2) ) * y_exp(i) / mag
       y(1:3,i) = y(1:3,i) + x0_der * f
    end do

    y_tot(1) =  sum(  y(1,:) ) * dx
    y_tot(2) =  sum(  y(2,:) ) * dx
    y_tot(3) =  sum(  y(3,:) ) * dx

  end subroutine broaden_spectrum_derivative


  subroutine get_xps_weights(weights, decay, positions)
    implicit none
    real*8, allocatable, intent(in) :: positions(:,:), decay
    real*8, allocatable, intent(out) :: weights(:)
    real*8 :: z_max, z
    integer :: n, i

    n = size(positions,2)
    allocate(weights(1:n))
    weights = 1.d0
    z_max = maxval(positions(3,1:n))

    if ( decay > 1e-5 )then
       do i = 1, n
          z = positions(3,i)
          weights(i) = exp( -(z_max - z) / decay )
       end do
    end if

  end subroutine get_xps_weights


  !#####################################!
  !###---   Experimental Forces   ---###!
  !#####################################!

  subroutine get_experimental_this_force(this_force, y, der_vec, y_der, norm, prefactor, n_samples, k)
    implicit none
    real*8, allocatable, intent(in) :: y(:)
    real*8, intent(in) :: norm
    real*8 :: sum_d1, sum_d2, sum_d3
    real*8, allocatable, intent(out) ::  y_der(:,:)
    real*8, allocatable, intent(in)  ::  der_vec(:,:,:), prefactor(:)
    integer, intent(in) :: n_samples, k
    integer :: l
    real*8, intent(out)  ::  this_force(1:3)


    sum_d1 = dot_product( y, der_vec(1,k,1:n_samples) )
    sum_d2 = dot_product( y, der_vec(2,k,1:n_samples) )
    sum_d3 = dot_product( y, der_vec(3,k,1:n_samples) )

    do l = 1, n_samples
       y_der(1, l) =  ( der_vec(1,k,l)  - y(l) * sum_d1) / norm
       y_der(2, l) =  ( der_vec(2,k,l)  - y(l) * sum_d2) / norm
       y_der(3, l) =  ( der_vec(3,k,l)  - y(l) * sum_d3) / norm
    end do

    this_force(1) =  dot_product( y_der(1, 1:n_samples), prefactor(1:n_samples) )
    this_force(2) =  dot_product( y_der(2, 1:n_samples), prefactor(1:n_samples) )
    this_force(3) =  dot_product( y_der(3, 1:n_samples), prefactor(1:n_samples) )

  end subroutine get_experimental_this_force


  subroutine get_experimental_forces(energy_scale, core_electron_be, core_electron_be_der,&
       & n_neigh, neighbors_list,  &
       & sigma, n_samples, norm, &
       & x, y_exp,  y,  do_forces,  xyz,&
       & energies_lp, forces0, virial, exp_similarity_type, quantity_name, rank, dx)
    implicit none
    integer, intent(in) :: n_neigh(:), neighbors_list(:)
    real*8, intent(in) :: sigma, core_electron_be(:),&
         & core_electron_be_der(:,:), xyz(:,:),&
         & energy_scale, norm
    real*8, allocatable, intent(inout) :: forces0(:,:)
    real*8, allocatable, intent(in) :: x(:), y_exp(:), y(:)
    real*8, intent(inout) :: energies_lp(:), dx
    real*8, intent(inout) :: virial(1:3,1:3)
    real*8 ::  this_force(1:3)
    real*8, allocatable ::  y_der(:,:), der_factor(:,:), der_vec(:,:,:), prefactor(:), dxa(:)
    integer, intent(in) :: n_samples, rank
    logical, intent(in) :: do_forces
    integer :: n_sites, n_pairs, n_sites0
    integer :: i, j, i2, j2, k, l,  k1, k2, mag_force_i
    logical :: vector_norm = .false.
    real*8 :: mag_force, max_mag_force, &
         sum_d1, sum_d2, sum_d3, yv
    character*32, intent(in) :: exp_similarity_type
    character*1024, intent(in) :: quantity_name


    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_sites0 = size(forces0, 2)

    if (exp_similarity_type == 'squared_diff')then

       energies_lp = + energy_scale / n_sites0 * ( dx * dot_product( (y - y_exp), (y - y_exp) ) )
       allocate(prefactor(1:n_samples))
       prefactor =  2.d0 * ( y - y_exp )
    end if


    ! Now, the full similarity is the overlap integral between the
    ! predicted, broadened spectra and experiment.

    ! The change of position of atom i, which has neighbours
    ! (j,k,l...,m) results in a change in the soap descriptors of
    ! (j,k,l,...m).
    ! Hence, the neighbors,
    ! number_of_neighbours = n_neigh(i)
    ! Index of the neighbour list for i => idx_i = sum_j=1^i( n_neigh(j) )
    ! Indexes of the neighbours => neighbor_list(idx_i: idx_i + n_neigh(i))

    ! We compute the expression
    !   Compute xps forces
    if( do_forces )then
       max_mag_force = 0.0
       mag_force_i  = 0
       forces0 = 0.d0
       virial = 0.d0

       allocate(y_der(1:3, 1:n_samples))
       allocate(der_factor(1:n_samples, 1:n_sites))

       ! This is the vector which has all the derivatives in for all the betas
       allocate(der_vec(1:3, 1:n_pairs, 1:n_samples))

       ! First get the derivative factors, to reduce multiplications

       der_factor = 0.d0

       do i = 1, n_samples
          do j = 1, n_sites
             yv = ( x(i) - core_electron_be(j) )
             der_factor(i,j) = yv  * exp( - yv**2 / (2.d0 * sigma**2)  ) / sigma**2
          end do
       end do


       !     First, we compute the forces acting on all the "SOAP-neighbors" of atom i
       !     (including i itself) due to the gradients of i's Hirshfeld volume wrt the
       !     positions of its neighbors

       ! First, have the derivatives of the normalisation factor in place

       der_vec = 0.d0


       do l = 1, n_samples
          k = 0
          do i = 1, n_sites
             i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
             do j = 1, n_neigh(i)
                k = k + 1
                j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
                if( .not. all(core_electron_be_der(1:3, k) == 0.d0) )then
                   if( exp_similarity_type == 'squared_diff' )then
                      ! (xi - cb(i)) * dxi_drka * exp( -(xi-cb(i))/2/sigma^2 )
                      der_vec(1:3, k, l) =  der_vec(1:3,k,l) + der_factor(l,i) * core_electron_be_der(1:3, k)
                      !                         if (l==1) print *, "i2=", i2, " j2=", j2, " cb(1:3,k)=", core_electron_be_der(1:3, k)
                   end if
                end if
             end do
          end do
       end do


       k = 0
       do i = 1, n_sites
          ! k is the index which keeps a number of the atom pairs
          ! i2 is the index of a particular atom
          i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
          do j = 1, n_neigh(i)
             ! Loop through the neighbors of atom i
             ! j2 is the index of a neighboring atom to i
             k = k + 1
             j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
             !             SOAP neighbors
             ! MAY NEED TO CHANGE THIS TO ACCOUNT FOR MACHINE PRECISIO

             if( .not. all(core_electron_be_der(1:3, k) == 0.d0) )then

                y_der = 0.d0

                if( exp_similarity_type == 'squared_diff' )then

                   if (vector_norm) then
                      sum_d1 = dot_product( y, der_vec(1,k,1:n_samples) )
                      sum_d2 = dot_product( y, der_vec(2,k,1:n_samples) )
                      sum_d3 = dot_product( y, der_vec(3,k,1:n_samples) )

                      do l = 1, n_samples
                         y_der(1, l) =  ( der_vec(1,k,l)  - y(l) * sum_d1) / norm
                         y_der(2, l) =  ( der_vec(2,k,l)  - y(l) * sum_d2) / norm
                         y_der(3, l) =  ( der_vec(3,k,l)  - y(l) * sum_d3) / norm
                      end do

                      this_force(1) =  dot_product( y_der(1, 1:n_samples), prefactor(1:n_samples) )
                      this_force(2) =  dot_product( y_der(2, 1:n_samples), prefactor(1:n_samples) )
                      this_force(3) =  dot_product( y_der(3, 1:n_samples), prefactor(1:n_samples) )

                   else
                      ! do sum
                      sum_d1 = dx * sum( der_vec(1,k,1:n_samples) )
                      sum_d2 = dx * sum( der_vec(2,k,1:n_samples) )
                      sum_d3 = dx * sum( der_vec(3,k,1:n_samples) )

                      do l = 1, n_samples
                         y_der(1, l) =  ( der_vec(1,k,l)  - y(l) * sum_d1 ) / norm
                         y_der(2, l) =  ( der_vec(2,k,l)  - y(l) * sum_d2 ) / norm
                         y_der(3, l) =  ( der_vec(3,k,l)  - y(l) * sum_d3 ) / norm
                      end do

                      this_force(1) =  dot_product( y_der(1, 1:n_samples), prefactor(1:n_samples) )
                      this_force(2) =  dot_product( y_der(2, 1:n_samples), prefactor(1:n_samples) )
                      this_force(3) =  dot_product( y_der(3, 1:n_samples), prefactor(1:n_samples) )

                   end if


                else
                   call broaden_spectrum_derivative(x(1:n_samples),&
                        & core_electron_be(i),&
                        & core_electron_be_der(1:3, k),&
                        & y_der(1:3,1:n_samples ), sigma, y_exp(1:n_samples), this_force(1:3), norm, dx)

                end if

                this_force(1:3) =  - energy_scale *  this_force(1:3)

                forces0(1:3, j2) = forces0(1:3, j2) + this_force(1:3)


                mag_force = norm2(forces0(1:3,j2))
                if (mag_force > max_mag_force)then
                   mag_force_i = j2
                   max_mag_force = mag_force
                end if

                !             Sign is plus because this force is acting on j2. Factor of one is because this is
                !             derived from a local energy
                !              virial = virial + dot_product(this_force(1:3), xyz(1:3,k))
                do k1 = 1, 3
                   do k2 =1, 3
                      virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
                   end do
                end do
             end if

             !           There is no net force acting on i2 from its periodic replicas...
             ! if( j2 /= i2 )then
             !    forces0(1:3,i2) = forces0(1:3,i2) + this_force(1:3)
             ! end if
             !           ... but the periodic replicas DO contribute to the virial
             !           Sign is minus because this force is acting on i2. Factor of 1/2 is because this is
             !           derived from a pair energy
             !            virial = virial - 0.5d0 * dot_product(this_force(1:3), xyz(1:3,k))
             ! do k1 = 1, 3
             !    do k2 =1, 3
             !       virial(k1, k2) = virial(k1, k2) - 0.25d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
             !    end do
             ! end do
          end do
       end do
       deallocate(y_der)
       if(allocated(prefactor)) deallocate(prefactor)
       deallocate(der_vec, der_factor)
       if(allocated(dxa)) deallocate(dxa)

    end if

  end subroutine get_experimental_forces



  ! subroutine get_xrd( n_neigh, neighbors_list, neighbor_species, rjs, xyz, &
  !      & positions, n_species, species_types, species, wavelength,&
  !      & damping, alpha, method, use_iwasa, x_min, x_max, &
  !      & n_samples, x_i_exp, y_i_pred, r_cut, forces0 )
  !   implicit none
  !   real*8, allocatable, intent(in) :: positions(:,:)
  !   integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
  !   real*8, intent(in) :: damping, wavelength, alpha, x_min, x_max, r_cut
  !   integer, intent(in) :: species(:)
  !   real*8, allocatable :: x(:), y(:), s(:), x_exp(:), y_exp(:), wfac(:), wfac_species(:), der_vec(:,:,:)
  !   real*8, allocatable, intent(inout) :: x_i_exp(:), y_i_pred(:)
  !   real*8, allocatable, intent(in) :: rjs(:), xyz(:)
  !   real*8, allocatable, intent(inout) :: forces0(:,:)
  !   logical, intent(in) :: use_iwasa
  !   integer, intent(in) :: n_samples, n_species
  !   character*8, allocatable :: species_types(:)
  !   character*32, intent(in) :: method
  !   real*8 :: prefactor, p, c, c2, rij, diff(1:3), intensity, mag, sth, wfaci, wfacj, kron
  !   integer :: i, j, k, l, n, n_sites, j2, i2, i3, j3, k3,  j4
  !   real*8 :: x_val, x_range, dx, pi=3.14159265359, ri(1:3), rj(1:3)

  !   ! Was going to use the rijs for this, but we want the /full/
  !   ! spectra, so we need the scattering contribution from all atoms



  !   n_sites = size(n_neigh)
  !   n_pairs = size(neighbors_list)
  !   n_sites0 = size(forces0, 2)

  !   dx = x(2) - x(1)

  !   ! This is essentially the same as the procedure for the hirshfeld gradients as in vdw.f90

  !   call linspace( x, x_min, x_max, n_samples, dx )

  !   ! We then broaden y to get the actual spectra
  !   allocate(y(1:n_samples))
  !   y = 0.d0
  !   n = size(x)

  !   allocate(s(1:n))
  !   allocate(wfac(1:n_sites))
  !   allocate(wfac_species(1:n_species))

  !   if (do_forces)then
  !      allocate(der_vec(1:3, 1:n_pairs, 1:n_samples))
  !      der_vec = 0.d0
  !   end if


  !   ! x is in units of 2θ for XRD and in units of q for SAXS
  !   if (trim(method) == "saxs")then
  !      do i = 1, n
  !         s(i) =  x(i) / 2.d0 / pi
  !      end do
  !   else
  !      ! Do XRD
  !      do i = 1, n
  !         s(i) = 2.d0 * sin( x(i) * pi / 180.d0 / 2.d0 ) / wavelength
  !      end do
  !   end if



  !   do l = 1, n

  !      prefactor = exp( - damping * s(l)**(2.d0) / 2.d0 )

  !      do i = 1, n_species
  !         call get_waasmaier(species_types(i), s(l), wfac_species(i))
  !      end do

  !      do i = 1, n_sites
  !         do j = 1, n_species
  !            if (species(i) == j)then
  !               wfac(i) = wfac_species(j)
  !            end if
  !         end do
  !      end do

  !      if (use_iwasa)then
  !         sth = ( wavelength * s(l) / 2.0d0 )
  !         p = 1.0d0 - sth * sth
  !         if (p < 0) p = 0.0d0
  !         c = sqrt(p)
  !         c2 = cos(2.0d0 * acos(c))
  !         prefactor = prefactor  *  c / (1.0d0 + alpha * c2**(2.d0))
  !      end if


  !      intensity = 0.d0

  !      k = 0
  !      do i = 1, n_sites
  !         i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
  !         ri = xyz(1:3, k+1)

  !         do j = 1, n_neigh(i)
  !            k = k + 1
  !            j2 = modulo(neighbors_list(k)-1, n_sites0) + 1

  !            rj = xyz(1:3, k)

  !            rij = rjs(k)

  !            if (rij > r_cut) cycle

  !            wfaci = wfac(i)
  !            ! Get the scattering factor of the neighbor species
  !            wfacj = wfac_species(neighbor_species(k))

  !            intensity = intensity + wfaci * wfacj * ( sinc( 2.d0 * s(l) * rij ) )

  !            if (do_forces)then
  !               if (rij < 1e-3) cycle

  !               ! Now we must cycle through the neighbours again to elucidate the derivative dependence
  !               do i3 = 1, n_sites
  !                  i4 = modulo(neighbors_list(k2+1)-1, n_sites0) + 1
  !                  do j3 = 1, n_neigh(i3)
  !                     k2 = k2 + 1
  !                     j4 = modulo(neighbors_list(k2)-1, n_sites0) + 1

  !                     ! we have a kronecker delta dependence
  !                     ! ( delta_jk - delta_ik )
  !                     ! Therefore, we must check the atomic indices to see what this factor is

  !                     if ( j4 /= i2 .and. j4 /= j2 ) cycle
  !                     if ( j4 == i2 .and. j4 == j2 ) cycle
  !                     if ( j4 /= i2 .and. j4 == j2 ) kron =  1.d0
  !                     if ( j4 == i2 .and. j4 /= j2 ) kron = -1.d0

  !                     ! the initial atom has index i2
  !                     ! the neighbour has index j2
  !                     ! j3 is the index of another atom
  !                     ! k is the index of another atom which is a neighbour
  !                     der_vec(1:3, k2, l) =  der_vec(1:3,k2,l) + &
  !                          & (1.d0 / rij) * ( 2.d0 * s(l) ) * wfaci * wfacj * kron * ( rj - ri ) * &
  !                          &   ( cosc(2.d0 * s(l) * rij) - sinc( 2.d0 * s(l) * rij ) / (2.d0 * s(l) * rij) )


  !                  end do
  !               end do

  !            end if



  !         end do
  !      end do

  !      y(l) = prefactor * intensity
  !      der_vec(1:3, 1:n_pairs, l) = prefactor * der_vec(1:3, 1:n_pairs, l)


  !   end do

  !   ! We wont actually calculate the forces here, we need to collect all the components together to obtain the forces
  !   if( do_forces )then
  !      forces0 = 0.d0
  !      virial = 0.d0

  !      k = 0
  !      do i = 1, n_sites
  !         ! k is the index which keeps a number of the atom pairs
  !         ! i2 is the index of a particular atom
  !         i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
  !         do j = 1, n_neigh(i)
  !            ! Loop through the neighbors of atom i
  !            ! j2 is the index of a neighboring atom to i
  !            k = k + 1
  !            j2 = modulo(neighbors_list(k)-1, n_sites0) + 1


  !            sum_d1 = dot_product( y, der_vec(1,k,1:n_samples) )
  !            sum_d2 = dot_product( y, der_vec(2,k,1:n_samples) )
  !            sum_d3 = dot_product( y, der_vec(3,k,1:n_samples) )

  !            do l = 1, n_samples
  !               y_der(1, l) =  ( der_vec(1,k,l)  - y(l) * sum_d1) / norm
  !               y_der(2, l) =  ( der_vec(2,k,l)  - y(l) * sum_d2) / norm
  !               y_der(3, l) =  ( der_vec(3,k,l)  - y(l) * sum_d3) / norm
  !            end do

  !            this_force(1) =  dot_product( y_der(1, 1:n_samples), prefactor(1:n_samples) )
  !            this_force(2) =  dot_product( y_der(2, 1:n_samples), prefactor(1:n_samples) )
  !            this_force(3) =  dot_product( y_der(3, 1:n_samples), prefactor(1:n_samples) )

  !               this_force(1:3) =  - energy_scale *  this_force(1:3)

  !               forces0(1:3, j2) = forces0(1:3, j2) + this_force(1:3)

  !               !             Sign is plus because this force is acting on j2. Factor of one is because this is
  !               !             derived from a local energy
  !               !              virial = virial + dot_product(this_force(1:3), xyz(1:3,k))
  !               do k1 = 1, 3
  !                  do k2 =1, 3
  !                     virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
  !                  end do
  !               end do
  !            end if

  !            !           There is no net force acting on i2 from its periodic replicas...
  !            ! if( j2 /= i2 )then
  !            !    forces0(1:3,i2) = forces0(1:3,i2) + this_force(1:3)
  !            ! end if
  !            !           ... but the periodic replicas DO contribute to the virial
  !            !           Sign is minus because this force is acting on i2. Factor of 1/2 is because this is
  !            !           derived from a pair energy
  !            !            virial = virial - 0.5d0 * dot_product(this_force(1:3), xyz(1:3,k))
  !            ! do k1 = 1, 3
  !            !    do k2 =1, 3
  !            !       virial(k1, k2) = virial(k1, k2) - 0.25d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
  !            !    end do
  !            ! end do
  !         end do
  !      end do
  !      deallocate(y_der)
  !      deallocate(der_vec)

  !   end if


  !   ! get the magnitude after
  !   ! mag = sqrt(dot_product(y,y))
  !   ! y = y / mag
  !   x_i_exp = x
  !   y_i_pred = y


  !   deallocate(x)
  !   deallocate(y)
  !   deallocate(s)
  !   deallocate(wfac_species)
  !   deallocate(wfac)

  ! end subroutine get_xrd


  !**************************************************************************
  subroutine get_xrd_single_process( positions, n_species,  species, wavelength, damping, alpha, &
       & method, use_iwasa, x_min, x_max,  n_samples, x_i_exp, y_i_pred )
    implicit none
    real*8, allocatable, intent(in) :: positions(:,:)
    real*8, intent(in) :: damping, wavelength, alpha, x_min, x_max
    integer, intent(in) :: species(:)
    real*8, allocatable :: x(:), y(:), s(:),  wfac(:), wfac_species(:)
    real*8, allocatable, intent(inout) :: x_i_exp(:), y_i_pred(:)
    logical, intent(in) :: use_iwasa
    integer, intent(in) :: n_samples, n_species
    character*32, intent(in) :: method
    real*8 :: prefactor, p, c, c2, rij, diff(1:3), intensity, mag, sth, wfac_n
    integer :: i, j, l, n, n_sites
    real*8 ::  dx, pi=3.14159265359

    ! Was going to use the rijs for this, but we want the /full/
    ! spectra, so we need the scattering contribution from all atoms

    n_sites = size(positions, 2)

    ! This is essentially the same as the procedure for the hirshfeld gradients as in vdw.f90

    call linspace( x, x_min, x_max, n_samples, dx )

    ! We then broaden y to get the actual spectra
    allocate(y(1:n_samples))
    y = 0.d0
    n = size(x)

    allocate(s(1:n))
    allocate(wfac(1:n_sites))
    allocate(wfac_species(1:n_species))


    ! x is in units of 2θ for XRD and in units of q for SAXS
    if (trim(method) == "saxs")then
       do i = 1, n
          s(i) =  x(i) / 2.d0 / pi
       end do
    else
       ! Do XRD
       do i = 1, n
          s(i) = 2.d0 * sin( x(i) * pi / 180.d0 / 2.d0 ) / wavelength
       end do
    end if

    do l = 1, n

       prefactor = exp( - damping * s(l)**(2.d0) / 2.d0 )

       do i = 1, n_species
          !            call get_scattering_factor(species_types(i), s(l)/2.d0, wfac_species(i))
       end do

       do i = 1, n_sites
          do j = 1, n_species
             if (species(i) == j)then
                wfac(i) = wfac_species(j)
             end if
          end do
       end do

       if (use_iwasa)then
          sth = ( wavelength * s(l) / 2.0d0 )
          p = 1.0d0 - sth * sth
          if (p < 0) p = 0.0d0
          c = sqrt(p)
          c2 = cos(2.0d0 * acos(c))
          prefactor = prefactor  *  c / (1.0d0 + alpha * c2**(2.d0))
       end if

       ! Accounting for rij == 0.0 with wfac_n
       wfac_n = 0.d0

       intensity = 0.d0
       do i = 1, n_sites
          wfac_n = wfac_n + wfac(i) * wfac(j)
          do j = i+1, n_sites
             !      if (i /= j)then
             diff(1:3) = positions(1:3,i) - positions(1:3,j)
             rij = sqrt( dot_product(diff, diff))
             ! Now should
             intensity = intensity + wfac(i) * wfac(j) * ( sinc( 2.d0 * s(l) * rij ) )
             !     end if
          end do
       end do

       y(l) = prefactor * ( wfac_n + 2.d0 * intensity )

    end do

    mag = sqrt(dot_product(y,y))
    y = y / mag

    x_i_exp = x
    y_i_pred = y


    deallocate(x)
    deallocate(y)
    deallocate(s)
    deallocate(wfac_species)
    deallocate(wfac)

  end subroutine get_xrd_single_process

  function sinc(x) !sinc function as used in DSP
    implicit none
    real*8 :: sinc
    real*8 :: x
    real*8, parameter :: pi = acos(-1.0)
    sinc = 1.0
    if (x /= 0.0) sinc = sin(x)/x
  end function sinc

  function cosc(x)
    implicit none
    real*8 :: cosc
    real*8 :: x
    real*8, parameter :: pi = acos(-1.0)
    cosc = 1.0
    if (x /= 0.0) cosc = cos(x)/x
  end function cosc


  function sincp(x) !sinc function as used in DSP
    implicit none
    real*8 :: sincp
    real*8 :: x
    real*8, parameter :: pi = acos(-1.0)
    x = x*pi
    sincp = 1.0
    if (x /= 0.0) sincp = sin(x)/x
  end function sincp

  function coscp(x)
    implicit none
    real*8 :: coscp
    real*8 :: x
    real*8, parameter :: pi = acos(-1.0)
    x = x*pi
    coscp = 1.0
    if (x /= 0.0) coscp = cos(x)/x
  end function coscp



  !************************************
!!!   XPS related Functions Below !!!
  !************************************

  subroutine get_moments_of_distribution(x, y, dx, moments, n_moments, mean_reference)
    implicit none
    real*8, allocatable, intent(in) ::  x(:), y(:)
    real*8, allocatable, intent(out) :: moments(:)
    real*8, allocatable :: xp(:)
    real*8, intent(in) :: dx, mean_reference
    integer, intent(in) :: n_moments
    integer :: i
    ! Moments are defined at mu^p = \int dx x^p f(x)

    allocate(moments(1:n_moments))
    allocate(xp(1:size(x)))
    moments = 0.d0

    ! Here, the first moment is the mean,
    ! was going to make the other moments central, but for now, just make them normal !

    ! could actually have the reference be the mean of the
    ! exp distribution which would simplify the mathematics

    xp = 1.d0
    do i = 1, n_moments

       if (i == 1)then
          xp = xp * x
       else
          xp = (x - mean_reference)**( real(i) )
       end if

       moments(i) = (dot_product( dx * xp, y )) !**(1.d0 / real(i))
       print *, moments(i)
    end do
  end subroutine get_moments_of_distribution





  subroutine get_exp_pred_spectra_energies_forces(energy_scale, core_electron_be, core_electron_be_der,&
       & n_neigh, neighbors_list, &
       & sigma, n_samples, norm, &
       & x, y_exp,  y, y_all, do_forces,  xyz,&
       &  energies_lp, forces0, virial, exp_similarity_type, rank)
    implicit none
    integer, intent(in) :: n_neigh(:), neighbors_list(:)
    real*8, intent(in) :: sigma, core_electron_be(:),&
         & core_electron_be_der(:,:), xyz(:,:),&
         & energy_scale, norm
    real*8, allocatable, intent(inout) :: forces0(:,:)
    real*8, allocatable, intent(in) :: x(:), y_exp(:), y(:)
    real*8, intent(in) :: y_all(:,:)
    real*8, intent(inout) :: energies_lp(:)
    real*8, intent(inout) :: virial(1:3,1:3)
    real*8 ::  this_force(1:3)
    real*8, allocatable ::  y_der(:,:), der_factor(:,:), der_vec(:,:,:), prefactor(:), dxa(:)
    integer, intent(in) :: n_samples, rank
    logical, intent(in) ::  do_forces
    integer :: n_sites, n_pairs, n_sites0
    integer :: i, j, i2, j2, k, l, k1, k2, mag_force_i
    logical :: vector_norm = .false.
    real*8 :: mag_force, max_mag_force,  dx, &
         sum_d1, sum_d2, sum_d3, yv
    character*32, intent(in) :: exp_similarity_type


    n_sites = size(n_neigh)
    n_pairs = size(neighbors_list)
    n_sites0 = size(forces0, 2)

    dx = x(2) - x(1)
    ! This is essentially the same as the procedure for the hirshfeld gradients as in vdw.f90


    ! The energetic contribution to this would then be, e_scale * ( 1 - \sum_i S_i ), this can be done after mpi_reduce
    ! S_i =  dot_product( y_all(:,i), y_exp )

    if( exp_similarity_type == 'similarity' .or. exp_similarity_type == 'overlap' )then
       do i = 1, n_sites
          energies_lp(i) = - energy_scale * ( dot_product( y_all(i,:), y_exp ))
       end do
    else if (exp_similarity_type == 'squared_diff')then
       do i = 1, n_sites
          energies_lp(i) = - energy_scale * ( 1.d0 / dfloat( n_sites0 ) -  dot_product( y_all(i,:), y_exp ))
       end do


       !energies_lp = + energy_scale / n_sites0 * ( dot_product( (y - y_exp), (y - y_exp) ) )
       ! Get the other terms for the squared_diff expression


       allocate(prefactor(1:n_samples))
       prefactor =  2.d0 * ( y - y_exp )
    end if


    ! Now, the full similarity is the overlap integral between the
    ! predicted, broadened spectra and experiment.

    ! The change of position of atom i, which has neighbours
    ! (j,k,l...,m) results in a change in the soap descriptors of
    ! (j,k,l,...m).
    ! Hence, the neighbors,
    ! number_of_neighbours = n_neigh(i)
    ! Index of the neighbour list for i => idx_i = sum_j=1^i( n_neigh(j) )
    ! Indexes of the neighbours => neighbor_list(idx_i: idx_i + n_neigh(i))

    ! We compute the expression



    !   Compute xps forces
    if( do_forces )then
       max_mag_force = 0.0
       mag_force_i  = 0
       forces0 = 0.d0
       virial = 0.d0

       allocate(y_der(1:3, 1:n_samples))
       allocate(der_factor(1:n_samples, 1:n_sites))

       ! This is the vector which has all the derivatives in for all the betas
       allocate(der_vec(1:3, 1:n_pairs, 1:n_samples))

       ! First get the derivative factors, to reduce multiplications

       der_factor = 0.d0

       do i = 1, n_samples
          do j = 1, n_sites
             yv = ( x(i) - core_electron_be(j) )
             der_factor(i,j) = yv  * exp( - yv**2 / (2.d0 * sigma**2)  ) / sigma**2
          end do
       end do


       !     First, we compute the forces acting on all the "SOAP-neighbors" of atom i
       !     (including i itself) due to the gradients of i's Hirshfeld volume wrt the
       !     positions of its neighbors

       ! First, have the derivatives of the normalisation factor in place

       der_vec = 0.d0


       do l = 1, n_samples
          k = 0
          do i = 1, n_sites
             i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
             do j = 1, n_neigh(i)
                k = k + 1
                j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
                if( .not. all(core_electron_be_der(1:3, k) == 0.d0) )then
                   if( exp_similarity_type == 'squared_diff' )then
                      ! (xi - cb(i)) * dxi_drka * exp( -(xi-cb(i))/2/sigma^2 )
                      der_vec(1:3, k, l) =  der_vec(1:3,k,l) + der_factor(l,i) * core_electron_be_der(1:3, k)
                      !                         if (l==1) print *, "i2=", i2, " j2=", j2, " cb(1:3,k)=", core_electron_be_der(1:3, k)
                   end if
                end if
             end do
          end do
       end do


       k = 0
       do i = 1, n_sites
          ! k is the index which keeps a number of the atom pairs
          ! i2 is the index of a particular atom
          i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
          do j = 1, n_neigh(i)
             ! Loop through the neighbors of atom i
             ! j2 is the index of a neighboring atom to i
             k = k + 1
             j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
             !             SOAP neighbors
             ! MAY NEED TO CHANGE THIS TO ACCOUNT FOR MACHINE PRECISIO

             if( .not. all(core_electron_be_der(1:3, k) == 0.d0) )then

                y_der = 0.d0

                if( exp_similarity_type == 'squared_diff' )then

                   if (vector_norm) then
                      sum_d1 = dot_product( y, der_vec(1,k,1:n_samples) )
                      sum_d2 = dot_product( y, der_vec(2,k,1:n_samples) )
                      sum_d3 = dot_product( y, der_vec(3,k,1:n_samples) )

                      do l = 1, n_samples
                         y_der(1, l) =  ( der_vec(1,k,l)  - y(l) * sum_d1) / norm
                         y_der(2, l) =  ( der_vec(2,k,l)  - y(l) * sum_d2) / norm
                         y_der(3, l) =  ( der_vec(3,k,l)  - y(l) * sum_d3) / norm
                      end do

                      this_force(1) =  dot_product( y_der(1, 1:n_samples), prefactor(1:n_samples) )
                      this_force(2) =  dot_product( y_der(2, 1:n_samples), prefactor(1:n_samples) )
                      this_force(3) =  dot_product( y_der(3, 1:n_samples), prefactor(1:n_samples) )

                   else
                      ! do sum
                      sum_d1 =  sum( der_vec(1,k,1:n_samples) )
                      sum_d2 =  sum( der_vec(2,k,1:n_samples) )
                      sum_d3 =  sum( der_vec(3,k,1:n_samples) )

                      do l = 1, n_samples
                         y_der(1, l) =  ( der_vec(1,k,l)  - y(l) * sum_d1 ) / norm
                         y_der(2, l) =  ( der_vec(2,k,l)  - y(l) * sum_d2 ) / norm
                         y_der(3, l) =  ( der_vec(3,k,l)  - y(l) * sum_d3 ) / norm
                      end do

                      this_force(1) =  dot_product( y_der(1, 1:n_samples), prefactor(1:n_samples) )
                      this_force(2) =  dot_product( y_der(2, 1:n_samples), prefactor(1:n_samples) )
                      this_force(3) =  dot_product( y_der(3, 1:n_samples), prefactor(1:n_samples) )

                   end if


                else
                   call broaden_spectrum_derivative(x(1:n_samples),&
                        & core_electron_be(i),&
                        & core_electron_be_der(1:3, k),&
                        & y_der(1:3,1:n_samples ), sigma, y_exp(1:n_samples), this_force(1:3), norm, dx)

                end if

                this_force(1:3) =  - energy_scale *  this_force(1:3)

                forces0(1:3, j2) = forces0(1:3, j2) + this_force(1:3)


                mag_force = norm2(forces0(1:3,j2))
                if (mag_force > max_mag_force)then
                   mag_force_i = j2
                   max_mag_force = mag_force
                end if

                !             Sign is plus because this force is acting on j2. Factor of one is because this is
                !             derived from a local energy
                !              virial = virial + dot_product(this_force(1:3), xyz(1:3,k))
                do k1 = 1, 3
                   do k2 =1, 3
                      virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
                   end do
                end do
             end if

             !           There is no net force acting on i2 from its periodic replicas...
             ! if( j2 /= i2 )then
             !    forces0(1:3,i2) = forces0(1:3,i2) + this_force(1:3)
             ! end if
             !           ... but the periodic replicas DO contribute to the virial
             !           Sign is minus because this force is acting on i2. Factor of 1/2 is because this is
             !           derived from a pair energy
             !            virial = virial - 0.5d0 * dot_product(this_force(1:3), xyz(1:3,k))
             ! do k1 = 1, 3
             !    do k2 =1, 3
             !       virial(k1, k2) = virial(k1, k2) - 0.25d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
             !    end do
             ! end do
          end do
       end do
       deallocate(y_der)
       if(allocated(prefactor)) deallocate(prefactor)
       deallocate(der_vec, der_factor)
       if(allocated(dxa)) deallocate(dxa)

    end if

  end subroutine get_exp_pred_spectra_energies_forces



  subroutine lerp(x, y, xi, yi)
    implicit none
    real*8, intent(inout) :: x(:), y(:)
    real*8, intent(in) :: xi(:), yi(:)
    real*8 :: x_min, x_max, x_range, t
    integer :: i, idx

    ! Here we will make x in the same range as xi
    x_min = xi(1)        !minval(xi)
    x_max = xi(size(xi)) !maxval(xi)
    x_range = x_max - x_min

    y = 0.d0

    do i = 1, size(x)-1
       idx = int( (((x(i)-x_min) / x_range) * dfloat(size(xi))) + 1 )

       !         print *, x(i), x(idx), x(idx+1)

       do while (x(i) < xi(idx))
          idx = idx-1
       end do

       do while (x(i) > xi(idx+1))
          idx = idx+1
       end do

       if (idx >= size(xi)) idx = size(xi)-1
       if (x(i) < x_min)then
          write(*,*) 'WARNING: Lerp value', xi(i), " is less than minimum of", x_min
          write(*,*) " Setting xi(i) to ", x_min
          x(i) = x_min
          idx=1
       else if (x(i) > x_max)then
          write(*,*) 'WARNING: Lerp value', xi(i), " is more than maximum of", x_max
          write(*,*) " Setting xi(i) to ", x_max
          x(i) = x_max
          idx = size(x)-1
       end if

       t = ( (x(i) - xi(idx)) / ( xi(idx+1) - xi(idx) )  )

       y(i) = (1.0 - t) * yi(idx) + t * yi(idx+1)

    end do

    x(size(x)) = x_max
    y(size(x)) = yi(size(yi))

  end subroutine lerp



  subroutine linspace( x, x_min, x_max, n_samples, dx )
    implicit none
    real*8, allocatable, intent(inout) :: x(:)
    real*8, intent(in) :: x_min, x_max
    real*8, intent(out) :: dx
    integer, intent(in) :: n_samples
    real*8 :: t, x_range
    integer :: i

    if(allocated(x)) deallocate(x)

    allocate(x(1:n_samples))

    x_range = x_max - x_min

    do i = 1, n_samples
       t = (dfloat(i-1) / dfloat(n_samples-1))
       x(i) = (1.d0 - t) * x_min  +  t * x_max !range
    end do

    dx = x(2) - x(1)

  end subroutine linspace



  subroutine interpolate_data( x, y, xi, yi, n_samples, dx )
    implicit none
    real*8,  intent(in) :: xi(:), yi(:)
    real*8, allocatable, intent(inout) :: x(:), y(:)
    real*8, intent(out) :: dx
    integer, intent(in) :: n_samples
    real*8 :: t, x_min, x_max, x_range
    integer :: i

    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)

    allocate(x(1:n_samples))
    allocate(y(1:n_samples))

    x_min = xi(1)
    x_max = xi(size(xi))
    x_range = x_max - x_min
    dx = x_range / real(n_samples)

    do i = 1, n_samples
       t = (real(i-1) / real(n_samples-1))
       x(i) = (1.d0 - t) * x_min  +  t * x_max !range
    end do
    ! Interpolate, this gets more x from xi
    call lerp(x, y, xi, yi)
  end subroutine interpolate_data


  subroutine get_waasmaier(element, s, f)
    implicit none
    real*8, intent(in) :: s ! scattering vector: s = q / 4pi [1/A]
    real*8 :: w(1:11,1:49)
    character*8 :: elements(1:49)
    !   Input variables
    character*8, intent(in) :: element
    !   Output variables
    logical :: is_in_database = .false.
    !   Internal variables
    real*8 ::  s_sq
    real*8, intent(out) :: f
    integer :: i, j




    ! D. Waasmaier and A. Kirfel, Acta Cryst. (1995). A51, 416-431

    ! |        a1 |        b1 |        a2 |         b2 |        a3 |         b3 |        a4 |         b4 |        a5 |         b5 |    C  ! | Atom |
    w(1:11,1:49) = reshape( &
         & (/   0.732354, 11.553918,  0.753896,   4.595831,  0.283819,   1.546299,  0.190003,&
         & 26.463964,  0.039139,   0.377523,   0.000487, &  ! He
         &      0.974637,  4.334946,  0.158472,   0.342451,  0.811855,  97.102969,  0.262416,&
         & 201.363824,  0.790108,   1.409234,   0.002542, &  ! Li
         &      1.533712, 42.662078,  0.638283,   0.595420,  0.601052,  99.106501,  0.106139,&
         & 0.151340,  1.118414,   1.843093,   0.002511, &  ! Be
         &      2.085185, 23.494069,  1.064580,   1.137894,  1.062788,  61.238975,  0.140515,&
         & 0.114886,  0.641784,   0.399036,   0.003823, &  ! B
         &      2.657506, 14.780758,  1.078079,   0.776775,  1.490909,  42.086843, -4.241070,&
         & -0.000294,  0.713791,   0.239535,   4.297983, &  ! C
         &     11.893780,  0.000158,  3.277479,  10.232723,  1.858092,  30.344690,  0.858927,&
         & 0.656065,  0.912985,   0.217287, -11.804902, &  ! N
         &      2.960427, 14.182259,  2.508818,   5.936858,  0.637853,   0.112726,  0.722838,&
         & 34.958481,  1.142756,   0.390240,   0.027014, &  ! O
         &      3.511943, 10.687859,  2.772244,   4.380466,  0.678385,   0.093982,  0.915159,&
         & 27.255203,  1.089261,   0.313066,   0.032557, &  ! F
         &      4.183749,  8.175457,  2.905726,   3.252536,  0.520513,   0.063295,  1.135641,&
         & 21.813909,  1.228065,   0.224952,   0.025576, &  ! Ne
         &      4.910127,  3.281434,  3.081783,   9.119178,  1.262067,   0.102763,  1.098938,&
         & 132.013942,  0.560991,   0.405878,   0.079712, &  ! Na
         &      4.708971,  4.875207,  1.194814, 108.506079,  1.558157,   0.111516,  1.170413,&
         & 48.292407,  3.239403,   1.928171,   0.126842, &  ! Mg
         &      4.730796,  3.628931,  2.313951,  43.051166,  1.541980,   0.095960,  1.117564,&
         & 108.932389,  3.154754,   1.555918,   0.139509, &  ! Al
         &      5.275329,  2.631338,  3.191038,  33.730728,  1.511514,   0.081119,  1.356849,&
         & 86.288640,  2.519114,   1.170087,   0.145073, &  ! Si
         &      1.950541,  0.908139,  4.146930,  27.044953,  1.494560,   0.071280,  1.522042,&
         & 67.520190,  5.729711,   1.981173,   0.155233, &  ! P
         &      6.372157,  1.514347,  5.154568,  22.092528,  1.473732,   0.061373,  1.635073,&
         & 55.445176,  1.209372,   0.646925,   0.154722, &  ! S
         &      1.446071,  0.052357,  6.870609,   1.193165,  6.151801,  18.343416,  1.750347,&
         & 46.398394,  0.634168,   0.401005,   0.146773, &  ! Cl
         &      7.188004,  0.956221,  6.638454,  15.339877,  0.454180,  15.339862,  1.929593,&
         & 39.043824,  1.523654,   0.062409,   0.265954, &  ! Ar
         &      8.163991, 12.816323,  7.146945,   0.808945,  1.070140, 210.327009,  0.877316,&
         & 39.597651,  1.486434,   0.052821,   0.253614, &  ! K
         &      8.593655, 10.460644,  1.477324,   0.041891,  1.436254,  81.390382,  1.182839,&
         & 169.847839,  7.113258,   0.688098,   0.196255, &  ! Ca
         &      1.476566, 53.131022,  1.487278,   0.035325,  1.600187, 137.319495,  9.177463,&
         & 9.098031,  7.099750,   0.602102,   0.157765, &  ! Sc
         &      9.818524,  8.001879,  1.522646,   0.029763,  1.703101,  39.885423,  1.768774,&
         & 120.158000,  7.082555,   0.532405,   0.102473, &  ! Ti
         &     10.473575,  7.081940,  1.547881,   0.026040,  1.986381,  31.909672,  1.865616,&
         & 108.022844,  7.056250,   0.474882,   0.067744, &  ! V
         &     11.007069,  6.366281,  1.555477,   0.023987,  2.985293,  23.244838,  1.347855,&
         & 105.774500,  7.034779,   0.429369,   0.065510, &  ! Cr
         &     11.709542,  5.597120,  1.733414,   0.017800,  2.673141,  21.788419,  2.023368,&
         & 89.517915,  7.003180,   0.383054,  -0.147293, &  ! Mn
         &     12.311098,  5.009415,  1.876623,   0.014461,  3.066177,  18.743041,  2.070451,&
         & 82.767874,  6.975185,   0.346506,  -0.304931, &  ! Fe
         &     12.914510,  4.507138,  2.481908,   0.009126,  3.466894,  16.438130,  2.106351,&
         & 76.987317,  6.960892,   0.314418,  -0.936572, &  ! Co
         &     13.521865,  4.077277,  6.947285,   0.286763,  3.866028,  14.622634,  2.135900,&
         & 71.966078,  4.284731,   0.004437,  -2.762697, &  ! Ni
         &     14.014192,  3.738280,  4.784577,   0.003744,  5.056806,  13.034982,  1.457971,&
         & 72.554793,  6.932996,   0.265666,  -3.254477, &  ! Cu
         &     14.741002,  3.388232,  6.907748,   0.243315,  4.642337,  11.903689,  2.191766,&
         & 63.312130, 38.424042,   0.000397,  36.915828, &  ! Zm
         &     15.758946,  3.121754,  6.841123,   0.226057,  4.121016,  12.482196,  2.714681,&
         & 66.203622,  2.395246,   0.007238,  -0.847395, &  ! Ga
         &     16.540614,  2.866618,  1.567900,   0.012198,  3.727829,  13.432163,  3.345098,&
         & 58.866046,  6.785079,   0.210974,   0.018726, &  ! Ge
         &     17.025643,  2.597739,  4.503441,   0.003012,  3.715904,  14.272119,  3.937200,&
         & 50.437997,  6.790175,   0.193015,  -2.984117, &  ! As
         &     17.354071,  2.349787,  4.653248,   0.002550,  4.259489,  15.579460,  4.136455,&
         & 45.181201,  6.749163,   0.177432,  -3.160982, &  ! Se
         &     17.550570,  2.119226,  5.411882,  16.557185,  3.937180,   0.002481,  3.880645,&
         & 42.164009,  6.707793,   0.162121,  -2.492088, &  ! Br
         &     17.655279,  1.908231,  6.848105,  16.606235,  4.171004,   0.001598,  3.446760,&
         & 39.917471,  6.685200,   0.146896,  -2.810592, &  ! Kr
         &      8.123134, 15.142385,  2.138042,  33.542666,  6.761702,   0.129372,  1.156051,&
         & 224.132506, 17.679547,   1.713368,   1.139548, &  ! Rb
         &     17.730219,  1.563060,  9.795867,  14.310868,  6.099763,   0.120574,  2.620025,&
         & 135.771318,  0.600053,   0.120574,   1.140251, &  ! Sr
         &     17.792040,  1.429691, 10.253252,  13.132816,  5.714949,   0.112173,  3.170516,&
         & 108.197029,  0.918251,   0.112173,   1.131787, &  ! Y
         &     17.859771,  1.310692, 10.911038,  12.319285,  5.821115,   0.104353,  3.512513,&
         & 91.777544,  0.746965,   0.104353,   1.124859, &  ! Zr
         &     17.958398,  1.211590, 12.063054,  12.246687,  5.007015,   0.098615,  3.287667,&
         & 75.011944,  1.531019,   0.098615,   1.123452, &  ! Nb
         &      6.236218,  0.090780, 17.987711,   1.108310, 12.973127,  11.468720,  3.451426,&
         & 66.684153,  0.210899,   0.090780,   1.108770, &  ! Mo
         &     17.840964,  1.005729,  3.428236,  41.901383,  1.373012, 119.320541, 12.947364,&
         & 9.781542,  6.335469,   0.083391,   1.074784, &  ! Tc
         &      6.271624,  0.077040, 17.906739,   0.928222, 14.123269,   9.555345,  3.746008,&
         & 35.860678,  0.908235, 123.552247,   1.043992, &  ! Ru
         &      6.216648,  0.070789, 17.919738,   0.856121,  3.854252,  33.889484,  0.840326,&
         & 121.686688, 15.173498,   9.029517,   0.995452, &  ! Rh
         &      6.121511,  0.062549,  4.784063,   0.784031, 16.631683,   8.751391,  4.318258,&
         & 34.489983, 13.246773,   0.784031,   0.883099, &  ! Pd
         &      6.073874,  0.055333, 17.155437,   7.896512,  4.173344,  28.443739,  0.852238,&
         & 110.376108, 17.988685,   0.716809,   0.756603, &  ! Ag
         &      6.080986,  0.048990, 18.019468,   7.273646,  4.018197,  29.119283,  1.303510,&
         & 95.831208, 17.974669,   0.661231,   0.603504, &  ! Cd
         &      6.196477,  0.042072, 18.816183,   6.695665,  4.050479,  31.009791,  1.638929,&
         & 103.284350, 17.962912,   0.610714,   0.333097, &  ! In
         &     19.325171,  6.118104,  6.281571,   0.036915,  4.498866,  32.529047,  1.856934,&
         & 95.037182, 17.917318,   0.565651,   0.119024 &  ! Sn
         & /), [11,49] )



    ! !       a1         b1         a2         b2         a3         b3          a4         b4         a5         b5         c
    !  w(1:11,1:12) = reshape( &
    !       & (/ 2.657506, 14.780758,  1.078079,  0.776775,  1.490909, 42.086843,&
    !       & -4.241070,  -0.000294,  0.713791, 0.239535,  4.297983, &  !   C(1:10)
    !       & 11.893780,  0.000158,  3.277479, 10.232723,  1.858092, 30.344690, &
    !       & 0.858927,   0.656065,  0.912985, 0.217287, -11.804902, &  !   N(1:10)
    !       &  2.960427, 14.182259,  2.5088111, 5.936858,  0.637053,  0.112726, &
    !       & 0.722838,  34.958481,  1.142756, 0.390240,   0.027014, & !   O(1:10)
    !       &  1.950541,  0.908139,  4.146930, 27.044953,  1.494560,  0.071280, &
    !       & 1.522042,  67.520190,  5.729711, 1.981173,   0.155233, & !   P(1:10)
    !       &  6.372157,  1.514347,  5.154568, 22.092528,  1.473732,  0.061373, &
    !       & 1.635073,  55.445176,  1.209372, 0.646925,   0.154722, & !   S(1:10)
    !       &  1.446071,  0.052357,  6.870609,  1.193165,  6.151801, 18.343416, &
    !       & 1.750347,  46.398394,  0.634168, 0.401005,   0.146773, & !  Cl(1:10)
    !       & 13.521865,  4.077277,  6.947285,  0.286763,  3.866028, 14.622634, &
    !       & 2.135900,  71.966078,  4.284731, 0.004437,  -2.762697, & !  Ni(1:10)
    !       & 14.014192,  3.738280,  4.784577,  0.003744,  5.056806, 13.034982, &
    !       & 1.457971,  72.554793,  6.932996, 0.265666,  -3.774477, & !  Cu(1:10)
    !       &  6.121511,  0.062549,  4.784063,  0.784031, 16.631683,  8.751391, &
    !       & 4.318258,  34.489983, 13.246773, 0.784031,   0.883099, & !  Pd(1:10)
    !       &  6.073874,  0.055333, 17.155437,  7.896512,  4.173344, 28.443739, &
    !       & 0.852238, 110.376108, 17.988685, 0.716809,   0.756603, & !  Ag(1:10)
    !       & 31.273891,  1.316992, 18.445441,  8.797154, 17.063745,  0.124741, &
    !       & 5.555933,  40.177994,  1.575270, 1.316997,   4.050394, & !  Pt(1:10)
    !       & 16.777389,  0.122737, 19.317156,  8.621570, 32.979682,  1.256902, &
    !       & 5.595453,  38.008821, 10.576854, 0.000601,  -6.279078 /),  [11,12] )  !  Au(1:10)


    !  elements = [ ' C', ' N', ' O', ' P', ' S', 'Cl', 'Na', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']

    elements = [ 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne',&
         & 'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', 'Ca',&
         & 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',&
         & 'Zm', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',&
         & ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',&
         & 'Cd', 'In', 'Sn']


    outer: do i = 1, size(elements, 1)
       if( trim(adjustl(elements(i))) == trim(adjustl(element)) )then
          is_in_database = .true.
          f = w(11, i)
          s_sq = s * s

          do j = 1, 5
             f = f + w((1 + 2*(j-1)), i) * exp( - w( (2 + 2*(j-1)), i ) * s_sq  )
          end do
          exit outer
       end if
    end do outer

    if ( .not. is_in_database )then
       write(*,'(A,1X,A,1X,A)') '!!! XRD Warning !!! The element ', element, ' is not in the database of Waasmaier '
    end if

  end subroutine get_waasmaier


  subroutine get_scattering_factor_params(element, wout)
    implicit none
    !      real*8, intent(in) :: s ! scattering vector: s = q / 2pi [1/A]
    real*8 :: w(1:9,1:215)
    character*20 :: elements(1:215)
    !   Input variables
    character*8, intent(in) :: element
    !   Output variables
    logical :: is_in_database = .false.
    !   Internal variables
    real*8, intent(out) :: wout(1:9)
    integer :: i

    elements = ['    H', "   H'", '    D', '  H1-', '   He', '  &
         & Li', ' Li1+', '   Be', ' Be2+', '    B', '    C', '&
         & Cval', '    N', '    O', '  O1-', '  O2-', '    F', ' &
         & F1-', '   Ne', '   Na', ' Na1+', '   Mg', ' Mg2+', '  &
         & Al', ' Al3+', '   Si', '  Siv', 'Sival', ' Si4+', '   &
         & P', '    S', '   Cl', ' Cl1-', '   Ar', '    K', '  K1+'&
         &, '   Ca', ' Ca2+', '   Sc', ' Sc3+', '   Ti', ' Ti2+', '&
         & Ti3+', ' Ti4+', '    V', '  V2+', '  V3+', '  V5+', '  &
         & Cr', ' Cr2+', ' Cr3+', '   Mn', ' Mn2+', ' Mn3+', ' Mn4&
         &+', '   Fe', ' Fe2+', ' Fe3+', '   Co', ' Co2+', ' Co3+',&
         & '   Ni', ' Ni2+', ' Ni3+', '   Cu', ' Cu1+', ' Cu2+', ' &
         &  Zn', ' Zn2+', '   Ga', ' Ga3+', '   Ge', ' Ge4+', '  &
         & As', '   Se', '   Br', ' Br1-', '   Kr', '   Rb', ' Rb1&
         &+', '   Sr', ' Sr2+', '    Y', '  Y3+', '   Zr', ' Zr4+',&
         & '   Nb', ' Nb3+', ' Nb5+', '   Mo', ' Mo3+', ' Mo5+', '&
         & Mo6+', '   Tc', '   Ru', ' Ru3+', ' Ru4+', '   Rh', '&
         & Rh3+', ' Rh4+', '   Pd', ' Pd2+', ' Pd4+', '   Ag', '&
         & Ag1+', ' Ag2+', '   Cd', ' Cd2+', '   In', ' In3+', '  &
         & Sn', ' Sn2+', ' Sn4+', '   Sb', ' Sb3+', ' Sb5+', '  &
         & Te', '    I', '  I1-', '   Xe', '   Cs', ' Cs1+', '  &
         & Ba', ' Ba2+', '   La', ' La3+', '   Ce', ' Ce3+', ' Ce4&
         &+', '   Pr', ' Pr3+', ' Pr4+', '   Nd', ' Nd3+', '   Pm',&
         & ' Pm3+', '   Sm', ' Sm3+', '   Eu', ' Eu2+', ' Eu3+', ' &
         &  Gd', ' Gd3+', '   Tb', ' Tb3+', '   Dy', ' Dy3+', '  &
         & Ho', ' Ho3+', '   Er', ' Er3+', '   Tm', ' Tm3+', '  &
         & Yb', ' Yb2+', ' Yb3+', '   Lu', ' Lu3+', '   Hf', ' Hf4&
         &+', '   Ta', ' Ta5+', '    W', '  W6+', '   Re', '   Os',&
         & ' Os4+', '   Ir', ' Ir3+', ' Ir4+', '   Pt', ' Pt2+', '&
         & Pt4+', '   Au', ' Au1+', ' Au3+', '   Hg', ' Hg1+', '&
         & Hg2+', '   Tl', ' Tl1+', ' Tl3+', '   Pb', ' Pb2+', '&
         & Pb4+', '   Bi', ' Bi3+', ' Bi5+', '   Po', '   At', '  &
         & Rn', '   Fr', '   Ra', ' Ra2+', '   Ac', ' Ac3+', '  &
         & Th', ' Th4+', '   Pa', '    U', '  U3+', '  U4+', '  U6&
         &+', '   Np', ' Np3+', ' Np4+', ' Np6+', '   Pu', ' Pu3+',&
         & ' Pu4+', ' Pu6+', '   Am', '   Cm', '   Bk', '   Cf']

    ! 9 rows with 215 columns. Each column corresponds to a given species.
    w(1:9,1:215) = reshape( &
         & (/ & ! a1, a2, a3, a4                  b1, b2, b3, b4,               c
         &  0.493002, 0.322912, 0.140191, 0.040810,   10.5109, 26.1257, 3.14236, 57.7997,  0.003038, & ! H
         &  0.489918, 0.262003, 0.196767, 0.049879,   20.6593, 7.74039, 49.5519, 2.20159,  0.001305, & ! H'
         &  0.489918, 0.262003, 0.196767, 0.049879,   20.6593, 7.74039, 49.5519, 2.20159,  0.001305, & ! D
         &  0.897661, 0.565616, 0.415815, 0.116973,   53.1368, 15.1870, 186.576, 3.56709,  0.002389, & ! H1-
         &  0.873400, 0.630900, 0.311200, 0.178000,   9.10370, 3.35680, 22.9276, 0.982100,  0.006400, & ! He
         &  1.12820, 0.750800, 0.617500, 0.465300,   3.95460, 1.05240, 85.3905, 168.261,  0.037700, & ! Li
         &  0.696800, 0.788800, 0.341400, 0.156300,   4.62370, 1.95570, 0.631600, 10.0953,  0.016700, & ! Li1+
         &  1.59190, 1.12780, 0.539100, 0.702900,   43.6427, 1.86230, 103.483, 0.542000,  0.038500, & ! Be
         &  6.26030, 0.884900, 0.799300, 0.164700,   0.002700, 0.831300, 2.27580, 5.11460,  -6.1092, & ! Be2+
         &  2.05450, 1.33260, 1.09790, 0.706800,   23.2185, 1.02100, 60.3498, 0.140300,  -0.19320, & ! B
         &  2.31000, 1.02000, 1.58860, 0.865000,   20.8439, 10.2075, 0.568700, 51.6512,  0.215600, & ! C
         &  2.26069, 1.56165, 1.05075, 0.839259,   22.6907, 0.656665, 9.75618, 55.5949,  0.286977, & ! Cval
         &  12.2126, 3.13220, 2.01250, 1.16630,   0.005700, 9.89330, 28.9975, 0.582600,  -11.529, & ! N
         &  3.04850, 2.28680, 1.54630, 0.867000,   13.2771, 5.70110, 0.323900, 32.9089,  0.250800, & ! O
         &  4.19160, 1.63969, 1.52673, -20.307,   12.8573, 4.17236, 47.0179, -0.01404,  21.9412, & ! O1-
         &  3.75040, 2.84294, 1.54298, 1.62091,   16.5151, 6.59203, 0.319201, 43.3486,  0.242060, & ! O2-
         &  3.53920, 2.64120, 1.51700, 1.02430,   10.2825, 4.29440, 0.261500, 26.1476,  0.277600, & ! F
         &  3.63220, 3.51057, 1.26064, 0.940706,   5.27756, 14.7353, 0.442258, 47.3437,  0.653396, & ! F1-
         &  3.95530, 3.11250, 1.45460, 1.12510,   8.40420, 3.42620, 0.230600, 21.7184,  0.351500, & ! Ne
         &  4.76260, 3.17360, 1.26740, 1.11280,   3.28500, 8.84220, 0.313600, 129.424,  0.676000, & ! Na
         &  3.25650, 3.93620, 1.39980, 1.00320,   2.66710, 6.11530, 0.200100, 14.0390,  0.404000, & ! Na1+
         &  5.42040, 2.17350, 1.22690, 2.30730,   2.82750, 79.2611, 0.380800, 7.19370,  0.858400, & ! Mg
         &  3.49880, 3.83780, 1.32840, 0.849700,   2.16760, 4.75420, 0.185000, 10.1411,  0.485300, & ! Mg2+
         &  6.42020, 1.90020, 1.59360, 1.96460,   3.03870, 0.742600, 31.5472, 85.0886,  1.11510, & ! Al
         &  4.17448, 3.38760, 1.20296, 0.528137,   1.93816, 4.14553, 0.228753, 8.28524,  0.706786, & ! Al3+
         &  6.29150, 3.03530, 1.98910, 1.54100,   2.43860, 32.3337, 0.678500, 81.6937,  1.14070, & ! Si
         &  6.29150, 3.03530, 1.98910, 1.54100,   2.43860, 32.3337, 0.678500, 81.6937,  1.14070, & ! Siv
         &  5.66269, 3.07164, 2.62446, 1.39320,   2.66520, 38.6634, 0.916946, 93.5458,  1.24707, & ! Sival
         &  4.43918, 3.20345, 1.19453, 0.416530,   1.64167, 3.43757, 0.214900, 6.65365,  0.746297, & ! Si4+
         &  6.43450, 4.17910, 1.78000, 1.49080,   1.90670, 27.1570, 0.526000, 68.1645,  1.11490, & ! P
         &  6.90530, 5.20340, 1.43790, 1.58630,   1.46790, 22.2151, 0.253600, 56.1720,  0.866900, & ! S
         &  11.4604, 7.19640, 6.25560, 1.64550,   0.010400, 1.16620, 18.5194, 47.7784,  -9.5574, & ! Cl
         &  18.2915, 7.20840, 6.53370, 2.33860,   0.006600, 1.17170, 19.5424, 60.4486,  -16.378, & ! Cl1-
         &  7.48450, 6.77230, 0.653900, 1.64420,   0.907200, 14.8407, 43.8983, 33.3929,  1.44450, & ! Ar
         &  8.21860, 7.43980, 1.05190, 0.865900,   12.7949, 0.774800, 213.187, 41.6841,  1.42280, & ! K
         &  7.95780, 7.49170, 6.35900, 1.19150,   12.6331, 0.767400, -0.00200, 31.9128,  -4.9978, & ! K1+
         &  8.62660, 7.38730, 1.58990, 1.02110,   10.4421, 0.659900, 85.7484, 178.437,  1.37510, & ! Ca
         &  15.6348, 7.95180, 8.43720, 0.853700,   -0.00740, 0.608900, 10.3116, 25.9905,  -14.875, & ! Ca2+
         &  9.18900, 7.36790, 1.64090, 1.46800,   9.02130, 0.572900, 136.108, 51.3531,  1.33290, & ! Sc
         &  13.4008, 8.02730, 1.65943, 1.57936,   0.298540, 7.96290, -0.28604, 16.0662,  -6.6667, & ! Sc3+
         &  9.75950, 7.35580, 1.69910, 1.90210,   7.85080, 0.500000, 35.6338, 116.105,  1.28070, & ! Ti
         &  9.11423, 7.62174, 2.27930, 0.087899,   7.52430, 0.457585, 19.5361, 61.6558,  0.897155, & ! Ti2+
         &  17.7344, 8.73816, 5.25691, 1.92134,   0.220610, 7.04716, -0.15762, 15.9768,  -14.652, & ! Ti3+
         &  19.5114, 8.23473, 2.01341, 1.52080,   0.178847, 6.67018, -0.29263, 12.9464,  -13.280, & ! Ti4+
         &  10.2971, 7.35110, 2.07030, 2.05710,   6.86570, 0.438500, 26.8938, 102.478,  1.21990, & ! V
         &  10.1060, 7.35410, 2.28840, 0.022300,   6.88180, 0.440900, 20.3004, 115.122,  1.22980, & ! V2+
         &  9.43141, 7.74190, 2.15343, 0.016865,   6.39535, 0.383349, 15.1908, 63.9690,  0.656565, & ! V3+
         &  15.6887, 8.14208, 2.03081, -9.5760,   0.679003, 5.40135, 9.97278, 0.940464,  1.71430, & ! V5+
         &  10.6406, 7.35370, 3.32400, 1.49220,   6.10380, 0.392000, 20.2626, 98.7399,  1.18320, & ! Cr
         &  9.54034, 7.75090, 3.58274, 0.509107,   5.66078, 0.344261, 13.3075, 32.4224,  0.616898, & ! Cr2+
         &  9.68090, 7.81136, 2.87603, 0.113575,   5.59463, 0.334393, 12.8288, 32.8761,  0.518275, & ! Cr3+
         &  11.2819, 7.35730, 3.01930, 2.24410,   5.34090, 0.343200, 17.8674, 83.7543,  1.08960, & ! Mn
         &  10.8061, 7.36200, 3.52680, 0.218400,   5.27960, 0.343500, 14.3430, 41.3235,  1.08740, & ! Mn2+
         &  9.84521, 7.87194, 3.56531, 0.323613,   4.91797, 0.294393, 10.8171, 24.1281,  0.393974, & ! Mn3+
         &  9.96253, 7.97057, 2.76067, 0.054447,   4.84850, 0.283303, 10.4852, 27.5730,  0.251877, & ! Mn4+
         &  11.7695, 7.35730, 3.52220, 2.30450,   4.76110, 0.307200, 15.3535, 76.8805,  1.03690, & ! Fe
         &  11.0424, 7.37400, 4.13460, 0.439900,   4.65380, 0.305300, 12.0546, 31.2809,  1.00970, & ! Fe2+
         &  11.1764, 7.38630, 3.39480, 0.072400,   4.61470, 0.300500, 11.6729, 38.5566,  0.970700, & ! Fe3+
         &  12.2841, 7.34090, 4.00340, 2.34880,   4.27910, 0.278400, 13.5359, 71.1692,  1.01180, & ! Co
         &  11.2296, 7.38830, 4.73930, 0.710800,   4.12310, 0.272600, 10.2443, 25.6466,  0.932400, & ! Co2+
         &  10.3380, 7.88173, 4.76795, 0.725591,   3.90969, 0.238668, 8.35583, 18.3491,  0.286667, & ! Co3+
         &  12.8376, 7.29200, 4.44380, 2.38000,   3.87850, 0.256500, 12.1763, 66.3421,  1.03410, & ! Ni
         &  11.4166, 7.40050, 5.34420, 0.977300,   3.67660, 0.244900, 8.87300, 22.1626,  0.861400, & ! Ni2+
         &  10.7806, 7.75868, 5.22746, 0.847114,   3.54770, 0.223140, 7.64468, 16.9673,  0.386044, & ! Ni3+
         &  13.3380, 7.16760, 5.61580, 1.67350,   3.58280, 0.247000, 11.3966, 64.8126,  1.19100, & ! Cu
         &  11.9475, 7.35730, 6.24550, 1.55780,   3.36690, 0.227400, 8.66250, 25.8487,  0.89000, & ! Cu1+
         &  11.8168, 7.11181, 5.78135, 1.14523,   3.37484, 0.244078, 7.98760, 19.8970,  1.14431, & ! Cu2+
         &  14.0743, 7.03180, 5.16520, 2.41000,   3.26550, 0.233300, 10.3163, 58.7097,  1.30410, & ! Zn
         &  11.9719, 7.38620, 6.46680, 1.39400,   2.99460, 0.203100, 7.08260, 18.0995,  0.780700, & ! Zn2+
         &  15.2354, 6.70060, 4.35910, 2.96230,   3.06690, 0.241200, 10.7805, 61.4135,  1.71890, & ! Ga
         &  12.6920, 6.69883, 6.06692, 1.00660,   2.81262, 0.227890, 6.36441, 14.4122,  1.53545, & ! Ga3+
         &  16.0816, 6.37470, 3.70680, 3.68300,   2.85090, 0.251600, 11.4468, 54.7625,  2.13130, & ! Ge
         &  12.9172, 6.70003, 6.06791, 0.859041,   2.53718, 0.205855, 5.47913, 11.6030,  1.45572, & ! Ge4+
         &  16.6723, 6.07010, 3.43130, 4.27790,   2.63450, 0.264700, 12.9479, 47.7972,  2.53100, & ! As
         &  17.0006, 5.81960, 3.97310, 4.35430,   2.40980, 0.272600, 15.2372, 43.8163,  2.84090, & ! Se
         &  17.1789, 5.23580, 5.63770, 3.98510,   2.17230, 16.5796, 0.260900, 41.4328,  2.95570, & ! Br
         &  17.1718, 6.33380, 5.57540, 3.72720,   2.20590, 19.3345, 0.287100, 58.1535,  3.17760, & ! Br1-
         &  17.3555, 6.72860, 5.54930, 3.53750,   1.93840, 16.5623, 0.226100, 39.3972,  2.82500, & ! Kr
         &  17.1784, 9.64350, 5.13990, 1.52920,   1.78880, 17.3151, 0.274800, 164.934,  3.48730, & ! Rb
         &  17.5816, 7.65980, 5.89810, 2.78170,   1.71390, 14.7957, 0.160300, 31.2087,  2.07820, & ! Rb1+
         &  17.5663, 9.81840, 5.42200, 2.66940,   1.55640, 14.0988, 0.166400, 132.376,  2.50640, & ! Sr
         &  18.0874, 8.13730, 2.56540, -34.193,   1.49070, 12.6963, 24.5651, -0.01380,  41.4025, & ! Sr2+
         &  17.7760, 10.2946, 5.72629, 3.26588,   1.40290, 12.8006, 0.125599, 104.354,  1.91213, & ! Y
         &  17.9268, 9.15310, 1.76795, -33.108,   1.35417, 11.2145, 22.6599, -0.01319,  40.2602, & ! Y3+
         &  17.8765, 10.9480, 5.41732, 3.65721,   1.27618, 11.9160, 0.117622, 87.6627,  2.06929, & ! Zr
         &  18.1668, 10.0562, 1.01118, -2.6479,   1.21480, 10.1483, 21.6054, -0.10276,  9.41454, & ! Zr4+
         &  17.6142, 12.0144, 4.04183, 3.53346,   1.18865, 11.7660, 0.204785, 69.7957,  3.75591, & ! Nb
         &  19.8812, 18.0653, 11.0177, 1.94715,   0.019175, 1.13305, 10.1621, 28.3389,  -12.912, & ! Nb3+
         &  17.9163, 13.3417, 10.7990, 0.337905,   1.12446, 0.028781, 9.28206, 25.7228,  -6.3934, & ! Nb5+
         &  3.70250, 17.2356, 12.8876, 3.74290,   0.277200, 1.09580, 11.0040, 61.6584,  4.38750, & ! Mo
         &  21.1664, 18.2017, 11.7423, 2.30951,   0.014734, 1.03031, 9.53659, 26.6307,  -14.421, & ! Mo3+
         &  21.0149, 18.0992, 11.4632, 0.740625,   0.014345, 1.02238, 8.78809, 23.3452,  -14.316, & ! Mo5+
         &  17.8871, 11.1750, 6.57891, 0.000000,   1.03649, 8.48061, 0.058881, 0.000000,  0.344941, & ! Mo6+
         &  19.1301, 11.0948, 4.64901, 2.71263,   0.864132, 8.14487, 21.5707, 86.8472,  5.40428, & ! Tc
         &  19.2674, 12.9182, 4.86337, 1.56756,   0.808520, 8.43467, 24.7997, 94.2928,  5.37874, & ! Ru
         &  18.5638, 13.2885, 9.32602, 3.00964,   0.847329, 8.37164, 0.017662, 22.8870,  -3.1892, & ! Ru3+
         &  18.5003, 13.1787, 4.71304, 2.18535,   0.844582, 8.12534, 0.36495, 20.8504,  1.42357, & ! Ru4+
         &  19.2957, 14.3501, 4.73425, 1.28918,   0.751536, 8.21758, 25.8749, 98.6062,  5.32800, & ! Rh
         &  18.8785, 14.1259, 3.32515, -6.1989,   0.764252, 7.84438, 21.2487, -0.01036,  11.8678, & ! Rh3+
         &  18.8545, 13.9806, 2.53464, -5.6526,   0.760825, 7.62436, 19.3317, -0.01020,  11.2835, & ! Rh4+
         &  19.3319, 15.5017, 5.29537, 0.605844,   0.698655, 7.98929, 25.2052, 76.8986,  5.26593, & ! Pd
         &  19.1701, 15.2096, 4.32234, 0.000000,   0.696219, 7.55573, 22.5057, 0.000000,  5.29160, & ! Pd2+
         &  19.2493, 14.7900, 2.89289, -7.9492,   0.683839, 7.14833, 17.9144, 0.005127,  13.0174, & ! Pd4+
         &  19.2808, 16.6885, 4.80450, 1.04630,   0.644600, 7.47260, 24.6605, 99.8156,  5.17900, & ! Ag
         &  19.1812, 15.9719, 5.27475, 0.357534,   0.646179, 7.19123, 21.7326, 66.1147,  5.21572, & ! Ag1+
         &  19.1643, 16.2456, 4.37090, 0.000000,   0.645643, 7.18544, 21.4072, 0.000000,  5.21404, & ! Ag2+
         &  19.2214, 17.6444, 4.46100, 1.60290,   0.594600, 6.90890, 24.7008, 87.4825,  5.06940, & ! Cd
         &  19.1514, 17.2535, 4.47128, 0.000000,   0.597922, 6.80639, 20.2521, 0.000000,  5.11937, & ! Cd2+
         &  19.1624, 18.5596, 4.29480, 2.03960,   0.547600, 6.37760, 25.8499, 92.8029,  4.93910, & ! In
         &  19.1045, 18.1108, 3.78897, 0.000000,   0.551522, 6.32470, 17.3595, 0.000000,  4.99635, & ! In3+
         &  19.1889, 19.1005, 4.45850, 2.46630,   5.83030, 0.503100, 26.8909, 83.9571,  4.78210, & ! Sn
         &  19.1094, 19.0548, 4.56480, 0.487000,   0.503600, 5.83780, 23.3752, 62.2061,  4.78610, & ! Sn2+
         &  18.9333, 19.7131, 3.41820, 0.019300,   5.76400, 0.465500, 14.0049, -0.75830,  3.91820, & ! Sn4+
         &  19.6418, 19.0455, 5.03710, 2.68270,   5.30340, 0.460700, 27.9074, 75.2825,  4.59090, & ! Sb
         &  18.9755, 18.9330, 5.10789, 0.288753,   0.467196, 5.22126, 19.5902, 55.5113,  4.69626, & ! Sb3+
         &  19.8685, 19.0302, 2.41253, 0.000000,   5.44853, 0.467973, 14.1259, 0.000000,  4.69263, & ! Sb5+
         &  19.9644, 19.0138, 6.14487, 2.52390,   4.81742, 0.420885, 28.5284, 70.8403,  4.35200, & ! Te
         &  20.1472, 18.9949, 7.51380, 2.27350,   4.34700, 0.381400, 27.7660, 66.8776,  4.07120, & ! I
         &  20.2332, 18.9970, 7.80690, 2.88680,   4.35790, 0.381500, 29.5259, 84.9304,  4.07140, & ! I1-
         &  20.2933, 19.0298, 8.97670, 1.99000,   3.92820, 0.344000, 26.4659, 64.2658,  3.71180, & ! Xe
         &  20.3892, 19.1062, 10.6620, 1.49530,   3.56900, 0.310700, 24.3879, 213.904,  3.33520, & ! Cs
         &  20.3524, 19.1278, 10.2821, 0.961500,   3.55200, 0.308600, 23.7128, 59.4565,  3.27910, & ! Cs1+
         &  20.3361, 19.2970, 10.8880, 2.69590,   3.21600, 0.275600, 20.2073, 167.202,  2.77310, & ! Ba
         &  20.1807, 19.1136, 10.9054, 0.77634,   3.21367, 0.283310, 20.0558, 51.7460,  3.02902, & ! Ba2+
         &  20.5780, 19.5990, 11.3727, 3.28719,   2.94817, 0.244475, 18.7726, 133.124,  2.14678, & ! La
         &  20.2489, 19.3763, 11.6323, 0.336048,   2.92070, 0.250698, 17.8211, 54.9453,  2.40860, & ! La3+
         &  21.1671, 19.7695, 11.8513, 3.33049,   2.81219, 0.226836, 17.6083, 127.113,  1.86264, & ! Ce
         &  20.8036, 19.5590, 11.9369, 0.612376,   2.77691, 0.231540, 16.5408, 43.1692,  2.09013, & ! Ce3+
         &  20.3235, 19.8186, 12.1233, 0.144583,   2.65941, 0.218850, 15.7992, 62.2355,  1.59180, & ! Ce4+
         &  22.0440, 19.6697, 12.3856, 2.82428,   2.77393, 0.222087, 16.7669, 143.644,  2.05830, & ! Pr
         &  21.3727, 19.7491, 12.1329, 0.975180,   2.64520, 0.214299, 15.3230, 36.4065,  1.77132, & ! Pr3+
         &  20.9413, 20.0539, 12.4668, 0.296689,   2.54467, 0.202481, 14.8137, 45.4643,  1.24285, & ! Pr4+
         &  22.6845, 19.6847, 12.7740, 2.85137,   2.66248, 0.210628, 15.8850, 137.903,  1.98486, & ! Nd
         &  21.9610, 19.9339, 12.1200, 1.51031,   2.52722, 0.199237, 14.1783, 30.8717,  1.47588, & ! Nd3+
         &  23.3405, 19.6095, 13.1235, 2.87516,   2.56270, 0.202088, 15.1009, 132.721,  2.02876, & ! Pm
         &  22.5527, 20.1108, 12.0671, 2.07492,   2.41740, 0.185769, 13.1275, 27.4491,  1.19499, & ! Pm3+
         &  24.0042, 19.4258, 13.4396, 2.89604,   2.47274, 0.196451, 14.3996, 128.007,  2.20963, & ! Sm
         &  23.1504, 20.2599, 11.9202, 2.71488,   2.31641, 0.174081, 12.1571, 24.8242,  0.954586, & ! Sm3+
         &  24.6274, 19.0886, 13.7603, 2.92270,   2.38790, 0.194200, 13.7546, 123.174,  2.57450, & ! Eu
         &  24.0063, 19.9504, 11.8034, 3.87243,   2.27783, 0.173530, 11.6096, 26.5156,  1.36389, & ! Eu2+
         &  23.7497, 20.3745, 11.8509, 3.26503,   2.22258, 0.163940, 11.3110, 22.9966,  0.759344, & ! Eu3+
         &  25.0709, 19.0798, 13.8518, 3.54545,   2.25341, 0.181951, 12.9331, 101.398,  2.41960, & ! Gd
         &  24.3466, 20.4208, 11.8708, 3.71490,   2.13553, 0.155525, 10.5782, 21.7029,  0.645089, & ! Gd3+
         &  25.8976, 18.2185, 14.3167, 2.95354,   2.24256, 0.196143, 12.6648, 115.362,  3.58324, & ! Tb
         &  24.9559, 20.3271, 12.2471, 3.77300,   2.05601, 0.149525, 10.0499, 21.2773,  0.691967, & ! Tb3+
         &  26.5070, 17.6383, 14.5596, 2.96577,   2.18020, 0.202172, 12.1899, 111.874,  4.29728, & ! Dy
         &  25.5395, 20.2861, 11.9812, 4.50073,   1.98040, 0.143384, 9.34972, 19.5810,  0.689690, & ! Dy3+
         &  26.9049, 17.2940, 14.5583, 3.63837,   2.07051, 0.197940, 11.4407, 92.6566,  4.56796, & ! Ho
         &  26.1296, 20.0994, 11.9788, 4.93676,   1.91072, 0.139358, 8.80018, 18.5908,  0.852795, & ! Ho3+
         &  27.6563, 16.4285, 14.9779, 2.98233,   2.07356, 0.223545, 11.3604, 105.703,  5.92046, & ! Er
         &  26.7220, 19.7748, 12.1506, 5.17379,   1.84659, 0.137290, 8.36225, 17.8974,  1.17613, & ! Er3+
         &  28.1819, 15.8851, 15.1542, 2.98706,   2.02859, 0.238849, 10.9975, 102.961,  6.75621, & ! Tm
         &  27.3083, 19.3320, 12.3339, 5.38348,   1.78711, 0.136974, 7.96778, 17.2922,  1.63929, & ! Tm3+
         &  28.6641, 15.4345, 15.3087, 2.98963,   1.98890, 0.257119, 10.6647, 100.417,  7.56672, & ! Yb
         &  28.1209, 17.6817, 13.3335, 5.14657,   1.78503, 0.159970, 8.18304, 20.3900,  3.70983, & ! Yb2+
         &  27.8917, 18.7614, 12.6072, 5.47647,   1.73272, 0.138790, 7.64412, 16.8153,  2.26001, & ! Yb3+
         &  28.9476, 15.2208, 15.1000, 3.71601,   1.90182, 9.98519, 0.261033, 84.3298,  7.97628, & ! Lu
         &  28.4628, 18.1210, 12.8429, 5.59415,   1.68216, 0.142292, 7.33727, 16.3535,  2.97573, & ! Lu3+
         &  29.1440, 15.1726, 14.7586, 4.30013,   1.83262, 9.59990, 0.275116, 72.0290,  8.58154, & ! Hf
         &  28.8131, 18.4601, 12.7285, 5.59927,   1.59136, 0.128903, 6.76232, 14.0366,  2.39699, & ! Hf4+
         &  29.2024, 15.2293, 14.5135, 4.76492,   1.77333, 9.37046, 0.295977, 63.3644,  9.24354, & ! Ta
         &  29.1587, 18.8407, 12.8268, 5.38695,   1.50711, 0.116741, 6.31524, 12.4244,  1.78555, & ! Ta5+
         &  29.0818, 15.4300, 14.4327, 5.11982,   1.72029, 9.22590, 0.321703, 57.0560,  9.88750, & ! W
         &  29.4936, 19.3763, 13.0544, 5.06412,   1.42755, 0.104621, 5.93667, 11.1972,  1.01074, & ! W6+
         &  28.7621, 15.7189, 14.5564, 5.44174,   1.67191, 9.09227, 0.350500, 52.0861,  10.4720, & ! Re
         &  28.1894, 16.1550, 14.9305, 5.67589,   1.62903, 8.97948, 0.382661, 48.1647,  11.0005, & ! Os
         &  30.4190, 15.2637, 14.7458, 5.06795,   1.37113, 6.84706, 0.165191, 18.0030,  6.49804, & ! Os4+
         &  27.3049, 16.7296, 15.6115, 5.83377,   1.59279, 8.86553, 0.417916, 45.0011,  11.4722, & ! Ir
         &  30.4156, 15.8620, 13.6145, 5.82008,   1.34323, 7.10909, 0.204633, 20.3254,  8.27903, & ! Ir3+
         &  30.7058, 15.5512, 14.2326, 5.53672,   1.30923, 6.71983, 0.167252, 17.4911,  6.96824, & ! Ir4+
         &  27.0059, 17.7639, 15.7131, 5.78370,   1.51293, 8.81174, 0.424593, 38.6103,  11.6883, & ! Pt
         &  29.8429, 16.7224, 13.2153, 6.35234,   1.32927, 7.38979, 0.263297, 22.9426,  9.85329, & ! Pt2+
         &  30.9612, 15.9829, 13.7348, 5.92034,   1.24813, 6.60834, 0.168640, 16.9392,  7.39534, & ! Pt4+
         &  16.8819, 18.5913, 25.5582, 5.86000,   0.461100, 8.62160, 1.48260, 36.3956,  12.0658, & ! Au
         &  28.0109, 17.8204, 14.3359, 6.58077,   1.35321, 7.73950, 0.356752, 26.4043,  11.2299, & ! Au1+
         &  30.6886, 16.9029, 12.7801, 6.52354,   1.21990, 6.82872, 0.212867, 18.6590,  9.09680, & ! Au3+
         &  20.6809, 19.0417, 21.6575, 5.96760,   0.545000, 8.44840, 1.57290, 38.3246,  12.6089, & ! Hg
         &  25.0853, 18.4973, 16.8883, 6.48216,   1.39507, 7.65105, 0.443378, 28.2262,  12.0205, & ! Hg1+
         &  29.5641, 18.0600, 12.8374, 6.89912,   1.21152, 7.05639, 0.284738, 20.7482,  10.6268, & ! Hg2+
         &  27.5446, 19.1584, 15.5380, 5.52593,   0.655150, 8.70751, 1.96347, 45.8149,  13.1746, & ! Tl
         &  21.3985, 20.4723, 18.7478, 6.82847,   1.47110, 0.517394, 7.43463, 28.8482,  12.5258, & ! Tl1+
         &  30.8695, 18.3841, 11.9328, 7.00574,   1.10080, 6.53852, 0.219074, 17.2114,  9.80270, & ! Tl3+
         &  31.0617, 13.0637, 18.4420, 5.96960,   0.690200, 2.35760, 8.61800, 47.2579,  13.4118, & ! Pb
         &  21.7886, 19.5682, 19.1406, 7.01107,   1.33660, 0.488383, 6.77270, 23.8132,  12.4734, & ! Pb2+
         &  32.1244, 18.8003, 12.0175, 6.96886,   1.00566, 6.10926, 0.147041, 14.7140,  8.08428, & ! Pb4+
         &  33.3689, 12.9510, 16.5877, 6.46920,   0.704000, 2.92380, 8.79370, 48.0093,  13.5782, & ! Bi
         &  21.8053, 19.5026, 19.1053, 7.10295,   1.23560, 6.24149, 0.469999, 20.3185,  12.4711, & ! Bi3+
         &  33.5364, 25.0946, 19.2497, 6.91555,   0.916540, 0.39042, 5.71414, 12.8285,  -6.7994, & ! Bi5+
         &  34.6726, 15.4733, 13.1138, 7.02588,   0.700999, 3.55078, 9.55642, 47.0045,  13.6770, & ! Po
         &  35.3163, 19.0211, 9.49887, 7.42518,   0.685870, 3.97458, 11.3824, 45.4715,  13.7108, & ! At
         &  35.5631, 21.2816, 8.00370, 7.44330,   0.663100, 4.06910, 14.0422, 44.2473,  13.6905, & ! Rn
         &  35.9299, 23.0547, 12.1439, 2.11253,   0.646453, 4.17619, 23.1052, 150.645,  13.7247, & ! Fr
         &  35.7630, 22.9064, 12.4739, 3.21097,   0.616341, 3.87135, 19.9887, 142.325,  13.6211, & ! Ra
         &  35.2150, 21.6700, 7.91342, 7.65078,   0.604909, 3.57670, 12.6010, 29.8436,  13.5431, & ! Ra2+
         &  35.6597, 23.1032, 12.5977, 4.08655,   0.589092, 3.65155, 18.5990, 117.020,  13.5266, & ! Ac
         &  35.1736, 22.1112, 8.19216, 7.05545,   0.579689, 3.41437, 12.9187, 25.9443,  13.4637, & ! Ac3+
         &  35.5645, 23.4219, 12.7473, 4.80703,   0.563359, 3.46204, 17.8309, 99.1722,  13.4314, & ! Th
         &  35.1007, 22.4418, 9.78554, 5.29444,   0.555054, 3.24498, 13.4661, 23.9533,  13.3760, & ! Th4+
         &  35.8847, 23.2948, 14.1891, 4.17287,   0.547751, 3.41519, 16.9235, 105.251,  13.4287, & ! Pa
         &  36.0228, 23.4128, 14.9491, 4.18800,   0.529300, 3.32530, 16.0927, 100.613,  13.3966, & ! U
         &  35.5747, 22.5259, 12.2165, 5.37073,   0.520480, 3.12293, 12.7148, 26.3394,  13.3092, & ! U3+
         &  35.3715, 22.5326, 12.0291, 4.79840,   0.516598, 3.05053, 12.5723, 23.4582,  13.2671, & ! U4+
         &  34.8509, 22.7584, 14.0099, 1.21457,   0.507079, 2.89030, 13.1767, 25.2017,  13.1665, & ! U6+
         &  36.1874, 23.5964, 15.6402, 4.18550,   0.511929, 3.25396, 15.3622, 97.4908,  13.3573, & ! Np
         &  35.7074, 22.6130, 12.9898, 5.43227,   0.502322, 3.03807, 12.1449, 25.4928,  13.2544, & ! Np3+
         &  35.5103, 22.5787, 12.7766, 4.92159,   0.498626, 2.96627, 11.9484, 22.7502,  13.2116, & ! Np4+
         &  35.0136, 22.7286, 14.3884, 1.75669,   0.489810, 2.81099, 12.3300, 22.6581,  13.1130, & ! Np6+
         &  36.5254, 23.8083, 16.7707, 3.47947,   0.499384, 3.26371, 14.9455, 105.980,  13.3812, & ! Pu
         &  35.8400, 22.7169, 13.5807, 5.66016,   0.484938, 2.96118, 11.5331, 24.3992,  13.1991, & ! Pu3+
         &  35.6493, 22.6460, 13.3595, 5.18831,   0.481422, 2.89020, 11.3160, 21.8301,  13.1555, & ! Pu4+
         &  35.1736, 22.7181, 14.7635, 2.28678,   0.473204, 2.73848, 11.5530, 20.9303,  13.0582, & ! Pu6+
         &  36.6706, 24.0992, 17.3415, 3.49331,   0.483629, 3.20647, 14.3136, 102.273,  13.3592, & ! Am
         &  36.6488, 24.4096, 17.3990, 4.21665,   0.465154, 3.08997, 13.4346, 88.4834,  13.2887, & ! Cm
         &  36.7881, 24.7736, 17.8919, 4.23284,   0.451018, 3.04619, 12.8946, 86.0030,  13.2754, & ! Bk
         &  36.9185, 25.1995, 18.3317, 4.24391,   0.437533, 3.00775, 12.4044, 83.7881,  13.2674 & ! Cf
         & /), [9,215] )




    outer: do i = 1, size(elements, 1)
       if( trim(adjustl(elements(i))) == trim(adjustl(element)) )then
          is_in_database = .true.
          wout(1:9) = w(1:9, i)
          ! f = w(9, i)
          ! s_sq = s * s

          ! do j = 1, 4
          !    f = f + w(j, i) * exp( - w( 4 + j, i ) * s_sq  )
          ! end do
          exit outer
       end if
    end do outer

    if ( .not. is_in_database )then
       write(*,'(A,1X,A,1X,A)') '!!! XRD Warning !!! The element ', element, ' is not in the scattering factor database! '
    end if

  end subroutine get_scattering_factor_params

  subroutine get_scattering_factor(f, w, s )
    implicit none
    real*8, intent(in) :: w(1:9), s
    real*8, intent(out) :: f
    real*8 :: s_sq
    integer :: j

    f = w(9)
    s_sq = s * s

    do j = 1, 4
       f = f + w(j) * exp( - w( j + 4 ) * s_sq  )
    end do
  end subroutine get_scattering_factor


end module exp_utils
