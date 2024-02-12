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

    subroutine get_partial_structure_factor( structure_factor_partial&
         &, pair_correlation_partial, q_list , rs, r_cut,&
         & n_samples_pc, n_samples_sf, n_species, n_atoms_of_species,&
         & n_sites, window)
      implicit none
      real*8 , intent(in), allocatable :: n_atoms_of_species(:)
      real*8, intent(out) :: structure_factor_partial(:,:,:),&
           & pair_correlation_partial(:,:,:), q_list(:), rs(:), r_cut
      integer, intent(in) :: n_samples_pc, n_samples_sf, n_species, n_sites
      real*8 :: r, q, dr, ca, cb, cabh, w
      integer :: i, j, k, l,  n
      real*8, parameter :: pi = acos(-1.0)
      logical, intent(in) :: window

      ! S_ab(q) = delta_ab + 4 pi rho (ca cb)^0.5
      !                  * int_0^r_cut dr r^2 [ g_ab(r) - 1 ] sin(qr)/(qr) * sin( pi r / R )/ (pi r /R)
      ! window corresponds to sin(pi r / R) / (pi r / R)
      ! rs has size n_samples_pc

      structure_factor_partial = 0.d0
      n = size(q_list)

      dr = rs(2) - rs(1)

      w = 1.d0

      do i = 1, n_species
         ca = n_atoms_of_species(i) / dfloat( n_sites )
         do j = 1, n_species

            if (i > j) then
               ! We'e calculated the structure factor
               structure_factor_partial(1:n, i, j ) = structure_factor_partial(1:n, j, i )
               cycle
            end if

            cb = n_atoms_of_species(j) / dfloat( n_sites )
            cabh = ( ca * cb )**(0.5)

            ! Iterate through each q
            do k = 1, n
               q = q_list(k)

               do l = 1, n_samples_pc
                  ! do the integral
                  if (window) w = sinc( pi * rs(l) / r_cut )

                  structure_factor_partial(k,i,j) = structure_factor_partial(k,i,j) + &
                       & dr * rs(l)**2 &
                       & * ( pair_correlation_partial(l, i, j) - 1.d0 ) &
                       & * sinc( q * rs(l) ) * w
               end do

               structure_factor_partial(k,i,j) = 4.d0 * pi * cabh * structure_factor_partial(k,i,j)

               if (i == j)then
                  structure_factor_partial(k,i,j)  = structure_factor_partial(k,i,j)  + 1.d0
               end if
            end do
         end do
      end do

    end subroutine get_partial_structure_factor



    subroutine get_pair_correlation(  n_sites0, &
         & neighbors_list, n_neigh, neighbor_species, rjs, r_min, r_max, n_samples, x,&
         & pair_correlation, r_cut, return_histogram,&
         & partial_rdf, species_1, species_2, kde_sigma )
      implicit none
      real*8,  intent(in) :: rjs(:), kde_sigma
      real*8 :: r_min, r_max, r_cut
      integer, intent(in) :: neighbors_list(:), n_neigh(:), neighbor_species(:)
      integer, intent(in) :: n_sites0, n_samples, species_1, species_2
      integer :: n_sites, n_pairs, count, count_species_1
      integer :: i, j, k, i2, j2, l, ind_bin_l, ind_bin_h, species_i, species_j
      real*8,  intent(inout) :: pair_correlation(1:n_samples), x(1:n_samples)
      real*8, allocatable :: bin_edges(:), dV(:), kde(:)
      real*8 :: r, n_pc
      real*8, parameter :: pi = acos(-1.0)
      logical, intent(in) :: partial_rdf
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

      x = 0.d0
      bin_edges = 0.d0
      pair_correlation = 0.d0

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

            if ( r < 1e-3 .or. r > r_cut)then
               cycle
            end if

            ! edge cases where r is not in range
            if (r < r_min)then
               !            print *, "pair_correlation_function: Given r is less than r_min! Continuing loop "
               cycle
            end if

            if (r > r_max)then
               !            print *, "pair_correlation_function: Given r is more than r_max! Continuing loop "
               cycle
            end if

            if (kde_sigma > 0.d0)then
               ! do a kernel density estimate
               count = count + 1
               kde = 0.d0
               do l = 1,n_samples
                  kde(l) = kde(l) +  exp( -( (x(l) - r) / kde_sigma )**2 / 2.d0 )
               end do
               pair_correlation(1:n_samples) = pair_correlation(1:n_samples) + &
                    & kde(1:n_samples)
            else

               call binary_search_real( r, bin_edges, 1, n_samples+1, ind_bin_l, ind_bin_h )
               pair_correlation(ind_bin_l) = pair_correlation(ind_bin_l) +  1.d0
            end if

         end do
      end do


      ! I don't know why the factor of 1 / n_samples / 0.1 works, but
      ! it seems to do so for all values of sigma and samples checked
      ! for...

      if (kde_sigma>0.d0)then
         pair_correlation(1:n_samples) =&
              & pair_correlation(1:n_samples) / sqrt(2.d0 * pi) /&
              & (kde_sigma) / dfloat(n_samples) / 0.1d0
      end if

      !( dfloat(count_species_1) )

      ! Remember that this is done for each process, so then the final
      ! result should be a reduction summation of all of these
      ! pair_correlation arrays multiplied by the 1/density factor ( pair_correlation *  (V / n_sites) )

      if ( .not. return_histogram )then
         ! Instead of giving the Radial distribution function (which goes as r^2), we give the pair distribution function
         pair_correlation =  pair_correlation / dV
         deallocate(dV)
      end if

      deallocate(bin_edges)
      if (kde_sigma>0.d0) deallocate(kde)

    end subroutine get_pair_correlation


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


    !**************************************************************************
    subroutine calculate_exp_interpolation(x, y, n_samples, data)
      implicit none
      real*8, allocatable, intent(in) ::  data(:,:)
      integer, intent(in) :: n_samples
      real*8, allocatable, intent(inout) :: x(:), y(:)
      real*8 :: dx

      allocate(x(1:n_samples))
      allocate(y(1:n_samples))

      x = 0.d0
      y = 0.d0

      call interpolate_data( x, y, data(1,:), data(2,:), n_samples, dx )

    end subroutine calculate_exp_interpolation


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
    

    subroutine get_experimental_forces(data, energy_scale, core_electron_be, core_electron_be_der,&
         & n_neigh, neighbors_list, neighbor_species, &
         & sigma, n_samples, norm, &
         & x, y_exp,  y, y_all, get_exp, do_forces, rjs, xyz,&
         & n_tot, energies_lp, forces0, virial, exp_similarity_type, quantity_name, rank)
      implicit none
      real*8, allocatable, intent(in) :: data(:,:)
      integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
      real*8, intent(in) :: sigma, core_electron_be(:),&
           & core_electron_be_der(:,:), rjs(:), xyz(:,:),&
           & energy_scale, norm
      real*8, allocatable, intent(inout) :: forces0(:,:)
      real*8, allocatable, intent(in) :: x(:), y_exp(:), y(:)
      real*8, intent(in) :: y_all(:,:)
      real*8, allocatable :: der_sum(:)
      real*8, intent(inout) :: energies_lp(:)
      real*8, intent(inout) :: virial(1:3,1:3)
      real*8 ::  this_force(1:3), t
      real*8, allocatable ::  y_der(:,:), der_factor(:,:), der_vec(:,:,:), prefactor(:), dxa(:)
      integer, intent(in) :: n_samples, n_tot, rank
      logical, intent(in) :: get_exp, do_forces
      integer :: n_sites, n_pairs, n_pairs_soap, n_species, n_sites0
      integer :: i, j, i2, j2, k, l,  n_in_buffer, k1, k2, mag_force_i
      logical :: vector_norm = .false.
      real*8 :: x_val, x_min, x_max, x_range, max_force, mag_force, max_mag_force, similarity, dx, &
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


    subroutine get_xrd_from_partial_structure_factors(&
         & structure_factor_partial, n_species, species_types,&
         & species, wavelength, damping, alpha, method, use_iwasa,&
         & x, y, n_atoms_of_species )
      implicit none
      real*8, intent(in) :: damping, wavelength, alpha, structure_factor_partial(:,:,:)
      integer, intent(in) :: species(:)
      logical, intent(in) :: use_iwasa
      integer, intent(in) :: n_species
      character*8, allocatable :: species_types(:)
      character*32, intent(in) :: method
      real*8 :: prefactor, p, c, c2, rij, diff(1:3), intensity, mag, sth, wfaci, wfacj, wfac_n, ntot
      integer :: i, j, k, l, n
      real*8 , intent(in), allocatable :: n_atoms_of_species(:)
      real*8, intent(in) :: x(:)
      real*8, intent(out) :: y(:)
      real*8, parameter :: pi = acos(-1.0)
      real*8, allocatable :: s(:)

      y = 0.d0
      n = size(x)

      allocate(s(1:n))

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

      ntot = sum(n_atoms_of_species)

      do l = 1, n
         prefactor = exp( - damping * s(l)**(2.d0) / 2.d0 )
         if (use_iwasa)then
            sth = ( wavelength * s(l) / 2.0d0 )
            p = 1.0d0 - sth * sth
            if (p < 0) p = 0.0d0
            c = sqrt(p)
            c2 = cos(2.0d0 * acos(c))
            prefactor = prefactor  *  c / (1.0d0 + alpha * c2**(2.d0))
         end if

         wfac_n = 0.d0
         do i = 1, n_species
            call get_waasmaier(species_types(i), s(l), wfaci)
            wfac_n = wfac_n + wfaci*wfaci*( n_atoms_of_species(i) / ntot )

            do j = 1, n_species
               call get_waasmaier(species_types(j), s(l), wfacj)


               y(l) = wfaci * wfacj * ( n_atoms_of_species(i) / ntot ) * ( n_atoms_of_species(j) / ntot ) &
                    & * structure_factor_partial(l,i,j)
            end do
         end do
         y(l) = y(l) / wfac_n
      end do


      
    end subroutine get_xrd_from_partial_structure_factors


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
    subroutine get_xrd_single_process( positions, n_species, species_types, species, wavelength, damping, alpha, &
                      & method, use_iwasa, x_min, x_max,  n_samples, x_i_exp, y_i_pred )
      implicit none
      real*8, allocatable, intent(in) :: positions(:,:)
      real*8, intent(in) :: damping, wavelength, alpha, x_min, x_max
      integer, intent(in) :: species(:)
      real*8, allocatable :: x(:), y(:), s(:), x_exp(:), y_exp(:), wfac(:), wfac_species(:)
      real*8, allocatable, intent(inout) :: x_i_exp(:), y_i_pred(:)
      logical, intent(in) :: use_iwasa
      integer, intent(in) :: n_samples, n_species
      character*8, allocatable :: species_types(:)
      character*32, intent(in) :: method
      real*8 :: prefactor, p, c, c2, rij, diff(1:3), intensity, mag, sth, wfac_n
      integer :: i, j, k, l, n, n_sites
      real*8 :: x_val, x_range, dx, pi=3.14159265359

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
            call get_waasmaier(species_types(i), s(l), wfac_species(i))
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

    
    subroutine get_waasmaier(element, s, f)
      implicit none
      real*8, intent(in) :: s ! scattering vector: s = q / 2pi [1/A]
      real*8 :: w(1:11,1:12)
      character*8 :: elements(1:12)
      !   Input variables
      character*8, intent(in) :: element
      !   Output variables
      logical :: is_in_database = .false.
      !   Internal variables
      real*8 ::  s_sq
      real*8, intent(out) :: f
      integer :: i, j

      ! D. Waasmaier and A. Kirfel, Acta Cryst. (1995). A51, 416-431
     !       a1         b1         a2         b2         a3         b3          a4         b4         a5         b5         c
      w(1:11,1:12) = reshape( &
           & (/ 2.657506, 14.780758,  1.078079,  0.776775,  1.490909, 42.086843,&
           & -4.241070,  -0.000294,  0.713791, 0.239535,  4.297983, &  !   C(1:10)
           & 11.893780,  0.000158,  3.277479, 10.232723,  1.858092, 30.344690, &
           & 0.858927,   0.656065,  0.912985, 0.217287, -11.804902, &  !   N(1:10)
           &  2.960427, 14.182259,  2.5088111, 5.936858,  0.637053,  0.112726, &
           & 0.722838,  34.958481,  1.142756, 0.390240,   0.027014, & !   O(1:10)
           &  1.950541,  0.908139,  4.146930, 27.044953,  1.494560,  0.071280, &
           & 1.522042,  67.520190,  5.729711, 1.981173,   0.155233, & !   P(1:10)
           &  6.372157,  1.514347,  5.154568, 22.092528,  1.473732,  0.061373, &
           & 1.635073,  55.445176,  1.209372, 0.646925,   0.154722, & !   S(1:10)
           &  1.446071,  0.052357,  6.870609,  1.193165,  6.151801, 18.343416, &
           & 1.750347,  46.398394,  0.634168, 0.401005,   0.146773, & !  Cl(1:10)
           & 13.521865,  4.077277,  6.947285,  0.286763,  3.866028, 14.622634, &
           & 2.135900,  71.966078,  4.284731, 0.004437,  -2.762697, & !  Ni(1:10)
           & 14.014192,  3.738280,  4.784577,  0.003744,  5.056806, 13.034982, &
           & 1.457971,  72.554793,  6.932996, 0.265666,  -3.774477, & !  Cu(1:10)
           &  6.121511,  0.062549,  4.784063,  0.784031, 16.631683,  8.751391, &
           & 4.318258,  34.489983, 13.246773, 0.784031,   0.883099, & !  Pd(1:10)
           &  6.073874,  0.055333, 17.155437,  7.896512,  4.173344, 28.443739, &
           & 0.852238, 110.376108, 17.988685, 0.716809,   0.756603, & !  Ag(1:10)
           & 31.273891,  1.316992, 18.445441,  8.797154, 17.063745,  0.124741, &
           & 5.555933,  40.177994,  1.575270, 1.316997,   4.050394, & !  Pt(1:10)
           & 16.777389,  0.122737, 19.317156,  8.621570, 32.979682,  1.256902, &
           & 5.595453,  38.008821, 10.576854, 0.000601,  -6.279078 /),  [11,12] )  !  Au(1:10)


      elements = [ ' C', ' N', ' O', ' P', ' S', 'Cl', 'Na', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']

      do i = 1, size(elements, 1)
         if( trim(adjustl(elements(i))) == trim(adjustl(element)) )then
            is_in_database = .true.
            f = w(11, i)
            s_sq = s * s

            do j = 1, 5
               f = f + w((1 + 2*(j-1)), i) * exp( - w( (2 + 2*(j-1)), i ) * s_sq  )
            end do
            exit
         end if
      end do

      if ( .not. is_in_database )then
         write(*,'(A,1X,A,1X,A)') '!!! XRD Warning !!! The element ', element, ' is not in the database of Waasmaier '
      end if

    end subroutine get_waasmaier


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
      integer :: i, m, n
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





    subroutine get_exp_pred_spectra_energies_forces(data, energy_scale, core_electron_be, core_electron_be_der,&
         & n_neigh, neighbors_list, neighbor_species, &
         & sigma, n_samples, norm, &
         & x, y_exp,  y, y_all, get_exp, do_forces, rjs, xyz,&
         & n_tot, energies_lp, forces0, virial, exp_similarity_type, rank)
      implicit none
      real*8, allocatable, intent(in) :: data(:,:)
      integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
      real*8, intent(in) :: sigma, core_electron_be(:),&
           & core_electron_be_der(:,:), rjs(:), xyz(:,:),&
           & energy_scale, norm
      real*8, allocatable, intent(inout) :: forces0(:,:)
      real*8, allocatable, intent(in) :: x(:), y_exp(:), y(:)
      real*8, intent(in) :: y_all(:,:)
      real*8, allocatable :: der_sum(:)
      real*8, intent(inout) :: energies_lp(:)
      real*8, intent(inout) :: virial(1:3,1:3)
      real*8 ::  this_force(1:3), t
      real*8, allocatable ::  y_der(:,:), der_factor(:,:), der_vec(:,:,:), prefactor(:), dxa(:)
      integer, intent(in) :: n_samples, n_tot, rank
      logical, intent(in) :: get_exp, do_forces
      integer :: n_sites, n_pairs, n_pairs_soap, n_species, n_sites0
      integer :: i, j, i2, j2, k, l,  n_in_buffer, k1, k2, mag_force_i
      logical :: vector_norm = .false.
      real*8 :: x_val, x_min, x_max, x_range, max_force, mag_force, max_mag_force, similarity, dx, &
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
            energies_lp(i) = - energy_scale * ( dot_product( y_all(:,i), y_exp ))
         end do
      else if (exp_similarity_type == 'squared_diff')then
         energies_lp = + energy_scale / n_sites0 * ( dx * dot_product( (y - y_exp), (y - y_exp) ) )
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

    end subroutine get_exp_pred_spectra_energies_forces



    subroutine get_data_similarity(x, y, y_pred, sim_exp_pred, exp_similarity_type)
      implicit none
      real*8, allocatable, intent(in) :: x(:), y(:), y_pred(:)
      character*32, intent(in) :: exp_similarity_type
      real*8, intent(out) :: sim_exp_pred

      if (exp_similarity_type == "squared_diff")then
         sim_exp_pred = - (x(2) - x(1) ) * dot_product(y - y_pred, y - y_pred)
      else
         sim_exp_pred =  dot_product(y, y_pred)
      end if
    end subroutine get_data_similarity



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

      call get_data_similarity(x_i_exp, y_i_exp, y_i_pred, sim_exp_pred, exp_similarity_type)

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
         y(i) = y(i) +  exp( -( x(i) - x0 )**2 / (2*sigma**2) )
         y_all(idx, i) = exp( -( x(i) - x0 )**2 / (2*sigma**2) )
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
      real*8 :: x_val, x_min, x_max, x_range, t, dx
      real*8, intent(out) :: mag
      logical, intent(in) :: broaden

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

      mag = sqrt(dot_product(y, y)) * dx !sum(y * dx) !
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
      real*8, allocatable :: g(:,:)
      real*8, intent(out) :: y_tot(1:3)
      integer :: i,idx
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
         idx = int( (((x(i)-x_min) / x_range) * real(size(xi))) + 1 )

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
      dx = x_range / real(n_samples)

      do i = 1, n_samples
         t = (real(i-1) / real(n_samples-1))
         x(i) = (1.d0 - t) * x_min  +  t * x_max !range
      end do

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


  end module exp_utils
