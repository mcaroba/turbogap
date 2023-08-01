! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, xps_utils.f90, is copyright (c) 2019-2023,
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

module xps_utils


  contains


    !**************************************************************************

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
      ! experimental distribution which would simplify the mathematics

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



    subroutine get_moment_spectra_energies_forces(data, energy_scales, core_electron_be, core_electron_be_der,&
         & n_neigh, neighbors_list, neighbor_species, &
         & sigma, n_samples, &
         & x_i_exp, y_i_exp,  y_i_pred, get_exp, do_forces, rjs, xyz, n_moments, moments, moments_exp,&
         & n_tot, energies_lp, forces0, virial)
      implicit none
      real*8, allocatable, intent(in) :: data(:,:)
      integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
      real*8, intent(in) :: sigma, core_electron_be(:), core_electron_be_der(:,:), rjs(:), xyz(:,:), energy_scales(:)
      real*8, allocatable, intent(inout) :: forces0(:,:)
      real*8, allocatable, intent(out) :: x_i_exp(:), y_i_exp(:), y_i_pred(:)
      real*8, intent(inout) :: energies_lp(:)
      real*8, intent(inout) :: virial(1:3,1:3)
      real*8, allocatable, intent(inout) :: moments_exp(:), moments(:)
      real*8 ::  this_force(1:3), mag, t, mean_reference
      real*8, allocatable :: x(:), y(:), x_exp(:), y_exp(:), y_all(:,:), y_der(:,:), int_der_moments(:,:), &
           x_int_der(:), der_diff_moments(:)
      integer, intent(in) :: n_samples, n_tot, n_moments
      logical, intent(in) :: get_exp, do_forces
      integer :: n_sites, n_pairs, n_pairs_soap, n_species, n_sites0
      integer :: i, j, i2, j2, k, k1, k2, mag_force_i, l
      character*32 :: mchar
      real*8 :: x_val, x_min, x_max, x_range, max_force, mag_force,&
           & max_mag_force, similarity, dx, e_tot, moment_force(1:3)

      n_sites = size(n_neigh)
      n_pairs = size(neighbors_list)
      n_sites0 = size(forces0, 2)


      ! This is essentially the same as the procedure for the hirshfeld gradients as in vdw.f90
      allocate(y_all(1:n_samples, 1:n_sites))
      allocate(y_der(1:3, 1:n_sites))
      allocate(x_exp(1:n_samples))
      allocate(y_exp(1:n_samples))

      call interpolate_data( x, y, data(1,:), data(2,:), n_samples, dx )

      x_exp = x
      y_exp = y
      y_der = 0.d0
      y_all = 0.d0

      ! We then broaden y to get the actual spectra
      y = 0.d0
      do i = 1, n_sites
         call broaden_spectrum(x, core_electron_be(i), y, y_all, i, sigma)
         !         y(1:n_samples) = y(1:n_samples) + y_all(1:n_samples,i)
      end do

      ! and then we normalise
      mag = sum(y * dx) ! real(n_tot) !sqrt( dot_product(y, y) )
      y = y / mag
      y_all = y_all / mag
      y_exp = y_exp / sum(y_exp * dx)  !sqrt( dot_product(y_exp, y_exp) )
      similarity = dot_product( y, y_exp * dx )

      print *, " Similarity: ", similarity

      ! This gets the interpolated spectra for both experiment and the prediction

      ! For the energy, just simply get the difference of the moments
      ! with that of the experiment and divide by the total number of
      ! atoms

      mean_reference =  dot_product( dx * x_exp, y_exp )
      if (get_exp)then
         call get_moments_of_distribution(x_exp, y_exp, dx, moments_exp, n_moments, mean_reference)
      end if
      call get_moments_of_distribution(x, y, dx, moments, n_moments, mean_reference)

      ! Now print out the moments for comparison
      e_tot = 0.d0
      write(*, *)  ' -- Moments of the distributions -- '
      do i = 1, n_moments
         write(mchar,'(A,I8,A)') 'mu_', i, " pred = "
         write(*, '(1X,A,1X,F20.8,1X,A,1X,F20.8)') mchar, moments(i),  " exp = ", moments_exp(i)
         e_tot = e_tot + energy_scales(i) * ( moments(i) - moments_exp(i) )**2
      end do

      energies_lp = e_tot / real(n_tot)


      ! Now get the forces !!
      ! First, we must build the function from the derivatives

      !   Compute xps forces
      if( do_forces )then
         max_mag_force = 0.0
         mag_force_i  = 0
         forces0 = 0.d0
         virial = 0.d0
         !     First, we compute the forces acting on all the "SOAP-neighbors" of atom i
         !     (including i itself) due to the gradients of i's Hirshfeld volume wrt the
         !     positions of its neighbors

         ! We have a contribution to the force from each of the moments
         ! We can precompute some of these quantities, such as the


         ! go over the n moments and get the corresponding integrals

         allocate(x_int_der(1:n_samples))
         allocate(int_der_moments(1:n_moments,1:n_sites))
         allocate(der_diff_moments(1:n_moments))

         x_int_der = 1.d0
         do i = 1, n_moments
            der_diff_moments(i) = 2.d0 * energy_scales(i) * ( moments(i) - moments_exp(i) )

            if (i == 1) then
               x_int_der = x_int_der * x
            else
               x_int_der = (x - mean_reference)**(real(i))
            end if

            do j = 1, size(core_electron_be)
               int_der_moments(i,j) = dx * dot_product( (x - core_electron_be(j)) / (sigma**2), &
                                                               & x_int_der * y_all(1:n_samples, j) )
               ! we multiply this by the gradient term when we calculate it, y_all is the individual spectrum
            end do
         end do

         ! Now we get the forces from this

         k = 0
         do i = 1, n_sites
            i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
            do j = 1, n_neigh(i)
               k = k + 1
               j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
               !             SOAP neighbors
               ! MAY NEED TO CHANGE THIS TO ACCOUNT FOR MACHINE PRECISIOv
               if( .not. all(core_electron_be_der(1:3, k) == 0.d0) )then
                  ! force comes from a sum of all the moment forces
                  this_force = 0.d0
                  do l = 1, n_moments
                     moment_force(1:3) = core_electron_be_der(1:3, k) * der_diff_moments(l) * int_der_moments(l,i)
                     this_force = this_force + moment_force
                  end do

                  forces0(1:3, j2) = forces0(1:3, j2) + - this_force(1:3)

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
               if( j2 /= i2 )then
                  forces0(1:3,i2) = forces0(1:3,i2) + this_force(1:3)
               end if
               !           ... but the periodic replicas DO contribute to the virial
               !           Sign is minus because this force is acting on i2. Factor of 1/2 is because this is
               !           derived from a pair energy
               !            virial = virial - 0.5d0 * dot_product(this_force(1:3), xyz(1:3,k))
               do k1 = 1, 3
                  do k2 =1, 3
                     virial(k1, k2) = virial(k1, k2) - 0.25d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
                  end do
               end do
            end do
         end do

         deallocate(x_int_der)
         deallocate(int_der_moments)
         deallocate(der_diff_moments)

         print *, "MAX FORCE FROM LOCAL PROP ", forces0(1:3,mag_force_i)
      end if

      allocate(x_i_exp(1:n_samples))
      allocate(y_i_exp(1:n_samples))
      allocate(y_i_pred(1:n_samples))

      x_i_exp = x
      y_i_exp = y_exp
      y_i_pred = y

      deallocate(x)
      deallocate(y)
      deallocate(y_all)
      deallocate(y_der)
      deallocate(x_exp)
      deallocate(y_exp)


    end subroutine get_moment_spectra_energies_forces


    

    subroutine get_exp_pred_spectra_energies_forces(data, energy_scale, core_electron_be, core_electron_be_der,&
         & n_neigh, neighbors_list, neighbor_species, &
         & sigma, n_samples, mag, &
         & x, y_exp,  y, y_all, get_exp, do_forces, rjs, xyz,&
         & n_tot, energies_lp, forces0, virial, similarity_type, rank)
      implicit none
      real*8, allocatable, intent(in) :: data(:,:)
      integer, intent(in) :: n_neigh(:), neighbors_list(:), neighbor_species(:)
      real*8, intent(in) :: sigma, core_electron_be(:),&
           & core_electron_be_der(:,:), rjs(:), xyz(:,:),&
           & energy_scale, mag
      real*8, allocatable, intent(inout) :: forces0(:,:)
      real*8, allocatable, intent(in) :: x(:), y_exp(:), y(:)
      real*8, intent(in) :: y_all(:,:)
      real*8, allocatable :: der_sum(:)
      real*8, intent(inout) :: energies_lp(:)
      real*8, intent(inout) :: virial(1:3,1:3)
      real*8 ::  this_force(1:3), t
      real*8, allocatable ::  y_der(:,:), prefactor(:)
      integer, intent(in) :: n_samples, n_tot, rank
      logical, intent(in) :: get_exp, do_forces
      integer :: n_sites, n_pairs, n_pairs_soap, n_species, n_sites0
      integer :: i, j, i2, j2, k, n_in_buffer, k1, k2, mag_force_i
      real*8 :: x_val, x_min, x_max, x_range, max_force, mag_force, max_mag_force, similarity, dx
      character*32, intent(in) :: similarity_type


      n_sites = size(n_neigh)
      n_pairs = size(neighbors_list)
      n_sites0 = size(forces0, 2)

      dx = x(2) - x(1)
      ! This is essentially the same as the procedure for the hirshfeld gradients as in vdw.f90
      allocate(y_der(1:3, 1:n_samples))


      ! The energetic contribution to this would then be, e_scale * ( 1 - \sum_i S_i ), this can be done after mpi_reduce
      ! S_i =  dot_product( y_all(:,i), y_exp )

      if( similarity_type == 'similarity' .or. similarity_type == 'overlap' )then
         do i = 1, n_sites
            energies_lp(i) = - energy_scale * ( dot_product( y_all(:,i), y_exp * dx))
         end do
      else if (similarity_type == 'lsquares')then
         energies_lp = + energy_scale / n_sites0 * ( dx * dot_product( (y - y_exp), (y - y_exp) ) )

         ! Get the other terms for the lsquares expression
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
         !     First, we compute the forces acting on all the "SOAP-neighbors" of atom i
         !     (including i itself) due to the gradients of i's Hirshfeld volume wrt the
         !     positions of its neighbors

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
                  call broaden_spectrum_derivative(x(1:n_samples),&
                       & core_electron_be(i),&
                       & core_electron_be_der(1:3, k),&
                       & y_der(1:3,1:n_samples ), sigma, y_exp(1:n_samples), this_force(1:3), mag, dx)

                  if( similarity_type == 'lsquares' )then

                     y_der(1, 1:n_samples) = y_der(1, 1:n_samples) * prefactor(1:n_samples)
                     y_der(2, 1:n_samples) = y_der(2, 1:n_samples) * prefactor(1:n_samples)
                     y_der(3, 1:n_samples) = y_der(3, 1:n_samples) * prefactor(1:n_samples)

                     this_force(1) =  - sum( y_der(1,:) ) * dx
                     this_force(2) =  - sum( y_der(2,:) ) * dx
                     this_force(3) =  - sum( y_der(3,:) ) * dx

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
               if( j2 /= i2 )then
                  forces0(1:3,i2) = forces0(1:3,i2) + this_force(1:3)
               end if
               !           ... but the periodic replicas DO contribute to the virial
               !           Sign is minus because this force is acting on i2. Factor of 1/2 is because this is
               !           derived from a pair energy
               !            virial = virial - 0.5d0 * dot_product(this_force(1:3), xyz(1:3,k))
               do k1 = 1, 3
                  do k2 =1, 3
                     virial(k1, k2) = virial(k1, k2) - 0.25d0 * (this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
                  end do
               end do
            end do
         end do

         print *, rank, "MAX FORCE FROM LOCAL PROP ", forces0(1:3,mag_force_i)
         deallocate(y_der)
         if(allocated(prefactor)) deallocate(prefactor)

      end if



    end subroutine get_exp_pred_spectra_energies_forces





    subroutine compare_exp_to_pred_spectra(data, core_electron_be, &
         & sigma, n_samples, mag, sim_exp_pred, &
         & x_i_exp, y_i_exp, x_i_pred, y_i_pred, y_i_pred_all, &
         & get_exp, similarity_type )
      implicit none
      real*8, allocatable, intent(in) :: data(:,:)
      real*8, intent(in) :: sigma, core_electron_be(:)
      real*8, allocatable, intent(inout) :: x_i_exp(:), y_i_exp(:), x_i_pred(:),&
           & y_i_pred(:), y_i_pred_all(:,:)
      real*8, intent(out) :: sim_exp_pred, mag
      integer, intent(in) :: n_samples
      logical, intent(in) :: get_exp
      character*32, intent(in) :: similarity_type

      if (get_exp)then
         ! Get interpolated experimental spectra
         call get_xps_spectra(data(1,:), data(2,:), sigma, n_samples, mag,&
              & x_i_exp,  y_i_exp, y_i_pred_all, core_electron_be,&
              &  .false.)
      end if

      ! Get interpolated predicted spectra
      call get_xps_spectra(data(1,:), data(2,:), sigma, n_samples, mag,&
           & x_i_pred, y_i_pred, y_i_pred_all, core_electron_be,&
           & .true.)

      ! Now get the similarity
      if (similarity_type == "lsquares")then
         sim_exp_pred = - (x_i_exp(2) - x_i_exp(1) ) * dot_product(y_i_exp - y_i_pred, y_i_exp - y_i_pred)
      else
         sim_exp_pred =  dot_product(y_i_exp, y_i_pred)
      end if

    end subroutine compare_exp_to_pred_spectra


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

      mag = sqrt(dot_product(y,y)) !sum(y * dx) !sqrt(dot_product(y, y))
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
         f = - ( ( x(i) - x0 ) / (sigma**2) ) * exp( -( x(i) - x0 )**2 / (2*sigma**2) ) * y_exp(i) / mag
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

      do i = 1, size(x)
         idx = floor( ((x(i)-x_min) / x_range) * real(size(xi)) ) + 1

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

    end subroutine lerp


  end module xps_utils
