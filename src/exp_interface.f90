! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, exp_interface.f90, is copyright (c) 2019-2023,
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

module exp_interface

   use kinds
   use types
   use read_files
#ifdef _MPIF90
   use mpi
#endif
   use exp_utils
   use soap_turbo_functions
#ifdef _GPU
   use F_B_C
   use iso_c_binding
#endif
!  The batched pair distribution traces every device allocation it makes. That
!  is per allocation, per batch, per step, and unbuffered, which costs more than
!  the work it describes; as a parameter the compiler removes it entirely. Set
!  it to .true. to get the trace back.
   logical, parameter :: debug_gpu_batches = .false.

contains

   ! This module implements the interfaces for the gradient of experimental functions

   subroutine get_write_condition(do_mc, do_md, mc_istep, md_istep, write_xyz, write_condition)
      implicit none
      logical, intent(in) :: do_mc
      logical, intent(in) :: do_md
      integer, intent(in) :: mc_istep
      integer, intent(in) :: md_istep
      integer, intent(in) :: write_xyz
      logical, intent(out) :: write_condition

      if (do_mc) then
         write_condition = mc_istep > -1 .and. (modulo(mc_istep, write_xyz) == 0)
      elseif (do_md) then
         write_condition = md_istep > -1 .and. (modulo(md_istep, write_xyz) == 0)
      else
         write_condition = .true.
      end if
   end subroutine get_write_condition

   subroutine get_overwrite_condition(do_mc, do_md, mc_istep, md_istep, write_xyz, write_condition)
      implicit none
      logical, intent(in) :: do_mc
      logical, intent(in) :: do_md
      integer, intent(in) :: mc_istep
      integer, intent(in) :: md_istep
      integer, intent(in) :: write_xyz
      logical, intent(out) :: write_condition

      if (do_mc) then
         write_condition = mc_istep == 0
      elseif (do_md) then
         write_condition = md_istep == 0
      else
         write_condition = .true.
      end if
   end subroutine get_overwrite_condition

   subroutine write_partial_exp(do_mc, do_md, mc_istep, md_istep,&
        & write_xyz, has_partials, n_species, n_samples, n_dim_partial&
        &, x, y, partials, species_types, name)
      implicit none
      logical, intent(in) :: do_mc
      logical, intent(in) :: do_md
      logical, intent(in) :: has_partials
      integer, intent(in) :: mc_istep
      integer, intent(in) :: md_istep
      integer, intent(in) :: write_xyz
      integer, intent(in) :: n_species
      integer, intent(in) :: n_samples
      integer, intent(in) :: n_dim_partial
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: partials(:, :)
      character*8, allocatable :: species_types(:)
      character(len=*), intent(in) :: name
      integer :: j
      integer :: k
      integer :: n_dim_idx
      character*1024 :: filename
      logical :: write_condition
      logical :: overwrite_condition

      call get_overwrite_condition(do_mc, do_md&
           &, mc_istep, md_istep, write_xyz,&
           & overwrite_condition)

      if (has_partials) then
         n_dim_idx = 1
         outer: do j = 1, n_species
            do k = 1, n_species

               if (j > k) cycle

               write (filename, '(A)')&
                    & name//'_'//trim(species_types(j))//'_'//trim(species_types(k))//&
                    & "_prediction.dat"
               call write_exp_datan(x(1:n_samples),&
                    & partials(1:n_samples, n_dim_idx),&
                    & overwrite_condition, filename, name)

               n_dim_idx = n_dim_idx + 1
               if (n_dim_idx > n_dim_partial) then
                  exit outer
               end if

            end do
         end do outer
      end if

      write (filename, '(A)')&
           & name//"_total.dat"
      call write_exp_datan(x(1:n_samples),&
           &y(1:n_samples),&
           & overwrite_condition, filename, name)
   end subroutine write_partial_exp

   subroutine get_pdf_sf_xrd_explicitly_kde(v_uc, n_sites0, n_species, species, species_types, &
        & neighbors_list, n_neigh, neighbor_species, rjs, xyz, r_cut, &
        & r_min, r_max, n_samples,  &!       & q_min, q_max, n_samples_sf, x_sf,  &
        & pair_distribution, pair_distribution_der, kde_sigma, do_derivatives, rank)
      ! & do_forces_pdf, do_forces_sf, do_forces_xrd, &
      ! &    forces_pdf,    forces_sf,    forces_xrd )
      implicit none
      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: kde_sigma
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(in) :: v_uc
      real(dp), intent(in) :: r_min
      real(dp), intent(in) :: r_max
      real(dp), intent(in) :: r_cut
      integer, intent(in) :: neighbors_list(:)
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: neighbor_species(:)
      integer, intent(in) :: species(:)
      integer, intent(in) :: n_sites0
      integer, intent(in) :: n_samples
      integer, intent(in) :: n_species
      integer, intent(in) :: rank
      logical, intent(in) :: do_derivatives !
      character*8, intent(in), allocatable :: species_types(:)
      real(dp), intent(out), allocatable :: pair_distribution(:, :)
      real(dp), intent(out), allocatable :: pair_distribution_der(:, :, :)!

      ! Internal  Variables
      integer :: n_sites
      integer :: n_pairs
      integer :: s
      integer :: species_i
      integer :: species_j
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: i2
      integer :: j2
      integer :: n_dim_partial
      integer :: n_dim_idx
      integer :: ierr
      real(dp) :: r
      real(dp) :: gauss
      real(dp) :: c
      real(dp), allocatable :: bin_edges(:)
      real(dp), allocatable :: dV(:)
      real(dp), allocatable :: factors(:)
      real(dp), allocatable :: pair_distribution_partial_temp(:, :)
      real(dp), allocatable :: pdf(:)
      real(dp), allocatable :: sf(:)
      real(dp), allocatable :: xrd(:)
      real(dp), allocatable :: forces_sf(:, :)
      real(dp), allocatable :: forces_xrd(:, :)
      real(dp), allocatable :: n_atoms_of_species(:)
      real(dp), allocatable :: prefactor_pdf(:)
      real(dp), allocatable :: prefactor_sf(:)
      real(dp), allocatable :: prefactor_xrd(:)
      real(dp), allocatable :: x_pdf(:)

      integer, allocatable :: species_to_ndim(:, :)
      integer, allocatable :: species1_partial(:)
      integer, allocatable :: species2_partial(:)
      character*1024 :: filename
      ! Parameters
      real(dp), parameter :: pi = acos(-1.0)

      n_sites = size(n_neigh)
      n_pairs = size(neighbors_list)

      allocate (n_atoms_of_species(1:n_species))
      n_atoms_of_species = 0.d0
      do i = 1, size(species, 1)
         s = species(i)
         n_atoms_of_species(s) = n_atoms_of_species(s) + 1
      end do

      allocate (bin_edges(1:n_samples + 1))
      allocate (dV(1:n_samples))
      allocate (x_pdf(1:n_samples))

      n_dim_partial = n_species*(n_species + 1)/2
      allocate (factors(1:n_dim_partial))
      allocate (species_to_ndim(1:n_species, 1:n_species))

      allocate (species1_partial(1:n_dim_partial))
      allocate (species2_partial(1:n_dim_partial))

      allocate (pair_distribution(1:n_samples, 1:n_dim_partial))
      pair_distribution = 0.d0

      if (do_derivatives) then
         allocate (pair_distribution_der(1:n_pairs, 1:n_samples, 1:n_dim_partial))
         pair_distribution_der = 0.d0
      end if

      n_dim_idx = 1
      outer: do i = 1, n_species
         do j = 1, n_species
            if (i > j) cycle

            if (i /= j) then
               factors(n_dim_idx) = 2.d0
            else
               factors(n_dim_idx) = 1.d0
            end if

            species_to_ndim(i, j) = n_dim_idx
            species_to_ndim(j, i) = n_dim_idx

            species1_partial(n_dim_idx) = i
            species2_partial(n_dim_idx) = j

            n_dim_idx = n_dim_idx + 1

            if (n_dim_idx > n_dim_partial) exit outer

         end do
      end do outer

      x_pdf = 0.d0
      bin_edges = 0.d0
      pair_distribution = 0.d0

      do i = 1, n_samples + 1
         bin_edges(i) = r_min + (real((i - 1))/real(n_samples))*(r_max - r_min)
      end do

      do i = 1, n_samples
         x_pdf(i) = (bin_edges(i) + bin_edges(i + 1))/2.d0
         dV(i) = 4.d0*pi*(bin_edges(i)**2*(bin_edges(i + 1) - bin_edges(i)))
      end do

      k = 0
      do i = 1, n_sites
         i2 = modulo(neighbors_list(k + 1) - 1, n_sites0) + 1
         species_i = neighbor_species(k + 1)
         k = k + 1
         do j = 2, n_neigh(i)
            k = k + 1
            j2 = modulo(neighbors_list(k) - 1, n_sites0) + 1
            species_j = neighbor_species(k)
            r = rjs(k) ! atom pair distance

            if (r > r_cut .or. r < r_min .or. r > r_max + kde_sigma*6.d0) cycle

            n_dim_idx = species_to_ndim(species_i, species_j)

            do l = 1, n_samples
               gauss = exp(-((x_pdf(l) - r)/kde_sigma)**2/2.d0)
               pair_distribution(l, n_dim_idx) = pair_distribution(l, n_dim_idx) + gauss

               ! Construct the initial derivatives here without the ri-rj term here
               if (do_derivatives) then
                  pair_distribution_der(k, l, n_dim_idx) = gauss*((x_pdf(l) - r)/kde_sigma**2)/r
               end if

            end do
         end do
      end do

      pair_distribution = pair_distribution*((r_max - r_min)/dfloat(n_samples)) &
           & /(sqrt(2.d0*pi)*kde_sigma)

      if (do_derivatives) pair_distribution_der = pair_distribution_der*((r_max - r_min)/dfloat(n_samples)) &
           & /(sqrt(2.d0*pi)*kde_sigma)

      do i = 1, n_dim_partial
         pair_distribution(1:n_samples, i) = pair_distribution(&
              & 1:n_samples, i)*v_uc/n_atoms_of_species(&
              & species1_partial(i))/n_atoms_of_species(&
              & species2_partial(i))/factors(i)

         pair_distribution(1:n_samples, i) = pair_distribution(1:n_samples, i)/dV

         if (do_derivatives) then
            pair_distribution_der(1:n_pairs, 1:n_samples, i) = pair_distribution_der(&
                 & 1:n_pairs, 1:n_samples, i)*v_uc/n_atoms_of_species(&
                 & species1_partial(i))/n_atoms_of_species(&
                 & species2_partial(i))/factors(i)

            do l = 1, n_pairs
               pair_distribution_der(l, 1:n_samples, i) = pair_distribution_der(l, 1:n_samples, i)/dV
            end do
         end if

      end do

#ifdef _MPIF90
      allocate (pair_distribution_partial_temp(1:n_samples, 1:n_dim_partial))
      pair_distribution_partial_temp = 0.d0

      call mpi_reduce(pair_distribution,&
           & pair_distribution_partial_temp, n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
           & 0, MPI_COMM_WORLD, ierr)

      pair_distribution = pair_distribution_partial_temp
      deallocate (pair_distribution_partial_temp)

      call mpi_bcast(pair_distribution, n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
           & MPI_COMM_WORLD, ierr)
#endif

      allocate (pdf(1:n_samples))
      pdf = 0.d0

      do i = 1, n_dim_partial
         c = (n_atoms_of_species(species1_partial(i))*&
              & n_atoms_of_species(species1_partial(i)))/&
              & dfloat(n_sites0)/dfloat(n_sites0)

         pdf(1:n_samples) = pdf(1:n_samples) + c*factors(i)*pair_distribution(1:n_samples, i)
      end do

      ! Write out the partial pair distribution functions
      if (rank == 0) then
         n_dim_idx = 1
         outer3: do j = 1, n_species
            do k = 1, n_species

               if (j > k) cycle

               write (filename, '(A)')&
                    & 'tpair_distribution_'//trim( &
                    & species_types(j))//'_'//trim( &
                    & species_types(k))//&
                    & "_prediction.dat"
               call write_exp_datan(x_pdf(1:n_samples),&
                    & pair_distribution(1:n_samples, n_dim_idx),&
                    & .true., filename, 'pair_distribution')

               n_dim_idx = n_dim_idx + 1
               if (n_dim_idx > n_dim_partial) then
                  exit outer3
               end if

            end do
         end do outer3

         write (filename, '(A)')&
              & "tpair_distribution_total.dat"
         call write_exp_datan(x_pdf(1:n_samples),&
              &pdf(1:n_samples),&
              & .true., filename, "pair_distribution")

      end if

      ! Now we have the pdf, we can construct the derivatives for the
      ! forces pdf only and we can construct the structure factor / xrd
      ! spectrum

      !    if (do_forces_pdf) allocate(prefactor_pdf(1:n_samples))
      !    if (do_forces_sf ) allocate(prefactor_sf(1:n_samples_sf))
      !    if (do_forces_xrd) allocate(prefactor_xrd(1:n_samples_sf))

      !    k = 0
      !    do i = 1, n_sites
      !       i2 = modulo(neighbors_list(k+1)-1, n_sites0) + 1
      !       species_i = neighbor_species(k+1)
      !       do j = 1, n_neigh(i)
      !          k = k + 1
      !          j2 = modulo(neighbors_list(k)-1, n_sites0) + 1
      !          species_j = neighbor_species(k)
      !          r = rjs(k) ! atom pair distance

      !          if ( r > r_cut .or. r < r_min .or. r > r_max + kde_sigma*6.d0  ) cycle

      !          n_dim_idx = species_to_ndim(species_i, species_j)

      !          if ( .not. all( xyz( 1:3, k ) == 0.d0 ) )then
      !             ! Actual derivative of the pair distribution function given here, no normalisation needed I think

      deallocate (species_to_ndim)
      deallocate (species1_partial)
      deallocate (species2_partial)
      deallocate (factors)
      deallocate (dV)
      deallocate (x_pdf)
      deallocate (bin_edges)
      deallocate (n_atoms_of_species)

   end subroutine get_pdf_sf_xrd_explicitly_kde

   ! subroutine get_this_exp_force(k, xyz, n_samples, n_dim_idx, pair_distribution_der, energy_scale, f, this_force)
   !   implicit none
   !   integer, intent(in) :: k, n_dim_idx, n_samples
   !   real, intent(in), allocatable  :: prefactor(:), pair_distribution_der(:,:,:)
   !   real, intent(in) :: rij(1:3), f, energy_scale
   !   real, intent(out) :: this_force(1:3)

   !   this_force(1) = dot_product( - 2.d0 * rij( 1 ) *&
   !        & pair_distribution_der(k, 1:n_samples,  n_dim_idx),&
   !        & prefactor(1:n_samples))

   !   this_force(2) = dot_product( - 2.d0 * rij( 2 ) *&
   !        & pair_distribution_der(k, 1:n_samples, n_dim_idx),&
   !        & prefactor(1:n_samples))

   !   this_force(3) = dot_product( - 2.d0 * rij( 3 ) *&
   !        & pair_distribution_der(k, 1:n_samples,  n_dim_idx),&
   !        & prefactor(1:n_samples))

   !   this_force(1:3) =  - f * energy_scale  * this_force(1:3)

   subroutine preprocess_exp_data(params, x, y, label, n_sites, V, input, output, exp)
      implicit none
      type(input_parameters), intent(in) :: params
      real(dp), intent(in), allocatable :: x(:)
      real(dp), intent(in) :: V
      real(dp), intent(inout), allocatable :: y(:)
      integer, intent(in) :: n_sites
      character*1024, intent(in) :: label
      real(dp), parameter :: pi = acos(-1.0)
      real(dp) :: mag
      real(dp) :: dx
      real(dp) :: rho
      logical, intent(in) :: exp
      character*32, intent(inout) :: output
      character*1024, intent(inout) :: input

      dx = x(2) - x(1)
      rho = dfloat(n_sites)/V

      if (trim(label) == "xps") then
         ! calculate the magnitude and normalize
         mag = sqrt(dot_product(y, y))
         y = y/mag
         output = "xps"

      elseif (trim(label) == "pair_distribution") then
         output = params%pair_distribution_output
         if (trim(params%pair_distribution_output) == "D(r)" .and. .not. (exp .and. trim(input) == "D(r)")) then
            ! D(r) = 4pi rho * r * G(r)
            ! G(r) = total pair distribution function
            y = 4.d0*pi*rho*x*(y - 1.d0)

         end if
      elseif (trim(label) == "xrd") then
         output = params%xrd_output

         if (trim(params%xrd_output) == "q*i(q)" .and. params%q_units &
              &== "q") then
            if (exp .and. (trim(input) == "i(q)" .or. trim(input) == "F(q)")) then
               y = x*(y - 1.d0)
            end if
         end if

!      The inverse: experimental data supplied as q*i(q) when i(q) is the
!      requested output.
         if (trim(params%xrd_output) == "i(q)" .and. params%q_units &
              &== "q") then
            if (exp .and. (trim(input) == "q*i(q)" .or. trim(input) == "q*F(q)")) then
               y = (y/x) + 1.d0
            end if
         end if
      elseif (trim(label) == "nd") then
         output = params%nd_output

         if (trim(params%nd_output) == "q*i(q)" .and. params%q_units &
              &== "q") then
            if (exp .and. (trim(input) == "i(q)" .or. trim(input) == "F(q)")) then
               y = x*(y - 1.d0)
            end if
         end if

!      The inverse: experimental data supplied as q*i(q) when i(q) is the
!      requested output.
         if (trim(params%nd_output) == "i(q)" .and. params%q_units &
              &== "q") then
            if (exp .and. (trim(input) == "q*i(q)" .or. trim(input) == "q*F(q)")) then
               y = (y/x) + 1.d0
            end if
         end if

      end if
   end subroutine preprocess_exp_data

   subroutine calculate_pair_distribution(params, x_pair_distribution&
        &, y_pair_distribution, y_pair_distribution_temp,&
        & pair_distribution_partial, pair_distribution_partial_temp, &
        & n_species, species_types, n_atoms_of_species, n_sites, a_box, b_box, c_box,&
        & indices, md_istep, mc_istep, i_beg, i_end, j_beg, j_end, ierr, rjs, xyz, &
        & neighbors_list, n_neigh, neighbor_species, species, rank,&
        & do_derivatives, pair_distribution_der, pair_distribution_partial_der,&
        & pair_distribution_partial_temp_der, energies_pair_distribution, forces_pair_distribution, virial)
      implicit none
      type(input_parameters), intent(inout) :: params
      real(dp), allocatable, intent(out) :: x_pair_distribution(:)
      real(dp), allocatable, intent(out) :: y_pair_distribution(:)
      real(dp), allocatable, intent(out) :: pair_distribution_partial(:, :)
      real(dp), allocatable, intent(out) :: n_atoms_of_species(:)
      real(dp), allocatable, intent(out) :: pair_distribution_partial_temp(:, :)
      real(dp), allocatable, intent(out) :: y_pair_distribution_temp(:)
      real(dp), allocatable, intent(out) :: pair_distribution_der(:, :)
      real(dp), allocatable, intent(out) :: pair_distribution_partial_der(:, :, :)
      real(dp), allocatable, intent(out) :: pair_distribution_partial_temp_der(:, :, :)
      real(dp), allocatable, intent(out) :: energies_pair_distribution(:)
      real(dp), allocatable, intent(out) :: forces_pair_distribution(:, :)
      character*8, allocatable, intent(in) :: species_types(:)
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: species(:)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(inout) :: virial(1:3, 1:3)
      real(dp) :: v_uc
      real(dp) :: f
!     -V dE/dV, the cell half of this observable's virial.
      real(dp) :: dedv
      integer, intent(in) :: n_species
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: rank
      integer, intent(inout) :: ierr
      real(dp), allocatable :: factors(:)
      real(dp), allocatable :: pair_distribution_der_temp(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: i2
      integer :: n_dim_partial
      integer :: n_dim_idx
      logical, intent(in) :: do_derivatives
      real(dp), parameter :: pi = acos(-1.0)
      logical :: write_condition
      logical :: overwrite_condition
      character*1024 :: filename

      ! Things that are allocated here:
      ! Always:
      !  > x_pair_distribution
      !  > y_pair_distribution
      ! if pair_distribution_partial == .true.
      !  > pair_distribution_partial( n_samples, n_spec * (n_spec + 1)/2 )
      !  if do_derivatives == .true.
      !    > pair_distribution_partial_der( n_samples, n_spec * (n_spec + 1)/2, j_beg : j_end )

      ! first allocate the necessary arrays for the
      ! calculation of the pair correlation function
      if (allocated(x_pair_distribution)) deallocate (x_pair_distribution)
      if (allocated(y_pair_distribution)) deallocate (y_pair_distribution)

      allocate (x_pair_distribution(1:params%pair_distribution_n_samples))
      allocate (y_pair_distribution(1:params%pair_distribution_n_samples))

      if (params%n_exp > 0) then
         do i = 1, params%n_exp
            if (trim(params%exp_data(i)%label) == 'pair_distribution') then
               x_pair_distribution = params%exp_data(i)%x
            end if
         end do
      end if

      if (params%pair_distribution_partial) then
         n_dim_partial = n_species*(n_species + 1)/2
         allocate (factors(1:n_dim_partial))

         n_dim_idx = 1
         outer: do i = 1, n_species
            do j = 1, n_species
               if (i > j) cycle

               if (i /= j) then
                  factors(n_dim_idx) = 2.d0
               else
                  factors(n_dim_idx) = 1.d0
               end if

               n_dim_idx = n_dim_idx + 1
               if (n_dim_idx > n_dim_partial) then
                  exit outer
               end if

            end do
         end do outer

         if (.not. allocated(pair_distribution_partial)) then   !deallocate(pair_distribution_partial)
            allocate (pair_distribution_partial(1:params%pair_distribution_n_samples,&
                 & 1:n_dim_partial))
         end if

         pair_distribution_partial = 0.d0

         if (params%do_forces .and. params%exp_forces) then
            allocate (pair_distribution_partial_der(1:params%pair_distribution_n_samples,&
              & 1:n_dim_partial, j_beg:j_end))
            pair_distribution_partial_der = 0.d0

            if (rank == 0 .and. md_istep == 0) write (*, '(A,1X,F7.4,1X,A)') "Gb/core: partial pdfder = ", dfloat(params&
                 &%pair_distribution_n_samples*n_dim_partial*j_end)&
                 & *8.d0/(dfloat(1024*1024*1024)), " Gb  |"
            if (rank == 0 .and. md_istep == 0) write (*, *) '                                       |'

         end if
      else
         if (params%do_forces .and. params%exp_forces) then
            allocate (pair_distribution_partial_der(1:params%pair_distribution_n_samples, 1:1, &
              &  j_beg:j_end))
            pair_distribution_partial_der = 0.d0
         end if

      end if

      if (allocated(n_atoms_of_species)) deallocate (n_atoms_of_species)
      allocate (n_atoms_of_species(1:n_species))

      do j = 1, n_species
         n_atoms_of_species(j) = 0.d0
         do i2 = 1, n_sites
            if (species(i2) == j) then
               n_atoms_of_species(j) = n_atoms_of_species(j) + 1.d0
            end if
         end do
      end do

#ifdef _MPIF90
      if (params%pair_distribution_partial) then
         allocate (pair_distribution_partial_temp(1:params%pair_distribution_n_samples, 1:n_dim_partial))

         pair_distribution_partial_temp = 0.0d0
      end if

      allocate (y_pair_distribution_temp(1:params%pair_distribution_n_samples))
      y_pair_distribution_temp = 0.d0

#endif
      v_uc = dot_product(cross_product(a_box,&
           & b_box), c_box)/(&
           & dfloat(indices(1)*indices(2)&
           &*indices(3)))

      !###---   Calculating the partial pair distribution functions   ---###!

      if (params%pair_distribution_partial) then
         n_dim_idx = 1
         outer1: do j = 1, n_species
            do k = 1, n_species

               if (j > k) cycle ! We have already calculated the pair correlation function!

               ! Note that with the calculation of the derivatives here,
               ! this is without the -2 * delta_ik (r_j^alpha - r_i^alpha)
               ! factor, which allows for some freeing of memory

               call get_pair_distribution(n_sites, &
                    & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                    & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                    &%r_range_min, params%r_range_max, params%pair_distribution_n_samples, x_pair_distribution,&
                    & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), params&
                    &%pair_distribution_rcut, .false.,&
                    & params%pair_distribution_partial, j, k,&
                    & params%pair_distribution_kde_sigma,&
                    & dfloat(n_sites)/v_uc, params%exp_forces,&
                    & pair_distribution_partial_der, n_dim_idx, &
                    & j_beg, j_end)

               n_dim_idx = n_dim_idx + 1

               if (n_dim_idx > n_dim_partial) then
                  exit outer1
               end if

            end do
         end do outer1

      else
         call get_pair_distribution(n_sites, &
              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3, j_beg:j_end),&
              & params%r_range_min, params%r_range_max, params &
              &%pair_distribution_n_samples, x_pair_distribution,&
              & y_pair_distribution, params &
              &%pair_distribution_rcut, .false., .false., 1, 1,&
              & params%pair_distribution_kde_sigma, dfloat(n_sites)&
              &/v_uc, params%do_forces .and. params%exp_forces, pair_distribution_partial_der, 1, &
              & j_beg, j_end)
      end if

      ! --- MPI communication is here  ---

      if (params%pair_distribution_partial) then
#ifdef _MPIF90
         call mpi_reduce(pair_distribution_partial,&
              & pair_distribution_partial_temp, params&
              &%pair_distribution_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
              & 0, MPI_COMM_WORLD, ierr)

         ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
         ! Note, we have only so far divided by 4 pi r^2 dr
         ! Therefore, we must scale by the density

         pair_distribution_partial = pair_distribution_partial_temp
         deallocate (pair_distribution_partial_temp)

         call mpi_bcast(pair_distribution_partial, params&
              &%pair_distribution_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
              & MPI_COMM_WORLD, ierr)

         ! Now, we have the derivatives of the partial pair distribution
         ! function with respect to the atom pairs in that rank
         !
         ! We can keep them in the rank and calculate forces

         ! call mpi_reduce(pair_distribution_partial_der,&
         !      & pair_distribution_partial_der_temp, params&
         !      &%pair_distribution_n_samples * n_species *&
         !      & n_species * 3 * n_pairs_tot, MPI_DOUBLE_PRECISION, MPI_SUM,&
         !      & 0, MPI_COMM_WORLD, ierr)

         ! ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
         ! ! Note, we have only so far divided by 4 pi r^2 dr
         ! ! Therefore, we must scale by the density

#endif

         if (params%valid_pdf) then
            allocate (energies_pair_distribution(1:n_sites))
            energies_pair_distribution = 0.d0

            if (params%do_forces .and. params%exp_forces) then
               allocate (forces_pair_distribution(1:3, 1:n_sites))
               forces_pair_distribution = 0.d0
            end if

         end if

         !###---   Accumulate the PDF   ---###!

         y_pair_distribution = 0.d0
         n_dim_idx = 1
         outer2: do j = 1, n_species
            do k = 1, n_species

               if (j > k) cycle

               pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) =&
                    & pair_distribution_partial(1:params&
                    &%pair_distribution_n_samples, n_dim_idx)*v_uc &
                    &  /n_atoms_of_species(j)/n_atoms_of_species(k)/factors(n_dim_idx) !real(n_sites)

               if (params%do_forces .and. params%exp_forces) then
                  pair_distribution_partial_der(1:params&
                       &%pair_distribution_n_samples, n_dim_idx, &
                       & j_beg:j_end) = pair_distribution_partial_der(1:params&
                       & %pair_distribution_n_samples, n_dim_idx, &
                       & j_beg:j_end)*v_uc/n_atoms_of_species(j)/ &
                       & n_atoms_of_species(k)/factors(n_dim_idx)!real(n_sites)

               end if

               y_pair_distribution(1:params%pair_distribution_n_samples) = &
                    & y_pair_distribution(1:params%pair_distribution_n_samples) +  &
                    &  factors(n_dim_idx)*(n_atoms_of_species(j)*n_atoms_of_species(k))* &
                    & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) &
                    &  /dfloat(n_sites)/dfloat(n_sites)

               n_dim_idx = n_dim_idx + 1

               if (n_dim_idx > n_dim_partial) then
                  exit outer2
               end if

            end do
         end do outer2

         ! --- Preprocess the pair distribution according to the output --- !
         if (trim(params%pair_distribution_output) == "D(r)") then
            y_pair_distribution = 4.d0*pi*(dfloat(n_sites)/v_uc)*x_pair_distribution*(y_pair_distribution - 1.d0)
         end if

         !###---   Calculate the forces   ---###!

         if (params%valid_pdf .and. allocated(params%exp_energy_scales)) then

            call get_energy_scale(params%do_md, params%do_mc,&
                 & md_istep, params%md_nsteps, mc_istep, params&
                 &%mc_nsteps, params &
                 &%exp_energy_scales_initial(params%pdf_idx), params &
                 &%exp_energy_scales_final(params%pdf_idx), params &
                 &%exp_energy_scales(params%pdf_idx))

            call get_exp_energies(params%exp_energy_scales(params&
                 &%pdf_idx), params%exp_data(params%pdf_idx)%y&
                 &, y_pair_distribution,&
                 & params%pair_distribution_n_samples, n_sites,&
                 & energies_pair_distribution(i_beg:i_end), params%exp_data(params%pdf_idx)%w)

            if (params%do_forces .and. params%exp_forces) then

               n_dim_idx = 1
               outerforces: do j = 1, n_species
                  do k = 1, n_species

                     if (j > k) cycle

                     call get_pair_distribution_forces(n_sites, params%exp_energy_scales(params%pdf_idx),&
                          & params%exp_data(params%pdf_idx)%x, params%exp_data(params%pdf_idx)%y,&
                          & forces_pair_distribution, virial,&
                          & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                          & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                          &%r_range_min, params%r_range_max, params&
                          &%pair_distribution_n_samples,&
                          & y_pair_distribution(1:params&
                          &%pair_distribution_n_samples), params%pair_distribution_rcut&
                          &, j, k, pair_distribution_partial_der(1:params &
                          &%pair_distribution_n_samples, n_dim_idx,&
                          & j_beg:j_end), params%pair_distribution_partial,&
                          & params%pair_distribution_kde_sigma,&
                          & ((n_atoms_of_species(j)*&
                          & n_atoms_of_species(k))/dfloat(n_sites)/&
                          & dfloat(n_sites)), (dfloat(n_sites)/v_uc), params%pair_distribution_output)

                     n_dim_idx = n_dim_idx + 1

                     if (n_dim_idx > n_dim_partial) then
                        exit outerforces
                     end if

                  end do
               end do outerforces

!              ---   The cell half of the virial   --- !
!
!              get_pair_distribution_forces differentiates the interatomic
!              distances; the pattern also depends on the cell directly,
!              because each partial is normalised by the number density, and a
!              homogeneous strain changes the volume as well as the distances.
!              That part is a property of the whole pattern rather than of one
!              (a,b) channel, so it is added here, once, and only by rank 0 --
!              the virial is summed over ranks afterwards.
!
!              g(r) is proportional to V outright, so V dy/dV is y itself.
!              D(r) = 4 pi rho r ( g(r) - 1 ) has the volume cancel in the g
!              term, leaving only the subtracted background, which goes as 1/V.
               if (rank == 0) then
                  if (trim(params%pair_distribution_output) == "D(r)") then
                     dedv = -params%exp_energy_scales(params%pdf_idx)*&
                          & sum((y_pair_distribution(1:params%pair_distribution_n_samples) - &
                          &      params%exp_data(params%pdf_idx)%y)*&
                          &     4.d0*pi*(dfloat(n_sites)/v_uc)*&
                          &     x_pair_distribution(1:params%pair_distribution_n_samples))
                  else
                     dedv = -params%exp_energy_scales(params%pdf_idx)*&
                          & sum((y_pair_distribution(1:params%pair_distribution_n_samples) - &
                          &      params%exp_data(params%pdf_idx)%y)*&
                          &     y_pair_distribution(1:params%pair_distribution_n_samples))
                  end if

                  do i = 1, 3
                     virial(i, i) = virial(i, i) + dedv
                  end do
               end if
            end if

            ! open(unit=1234, file="grad", status="unknown")
            ! do i = 100, 110
            !    write(1234,  '(A,1X,I8,1X,F20.8)'), "dg_dr_0^1 ", i, pair_distribution_der_temp( i )
            ! end do
            ! close(unit=1234)

         end if

         !###---   If not doing partial pair distribution functions   ---###!

      else
#ifdef _MPIF90
         call mpi_reduce(y_pair_distribution,&
              & y_pair_distribution_temp, params&
              &%pair_distribution_n_samples,&
              & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
              & MPI_COMM_WORLD, ierr)

         y_pair_distribution = y_pair_distribution_temp
         deallocate (y_pair_distribution_temp)

         call mpi_bcast(y_pair_distribution, params&
              &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
              & MPI_COMM_WORLD, ierr)

         !    call mpi_reduce(forces_pair_distribution,&
         !         & forces_pair_distribution_temp, 3 * n_sites,&
         !         & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         !         & MPI_COMM_WORLD, ierr)

         !    call mpi_bcast(forces_pair_distribution, 3*n_sites, MPI_DOUBLE_PRECISION, 0,&
         !         & MPI_COMM_WORLD, ierr)

#endif
         y_pair_distribution = y_pair_distribution* &
              & dot_product(cross_product(a_box, b_box),&
              & c_box)/(dfloat(indices(1)*indices(2)&
              &*indices(3)))/dfloat(n_sites)/dfloat(n_sites)

      end if

      if (params%pair_distribution_partial .and. allocated(factors)) deallocate (factors)

      ! Write out the partial pair distribution functions
      call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

      if (rank == 0 .and. params%write_pair_distribution .and. write_condition) then
         ! call write_partial_exp(params%do_mc, params%do_md, mc_istep, md_istep,&
         !      & params%write_xyz, params%pair_distribution_partial,&
         !      & n_species, params%pair_distribution_n_samples,&
         !      & n_dim_partial , x_pair_distribution(1:params&
         !      &%pair_distribution_n_samples), y_pair_distribution(1:params&
         !      &%pair_distribution_n_samples), pair_distribution_partial(1:params &
         !      &%pair_distribution_n_samples, 1:n_dim_partial),&
         !      & species_types , 'pair_distribution')

         call get_overwrite_condition(params%do_mc, params%do_md,&
              & mc_istep, md_istep, params%write_xyz,&
              & overwrite_condition)

         if (params%pair_distribution_partial) then
            n_dim_idx = 1
            outer3: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle

                  write (filename, '(A)')&
                       & 'pair_distribution_'//trim(params&
                       &%species_types(j))//'_'//trim(params&
                       &%species_types(k))//&
                       & "_prediction.dat"
                  call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
                       & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx),&
                       & overwrite_condition, filename, 'pair_distribution')

                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) then
                     exit outer3
                  end if

               end do
            end do outer3
         end if

         write (filename, '(A)')&
              & "pair_distribution_total.dat"
         call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
              &y_pair_distribution(1:params%pair_distribution_n_samples),&
              & overwrite_condition, filename, "pair_distribution  output: "//trim(params&
              &%pair_distribution_output))

      end if

   end subroutine calculate_pair_distribution

   subroutine finalize_pair_distribution(params, x_pair_distribution&
        &, y_pair_distribution, y_pair_distribution_temp,&
        & pair_distribution_partial, pair_distribution_partial_temp,&
        & do_derivatives, pair_distribution_der, pair_distribution_partial_der,&
        & pair_distribution_partial_temp_der, n_atoms_of_species, rank)
      implicit none
      type(input_parameters), intent(in) :: params
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: x_pair_distribution(:)
      real(dp), allocatable, intent(inout) :: y_pair_distribution(:)
      real(dp), allocatable, intent(inout) :: pair_distribution_partial(:, :)
      real(dp), allocatable, intent(inout) :: n_atoms_of_species(:)
      real(dp), allocatable, intent(inout) :: pair_distribution_partial_temp(:, :)
      real(dp), allocatable, intent(inout) :: y_pair_distribution_temp(:)
      real(dp), allocatable, intent(inout) :: pair_distribution_der(:, :)
      real(dp), allocatable, intent(inout) :: pair_distribution_partial_der(:, :, :)
      real(dp), allocatable, intent(inout) :: pair_distribution_partial_temp_der(:, :, :)
      logical, intent(in) :: do_derivatives

      ! Naive finalization, include the logic of how things are actually
      ! allocated above rather then allocating and deallocating

      if (allocated(x_pair_distribution)) deallocate (x_pair_distribution)
      if (allocated(y_pair_distribution)) deallocate (y_pair_distribution)
      if (allocated(y_pair_distribution_temp)) deallocate (y_pair_distribution_temp)
      if (allocated(pair_distribution_partial)) deallocate (pair_distribution_partial)
      if (allocated(pair_distribution_partial_temp)) deallocate (pair_distribution_partial_temp)
      if (allocated(pair_distribution_der)) deallocate (pair_distribution_der)
      if (allocated(pair_distribution_partial_der)) deallocate (pair_distribution_partial_der)
      if (allocated(pair_distribution_partial_temp_der)) deallocate (pair_distribution_partial_temp_der)
      !    if ( allocated( forces_pair_distribution )          ) deallocate(forces_pair_distribution)
      if (allocated(n_atoms_of_species)) deallocate (n_atoms_of_species)

   end subroutine finalize_pair_distribution

#ifdef _GPU
   subroutine calculate_structure_factor(params, x_structure_factor, x_structure_factor_temp,&
        & y_structure_factor, y_structure_factor_temp,&
        & structure_factor_partial, structure_factor_partial_temp,&
        & x_pair_distribution, y_pair_distribution, &
        & pair_distribution_partial, n_species, species_types, n_atoms_of_species,&
        & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
        & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
        & neighbor_species, species, rank, q_beg, q_end, ntasks,&
        & sinc_factor_matrix, do_derivatives, &
        & nk, nk_d, k_index_d, j2_index_d, xyz_k_d, pair_distribution_partial_d, pair_distribution_partial_der_d, &
        & st_nk_d, st_k_index_d, st_j2_index_d, st_pair_distribution_partial_d, st_pair_distribution_partial_der_d, &
        pair_distribution_partial_der,&
        & energies_sf, forces_sf, virial_sf, use_matrix_forces, cublas_handle, gpu_stream, gpu_host_storage, gpu_low_memory)
      implicit none
      type(input_parameters), intent(inout) :: params
      real(dp), allocatable, intent(out) :: x_structure_factor(:)
      real(dp), allocatable, intent(out) :: x_structure_factor_temp(:)
      real(dp), allocatable, intent(out) :: y_structure_factor(:)
      real(dp), allocatable, intent(out) :: structure_factor_partial(:, :)
      real(dp), allocatable, intent(out) :: structure_factor_partial_temp(:, :)
      real(dp), allocatable, intent(out) :: y_structure_factor_temp(:)
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      real(dp), intent(in), allocatable :: x_pair_distribution(:)
      real(dp), intent(in), allocatable :: y_pair_distribution(:)
      real(dp), intent(in), allocatable :: pair_distribution_partial(:, :)
      real(dp), intent(in), allocatable :: n_atoms_of_species(:)
      real(dp), intent(in), allocatable :: pair_distribution_partial_der(:, :, :)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: species(:)
      integer, intent(in), allocatable :: nk(:)
      character*8, allocatable, intent(in) :: species_types(:)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), allocatable, intent(inout) :: sinc_factor_matrix(:, :)
      real(dp), allocatable, intent(inout) :: energies_sf(:)
      real(dp), allocatable, intent(inout) :: forces_sf(:, :)
      real(dp), intent(inout) :: virial_sf(1:3, 1:3)
      real(dp), allocatable :: sinc_factor_matrix_temp(:, :)
      real(dp), allocatable :: temp_pdf(:, :)
      real(dp) :: v_uc
      integer, intent(in) :: n_species
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: ntasks
      integer, intent(out) :: q_beg
      integer, intent(out) :: q_end
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: rank
      integer, intent(out) :: ierr
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: i2
      integer :: n_dim_partial
      integer :: n_dim_idx
      integer :: n
      integer :: m
      real(dp) :: dq
      real(dp) :: f
      real(dp) :: cabh
      real(dp) :: delta
      real(dp), parameter :: pi = acos(-1.0)
      character*1024 :: filename
      logical :: overwrite_condition
      logical :: write_condition
      logical, intent(in) :: do_derivatives
      logical, intent(in) :: use_matrix_forces

      type(c_ptr) :: cublas_handle
      type(c_ptr) :: gpu_stream
      type(c_ptr), allocatable :: nk_d(:)
      type(c_ptr), allocatable :: nk_flags_d(:)
      type(c_ptr), allocatable :: nk_flags_sum_d(:)
      type(c_ptr), allocatable :: k_index_d(:)
      type(c_ptr), allocatable :: j2_index_d(:)
      type(c_ptr), allocatable :: rjs_index_d(:)
      type(c_ptr), allocatable :: xyz_k_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_der_d(:)
      integer(c_size_t), allocatable :: st_nk_d(:)
      integer(c_size_t), allocatable :: st_k_index_d(:)
      integer(c_size_t), allocatable :: st_j2_index_d(:)
      integer(c_size_t), allocatable :: st_rjs(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_d(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_der_d(:)
      type(gpu_host_storage_type), intent(inout), allocatable, target :: gpu_host_storage(:)
      logical, intent(in) :: gpu_low_memory

      v_uc = dot_product(cross_product(a_box,&
              & b_box), c_box)/(&
              & dfloat(indices(1)*indices(2)&
              &*indices(3)))

      if (params%structure_factor_from_pdf) then
         q_beg = 1
         q_end = params%structure_factor_n_samples
#ifdef _MPIF90
         ! We need to do an integral for each q value
         ! This can be split among the processes

         ! Split each q the integrals among each of the
         ! processes, just like we do with the atoms,
         ! and then collect accordingly

         if (rank < mod(params%structure_factor_n_samples, ntasks)) then
            q_beg = 1 + rank*(params%structure_factor_n_samples/ntasks + 1)
         else
            q_beg = 1 + mod(params%structure_factor_n_samples, ntasks)*(params&
                 &%structure_factor_n_samples/ntasks + 1) &
                 &+ (rank - mod(params%structure_factor_n_samples, ntasks))*(params&
                 &%structure_factor_n_samples/ntasks)
         end if
         if (rank < mod(params%structure_factor_n_samples, ntasks)) then
            q_end = (rank + 1)*(params%structure_factor_n_samples/ntasks + 1)
         else
            q_end = q_beg + params%structure_factor_n_samples/ntasks - 1
         end if

#endif
         if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then

            n_dim_partial = n_species*(n_species + 1)/2

            if (allocated(structure_factor_partial)) deallocate (structure_factor_partial)
            allocate (structure_factor_partial(1:params&
                 &%structure_factor_n_samples, 1:n_dim_partial))
            structure_factor_partial = 0.d0
         end if

      end if

      if (allocated(y_structure_factor)) deallocate (y_structure_factor)
      allocate (y_structure_factor(1:params%structure_factor_n_samples))
      y_structure_factor = 0.d0

      if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
#ifdef _MPIF90
         allocate (structure_factor_partial_temp(1:params&
              &%structure_factor_n_samples, 1:n_dim_partial))

         structure_factor_partial_temp = 0.0d0
#endif

      else

#ifdef _MPIF90
         allocate (y_structure_factor_temp(1:params&
              &%structure_factor_n_samples))

         y_structure_factor_temp = 0.0d0
#endif

      end if

      if (allocated(x_structure_factor)) deallocate (x_structure_factor)
      if (allocated(x_structure_factor_temp)) deallocate (x_structure_factor_temp)
      call linspace(x_structure_factor, params&
           &%q_range_min, params%q_range_max, params&
           &%structure_factor_n_samples, dq)

      call linspace(x_structure_factor_temp, params&
           &%q_range_min, params%q_range_max, params&
           &%structure_factor_n_samples, dq)

      if (trim(params%q_units) == "xrd" .or. params%q_units == "twotheta") then
         ! assume that theta is given for the Q range
         do i = 1, params%structure_factor_n_samples
            ! This gives s, but Q  = 2 pi * s = 4pi sin(theta) / lambda
            x_structure_factor(i) = 2.d0*sin(pi*&
                 & x_structure_factor(i)/180.d0/2.d0)/&
                 & params%xrd_wavelength
         end do
      elseif (trim(params%q_units) == "saxs" .or. params%q_units == "q") then
         do i = 1, params%structure_factor_n_samples
            x_structure_factor(i) = x_structure_factor(i)/2.d0/pi
         end do
      end if

      ! Here the units of q range are that called Q (1/A) in
      ! the literature, small q in the literature are =
      ! 2sin(theta)/lambda = Q/2/pi

      if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then

         v_uc = dot_product(cross_product(a_box,&
              & b_box), c_box)/(&
              & dfloat(indices(1)*indices(2)&
              &*indices(3)))

         if (.not. params%structure_factor_matrix) then
            call get_partial_structure_factor(q_beg, q_end, &
                 & pair_distribution_partial,&
                 & x_structure_factor(1:params%structure_factor_n_samples),&
                 & x_pair_distribution, params%pair_distribution_rcut,&
                 & params%pair_distribution_n_samples, params&
                 &%structure_factor_n_samples, n_species, n_dim_partial,&
                 & n_atoms_of_species, n_sites, dfloat(n_sites)/v_uc, &
                 & params%structure_factor_window,&
                 & structure_factor_partial)

#ifdef _MPIF90
            call mpi_reduce(structure_factor_partial,&
                 & structure_factor_partial_temp, params&
                 &%structure_factor_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                 & MPI_COMM_WORLD, ierr)
            structure_factor_partial = structure_factor_partial_temp
            deallocate (structure_factor_partial_temp)

            call mpi_bcast(structure_factor_partial, params&
                 &%structure_factor_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
                 & MPI_COMM_WORLD, ierr)

#endif
         else
            if (.not. allocated(sinc_factor_matrix)) then
               call get_sinc_factor_matrix(q_beg, q_end,&
                    & x_structure_factor, x_pair_distribution,&
                    & params%pair_distribution_rcut, params&
                    &%pair_distribution_n_samples, params&
                    &%structure_factor_n_samples, params%structure_factor_window,&
                    & sinc_factor_matrix)

#ifdef _MPIF90
               allocate (sinc_factor_matrix_temp(1:params%structure_factor_n_samples, 1:params&
                    &%pair_distribution_n_samples))
               sinc_factor_matrix_temp = 0.d0

               call mpi_reduce(sinc_factor_matrix, sinc_factor_matrix_temp&
                    &, params%structure_factor_n_samples*params&
                    &%pair_distribution_n_samples,&
                    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,&
                    & ierr)

               sinc_factor_matrix = sinc_factor_matrix_temp

               deallocate (sinc_factor_matrix_temp)
               call mpi_bcast(sinc_factor_matrix, params&
                    &%structure_factor_n_samples*params&
                    &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
                    & MPI_COMM_WORLD, ierr)

#endif
            end if

            ! Now make a blas call to perform the matrix multiplication necessary to obtain the structure factor.

            n = params%structure_factor_n_samples
            m = params%pair_distribution_n_samples
            k = n_dim_partial

            ! [sinc_factor_matrix] = n_samples_sf * n_samples_pc  ( N x M )
            ! [g_ab]               = n_samples_pc * n_dim_partial ( M x K )
            ! [S_ab]               = n_samples_sf * n_dim_partial ( N x K )
            ! S_ab = [sinc_factor_matrix] x [g_ab - 1]
            !      = ( N x M ) . ( M x K ) -> ( N x K )
            ! alpha * (A * B) + beta * C
            ! A === sinc_factor_partial
            ! B === [g_ab - 1]
            ! C === S_ab
            !   TRANS_A, TRANS_B  N_ROWS_A  N_COLS_B  N_COLS_A ( == N_ROWS_B ) alpha,
            !                                                   A,  first_dim_A,    B,      first_dim_B,  beta,  C, first_dim_C

            call dgemm("N", "N", n, k, m, 1.d0, sinc_factor_matrix, n,&
                 & pair_distribution_partial - 1.d0, m, 0.d0,&
                 & structure_factor_partial, n)

            !###---   Do derivatives!    ---###!

            ! Should do this smartly with memory allocation

            ! All processes do this calculation, so they have the whole of the partial structure factors to work with.
            n_dim_idx = 1
            outersfm: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle

                  if (j == k) f = 1.d0
                  if (j /= k) f = 0.d0

                  cabh = ((n_atoms_of_species(j)/dfloat(n_sites)) &
                       &*(n_atoms_of_species(k)/dfloat(n_sites)) &
                       & )**(0.5)

                  structure_factor_partial(1:params&
                       &%structure_factor_n_samples, n_dim_idx) = f +&
                       & 4.d0*pi*cabh*(dfloat(n_sites)/v_uc)* &
                       & structure_factor_partial(1:params&
                       &%structure_factor_n_samples, n_dim_idx)

                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) exit outersfm
               end do
            end do outersfm

         end if

! #ifdef _MPIF90
!           call mpi_reduce(structure_factor_partial,&
!                & structure_factor_partial_temp, params&
!                &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
!                & MPI_COMM_WORLD, ierr)
!           structure_factor_partial = structure_factor_partial_temp
!           deallocate(structure_factor_partial_temp)

!           call mpi_bcast(structure_factor_partial, params&
!                &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
!                & MPI_COMM_WORLD, ierr)

! #endif

         ! using dq as a temp variable
         y_structure_factor = 0.d0
         dq = 0.d0
         n_dim_idx = 1
         outer2: do j = 1, n_species
            dq = dq + (n_atoms_of_species(j)/dfloat(n_sites))
            do k = 1, n_species

               if (j > k) cycle

               if (j == k) f = 1.d0
               if (j /= k) f = 2.d0

               if (j == k) delta = 1.d0
               if (j /= k) delta = 0.d0

               y_structure_factor(1:params%structure_factor_n_samples) = &
                    & y_structure_factor(1:params%structure_factor_n_samples) +  &
                    & f*(n_atoms_of_species(j)*n_atoms_of_species(k))**0.5* &
                    &  (structure_factor_partial(1:params%structure_factor_n_samples, n_dim_idx) &
                    &  - delta)/dfloat(n_sites) !/ dfloat(n_sites)

               n_dim_idx = n_dim_idx + 1

               if (n_dim_idx > n_dim_partial) exit outer2

            end do
         end do outer2

         y_structure_factor = y_structure_factor + 1.d0
         ! --- Preprocess the structure factor according to the output --- !
         if (trim(params%sf_output) == "q*i(q)") then
            y_structure_factor = x_structure_factor*(y_structure_factor - 1.d0)
         end if

         if (params%valid_sf) then

            allocate (energies_sf(1:n_sites))
            energies_sf = 0.d0

            if (params%valid_sf .and. allocated(params%exp_energy_scales)) then

               call get_energy_scale(params%do_md, params%do_mc,&
                    & md_istep, params%md_nsteps, mc_istep, params&
                    &%mc_nsteps, params%exp_energy_scales_initial(params%sf_idx), &
                    & params%exp_energy_scales_final(params%sf_idx), &
                    & params%exp_energy_scales(params%sf_idx))

               call get_exp_energies(params%exp_energy_scales(params&
                    &%sf_idx), params%exp_data(params%sf_idx)%y&
                    &, y_structure_factor,&
                    & params%structure_factor_n_samples, n_sites,&
                    & energies_sf(i_beg:i_end), params%exp_data(params%sf_idx)%w)

               if (params%do_forces .and. params%exp_forces) then

                  allocate (forces_sf(1:3, 1:n_sites))
                  forces_sf = 0.d0
                  virial_sf = 0.d0

                  n_dim_idx = 1
                  outerf: do j = 1, n_species
                     do k = 1, n_species

                        if (j > k) cycle

                        if (j == k) f = 1.d0
                        if (j /= k) f = 2.d0

                        if (use_matrix_forces) then
                           call get_structure_factor_forces_matrix(i_beg, i_end, n_sites, params%exp_energy_scales(params%sf_idx),&
                                & params%exp_data(params%sf_idx)%x, params%exp_data(params%sf_idx)%y,&
                                & forces_sf, virial_sf,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), species_types, species(i_beg:i_end), rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_structure_factor(1:params&
                                &%structure_factor_n_samples), y_structure_factor(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .false., params%xrd_output,&
                                & n_atoms_of_species, .false., rank, cublas_handle, gpu_stream, &
                & nk(n_dim_idx), nk_d(n_dim_idx), k_index_d(n_dim_idx), j2_index_d(n_dim_idx), xyz_k_d(n_dim_idx), pair_distribution_partial_d(n_dim_idx), pair_distribution_partial_der_d(n_dim_idx), &
                & st_nk_d(n_dim_idx), st_k_index_d(n_dim_idx), st_j2_index_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), &
                                gpu_host_storage(n_dim_idx), gpu_low_memory)
                        else
                           call get_structure_factor_forces(n_sites, params%exp_energy_scales(params%sf_idx),&
                                & params%exp_data(params%sf_idx)%x, params%exp_data(params%sf_idx)%y,&
                                & forces_sf, virial_sf,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                                & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_structure_factor(1:params&
                                &%structure_factor_n_samples), y_structure_factor(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .false., params%xrd_output,&
                                & n_atoms_of_species, .false., rank)
                        end if

                        n_dim_idx = n_dim_idx + 1

                        if (n_dim_idx > n_dim_partial) exit outerf

                     end do
                  end do outerf
               end if
            end if
         end if

         ! y_structure_factor(1:params%structure_factor_n_samples)&
         !      & = y_structure_factor(1:params&
         !      &%structure_factor_n_samples) !/ dq

      elseif (params%structure_factor_from_pdf) then

         call get_structure_factor_from_pdf(q_beg, q_end, &
              & y_structure_factor(1:params%structure_factor_n_samples),&
              & y_pair_distribution,&
              & x_structure_factor(1:params%structure_factor_n_samples),&
              & x_pair_distribution, params%pair_distribution_rcut,&
              & params%pair_distribution_n_samples, params&
              &%structure_factor_n_samples, n_species,&
              & n_atoms_of_species, n_sites, dfloat(n_sites)/v_uc, params%structure_factor_window)

      else

         call get_structure_factor_explicit(n_sites, &
              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
              & params%structure_factor_n_samples,&
              & x_structure_factor(1:params&
              &%structure_factor_n_samples), y_structure_factor(1:params&
              &%structure_factor_n_samples), params%r_range_min, params%r_range_max,&
              & params%pair_distribution_rcut, .false., 1&
              &, 1, params%structure_factor_window)
      end if

      if (.not. (params%structure_factor_from_pdf .and. params%pair_distribution_partial)) then
#ifdef _MPIF90

         call mpi_reduce(y_structure_factor,&
              & y_structure_factor_temp, params&
              &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
              & MPI_COMM_WORLD, ierr)

         y_structure_factor = y_structure_factor_temp
         deallocate (y_structure_factor_temp)
#endif
         y_structure_factor = y_structure_factor + 1.d0

#ifdef _MPIF90
         call mpi_bcast(y_structure_factor, params &
              &%structure_factor_n_samples,&
              & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

      end if
      ! Write out the partial structure functions
      call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

      if (rank == 0 .and. params%write_structure_factor .and. write_condition) then

         call get_overwrite_condition(params%do_mc, params%do_md&
              &, mc_istep, md_istep, params%write_xyz,&
              & overwrite_condition)

         if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
            n_dim_idx = 1
            outer3: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle

                  ! write with the temp data
                  write (filename, '(A)')&
                       & 'structure_factor_'//trim(species_types(j))//'_'//trim(species_types(k))//&
                       & "_prediction.dat"
                  call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples)&
                       &,&
                       & structure_factor_partial(1:params&
                       &%structure_factor_n_samples, n_dim_idx),&
                       & overwrite_condition, filename, "structure_factor: units of "//trim(params%q_units))

                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) exit outer3

               end do
            end do outer3
         end if

         write (filename, '(A)')&
              & 'structure_factor_total.dat'
         call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples),&
              & y_structure_factor(1:params&
              &%structure_factor_n_samples),&
              & overwrite_condition, filename, "structure_factor: units&
              & of "//trim(params%q_units)//" output: "//trim(params&
              &%xrd_output))

      end if
   end subroutine calculate_structure_factor
#else
   subroutine calculate_structure_factor(params, x_structure_factor, x_structure_factor_temp,&
        & y_structure_factor, y_structure_factor_temp,&
        & structure_factor_partial, structure_factor_partial_temp,&
        & x_pair_distribution, y_pair_distribution, &
        & pair_distribution_partial, n_species, species_types, n_atoms_of_species,&
        & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
        & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
        & neighbor_species, species, rank, q_beg, q_end, ntasks,&
        & sinc_factor_matrix, do_derivatives, pair_distribution_partial_der,&
        & energies_sf, forces_sf, virial_sf, use_matrix_forces)
      implicit none
      type(input_parameters), intent(inout) :: params
      real(dp), allocatable, intent(out) :: x_structure_factor(:)
      real(dp), allocatable, intent(out) :: x_structure_factor_temp(:)
      real(dp), allocatable, intent(out) :: y_structure_factor(:)
      real(dp), allocatable, intent(out) :: structure_factor_partial(:, :)
      real(dp), allocatable, intent(out) :: structure_factor_partial_temp(:, :)
      real(dp), allocatable, intent(out) :: y_structure_factor_temp(:)
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      real(dp), intent(in), allocatable :: x_pair_distribution(:)
      real(dp), intent(in), allocatable :: y_pair_distribution(:)
      real(dp), intent(in), allocatable :: pair_distribution_partial(:, :)
      real(dp), intent(in), allocatable :: n_atoms_of_species(:)
      real(dp), intent(in), allocatable :: pair_distribution_partial_der(:, :, :)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: species(:)
      character*8, allocatable, intent(in) :: species_types(:)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), allocatable, intent(inout) :: sinc_factor_matrix(:, :)
      real(dp), allocatable, intent(inout) :: energies_sf(:)
      real(dp), allocatable, intent(inout) :: forces_sf(:, :)
      real(dp), intent(inout) :: virial_sf(1:3, 1:3)
      real(dp), allocatable :: sinc_factor_matrix_temp(:, :)
      real(dp), allocatable :: temp_pdf(:, :)
      real(dp) :: v_uc
      integer, intent(in) :: n_species
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: ntasks
      integer, intent(out) :: q_beg
      integer, intent(out) :: q_end
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: rank
      integer, intent(out) :: ierr
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: i2
      integer :: n_dim_partial
      integer :: n_dim_idx
      integer :: n
      integer :: m
      real(dp) :: dq
      real(dp) :: f
      real(dp) :: cabh
      real(dp) :: delta
      real(dp), parameter :: pi = acos(-1.0)
      character*1024 :: filename
      logical :: overwrite_condition
      logical :: write_condition
      logical, intent(in) :: do_derivatives
      logical, intent(in) :: use_matrix_forces

      v_uc = dot_product(cross_product(a_box,&
              & b_box), c_box)/(&
              & dfloat(indices(1)*indices(2)&
              &*indices(3)))

      if (params%structure_factor_from_pdf) then
         q_beg = 1
         q_end = params%structure_factor_n_samples
#ifdef _MPIF90
         ! We need to do an integral for each q value
         ! This can be split among the processes

         ! Split each q the integrals among each of the
         ! processes, just like we do with the atoms,
         ! and then collect accordingly

         if (rank < mod(params%structure_factor_n_samples, ntasks)) then
            q_beg = 1 + rank*(params%structure_factor_n_samples/ntasks + 1)
         else
            q_beg = 1 + mod(params%structure_factor_n_samples, ntasks)*(params&
                 &%structure_factor_n_samples/ntasks + 1) &
                 &+ (rank - mod(params%structure_factor_n_samples, ntasks))*(params&
                 &%structure_factor_n_samples/ntasks)
         end if
         if (rank < mod(params%structure_factor_n_samples, ntasks)) then
            q_end = (rank + 1)*(params%structure_factor_n_samples/ntasks + 1)
         else
            q_end = q_beg + params%structure_factor_n_samples/ntasks - 1
         end if

#endif
         if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then

            n_dim_partial = n_species*(n_species + 1)/2

            if (allocated(structure_factor_partial)) deallocate (structure_factor_partial)
            allocate (structure_factor_partial(1:params&
                 &%structure_factor_n_samples, 1:n_dim_partial))
            structure_factor_partial = 0.d0
         end if

      end if

      if (allocated(y_structure_factor)) deallocate (y_structure_factor)
      allocate (y_structure_factor(1:params%structure_factor_n_samples))
      y_structure_factor = 0.d0

      if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
#ifdef _MPIF90
         allocate (structure_factor_partial_temp(1:params&
              &%structure_factor_n_samples, 1:n_dim_partial))

         structure_factor_partial_temp = 0.0d0
#endif

      else

#ifdef _MPIF90
         allocate (y_structure_factor_temp(1:params&
              &%structure_factor_n_samples))

         y_structure_factor_temp = 0.0d0
#endif

      end if

      if (allocated(x_structure_factor)) deallocate (x_structure_factor)
      if (allocated(x_structure_factor_temp)) deallocate (x_structure_factor_temp)
      call linspace(x_structure_factor, params&
           &%q_range_min, params%q_range_max, params&
           &%structure_factor_n_samples, dq)

      call linspace(x_structure_factor_temp, params&
           &%q_range_min, params%q_range_max, params&
           &%structure_factor_n_samples, dq)

      if (trim(params%q_units) == "xrd" .or. params%q_units == "twotheta") then
         ! assume that theta is given for the Q range
         do i = 1, params%structure_factor_n_samples
            ! This gives s, but Q  = 2 pi * s = 4pi sin(theta) / lambda
            x_structure_factor(i) = 2.d0*sin(pi*&
                 & x_structure_factor(i)/180.d0/2.d0)/&
                 & params%xrd_wavelength
         end do
      elseif (trim(params%q_units) == "saxs" .or. params%q_units == "q") then
         do i = 1, params%structure_factor_n_samples
            x_structure_factor(i) = x_structure_factor(i)/2.d0/pi
         end do
      end if

      ! Here the units of q range are that called Q (1/A) in
      ! the literature, small q in the literature are =
      ! 2sin(theta)/lambda = Q/2/pi

      if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then

         v_uc = dot_product(cross_product(a_box,&
              & b_box), c_box)/(&
              & dfloat(indices(1)*indices(2)&
              &*indices(3)))

         if (.not. params%structure_factor_matrix) then
            call get_partial_structure_factor(q_beg, q_end, &
                 & pair_distribution_partial,&
                 & x_structure_factor(1:params%structure_factor_n_samples),&
                 & x_pair_distribution, params%pair_distribution_rcut,&
                 & params%pair_distribution_n_samples, params&
                 &%structure_factor_n_samples, n_species, n_dim_partial,&
                 & n_atoms_of_species, n_sites, dfloat(n_sites)/v_uc, &
                 & params%structure_factor_window,&
                 & structure_factor_partial)

#ifdef _MPIF90
            call mpi_reduce(structure_factor_partial,&
                 & structure_factor_partial_temp, params&
                 &%structure_factor_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
                 & MPI_COMM_WORLD, ierr)
            structure_factor_partial = structure_factor_partial_temp
            deallocate (structure_factor_partial_temp)

            call mpi_bcast(structure_factor_partial, params&
                 &%structure_factor_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
                 & MPI_COMM_WORLD, ierr)

#endif
         else
            if (.not. allocated(sinc_factor_matrix)) then
               call get_sinc_factor_matrix(q_beg, q_end,&
                    & x_structure_factor, x_pair_distribution,&
                    & params%pair_distribution_rcut, params&
                    &%pair_distribution_n_samples, params&
                    &%structure_factor_n_samples, params%structure_factor_window,&
                    & sinc_factor_matrix)

#ifdef _MPIF90
               allocate (sinc_factor_matrix_temp(1:params%structure_factor_n_samples, 1:params&
                    &%pair_distribution_n_samples))
               sinc_factor_matrix_temp = 0.d0

               call mpi_reduce(sinc_factor_matrix, sinc_factor_matrix_temp&
                    &, params%structure_factor_n_samples*params&
                    &%pair_distribution_n_samples,&
                    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,&
                    & ierr)

               sinc_factor_matrix = sinc_factor_matrix_temp

               deallocate (sinc_factor_matrix_temp)
               call mpi_bcast(sinc_factor_matrix, params&
                    &%structure_factor_n_samples*params&
                    &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
                    & MPI_COMM_WORLD, ierr)

#endif
            end if

            ! Now make a blas call to perform the matrix multiplication necessary to obtain the structure factor.

            n = params%structure_factor_n_samples
            m = params%pair_distribution_n_samples
            k = n_dim_partial

            ! [sinc_factor_matrix] = n_samples_sf * n_samples_pc  ( N x M )
            ! [g_ab]               = n_samples_pc * n_dim_partial ( M x K )
            ! [S_ab]               = n_samples_sf * n_dim_partial ( N x K )
            ! S_ab = [sinc_factor_matrix] x [g_ab - 1]
            !      = ( N x M ) . ( M x K ) -> ( N x K )
            ! alpha * (A * B) + beta * C
            ! A === sinc_factor_partial
            ! B === [g_ab - 1]
            ! C === S_ab
            !   TRANS_A, TRANS_B  N_ROWS_A  N_COLS_B  N_COLS_A ( == N_ROWS_B ) alpha,
            !                                                   A,  first_dim_A,    B,      first_dim_B,  beta,  C, first_dim_C

            call dgemm("N", "N", n, k, m, 1.d0, sinc_factor_matrix, n,&
                 & pair_distribution_partial - 1.d0, m, 0.d0,&
                 & structure_factor_partial, n)

            !###---   Do derivatives!    ---###!

            ! Should do this smartly with memory allocation

            ! All processes do this calculation, so they have the whole of the partial structure factors to work with.
            n_dim_idx = 1
            outersfm: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle

                  if (j == k) f = 1.d0
                  if (j /= k) f = 0.d0

                  cabh = ((n_atoms_of_species(j)/dfloat(n_sites)) &
                       &*(n_atoms_of_species(k)/dfloat(n_sites)) &
                       & )**(0.5)

                  structure_factor_partial(1:params&
                       &%structure_factor_n_samples, n_dim_idx) = f +&
                       & 4.d0*pi*cabh*(dfloat(n_sites)/v_uc)* &
                       & structure_factor_partial(1:params&
                       &%structure_factor_n_samples, n_dim_idx)

                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) exit outersfm
               end do
            end do outersfm

         end if

! #ifdef _MPIF90
!           call mpi_reduce(structure_factor_partial,&
!                & structure_factor_partial_temp, params&
!                &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
!                & MPI_COMM_WORLD, ierr)
!           structure_factor_partial = structure_factor_partial_temp
!           deallocate(structure_factor_partial_temp)

!           call mpi_bcast(structure_factor_partial, params&
!                &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
!                & MPI_COMM_WORLD, ierr)

! #endif

         ! using dq as a temp variable
         y_structure_factor = 0.d0
         dq = 0.d0
         n_dim_idx = 1
         outer2: do j = 1, n_species
            dq = dq + (n_atoms_of_species(j)/dfloat(n_sites))
            do k = 1, n_species

               if (j > k) cycle

               if (j == k) f = 1.d0
               if (j /= k) f = 2.d0

               if (j == k) delta = 1.d0
               if (j /= k) delta = 0.d0

               y_structure_factor(1:params%structure_factor_n_samples) = &
                    & y_structure_factor(1:params%structure_factor_n_samples) +  &
                    & f*(n_atoms_of_species(j)*n_atoms_of_species(k))**0.5* &
                    &  (structure_factor_partial(1:params%structure_factor_n_samples, n_dim_idx) &
                    &  - delta)/dfloat(n_sites) !/ dfloat(n_sites)

               n_dim_idx = n_dim_idx + 1

               if (n_dim_idx > n_dim_partial) exit outer2

            end do
         end do outer2

         y_structure_factor = y_structure_factor + 1.d0
         ! --- Preprocess the structure factor according to the output --- !
         if (trim(params%sf_output) == "q*i(q)") then
            y_structure_factor = x_structure_factor*(y_structure_factor - 1.d0)
         end if

         if (params%valid_sf) then

            allocate (energies_sf(1:n_sites))
            energies_sf = 0.d0

            if (params%valid_sf .and. allocated(params%exp_energy_scales)) then

               call get_energy_scale(params%do_md, params%do_mc,&
                    & md_istep, params%md_nsteps, mc_istep, params&
                    &%mc_nsteps, params%exp_energy_scales_initial(params%sf_idx), &
                    & params%exp_energy_scales_final(params%sf_idx), &
                    & params%exp_energy_scales(params%sf_idx))

               call get_exp_energies(params%exp_energy_scales(params&
                    &%sf_idx), params%exp_data(params%sf_idx)%y&
                    &, y_structure_factor,&
                    & params%structure_factor_n_samples, n_sites,&
                    & energies_sf(i_beg:i_end), params%exp_data(params%sf_idx)%w)

               if (params%do_forces .and. params%exp_forces) then

                  allocate (forces_sf(1:3, 1:n_sites))
                  forces_sf = 0.d0

                  n_dim_idx = 1
                  outerf: do j = 1, n_species
                     do k = 1, n_species

                        if (j > k) cycle

                        if (j == k) f = 1.d0
                        if (j /= k) f = 2.d0

                        if (use_matrix_forces) then
                           call get_structure_factor_forces_matrix(n_sites, params%exp_energy_scales(params%sf_idx),&
                                & params%exp_data(params%sf_idx)%x, params%exp_data(params%sf_idx)%y,&
                                & forces_sf, virial_sf,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                                & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_structure_factor(1:params&
                                &%structure_factor_n_samples), y_structure_factor(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .false., params%xrd_output,&
                                & n_atoms_of_species, .false., rank, w=params%exp_data(params%sf_idx)%w)
                        else
                           call get_structure_factor_forces(n_sites, params%exp_energy_scales(params%sf_idx),&
                                & params%exp_data(params%sf_idx)%x, params%exp_data(params%sf_idx)%y,&
                                & forces_sf, virial_sf,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                                & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_structure_factor(1:params&
                                &%structure_factor_n_samples), y_structure_factor(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .false., params%xrd_output,&
                                & n_atoms_of_species, .false., rank)
                        end if

                        n_dim_idx = n_dim_idx + 1

                        if (n_dim_idx > n_dim_partial) exit outerf

                     end do
                  end do outerf
               end if
            end if
         end if

         ! y_structure_factor(1:params%structure_factor_n_samples)&
         !      & = y_structure_factor(1:params&
         !      &%structure_factor_n_samples) !/ dq

      elseif (params%structure_factor_from_pdf) then

         call get_structure_factor_from_pdf(q_beg, q_end, &
              & y_structure_factor(1:params%structure_factor_n_samples),&
              & y_pair_distribution,&
              & x_structure_factor(1:params%structure_factor_n_samples),&
              & x_pair_distribution, params%pair_distribution_rcut,&
              & params%pair_distribution_n_samples, params&
              &%structure_factor_n_samples, n_species,&
              & n_atoms_of_species, n_sites, dfloat(n_sites)/v_uc, params%structure_factor_window)

      else

         call get_structure_factor_explicit(n_sites, &
              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
              & params%structure_factor_n_samples,&
              & x_structure_factor(1:params&
              &%structure_factor_n_samples), y_structure_factor(1:params&
              &%structure_factor_n_samples), params%r_range_min, params%r_range_max,&
              & params%pair_distribution_rcut, .false., 1&
              &, 1, params%structure_factor_window)
      end if

      if (.not. (params%structure_factor_from_pdf .and. params%pair_distribution_partial)) then
#ifdef _MPIF90

         call mpi_reduce(y_structure_factor,&
              & y_structure_factor_temp, params&
              &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
              & MPI_COMM_WORLD, ierr)

         y_structure_factor = y_structure_factor_temp
         deallocate (y_structure_factor_temp)
#endif
         y_structure_factor = y_structure_factor + 1.d0

#ifdef _MPIF90
         call mpi_bcast(y_structure_factor, params &
              &%structure_factor_n_samples,&
              & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

      end if
      ! Write out the partial structure functions
      call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

      if (rank == 0 .and. params%write_structure_factor .and. write_condition) then

         call get_overwrite_condition(params%do_mc, params%do_md&
              &, mc_istep, md_istep, params%write_xyz,&
              & overwrite_condition)

         if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
            n_dim_idx = 1
            outer3: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle

                  ! write with the temp data
                  write (filename, '(A)')&
                       & 'structure_factor_'//trim(species_types(j))//'_'//trim(species_types(k))//&
                       & "_prediction.dat"
                  call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples)&
                       &,&
                       & structure_factor_partial(1:params&
                       &%structure_factor_n_samples, n_dim_idx),&
                       & overwrite_condition, filename, "structure_factor: units of "//trim(params%q_units))

                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) exit outer3

               end do
            end do outer3
         end if

         write (filename, '(A)')&
              & 'structure_factor_total.dat'
         call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples),&
              & y_structure_factor(1:params&
              &%structure_factor_n_samples),&
              & overwrite_condition, filename, "structure_factor: units&
              & of "//trim(params%q_units)//" output: "//trim(params&
              &%xrd_output))

      end if
   end subroutine calculate_structure_factor
#endif

   subroutine finalize_structure_factor(params, x_structure_factor, x_structure_factor_temp,&
        & y_structure_factor, y_structure_factor_temp,&
        & structure_factor_partial, structure_factor_partial_temp,&
        & x_pair_distribution, y_pair_distribution, &
        & pair_distribution_partial, sinc_factor_matrix)
      implicit none
      type(input_parameters), intent(in) :: params
      real(dp), allocatable, intent(inout) :: x_structure_factor(:)
      real(dp), allocatable, intent(inout) :: x_structure_factor_temp(:)
      real(dp), allocatable, intent(inout) :: y_structure_factor(:)
      real(dp), allocatable, intent(inout) :: structure_factor_partial(:, :)
      real(dp), allocatable, intent(inout) :: structure_factor_partial_temp(:, :)
      real(dp), allocatable, intent(inout) :: y_structure_factor_temp(:)
      real(dp), intent(inout), allocatable :: x_pair_distribution(:)
      real(dp), intent(inout), allocatable :: y_pair_distribution(:)
      real(dp), intent(inout), allocatable :: pair_distribution_partial(:, :)
      real(dp), allocatable, intent(inout) :: sinc_factor_matrix(:, :)

      if (allocated(x_structure_factor)) deallocate (x_structure_factor)
      if (allocated(x_structure_factor_temp)) deallocate (x_structure_factor_temp)
      if (allocated(y_structure_factor)) deallocate (y_structure_factor)
      if (allocated(structure_factor_partial)) deallocate (structure_factor_partial)
      if (allocated(structure_factor_partial_temp)) deallocate (structure_factor_partial_temp)
      if (allocated(y_structure_factor_temp)) deallocate (y_structure_factor_temp)
      if (allocated(x_pair_distribution)) deallocate (x_pair_distribution)
      if (allocated(y_pair_distribution)) deallocate (y_pair_distribution)
      if (allocated(pair_distribution_partial)) deallocate (pair_distribution_partial)
      if (allocated(sinc_factor_matrix)) deallocate (sinc_factor_matrix)

   end subroutine finalize_structure_factor

#ifdef _GPU
   subroutine calculate_xrd(params, x_xrd, x_xrd_temp,&
        & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
        & structure_factor_partial, structure_factor_partial_temp,&
        & n_species, species_types, n_atoms_of_species,&
        & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
        & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
        & neighbor_species, species, rank, q_beg, q_end, ntasks,&
        & sinc_factor_matrix, do_derivatives, &
        & nk, nk_d, k_index_d, j2_index_d, xyz_k_d, pair_distribution_partial_d, pair_distribution_partial_der_d, &
        & st_nk_d, st_k_index_d, st_j2_index_d, st_pair_distribution_partial_d, st_pair_distribution_partial_der_d, &
         pair_distribution_partial_der,&
    & energies_xrd, forces_xrd, virial_xrd, neutron, use_matrix_forces, cublas_handle, gpu_stream, gpu_host_storage, gpu_low_memory)
      implicit none
      type(input_parameters), intent(inout) :: params
      real(dp), allocatable, intent(out) :: x_xrd(:)
      real(dp), allocatable, intent(out) :: x_xrd_temp(:)
      real(dp), allocatable, intent(out) :: y_xrd(:)
      real(dp), allocatable, intent(out) :: y_xrd_temp(:)
      real(dp), allocatable, intent(out) :: energies_xrd(:)
      real(dp), allocatable, intent(out) :: forces_xrd(:, :)
      real(dp), allocatable, intent(in) :: structure_factor_partial(:, :)
      real(dp), allocatable, intent(in) :: structure_factor_partial_temp(:, :)
      real(dp), allocatable, intent(in) :: x_structure_factor(:)
      real(dp), allocatable, intent(in) :: x_structure_factor_temp(:)
      real(dp), allocatable, intent(in) :: sinc_factor_matrix(:, :)
      real(dp), allocatable, intent(in) :: pair_distribution_partial_der(:, :, :)
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      real(dp), intent(in), allocatable :: n_atoms_of_species(:)
      real(dp), intent(out) :: virial_xrd(1:3, 1:3)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: species(:)
      integer, intent(in), allocatable :: nk(:)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      character*8, allocatable, intent(in) :: species_types(:)
      real(dp) :: v_uc
      integer, intent(in) :: n_species
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: ntasks
      integer, intent(out) :: q_beg
      integer, intent(out) :: q_end
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: rank
      integer, intent(out) :: ierr
      real(dp), allocatable :: y_sub(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: i2
      integer :: n_dim_idx
      integer :: n_dim_partial
      real(dp) :: dq
      real(dp) :: f
      real(dp), parameter :: pi = acos(-1.0)
      character*1024 :: filename
      logical :: write_condition
      logical :: overwrite_condition
      logical :: valid_xrd
      logical, intent(in) :: do_derivatives
      logical, intent(in) :: neutron
      logical, intent(in) :: use_matrix_forces
      integer :: xrd_idx
      character*32 :: xrd_output

      type(c_ptr) :: cublas_handle
      type(c_ptr) :: gpu_stream
      type(c_ptr), allocatable :: nk_d(:)
      type(c_ptr), allocatable :: nk_flags_d(:)
      type(c_ptr), allocatable :: nk_flags_sum_d(:)
      type(c_ptr), allocatable :: k_index_d(:)
      type(c_ptr), allocatable :: j2_index_d(:)
      type(c_ptr), allocatable :: rjs_index_d(:)
      type(c_ptr), allocatable :: xyz_k_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_der_d(:)
      integer(c_size_t), allocatable :: st_nk_d(:)
      integer(c_size_t), allocatable :: st_k_index_d(:)
      integer(c_size_t), allocatable :: st_j2_index_d(:)
      integer(c_size_t), allocatable :: st_rjs(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_d(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_der_d(:)
      type(gpu_host_storage_type), intent(inout), allocatable, target :: gpu_host_storage(:)
      logical, intent(in) :: gpu_low_memory

      if (neutron) xrd_idx = params%nd_idx
      if (.not. neutron) xrd_idx = params%xrd_idx

      if (neutron) xrd_output = params%nd_output
      if (.not. neutron) xrd_output = params%xrd_output

      if (neutron) valid_xrd = params%valid_nd
      if (.not. neutron) valid_xrd = params%valid_xrd

      v_uc = dot_product(cross_product(a_box,&
           & b_box), c_box)/(&
           & dfloat(indices(1)*indices(2)&
           &*indices(3)))

      n_dim_partial = n_species*(n_species + 1)/2

      ! Get the XRD from the partial structure factors!
      if (allocated(x_xrd)) deallocate (x_xrd)
      allocate (x_xrd(1:params%structure_factor_n_samples))
      if (allocated(x_xrd_temp)) deallocate (x_xrd_temp)
      allocate (x_xrd_temp(1:params%structure_factor_n_samples))

      x_xrd = x_structure_factor
      x_xrd_temp = x_structure_factor_temp

      if (allocated(y_xrd)) deallocate (y_xrd)
      allocate (y_xrd(1:params%structure_factor_n_samples))
      y_xrd = 0.d0

      ! allocate for the structure factor parameters, so we don't have to look through the horrible list

#ifdef _MPIF90
      allocate (y_xrd_temp(1:params&
           &%structure_factor_n_samples))

      y_xrd_temp = 0.0d0
#endif

      ! Have the same range for the xrd as the structure factor

      allocate (y_sub(1:params&
           &%structure_factor_n_samples))

      if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
         call get_xrd_from_partial_structure_factors(q_beg, q_end, &
              & structure_factor_partial(1:params&
              &%structure_factor_n_samples, 1:n_dim_partial), n_species, species_types&
              &, species, params%xrd_wavelength, params &
              &%xrd_damping, params%xrd_alpha, params &
              &%xrd_method, params%xrd_iwasa, xrd_output, x_xrd(1:params&
              &%structure_factor_n_samples), y_xrd(1:params&
              &%structure_factor_n_samples),&
              & n_atoms_of_species, y_sub, neutron)

      else

         call get_xrd_explicit(n_sites, species_types, n_species,  &
              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
              & params%structure_factor_n_samples,&
              & x_xrd(1:params&
              &%structure_factor_n_samples), y_xrd(1:params&
              &%structure_factor_n_samples), params%r_range_min, params%r_range_max,&
              & params%pair_distribution_rcut, .false., 1&
              &, 1, params%structure_factor_window)

      end if

      !###---   Can calculate the Structure factors related to XRD / Neutron here   ---###!

#ifdef _MPIF90

      call mpi_reduce(y_xrd,&
           & y_xrd_temp, params&
           &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
           & MPI_COMM_WORLD, ierr)
      y_xrd = y_xrd_temp
      deallocate (y_xrd_temp)

      call mpi_bcast(y_xrd, params&
           &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, 0,&
           & MPI_COMM_WORLD, ierr)

#endif

      ! Already preprocessed the xrd !
      ! --- Preprocess the structure factor according to the output --- !
      ! if ( trim( params%xrd_output ) == "q*i(q)" )then
      !    y_xrd = 2.d0 * pi * x_xrd * ( y_xrd )
      ! end if

      if (allocated(sinc_factor_matrix)) then
         if (valid_xrd) then

            allocate (energies_xrd(1:n_sites))
            energies_xrd = 0.d0

            if (valid_xrd .and. allocated(params%exp_energy_scales)) then
               call get_energy_scale(params%do_md, params%do_mc,&
                    & md_istep, params%md_nsteps, mc_istep, params &
                    &%mc_nsteps, params&
                    &%exp_energy_scales_initial(xrd_idx), params&
                    &%exp_energy_scales_final(xrd_idx), params&
                    &%exp_energy_scales(xrd_idx))

               call get_exp_energies(params%exp_energy_scales(xrd_idx), params%exp_data(xrd_idx)%y&
                    &, y_xrd,&
                    & params%structure_factor_n_samples, n_sites,&
                    & energies_xrd(i_beg:i_end), params%exp_data(xrd_idx)%w)

               if (params%do_forces .and. params%exp_forces) then

                  allocate (forces_xrd(1:3, 1:n_sites))
                  forces_xrd = 0.d0
                  virial_xrd = 0.d0
                  n_dim_idx = 1
                  outerf: do j = 1, n_species
                     do k = 1, n_species

                        if (j > k) cycle

                        if (j == k) f = 1.d0
                        if (j /= k) f = 2.d0

                        if (use_matrix_forces) then

                           call get_structure_factor_forces_matrix(i_beg, i_end, n_sites, params%exp_energy_scales(xrd_idx),&
                                & x_xrd, params%exp_data(xrd_idx)%y,&
                                & forces_xrd, virial_xrd,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), species_types, species(i_beg:i_end), rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_xrd(1:params&
                                &%structure_factor_n_samples), y_xrd(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .true., xrd_output, n_atoms_of_species, neutron, rank, cublas_handle, gpu_stream,&
                & nk(n_dim_idx), nk_d(n_dim_idx), k_index_d(n_dim_idx), j2_index_d(n_dim_idx), xyz_k_d(n_dim_idx), pair_distribution_partial_d(n_dim_idx), pair_distribution_partial_der_d(n_dim_idx), &
                & st_nk_d(n_dim_idx), st_k_index_d(n_dim_idx), st_j2_index_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), &
                                gpu_host_storage(n_dim_idx), gpu_low_memory, params%exp_data(xrd_idx)%w)
                        else
                           call get_structure_factor_forces(n_sites, params%exp_energy_scales(xrd_idx),&
                                & x_xrd, params%exp_data(xrd_idx)%y,&
                                & forces_xrd, virial_xrd,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                                & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_xrd(1:params&
                                &%structure_factor_n_samples), y_xrd(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .true., xrd_output, n_atoms_of_species, neutron, rank, params%exp_data(xrd_idx)%w)
                        end if

                        n_dim_idx = n_dim_idx + 1

                        if (n_dim_idx > n_dim_partial) exit outerf

                     end do
                  end do outerf
               end if
            end if
         end if
      end if

      call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

      if (rank == 0 .and. (params%write_xrd .or. params%write_nd) .and. write_condition) then

         call get_overwrite_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & overwrite_condition)

         if (.not. neutron) then
            write (filename, '(A)')&
                 & 'xrd_prediction.dat'
            call write_exp_datan(x_xrd_temp(1:params%structure_factor_n_samples),&
                 & y_xrd(1:params&
                 &%structure_factor_n_samples),&
                 & overwrite_condition, filename, "xrd: units of "//&
                 & trim(params%q_units)//" output: "//trim(xrd_output))
         else
            write (filename, '(A)')&
                 & 'nd_prediction.dat'
            call write_exp_datan(x_xrd_temp(1:params%structure_factor_n_samples),&
                 & y_xrd(1:params&
                 &%structure_factor_n_samples),&
                 & overwrite_condition, filename, "nd: units of "//&
                 & trim(params%q_units)//" output: "//trim(xrd_output))
         end if

      end if

      ! if ( trim(params%xrd_output) == "q*i(q)" .or. trim(params%xrd_output) == "q*F(q)")then
      !    ! output q * i(q) === q * F_x(q)
      !    do l = q_beg, q_end
      !       y_xrd(l) = x_xrd_temp(l) * ( y_xrd(l) - y_sub(l) )
      !    end do

      !    elseif( trim(output) == "F(q)" .or. trim(output) == "i(q)")
      !       ! do nothing,
      !       ! Output the total scattering functon, i(q) === F_x(q)

      if (allocated(y_sub)) deallocate (y_sub)

   end subroutine calculate_xrd
#else
   subroutine calculate_xrd(params, x_xrd, x_xrd_temp,&
        & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
        & structure_factor_partial, structure_factor_partial_temp,&
        & n_species, species_types, n_atoms_of_species,&
        & n_sites, a_box, b_box, c_box, indices, md_istep, mc_istep, i_beg,&
        & i_end, j_beg, j_end, ierr, rjs, xyz, neighbors_list, n_neigh,&
        & neighbor_species, species, rank, q_beg, q_end, ntasks,&
        & sinc_factor_matrix, do_derivatives, pair_distribution_partial_der,&
        & energies_xrd, forces_xrd, virial_xrd, neutron, use_matrix_forces)
      implicit none
      type(input_parameters), intent(inout) :: params
      real(dp), allocatable, intent(out) :: x_xrd(:)
      real(dp), allocatable, intent(out) :: x_xrd_temp(:)
      real(dp), allocatable, intent(out) :: y_xrd(:)
      real(dp), allocatable, intent(out) :: y_xrd_temp(:)
      real(dp), allocatable, intent(out) :: energies_xrd(:)
      real(dp), allocatable, intent(out) :: forces_xrd(:, :)
      real(dp), allocatable, intent(in) :: structure_factor_partial(:, :)
      real(dp), allocatable, intent(in) :: structure_factor_partial_temp(:, :)
      real(dp), allocatable, intent(in) :: x_structure_factor(:)
      real(dp), allocatable, intent(in) :: x_structure_factor_temp(:)
      real(dp), allocatable, intent(in) :: sinc_factor_matrix(:, :)
      real(dp), allocatable, intent(in) :: pair_distribution_partial_der(:, :, :)
      real(dp), intent(in), allocatable :: rjs(:)
      real(dp), intent(in), allocatable :: xyz(:, :)
      real(dp), intent(in), allocatable :: n_atoms_of_species(:)
      real(dp), intent(out) :: virial_xrd(1:3, 1:3)
      integer, intent(in), allocatable :: neighbors_list(:)
      integer, intent(in), allocatable :: n_neigh(:)
      integer, intent(in), allocatable :: neighbor_species(:)
      integer, intent(in), allocatable :: species(:)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      character*8, allocatable, intent(in) :: species_types(:)
      real(dp) :: v_uc
      integer, intent(in) :: n_species
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: ntasks
      integer, intent(out) :: q_beg
      integer, intent(out) :: q_end
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: rank
      integer, intent(out) :: ierr
      real(dp), allocatable :: y_sub(:)
      real(dp), allocatable :: lp(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: i2
      integer :: n_dim_idx
      integer :: n_dim_partial
      real(dp) :: dq
      real(dp) :: f
      real(dp), parameter :: pi = acos(-1.0)
      character*1024 :: filename
      logical :: write_condition
      logical :: overwrite_condition
      logical :: valid_xrd
      logical, intent(in) :: do_derivatives
      logical, intent(in) :: neutron
      logical, intent(in) :: use_matrix_forces
      integer :: xrd_idx
      character*32 :: xrd_output

      if (neutron) xrd_idx = params%nd_idx
      if (.not. neutron) xrd_idx = params%xrd_idx

      if (neutron) xrd_output = params%nd_output
      if (.not. neutron) xrd_output = params%xrd_output

      if (neutron) valid_xrd = params%valid_nd
      if (.not. neutron) valid_xrd = params%valid_xrd

      v_uc = dot_product(cross_product(a_box,&
           & b_box), c_box)/(&
           & dfloat(indices(1)*indices(2)&
           &*indices(3)))

      n_dim_partial = n_species*(n_species + 1)/2

      ! Get the XRD from the partial structure factors!
      if (allocated(x_xrd)) deallocate (x_xrd)
      allocate (x_xrd(1:params%structure_factor_n_samples))
      if (allocated(x_xrd_temp)) deallocate (x_xrd_temp)
      allocate (x_xrd_temp(1:params%structure_factor_n_samples))

      x_xrd = x_structure_factor
      x_xrd_temp = x_structure_factor_temp

      if (allocated(y_xrd)) deallocate (y_xrd)
      allocate (y_xrd(1:params%structure_factor_n_samples))
      y_xrd = 0.d0

      ! allocate for the structure factor parameters, so we don't have to look through the horrible list

#ifdef _MPIF90
      allocate (y_xrd_temp(1:params&
           &%structure_factor_n_samples))

      y_xrd_temp = 0.0d0
#endif

      ! Have the same range for the xrd as the structure factor

      allocate (y_sub(1:params&
           &%structure_factor_n_samples))

      if (params%structure_factor_from_pdf .and. params%pair_distribution_partial) then
         call get_xrd_from_partial_structure_factors(q_beg, q_end, &
              & structure_factor_partial(1:params&
              &%structure_factor_n_samples, 1:n_dim_partial), n_species, species_types&
              &, species, params%xrd_wavelength, params &
              &%xrd_damping, params%xrd_alpha, params &
              &%xrd_method, params%xrd_iwasa, xrd_output, x_xrd(1:params&
              &%structure_factor_n_samples), y_xrd(1:params&
              &%structure_factor_n_samples),&
              & n_atoms_of_species, y_sub, neutron)

      else

         call get_xrd_explicit(n_sites, species_types, n_species,  &
              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
              & params%structure_factor_n_samples,&
              & x_xrd(1:params&
              &%structure_factor_n_samples), y_xrd(1:params&
              &%structure_factor_n_samples), params%r_range_min, params%r_range_max,&
              & params%pair_distribution_rcut, .false., 1&
              &, 1, params%structure_factor_window)

      end if

      !###---   Can calculate the Structure factors related to XRD / Neutron here   ---###!

#ifdef _MPIF90

      call mpi_reduce(y_xrd,&
           & y_xrd_temp, params&
           &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
           & MPI_COMM_WORLD, ierr)
      y_xrd = y_xrd_temp
      deallocate (y_xrd_temp)

      call mpi_bcast(y_xrd, params&
           &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, 0,&
           & MPI_COMM_WORLD, ierr)

#endif

!     ---   The Lorentz-polarization factor   --- !
!
!     A weight per q sample, applied here rather than inside
!     get_xrd_from_partial_structure_factors so that it lands once on the
!     whole pattern -- each rank builds only its own q slice, and this is
!     after the reduce. It multiplies whatever xrd_output asked for, and so
!     covers the get_xrd_explicit branch above too.
!
!     Everything downstream sees the weighted pattern: the energy, the
!     residual (y - y_exp) the forces are built from, and the file written at
!     the end. What the forces additionally need is the same weight on
!     dy/dr, which is why lp is handed to the two force routines below.
!
!     x_xrd holds 2 sin(theta)/lambda, the grid get_lorentz_polarization
!     expects. lp stays allocated and equal to one when the factor is off, so
!     the calls below need no second form; multiplying by 1.0d0 is exact.
      allocate (lp(1:params%structure_factor_n_samples))
      lp = 1.d0

      if (params%xrd_lorentz_polarization) then
         call get_lorentz_polarization(x_xrd(1:params%structure_factor_n_samples), &
              & params%xrd_wavelength, params%xrd_lp_polarization, &
              & params%xrd_lp_sin_theta_min, lp)
         y_xrd = lp*y_xrd
      end if

      ! Already preprocessed the xrd !
      ! --- Preprocess the structure factor according to the output --- !
      ! if ( trim( params%xrd_output ) == "q*i(q)" )then
      !    y_xrd = 2.d0 * pi * x_xrd * ( y_xrd )
      ! end if

      if (allocated(sinc_factor_matrix)) then
         if (valid_xrd) then

            allocate (energies_xrd(1:n_sites))
            energies_xrd = 0.d0

            if (valid_xrd .and. allocated(params%exp_energy_scales)) then
               call get_energy_scale(params%do_md, params%do_mc,&
                    & md_istep, params%md_nsteps, mc_istep, params&
                    &%mc_nsteps, params%exp_energy_scales_initial(xrd_idx), &
                    & params%exp_energy_scales_final(xrd_idx), &
                    & params%exp_energy_scales(xrd_idx))

               call get_exp_energies(params%exp_energy_scales(xrd_idx), params%exp_data(xrd_idx)%y&
                    &, y_xrd,&
                    & params%structure_factor_n_samples, n_sites,&
                    & energies_xrd(i_beg:i_end), params%exp_data(xrd_idx)%w)

               if (params%do_forces .and. params%exp_forces) then

                  allocate (forces_xrd(1:3, 1:n_sites))
                  forces_xrd = 0.d0

                  n_dim_idx = 1
                  outerf: do j = 1, n_species
                     do k = 1, n_species

                        if (j > k) cycle

                        if (j == k) f = 1.d0
                        if (j /= k) f = 2.d0

                        if (use_matrix_forces) then

                           call get_structure_factor_forces_matrix(n_sites, params%exp_energy_scales(xrd_idx),&
                                & x_xrd, params%exp_data(xrd_idx)%y,&
                                & forces_xrd, virial_xrd,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                                & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_xrd(1:params&
                                &%structure_factor_n_samples), y_xrd(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .true., xrd_output, n_atoms_of_species, neutron, rank, lp, w=params%exp_data(xrd_idx)%w)
                        else
                           call get_structure_factor_forces(n_sites, params%exp_energy_scales(xrd_idx),&
                                & x_xrd, params%exp_data(xrd_idx)%y,&
                                & forces_xrd, virial_xrd,&
                                & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                                & neighbor_species(j_beg:j_end), species_types, rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                                &%r_range_min, params%r_range_max, params&
                                &%pair_distribution_n_samples, params&
                                &%structure_factor_n_samples, n_species,&
                                & x_xrd(1:params&
                                &%structure_factor_n_samples), y_xrd(1:params&
                                &%structure_factor_n_samples), params%pair_distribution_rcut&
                                &, j, k, pair_distribution_partial_der(1:params &
                                &%pair_distribution_n_samples, n_dim_idx,&
                                & j_beg:j_end), params%pair_distribution_partial,&
                                & params%pair_distribution_kde_sigma,&
                                & 4.d0*pi*f*((n_atoms_of_species(j)*&
                                & n_atoms_of_species(k))/dfloat(n_sites)/&
                                & dfloat(n_sites))*(dfloat(n_sites)/&
                                & v_uc), sinc_factor_matrix, n_dim_idx,&
                                & .true., xrd_output, n_atoms_of_species, neutron, rank, lp)
                        end if

                        n_dim_idx = n_dim_idx + 1

                        if (n_dim_idx > n_dim_partial) exit outerf

                     end do
                  end do outerf
               end if
            end if
         end if
      end if

      call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

      if (rank == 0 .and. (params%write_xrd .or. params%write_nd) .and. write_condition) then

         call get_overwrite_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & overwrite_condition)

         if (.not. neutron) then
            write (filename, '(A)')&
                 & 'xrd_prediction.dat'
            call write_exp_datan(x_xrd_temp(1:params%structure_factor_n_samples),&
                 & y_xrd(1:params&
                 &%structure_factor_n_samples),&
                 & overwrite_condition, filename, "xrd: units of "//&
                 & trim(params%q_units)//" output: "//trim(xrd_output))
         else
            write (filename, '(A)')&
                 & 'nd_prediction.dat'
            call write_exp_datan(x_xrd_temp(1:params%structure_factor_n_samples),&
                 & y_xrd(1:params&
                 &%structure_factor_n_samples),&
                 & overwrite_condition, filename, "nd: units of "//&
                 & trim(params%q_units)//" output: "//trim(xrd_output))
         end if

      end if

      ! if ( trim(params%xrd_output) == "q*i(q)" .or. trim(params%xrd_output) == "q*F(q)")then
      !    ! output q * i(q) === q * F_x(q)
      !    do l = q_beg, q_end
      !       y_xrd(l) = x_xrd_temp(l) * ( y_xrd(l) - y_sub(l) )
      !    end do

      !    elseif( trim(output) == "F(q)" .or. trim(output) == "i(q)")
      !       ! do nothing,
      !       ! Output the total scattering functon, i(q) === F_x(q)

      if (allocated(y_sub)) deallocate (y_sub)
      if (allocated(lp)) deallocate (lp)

   end subroutine calculate_xrd
#endif

   subroutine finalize_xrd(params, x_xrd, x_xrd_temp,&
        & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
        & structure_factor_partial, structure_factor_partial_temp)
      implicit none
      type(input_parameters), intent(in) :: params
      real(dp), allocatable, intent(inout) :: x_xrd(:)
      real(dp), allocatable, intent(inout) :: x_xrd_temp(:)
      real(dp), allocatable, intent(inout) :: y_xrd(:)
      real(dp), allocatable, intent(inout) :: y_xrd_temp(:)
      real(dp), allocatable, intent(inout) :: structure_factor_partial(:, :)
      real(dp), allocatable, intent(inout) :: structure_factor_partial_temp(:, :)
      real(dp), allocatable, intent(inout) :: x_structure_factor(:)
      real(dp), allocatable, intent(inout) :: x_structure_factor_temp(:)

      if (allocated(x_xrd)) deallocate (x_xrd)
      if (allocated(x_xrd_temp)) deallocate (x_xrd_temp)
      if (allocated(y_xrd)) deallocate (y_xrd)
      if (allocated(y_xrd_temp)) deallocate (y_xrd_temp)
      if (allocated(x_structure_factor)) deallocate (x_structure_factor)
      if (allocated(x_structure_factor_temp)) deallocate (x_structure_factor_temp)
      if (allocated(structure_factor_partial)) deallocate (structure_factor_partial)
      if (allocated(structure_factor_partial_temp)) deallocate (structure_factor_partial_temp)
   end subroutine finalize_xrd

!  The XRD (or ND) pattern from the Debye scattering equation, with the
!  energies, forces and virial that fitting it against experiment produces.
!
!  This is the counterpart of calculate_xrd, and takes the other of the two
!  routes to the same observable. calculate_xrd needs a pair distribution to
!  have been binned and Fourier transformed first; this one goes straight from
!  the positions, so params%do_pair_distribution and params%do_structure_factor
!  can both be off and nothing above this routine has to have run.
!
!  What arrives from get_xrd_debye is the raw intensity per atom,
!
!     I(q) = 1/N sum_ij f_i f_j sinc(q r_ij) w(r_ij),
!
!  self term included. Everything after the MPI reduction is one affine map
!  per q sample,
!
!     y(l) = lp(l) * ( a(l) * ( I(l) - b(l) ) + c(l) ),
!
!  chosen so that the four xrd_output conventions here mean exactly what they
!  mean in get_xrd_from_partial_structure_factors, whose y is the interference
!  part I - sum_i c_i f_i^2. Because the map is affine in I with coefficients
!  that do not depend on the positions, the gradient is the same map without
!  the constants: dy(l)/dr = lp(l) a(l) dI(l)/dr. That is what keeps the
!  forces below consistent with the energy no matter which output, and which
!  the pdf/sf route cannot say as cheaply.
   subroutine calculate_xrd_debye(params, x_xrd, x_xrd_temp, y_xrd, y_xrd_temp, &
        & n_sites, positions, species, md_istep, mc_istep, i_beg, i_end, &
        & ierr, rank, do_derivatives, energies_xrd, forces_xrd, virial_xrd, neutron)
      implicit none
      type(input_parameters), intent(inout) :: params
      real(dp), allocatable, intent(out) :: x_xrd(:)
      real(dp), allocatable, intent(out) :: x_xrd_temp(:)
      real(dp), allocatable, intent(out) :: y_xrd(:)
      real(dp), allocatable, intent(out) :: y_xrd_temp(:)
      real(dp), allocatable, intent(out) :: energies_xrd(:)
      real(dp), allocatable, intent(out) :: forces_xrd(:, :)
      real(dp), intent(out) :: virial_xrd(1:3, 1:3)
      real(dp), intent(in), allocatable :: positions(:, :)
      integer, intent(in), allocatable :: species(:)
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: rank
      integer, intent(out) :: ierr
      logical, intent(in) :: do_derivatives
      logical, intent(in) :: neutron
      real(dp), allocatable :: y_debye(:)
      real(dp), allocatable :: y_debye_der(:, :, :)
      real(dp), allocatable :: f_table(:, :)
      real(dp), allocatable :: channel(:)
      real(dp), allocatable :: self_term(:)
      real(dp), allocatable :: lp(:)
      real(dp), allocatable :: prefactor(:)
      real(dp), allocatable :: concentration(:)
      real(dp) :: this_force(1:3)
      real(dp) :: sth
      real(dp) :: dq
      real(dp) :: energy_scale
!     Double-precision pi, unlike the acos(-1.0) the rest of this module uses
!     -- see the note in get_xrd_debye.
      real(dp), parameter :: pi = acos(-1.0d0)
      integer :: i
      integer :: l
      integer :: k1
      integer :: k2
      integer :: n_species
      integer :: n_samples
      integer :: xrd_idx
      character*32 :: xrd_output
      character*1024 :: filename
      logical :: write_condition
      logical :: overwrite_condition
      logical :: valid_xrd
      logical :: want_forces

      if (neutron) then
         xrd_idx = params%nd_idx
         xrd_output = params%nd_output
         valid_xrd = params%valid_nd
      else
         xrd_idx = params%xrd_idx
         xrd_output = params%xrd_output
         valid_xrd = params%valid_xrd
      end if

      n_species = size(params%species_types)
      n_samples = params%structure_factor_n_samples
      virial_xrd = 0.d0
      ierr = 0

!     The q grid, built exactly as calculate_structure_factor builds it so
!     that a deck can be switched between the two routes without touching
!     q_range_min/max or q_units. x_xrd is the working grid,
!     x = 2 sin(theta)/lambda; x_xrd_temp keeps the user's units for output.
      call linspace(x_xrd, params%q_range_min, params%q_range_max, n_samples, dq)
      call linspace(x_xrd_temp, params%q_range_min, params%q_range_max, n_samples, dq)

      if (trim(params%q_units) == "xrd" .or. params%q_units == "twotheta") then
         do i = 1, n_samples
            x_xrd(i) = 2.d0*sin(pi*x_xrd(i)/180.d0/2.d0)/params%xrd_wavelength
         end do
      elseif (trim(params%q_units) == "saxs" .or. params%q_units == "q") then
         do i = 1, n_samples
            x_xrd(i) = x_xrd(i)/2.d0/pi
         end do
      end if

      want_forces = do_derivatives .and. params%do_forces .and. valid_xrd &
           & .and. allocated(params%exp_energy_scales)

      call get_xrd_debye(n_sites, positions, i_beg, i_end, n_species, species, &
           & params%species_types, x_xrd, params%structure_factor_window, &
           & params%xrd_rcut, want_forces, neutron, y_debye, y_debye_der)

#ifdef _MPIF90
!     Each rank summed the pairs of the atoms it owns; the pattern is the sum
!     over ranks. The derivative slices are already complete per atom and are
!     deliberately not reduced here -- the driver reduces the forces built
!     from them.
      allocate (y_xrd_temp(1:n_samples))
      y_xrd_temp = 0.d0
      call mpi_reduce(y_debye, y_xrd_temp, n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
           & MPI_COMM_WORLD, ierr)
      y_debye = y_xrd_temp
      deallocate (y_xrd_temp)
      call mpi_bcast(y_debye, n_samples, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

!     Concentrations and the two species averages the output conventions need:
!     the self term sum_i c_i f_i(q)^2 that the interference part excludes,
!     and sth = sum_i c_i f_i(q), whose square normalises F(q).
      allocate (concentration(1:n_species))
      concentration = 0.d0
      do i = 1, n_sites
         concentration(species(i)) = concentration(species(i)) + 1.d0
      end do
      concentration = concentration/dfloat(n_sites)

      call get_form_factor_table(n_species, params%species_types, x_xrd, neutron, f_table)

      allocate (self_term(1:n_samples))
      allocate (channel(1:n_samples))
      self_term = 0.d0
      channel = 1.d0

      allocate (y_xrd(1:n_samples))

      if (trim(xrd_output) == "xrd") then
!        The raw intensity per atom, self term included. Nothing to do.
         y_xrd = y_debye
      else
         do l = 1, n_samples
            sth = 0.d0
            do i = 1, n_species
               self_term(l) = self_term(l) + concentration(i)*f_table(l, i)*f_table(l, i)
               sth = sth + concentration(i)*f_table(l, i)
            end do
            channel(l) = 1.d0/sth**2
         end do

         if (trim(xrd_output) == "q*i(q)" .or. trim(xrd_output) == "q*F(q)") then
!           q here is the wavevector 2 pi x, matching what the pdf/sf route
!           multiplies by.
            channel = 2.d0*pi*x_xrd*channel
            y_xrd = channel*(y_debye - self_term)
         elseif (trim(xrd_output) == "F(q)" .or. trim(xrd_output) == "i(q)") then
            y_xrd = channel*(y_debye - self_term) + 1.d0
         else
            if (rank == 0) write (*, *) "WARNING: unrecognised xrd_output '"//trim(xrd_output)// &
                 & "' for the Debye route; falling back to the raw intensity"
            channel = 1.d0
            y_xrd = y_debye
         end if
      end if

!     The Lorentz-polarization factor multiplies whatever pattern was asked
!     for. It is a per-q weight, so it scales the gradient channel by the same
!     number and the forces stay exact. It is meant for xrd_output = 'xrd',
!     where the prediction really is a raw powder intensity.
      if (params%xrd_lorentz_polarization) then
         allocate (lp(1:n_samples))
         call get_lorentz_polarization(x_xrd, params%xrd_wavelength, &
              & params%xrd_lp_polarization, params%xrd_lp_sin_theta_min, lp)
         y_xrd = lp*y_xrd
         channel = lp*channel
         deallocate (lp)
      end if

!     ---   Energies and forces against the experimental pattern   --- !
      if (valid_xrd .and. allocated(params%exp_energy_scales)) then
         allocate (energies_xrd(1:n_sites))
         energies_xrd = 0.d0

         call get_energy_scale(params%do_md, params%do_mc, md_istep, params%md_nsteps, &
              & mc_istep, params%mc_nsteps, params%exp_energy_scales_initial(xrd_idx), &
              & params%exp_energy_scales_final(xrd_idx), params%exp_energy_scales(xrd_idx))

         energy_scale = params%exp_energy_scales(xrd_idx)

         call get_exp_energies(energy_scale, params%exp_data(xrd_idx)%y, y_xrd, &
              & n_samples, n_sites, energies_xrd(i_beg:i_end), params%exp_data(xrd_idx)%w)

         if (want_forces) then
            allocate (forces_xrd(1:3, 1:n_sites))
            forces_xrd = 0.d0

!           E = 1/2 s sum_l ( y_l - y_exp_l )^2, so the force on atom a is
!           -dE/dr_a = -s sum_l ( y_l - y_exp_l ) dy_l/dr_a, with
!           dy_l/dr_a = channel(l) * y_debye_der(:,l,a). Only the atoms this
!           rank owns have a complete derivative, and only those are written.
            allocate (prefactor(1:n_samples))
            prefactor = energy_scale*channel &
                       & *(y_xrd - params%exp_data(xrd_idx)%y(1:n_samples))

            do i = i_beg, i_end
               do k1 = 1, 3
                  this_force(k1) = -dot_product(y_debye_der(k1, 1:n_samples, i), prefactor)
               end do
               forces_xrd(1:3, i) = this_force(1:3)

!              The Debye energy depends on the positions only through the
!              interatomic distances of atoms in the cell -- no periodic
!              images enter the sum -- so sum_i F_i (x) r_i telescopes into
!              the pair form sum_{i<j} F_ij (x) r_ij and is the virial.
               do k1 = 1, 3
                  do k2 = 1, 3
                     virial_xrd(k1, k2) = virial_xrd(k1, k2) + &
                          & 0.5d0*(this_force(k1)*positions(k2, i) + this_force(k2)*positions(k1, i))
                  end do
               end do
            end do

            deallocate (prefactor)
         end if
      end if

!     ---   Write the prediction   --- !
      call get_write_condition(params%do_mc, params%do_md, mc_istep, md_istep, &
           & params%write_xyz, write_condition)

      if (rank == 0 .and. (params%write_xrd .or. params%write_nd) .and. write_condition) then
         call get_overwrite_condition(params%do_mc, params%do_md, mc_istep, md_istep, &
              & params%write_xyz, overwrite_condition)

         if (.not. neutron) then
            write (filename, '(A)') 'xrd_prediction.dat'
            call write_exp_datan(x_xrd_temp, y_xrd, overwrite_condition, filename, &
                 & "xrd (debye): units of "//trim(params%q_units)//" output: "//trim(xrd_output))
         else
            write (filename, '(A)') 'nd_prediction.dat'
            call write_exp_datan(x_xrd_temp, y_xrd, overwrite_condition, filename, &
                 & "nd (debye): units of "//trim(params%q_units)//" output: "//trim(xrd_output))
         end if
      end if

      deallocate (y_debye)
      if (allocated(y_debye_der)) deallocate (y_debye_der)
      deallocate (f_table)
      deallocate (channel)
      deallocate (self_term)
      deallocate (concentration)

   end subroutine calculate_xrd_debye

   subroutine finalize_xrd_debye(x_xrd, x_xrd_temp, y_xrd, y_xrd_temp)
      implicit none
      real(dp), allocatable, intent(inout) :: x_xrd(:)
      real(dp), allocatable, intent(inout) :: x_xrd_temp(:)
      real(dp), allocatable, intent(inout) :: y_xrd(:)
      real(dp), allocatable, intent(inout) :: y_xrd_temp(:)

      if (allocated(x_xrd)) deallocate (x_xrd)
      if (allocated(x_xrd_temp)) deallocate (x_xrd_temp)
      if (allocated(y_xrd)) deallocate (y_xrd)
      if (allocated(y_xrd_temp)) deallocate (y_xrd_temp)
   end subroutine finalize_xrd_debye

#ifdef _GPU
   subroutine gpu_copy_pdf(n_samples, pdf_d, pdf, gpu_stream)
      implicit none
      integer :: n_samples
      real(dp), intent(inout), target :: pdf(1:n_samples)
      type(c_ptr) :: pdf_d
      type(c_ptr) :: gpu_stream
      integer(c_size_t) :: size

      size = int(n_samples, c_size_t)*c_double

      call cpy_dtoh(pdf_d, c_loc(pdf), size, gpu_stream)

   end subroutine gpu_copy_pdf
#endif

   subroutine estimate_device_memory_usage(n_sites, n_pairs, nk, n_samples, n_samples_sf, total, standard)
      implicit none
      integer :: n_sites
      integer :: nk
      integer :: n_samples
      integer :: n_samples_sf
      integer :: n_pairs
      real(dp), intent(inout) :: total
      real(dp) :: total_exp = 0.d0
      real(dp) :: total_standard = 0.d0
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
      logical :: standard

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

      if (standard) then
         write (*, '(A)') "\n > Estimating memory for normal allocations \n"
         total_standard = total_standard + neigh_list
         write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device neigh_list_d  ", neigh_list, " Gb\n"

         total_standard = total_standard + neigh_list
         write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device neigh_spec_d  ", neigh_list, " Gb\n"

         total_standard = total_standard + neigh_list*2.d0*3.d0
         write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device xyz_d  ", neigh_list*2.d0*3.d0, " Gb\n"

         total_standard = total_standard + neigh_list*2.d0
         write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device rjs_d  ", neigh_list*2.d0, " Gb\n"

         total_standard = total_standard + neigh_list*2.d0
         write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device rjs_d  ", neigh_list*2.d0, " Gb\n"

         total_standard = total_standard + forces/3.d0
         write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device species_d  ", forces/3.d0, " Gb\n"
      end if

      write (*, '(A)') "\n > Estimating memory xrd and pdf calculation\n"
      total_exp = total_exp + nk_int
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device k_index_d  ", nk_int, " Gb\n"

      total_exp = total_exp + 3.d0*Gk
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device      Gk_d  ", Gk, " Gb\n"
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device     Gka_d  ", Gk, " Gb\n"

      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device par_pdf_d  ", Gk, " Gb\n"

      total_exp = total_exp + dermat
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device  dermat_d  ", dermat, " Gb\n"

      total_exp = total_exp + fi
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device      fi_d  ", fi, " Gb\n"

      total_exp = total_exp + xyz
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device     xyz_d  ", xyz, " Gb\n"

      total_exp = total_exp + forces
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device   forces_d  ", forces, " Gb\n"

      total_exp = total_exp + pref
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device    pref_d  ", pref, " Gb\n"

      total_exp = total_exp + pref
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device  scat_f_d  ", pref, " Gb\n"

      total_exp = total_exp + sf
      write (*, '(A,1X,F7.4,1X,A)') "Gb/core: device  sinc_f_d  ", sf, " Gb\n"

      total = total + total_exp + total_standard
      write (*, '(A,1X,F7.4,1X,A)') "\nTotal device memory usage in block:", total_standard + total_exp, " Gb\n"
      write (*, '(A,1X,F7.4,1X,A)') "\nTotal device memory usage:", total, " Gb\n"

   end subroutine estimate_device_memory_usage

   subroutine get_n_atoms_of_species(n_atoms_of_species, n_sites, species, n_species)
      implicit none
      integer, intent(in) :: n_sites
      integer, intent(in) :: n_species
      integer, allocatable, intent(in) :: species(:)
      real(dp), allocatable, intent(out) :: n_atoms_of_species(:)
      integer :: j
      integer :: i2

      if (allocated(n_atoms_of_species)) deallocate (n_atoms_of_species)
      allocate (n_atoms_of_species(1:n_species))

      do j = 1, n_species
         n_atoms_of_species(j) = 0.d0
         do i2 = 1, n_sites
            if (species(i2) == j) then
               n_atoms_of_species(j) = n_atoms_of_species(j) + 1.d0
            end if
         end do
      end do
   end subroutine get_n_atoms_of_species

   subroutine setup_batched_pair_distribution(r_min, r_max, r_cut, n_samples, x, dV)
      implicit none
      real(dp), intent(in) :: r_min
      real(dp), intent(in) :: r_max
      real(dp), intent(in) :: r_cut
      integer, intent(in) :: n_samples
      real(dp), allocatable, intent(out) :: x(:)
      real(dp), allocatable, intent(out) :: dV(:)

      if (allocated(x)) deallocate (x)
      if (allocated(dV)) deallocate (dV)

      allocate (x(1:n_samples))
      allocate (dV(1:n_samples))

      call setup_pdf_arrays(r_min, r_max, r_cut, n_samples, x, dV)

   end subroutine setup_batched_pair_distribution

#ifdef _GPU
   subroutine gpu_malloc_neighbors(gpu_neigh, n_sites, n_pairs,&
        & n_neigh, species, neighbor_species, neighbors_list, rjs, &
        & xyz, gpu_stream, rank)
      implicit none
      type(gpu_neigh_storage_type) :: gpu_neigh
      integer, intent(in):: n_sites
      integer, intent(in):: n_pairs
      integer, intent(in):: rank
      integer, intent(in), target :: n_neigh(:)
      integer, intent(in), target :: species(:)
      integer, intent(in), target :: neighbor_species(:)
      integer, intent(in), target :: neighbors_list(:)
      real(dp), intent(in), target :: rjs(:)
      real(dp), intent(in), target :: xyz(:, :)

      integer, allocatable, target :: n_neigh_temp(:)
      integer, allocatable, target :: species_temp(:)
      integer, allocatable, target :: neighbor_species_temp(:)
      integer, allocatable, target :: neighbors_list_temp(:)
      real(dp), allocatable, target :: rjs_temp(:)
      real(dp), allocatable, target :: xyz_temp(:, :)

      integer(c_size_t) :: st_n_sites_int
      integer(c_size_t) :: st_n_atom_pairs_int
      integer(c_size_t) :: st_n_atom_pairs_double
      type(c_ptr) :: gpu_stream
      integer :: i

      st_n_sites_int = int(n_sites, c_size_t)*c_int
      gpu_neigh%st_n_neigh_d = st_n_sites_int
      gpu_neigh%st_species_d = st_n_sites_int
      call gpu_malloc_async(gpu_neigh%n_neigh_d, st_n_sites_int, gpu_stream)
      call cpy_htod(c_loc(n_neigh), gpu_neigh%n_neigh_d, st_n_sites_int, gpu_stream)
      call gpu_malloc_async(gpu_neigh%species_d, st_n_sites_int, gpu_stream)
      call cpy_htod(c_loc(species), gpu_neigh%species_d, st_n_sites_int, gpu_stream)

      st_n_atom_pairs_int = int(n_pairs, c_size_t)*c_int
      gpu_neigh%st_neighbor_species_d = st_n_atom_pairs_int
      gpu_neigh%st_neighbors_list_d = st_n_atom_pairs_int
      call gpu_malloc_async(gpu_neigh%neighbor_species_d, st_n_atom_pairs_int, gpu_stream)
      call cpy_htod(c_loc(neighbor_species), gpu_neigh%neighbor_species_d, st_n_atom_pairs_int, gpu_stream)
      call gpu_malloc_async(gpu_neigh%neighbors_list_d, st_n_atom_pairs_int, gpu_stream)
      call cpy_htod(c_loc(neighbors_list), gpu_neigh%neighbors_list_d, st_n_atom_pairs_int, gpu_stream)

      st_n_atom_pairs_double = int(n_pairs, c_size_t)*c_double
      gpu_neigh%st_rjs_d = st_n_atom_pairs_double
      gpu_neigh%st_xyz_d = 3*st_n_atom_pairs_double

      call gpu_malloc_async(gpu_neigh%rjs_d, st_n_atom_pairs_double, gpu_stream)
      call cpy_htod(c_loc(rjs), gpu_neigh%rjs_d, st_n_atom_pairs_double, gpu_stream)
      call gpu_malloc_async(gpu_neigh%xyz_d, 3*st_n_atom_pairs_double, gpu_stream)
      call cpy_htod(c_loc(xyz), gpu_neigh%xyz_d, 3*st_n_atom_pairs_double, gpu_stream)

   if (debug_gpu_batches) print *, "-- Rank ", rank, " ", " malloc neighbors: n_sites_temp = ", n_sites, " n_pairs_temp = ", n_pairs

   end subroutine gpu_malloc_neighbors
#endif

#ifdef _GPU
   subroutine gpu_free_neighbors(gpu_neigh, gpu_stream)
      implicit none
      type(gpu_neigh_storage_type) :: gpu_neigh
      type(c_ptr) :: gpu_stream

      call gpu_free_async(gpu_neigh%n_neigh_d, gpu_stream)
      call gpu_free_async(gpu_neigh%species_d, gpu_stream)
      call gpu_free_async(gpu_neigh%neighbor_species_d, gpu_stream)
      call gpu_free_async(gpu_neigh%neighbors_list_d, gpu_stream)
      call gpu_free_async(gpu_neigh%rjs_d, gpu_stream)
      call gpu_free_async(gpu_neigh%xyz_d, gpu_stream)
      call gpu_stream_sync(gpu_stream)
   end subroutine gpu_free_neighbors
#endif

#ifdef _GPU
   subroutine collect_batched_pair_distribution(n_batches, gpu_host, n_dim_partial, n_samples, &
                                                pair_distribution_partial, n_species, n_atoms_of_species, v_uc)
      implicit none
      integer, intent(in) :: n_dim_partial
      integer, intent(in) :: n_samples
      integer, intent(in) :: n_batches
      integer, intent(in) :: n_species
      type(gpu_host_batch_storage_type), intent(in), allocatable :: gpu_host(:)
      real(dp), allocatable, intent(out) :: pair_distribution_partial(:, :)
      real(dp), allocatable, intent(in) :: n_atoms_of_species(:)
      real(dp), intent(in) :: v_uc
      real(dp), allocatable :: pair_distribution_partial_temp(:, :)
      real(dp), allocatable :: factors(:)
      real(dp) :: f
      integer :: i
      integer :: j
      integer :: k
      integer :: n_dim_idx
      integer :: ierr

      allocate (factors(1:n_dim_partial))

      n_dim_idx = 1
      outer: do j = 1, n_species
         do k = 1, n_species

            if (j > k) cycle

            if (j == k) f = v_uc/n_atoms_of_species(j)/n_atoms_of_species(k)
            if (j /= k) f = v_uc/n_atoms_of_species(j)/n_atoms_of_species(k)/2.d0

            factors(n_dim_idx) = f

            n_dim_idx = n_dim_idx + 1
            if (n_dim_idx > n_dim_partial) exit outer
         end do
      end do outer

      allocate (pair_distribution_partial(1:n_samples, 1:n_dim_partial))
      pair_distribution_partial = 0.d0

      do j = 1, n_dim_partial
         f = factors(j)
         do i = 1, n_batches

            pair_distribution_partial(1:n_samples, j) = pair_distribution_partial(1:n_samples, j) &
                                                        + gpu_host(i)%host(j)%pair_distribution_partial_h(1:n_samples)*f

         end do
      end do

      deallocate (factors)

#ifdef _MPIF90
      allocate (pair_distribution_partial_temp(1:n_samples, 1:n_dim_partial))
      pair_distribution_partial_temp = 0.d0

      call mpi_reduce(pair_distribution_partial,&
           & pair_distribution_partial_temp, n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
           & 0, MPI_COMM_WORLD, ierr)

      pair_distribution_partial = pair_distribution_partial_temp
      deallocate (pair_distribution_partial_temp)

      call mpi_bcast(pair_distribution_partial, n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
           & MPI_COMM_WORLD, ierr)

#endif

   end subroutine collect_batched_pair_distribution
#endif

#ifdef _GPU
   subroutine collect_batched_forces(n_batches, gpu_host, n_dim_partial, &
                                     forces_out, virial_out, n_sites)
      implicit none
      integer, intent(in) :: n_batches
      integer, intent(in) :: n_sites
      integer, intent(in) :: n_dim_partial
      type(gpu_host_batch_storage_type), intent(in), allocatable :: gpu_host(:)
      real(dp), allocatable :: forces_out(:, :)
      real(dp) :: virial_out(1:3, 1:3)
      integer :: i
      integer :: j
      integer :: k
      integer :: n_dim_idx

      do i = 1, n_batches
         do j = 1, n_dim_partial

            forces_out = forces_out + gpu_host(i)%host(j)%forces_h
            virial_out = virial_out + gpu_host(i)%host(j)%virial_h

         end do
      end do

   end subroutine collect_batched_forces
#endif

#ifdef _GPU
   subroutine free_host_batches(gpu_host, n_batches, n_dim_partial)
      implicit none
      type(gpu_host_batch_storage_type), allocatable :: gpu_host(:)
      integer, intent(in) :: n_batches
      integer, intent(in) :: n_dim_partial
      integer :: i
      integer :: j

      do i = 1, n_batches
         do j = 1, n_dim_partial

            if (allocated(gpu_host(i)%host(j)%xyz_k_h)) &
               deallocate (gpu_host(i)%host(j)%xyz_k_h)

            if (allocated(gpu_host(i)%host(j)%pair_distribution_partial_h)) &
               deallocate (gpu_host(i)%host(j)%pair_distribution_partial_h)

            if (allocated(gpu_host(i)%host(j)%pair_distribution_partial_der_h)) &
               deallocate (gpu_host(i)%host(j)%pair_distribution_partial_der_h)

            if (allocated(gpu_host(i)%host(j)%forces_h)) &
               deallocate (gpu_host(i)%host(j)%forces_h)

            if (allocated(gpu_host(i)%host(j)%rjs_index_h)) &
               deallocate (gpu_host(i)%host(j)%rjs_index_h)

            if (allocated(gpu_host(i)%host(j)%k_index_h)) &
               deallocate (gpu_host(i)%host(j)%k_index_h)

         end do
         deallocate (gpu_host(i)%host)
      end do
      deallocate (gpu_host)

   end subroutine free_host_batches
#endif

#ifdef _GPU
   subroutine free_exp_batches(gpu_exp, n_batches)
      implicit none
      type(gpu_storage_type), allocatable :: gpu_exp(:)
      integer, intent(in) :: n_batches
      integer :: i
      integer :: j

      do i = 1, n_batches

         if (allocated(gpu_exp(i)%nk)) deallocate (gpu_exp(i)%nk)

         if (allocated(gpu_exp(i)%nk_d)) deallocate (gpu_exp(i)%nk_d)

         if (allocated(gpu_exp(i)%k_index_d)) deallocate (gpu_exp(i)%k_index_d)

         if (allocated(gpu_exp(i)%j2_index_d)) deallocate (gpu_exp(i)%j2_index_d)

         if (allocated(gpu_exp(i)%xyz_k_d)) deallocate (gpu_exp(i)%xyz_k_d)

         if (allocated(gpu_exp(i)%pair_distribution_partial_d)) deallocate (gpu_exp(i)%pair_distribution_partial_d)

         if (allocated(gpu_exp(i)%pair_distribution_partial_der_d)) deallocate (gpu_exp(i)%pair_distribution_partial_der_d)

         if (allocated(gpu_exp(i)%st_nk_d)) deallocate (gpu_exp(i)%st_nk_d)

         if (allocated(gpu_exp(i)%st_k_index_d)) deallocate (gpu_exp(i)%st_k_index_d)

         if (allocated(gpu_exp(i)%st_j2_index_d)) deallocate (gpu_exp(i)%st_j2_index_d)

         if (allocated(gpu_exp(i)%rjs_index_d)) deallocate (gpu_exp(i)%rjs_index_d)

         if (allocated(gpu_exp(i)%st_pair_distribution_partial_d)) deallocate (gpu_exp(i)%st_pair_distribution_partial_d)

         if (allocated(gpu_exp(i)%st_pair_distribution_partial_der_d)) deallocate (gpu_exp(i)%st_pair_distribution_partial_der_d)

         if (allocated(gpu_exp(i)%nk_flags_d)) deallocate (gpu_exp(i)%nk_flags_d)
         if (allocated(gpu_exp(i)%nk_flags_sum_d)) deallocate (gpu_exp(i)%nk_flags_sum_d)

      end do
      deallocate (gpu_exp)

   end subroutine free_exp_batches
#endif

#ifdef _GPU
   subroutine total_gpu_memory(add)
      real(dp), intent(in) :: add
      real(dp), save :: total
      logical, save :: first_call = .true.

      if (first_call) then
         total = 0.0d0
         first_call = .false.
      end if

      total = total + add

      if (debug_gpu_batches) print *, " GPU mem = ", total/1024.d0/1024.d0/1024.d0, " Gb"
      call flush (101)

   end subroutine
#endif

#ifdef _GPU
   subroutine calculate_batched_pair_distribution( &
      gpu_exp, &
      gpu_host, &
      gpu_neigh, &
      x, &
      dV, &
      n_atoms_of_species, &
      n_species, &
      n_sites, &
      i_beg, &
      i_end, &
      j_beg, &
      j_end, &
      n_samples, &
      r_min, &
      r_max, &
      r_cut, &
      kde_sigma, &
      gpu_stream, &
      x_d, &
      dV_d, &
      v_uc, &
      rank)
      implicit none
      integer, intent(in) :: n_species
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: n_sites
      integer, intent(in) :: n_samples
      integer, intent(in) :: rank
      real(dp), allocatable, intent(in) :: n_atoms_of_species(:)
      real(dp), allocatable, intent(in) :: x(:)
      real(dp), allocatable, intent(in) :: dV(:)
      real(dp), intent(in) :: r_min
      real(dp), intent(in) :: r_max
      real(dp), intent(in) :: r_cut
      real(dp), intent(in) :: kde_sigma
      real(dp), intent(in) :: v_uc
      type(gpu_storage_type), intent(inout) :: gpu_exp
      type(gpu_host_batch_storage_type), intent(inout), target :: gpu_host
      type(gpu_neigh_storage_type), intent(in) :: gpu_neigh

      real(dp), parameter :: pi = acos(-1.0)
      integer :: n_dim_partial
      integer :: n_dim_idx
      type(c_ptr), intent(in) :: x_d
      type(c_ptr), intent(in) :: dV_d
      integer(c_size_t) :: st_x_d
      integer, target :: nk_temp(1)
      type(c_ptr) :: nk_flags_d
      type(c_ptr) :: nk_flags_sum_d
      integer(c_size_t) :: st_nk_flags
      integer(c_size_t) :: st_nk_temp
      type(c_ptr) :: rjs_index_d
      type(c_ptr) :: pdf_to_reduce_d
      integer(c_size_t) :: st_rjs_index_d
      integer(c_size_t) :: st_k_index_d
      integer(c_size_t) :: st_pdf_to_reduce_d

      real(dp) :: pdf_factor
      real(dp) :: der_factor = 0.d0
      real(dp) :: f
      integer :: i
      integer :: j
      integer :: k
      integer :: l

      type(c_ptr) :: gpu_stream

      real(dp), allocatable, target :: x_check(:)

      n_dim_partial = n_species*(n_species + 1)/2

      allocate (gpu_host%host(1:n_dim_partial))

      allocate (gpu_exp%nk(1:n_dim_partial))
      allocate (gpu_exp%nk_d(1:n_dim_partial))
      allocate (gpu_exp%k_index_d(1:n_dim_partial))
      allocate (gpu_exp%j2_index_d(1:n_dim_partial))
      allocate (gpu_exp%xyz_k_d(1:n_dim_partial))
      allocate (gpu_exp%pair_distribution_partial_d(1:n_dim_partial))
      allocate (gpu_exp%nk_flags_sum_d(1:n_dim_partial))
      allocate (gpu_exp%nk_flags_d(1:n_dim_partial))
      allocate (gpu_exp%rjs_index_d(1:n_dim_partial))

      allocate (gpu_exp%st_nk_d(1:n_dim_partial))
      allocate (gpu_exp%st_k_index_d(1:n_dim_partial))
      allocate (gpu_exp%st_j2_index_d(1:n_dim_partial))
      allocate (gpu_exp%st_pair_distribution_partial_d(1:n_dim_partial))

      n_dim_idx = 1
      outer1: do j = 1, n_species
         do k = 1, n_species

            if (j > k) cycle ! We have already calculated the pair correlation function!

            st_nk_temp = int(1, c_size_t)*c_int
            call gpu_malloc_async(gpu_exp%nk_d(n_dim_idx), st_nk_temp, gpu_stream)
            st_nk_flags = int((j_end - j_beg + 1), c_size_t)*c_int
            call gpu_malloc_async(gpu_exp%nk_flags_d(n_dim_idx), st_nk_flags, gpu_stream)
            call gpu_memset_async(gpu_exp%nk_flags_d(n_dim_idx), 0, st_nk_flags, gpu_stream)
            call gpu_malloc_async(gpu_exp%nk_flags_sum_d(n_dim_idx), st_nk_flags, gpu_stream)

            call total_gpu_memory(dfloat(int((j_end - j_beg + 1), c_size_t)*2*4))

            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "int(gpu_exp % nk_d(n_dim_idx))"
            if (debug_gpu_batches) call gpu_print_pointer_int(gpu_exp%nk_d(n_dim_idx))
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "int(gpu_exp % nk_flags_d(n_dim_idx))"
            if (debug_gpu_batches) call gpu_print_pointer_int(gpu_exp%nk_flags_d(n_dim_idx))
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "int(gpu_exp % nk_flags_sum_d(n_dim_idx))"
            if (debug_gpu_batches) call gpu_print_pointer_int(gpu_exp%nk_flags_sum_d(n_dim_idx))
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "int(gpu_neigh % neighbors_list_d  )"
            if (debug_gpu_batches) call gpu_print_pointer_int(gpu_neigh%neighbors_list_d)
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "int(gpu_neigh % n_neigh_d         )"
            if (debug_gpu_batches) call gpu_print_pointer_int(gpu_neigh%n_neigh_d)
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "int(gpu_neigh % neighbor_species_d)"
            if (debug_gpu_batches) call gpu_print_pointer_int(gpu_neigh%neighbor_species_d)
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "int(gpu_neigh % species_d         )"
            if (debug_gpu_batches) call gpu_print_pointer_int(gpu_neigh%species_d)
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "double(gpu_neigh % rjs_d             )"
            if (debug_gpu_batches) call gpu_print_pointer_double(gpu_neigh%rjs_d)
            if (debug_gpu_batches) print *, "-- Rank ", rank, " ", "double(gpu_neigh % xyz_d             )"
            if (debug_gpu_batches) call gpu_print_pointer_double(gpu_neigh%xyz_d)
            call flush (101)

            call gpu_get_pair_distribution_nk( &
               1, &
               i_end - i_beg + 1, &
               j_end - j_beg + 1, &
               n_sites, &
               gpu_neigh%neighbors_list_d, &
               gpu_neigh%n_neigh_d, &
               gpu_neigh%neighbor_species_d, &
               gpu_neigh%species_d, &
               gpu_neigh%rjs_d, &
               gpu_neigh%xyz_d, &
               r_min, &
               r_max, &
               r_cut, &
               6.d0*kde_sigma, &
               gpu_exp%nk_d(n_dim_idx), &
               gpu_exp%nk_flags_d(n_dim_idx), &
               gpu_exp%nk_flags_sum_d(n_dim_idx), &
               j, &
               k, &
               gpu_stream)

            call gpu_stream_sync(gpu_stream)

            call gpu_free_async(gpu_exp%nk_flags_d(n_dim_idx), gpu_stream)

            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            call total_gpu_memory(dfloat(-int((j_end - j_beg + 1), c_size_t)*4))

            ! Now copy the value of nk from the gpu
!          print *, "out of pdf nk kernel "
            st_nk_temp = int(1, c_size_t)*c_int
            call cpy_dtoh(gpu_exp%nk_d(n_dim_idx), c_loc(nk_temp), st_nk_temp, gpu_stream)
            call gpu_stream_sync(gpu_stream)

            gpu_exp%nk(n_dim_idx) = nk_temp(1)
            call gpu_free_async(gpu_exp%nk_d(n_dim_idx), gpu_stream)

            ! Now we create temporary arrays for the k indices

            st_rjs_index_d = int(gpu_exp%nk(n_dim_idx), c_size_t)*c_double
            if (debug_gpu_batches) print *, " allocating rjs "
            call total_gpu_memory(dfloat(int(gpu_exp%nk(n_dim_idx), c_size_t)*8))
            call gpu_malloc_async(gpu_exp%rjs_index_d(n_dim_idx), st_rjs_index_d, gpu_stream)
            call gpu_memset_async(gpu_exp%rjs_index_d(n_dim_idx), 0, st_rjs_index_d, gpu_stream)

            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            gpu_exp%st_k_index_d(n_dim_idx) = int(gpu_exp%nk(n_dim_idx), c_size_t)*c_int
            if (debug_gpu_batches) print *, " allocating k index "
            call total_gpu_memory(dfloat(int(gpu_exp%nk(n_dim_idx), c_size_t)*4))
            call gpu_malloc_async(gpu_exp%k_index_d(n_dim_idx), gpu_exp%st_k_index_d(n_dim_idx), gpu_stream)
            call gpu_memset_async(gpu_exp%k_index_d(n_dim_idx), 0, gpu_exp%st_k_index_d(n_dim_idx), gpu_stream)

            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            if (debug_gpu_batches) print *, " allocating j2 index "
            call total_gpu_memory(dfloat(int(gpu_exp%nk(n_dim_idx), c_size_t)*4))
            call gpu_malloc_async(gpu_exp%j2_index_d(n_dim_idx), gpu_exp%st_k_index_d(n_dim_idx), gpu_stream)
            call gpu_memset_async(gpu_exp%j2_index_d(n_dim_idx), 0, gpu_exp%st_k_index_d(n_dim_idx), gpu_stream)

            if (debug_gpu_batches) print *, " allocating j2 index "
            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            call total_gpu_memory(dfloat(int(gpu_exp%nk(n_dim_idx), c_size_t)*8*3))
            call gpu_malloc_async(gpu_exp%xyz_k_d(n_dim_idx), 3*st_rjs_index_d, gpu_stream)
            call gpu_memset_async(gpu_exp%xyz_k_d(n_dim_idx), 0, 3*st_rjs_index_d, gpu_stream)

            call gpu_set_pair_distribution_k_index(1, i_end - i_beg + 1, j_end - j_beg + 1, n_sites, & ! i_beg, i_end, j_end, n_sites,&
                                                   gpu_neigh%neighbors_list_d, &
                                                   gpu_neigh%rjs_d, &
                                                   gpu_neigh%xyz_d, &
                                                   gpu_exp%k_index_d(n_dim_idx), &
                                                   gpu_exp%j2_index_d(n_dim_idx), &
                                                   gpu_exp%rjs_index_d(n_dim_idx), &
                                                   gpu_exp%xyz_k_d(n_dim_idx), &
                                                   gpu_exp%nk_flags_d(n_dim_idx), gpu_exp%nk_flags_sum_d(n_dim_idx), &
                                                   gpu_stream)

            call gpu_free_async(gpu_exp%nk_flags_sum_d(n_dim_idx), gpu_stream)
            if (debug_gpu_batches) print *, "deallocing flags sum"
            call total_gpu_memory(dfloat(-int((j_end - j_beg + 1), c_size_t)*4))
            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            allocate (gpu_host%host(n_dim_idx)%k_index_h(1:gpu_exp%nk(n_dim_idx)))

            if (debug_gpu_batches) print *, "deallocing k index"
            call total_gpu_memory(dfloat(-int(gpu_exp%nk(n_dim_idx), c_size_t)*4))

            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            call cpy_dtoh( &
               gpu_exp%k_index_d(n_dim_idx), &
               c_loc(gpu_host%host(n_dim_idx)%k_index_h), &
               gpu_exp%st_k_index_d(n_dim_idx), &
               gpu_stream)
            call gpu_free_async(gpu_exp%k_index_d(n_dim_idx), gpu_stream)

            allocate (gpu_host%host(n_dim_idx)%j2_index_h(1:gpu_exp%nk(n_dim_idx)))

            if (debug_gpu_batches) print *, "deallocing j2 index"
            call total_gpu_memory(dfloat(-int(gpu_exp%nk(n_dim_idx), c_size_t)*4))
            call cpy_dtoh( &
               gpu_exp%j2_index_d(n_dim_idx), &
               c_loc(gpu_host%host(n_dim_idx)%j2_index_h), &
               gpu_exp%st_k_index_d(n_dim_idx), &
               gpu_stream)
            call gpu_free_async(gpu_exp%j2_index_d(n_dim_idx), gpu_stream)
            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            allocate (gpu_host%host(n_dim_idx)%rjs_index_h(1:gpu_exp%nk(n_dim_idx)))

            call cpy_dtoh( &
               gpu_exp%rjs_index_d(n_dim_idx), &
               c_loc(gpu_host%host(n_dim_idx)%rjs_index_h), &
               st_rjs_index_d, &
               gpu_stream)

            allocate (gpu_host%host(n_dim_idx)%xyz_k_h(1:3, 1:gpu_exp%nk(n_dim_idx)))
            if (debug_gpu_batches) print *, "deallocing xyz_k"
            call total_gpu_memory(dfloat(-int(gpu_exp%nk(n_dim_idx), c_size_t)*3*8))
            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            call cpy_dtoh( &
               gpu_exp%xyz_k_d(n_dim_idx), &
               c_loc(gpu_host%host(n_dim_idx)%xyz_k_h), &
               3*st_rjs_index_d, &
               gpu_stream)
            call gpu_free_async(gpu_exp%xyz_k_d(n_dim_idx), gpu_stream)

            if (debug_gpu_batches) print *, "allocing pdf"
            call total_gpu_memory(dfloat(int(n_samples, c_size_t)*8))
            call gpu_stream_sync(gpu_stream)
            if (debug_gpu_batches) call gpu_check_error()
            gpu_exp%st_pair_distribution_partial_d(n_dim_idx) = int(n_samples, c_size_t)*c_double
            call gpu_malloc_async(gpu_exp%pair_distribution_partial_d(n_dim_idx), &
                                  gpu_exp%st_pair_distribution_partial_d(n_dim_idx), gpu_stream)
            call gpu_memset_async(gpu_exp%pair_distribution_partial_d(n_dim_idx), 0, &
                                  gpu_exp%st_pair_distribution_partial_d(n_dim_idx), gpu_stream)

            if (debug_gpu_batches) print *, "allocing pdf to reduce "
            call total_gpu_memory(dfloat(int(gpu_exp%nk(n_dim_idx), c_size_t)*n_samples*8))

            st_pdf_to_reduce_d = int(gpu_exp%nk(n_dim_idx), c_size_t)*n_samples*c_double
            call gpu_malloc_async(pdf_to_reduce_d, st_pdf_to_reduce_d, gpu_stream)
            call gpu_memset_async(pdf_to_reduce_d, 0, st_pdf_to_reduce_d, gpu_stream)

            pdf_factor = ((r_max - r_min)/dfloat(n_samples))/(sqrt(2.d0*pi)*kde_sigma)

            !        if ( j == k ) f = 1.d0
!           if ( j /= k ) f = 2.d0

!            pdf_factor = pdf_factor * der_factor
! ! !
!          call gpu_stream_sync(gpu_stream)
            !          print *, " >> Getting pdf batch"

            if (debug_gpu_batches) print *, "gpu_exp%pair_distribution_partial_d(n_dim_idx)"
            if (debug_gpu_batches) call gpu_print_pointer_double(gpu_exp%pair_distribution_partial_d(n_dim_idx))
            if (debug_gpu_batches) print *, "pdf_to_reduce_d "
            if (debug_gpu_batches) call gpu_print_pointer_double(pdf_to_reduce_d)
            if (debug_gpu_batches) print *, "x_d "
            if (debug_gpu_batches) call gpu_print_pointer_double(x_d)
            if (debug_gpu_batches) print *, "dV_d "
            if (debug_gpu_batches) call gpu_print_pointer_double(dV_d)
            if (debug_gpu_batches) print *, "gpu_exp%rjs_index_d(n_dim_idx) "
            if (debug_gpu_batches) call gpu_print_pointer_double(gpu_exp%rjs_index_d(n_dim_idx))
            call flush (101)

            der_factor = 0.0d0
            call gpu_get_pair_distribution_only_falloc( &
               gpu_exp%pair_distribution_partial_d(n_dim_idx), &
               pdf_to_reduce_d, &
               gpu_exp%nk(n_dim_idx), &
               n_samples, &
               kde_sigma, &
               x_d, &
               dV_d, &
               gpu_exp%rjs_index_d(n_dim_idx), &
               pdf_factor, &
               der_factor, &
               gpu_stream)

            !--- check x ---!
            ! allocate( x_check(1:n_samples) )
            ! st_x_d = n_samples * c_double
            !          ! call cpy_dtoh_event(&
!          call cpy_dtoh(&
            !      x_d, &
            !      c_loc(x_check), &
            !      st_x_d, &
            !      gpu_stream)

            ! call cpy_dtoh_event(&
            !      dV_d, &
            !      c_loc(x_check), &
            !      st_x_d, &
            !      gpu_stream)

            call gpu_free_async(pdf_to_reduce_d, gpu_stream)
            call gpu_free_async(gpu_exp%rjs_index_d(n_dim_idx), gpu_stream)

            allocate (gpu_host%host(n_dim_idx)%pair_distribution_partial_h(1:n_samples))
            call cpy_dtoh( &
               gpu_exp%pair_distribution_partial_d(n_dim_idx), &
               c_loc(gpu_host%host(n_dim_idx)%pair_distribution_partial_h), &
               gpu_exp%st_pair_distribution_partial_d(n_dim_idx), &
               gpu_stream)

            call gpu_free_async(gpu_exp%pair_distribution_partial_d(n_dim_idx), gpu_stream)

!          call gpu_stream_sync( gpu_stream )
            ! call gpu_copy_pdf( n_samples, gpu_exp % pair_distribution_partial_d(n_dim_idx), &
            !      gpu_host % pair_distribution_partial_h(n_dim)pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
            !      gpu_stream )

            n_dim_idx = n_dim_idx + 1

            if (n_dim_idx > n_dim_partial) then
               exit outer1
            end if

         end do
      end do outer1

   end subroutine calculate_batched_pair_distribution
#endif

#ifdef _GPU
   subroutine calculate_batched_pair_distribution_der(gpu_exp, gpu_host,&
        &  x, dV, n_atoms_of_species, n_species, n_sites,&
        & i_beg, i_end, j_beg, j_end, n_samples, r_min, r_max, r_cut, kde_sigma, &
        & gpu_stream, x_d, dV_d, j, k, n_dim_idx, v_uc)
      implicit none
      integer, intent(in) :: n_species
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: n_sites
      integer, intent(in) :: n_samples
      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: n_dim_idx
      real(dp), allocatable, intent(in) :: n_atoms_of_species(:)
      real(dp), allocatable, intent(in) :: x(:)
      real(dp), allocatable, intent(in) :: dV(:)
      real(dp), intent(in) :: r_min
      real(dp), intent(in) :: r_max
      real(dp), intent(in) :: r_cut
      real(dp), intent(in) :: kde_sigma
      real(dp), intent(in) :: v_uc
      type(gpu_storage_type), intent(inout) :: gpu_exp
      type(gpu_host_batch_storage_type), intent(inout), target :: gpu_host

      real(dp), parameter :: pi = acos(-1.0)
      integer :: n_dim_partial
      type(c_ptr), intent(in) :: x_d
      type(c_ptr), intent(in) :: dV_d
      integer(c_size_t) :: st_x_d
      integer, target :: nk_temp(1)
      type(c_ptr) :: nk_flags_d
      type(c_ptr) :: nk_sum_flags_d
      integer(c_size_t) :: st_nk_flags
      integer(c_size_t) :: st_nk_temp
      type(c_ptr) :: rjs_index_d
      integer(c_size_t) :: st_rjs_index_d

      real(dp) :: pdf_factor
      real(dp) :: der_factor = 0.d0
      real(dp) :: f
      integer :: i

      type(c_ptr) :: gpu_stream

      n_dim_partial = n_species*(n_species + 1)/2

    if( .not. allocated( gpu_exp % st_pair_distribution_partial_der_d ) ) allocate( gpu_exp % st_pair_distribution_partial_der_d(1:n_dim_partial) )
   if (.not. allocated(gpu_exp%pair_distribution_partial_der_d)) allocate (gpu_exp%pair_distribution_partial_der_d(1:n_dim_partial))

      st_rjs_index_d = int(gpu_exp%nk(n_dim_idx), c_size_t)*c_double
      call gpu_malloc_async(rjs_index_d, st_rjs_index_d, gpu_stream)
      call gpu_memset_async(rjs_index_d, 0, st_rjs_index_d, gpu_stream)

      call cpy_htod( &
         c_loc(gpu_host%host(n_dim_idx)%rjs_index_h), &
         rjs_index_d, &
         st_rjs_index_d, &
         gpu_stream)

      gpu_exp%st_pair_distribution_partial_der_d(n_dim_idx) = int(n_samples, c_size_t)*gpu_exp%nk(n_dim_idx)*c_double
    call gpu_malloc_async(gpu_exp %pair_distribution_partial_der_d(n_dim_idx), gpu_exp %st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)
    call gpu_memset_async(gpu_exp %pair_distribution_partial_der_d(n_dim_idx), 0, gpu_exp %st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

      pdf_factor = ((r_max - r_min)/ &
                    dfloat(n_samples))/(sqrt(2.d0*pi)*kde_sigma)

      if (j == k) f = 1.d0
      if (j /= k) f = 2.d0

      der_factor = v_uc/n_atoms_of_species(j)/ &
           & n_atoms_of_species(k)/f

      call gpu_get_pair_distribution_der_only( &
         gpu_exp%pair_distribution_partial_der_d(n_dim_idx), &
         gpu_exp%nk(n_dim_idx), &
         n_samples, &
         kde_sigma, &
         x_d, dV_d, &
         rjs_index_d, pdf_factor, der_factor, gpu_stream)

      call gpu_free_async(rjs_index_d, gpu_stream)

   end subroutine calculate_batched_pair_distribution_der
#endif

#ifdef _GPU
   subroutine setup_gpu_xrd_forces(gpu_exp, gpu_host, n_dim_idx, gpu_stream)
      implicit none
      integer, intent(in) :: n_dim_idx
      type(gpu_storage_type), intent(inout) :: gpu_exp
      type(gpu_host_batch_storage_type), intent(inout), target :: gpu_host
      integer(c_size_t) :: st_rjs_index_d
      type(c_ptr) :: gpu_stream
      ! copy the xyz, j2 and k_index_d arrays

      if (debug_gpu_batches) print *, "> nk = ", gpu_exp%nk(n_dim_idx)
      st_rjs_index_d = int(gpu_exp%nk(n_dim_idx), c_size_t)*c_double
      call gpu_malloc_async(gpu_exp%k_index_d(n_dim_idx), gpu_exp%st_k_index_d(n_dim_idx), gpu_stream)
      call gpu_malloc_async(gpu_exp%j2_index_d(n_dim_idx), gpu_exp%st_k_index_d(n_dim_idx), gpu_stream)
      call gpu_malloc_async(gpu_exp%xyz_k_d(n_dim_idx), 3*st_rjs_index_d, gpu_stream)

      call cpy_htod( &
         c_loc(gpu_host%host(n_dim_idx)%k_index_h), &
         gpu_exp%k_index_d(n_dim_idx), &
         gpu_exp%st_k_index_d(n_dim_idx), &
         gpu_stream)

      call cpy_htod( &
         c_loc(gpu_host%host(n_dim_idx)%j2_index_h), &
         gpu_exp%j2_index_d(n_dim_idx), &
         gpu_exp%st_k_index_d(n_dim_idx), &
         gpu_stream)

!    st_rjs_index_d = gpu_exp % nk( n_dim_idx ) * c_double

      call cpy_htod( &
         c_loc(gpu_host%host(n_dim_idx)%xyz_k_h), &
         gpu_exp%xyz_k_d(n_dim_idx), &
         3*st_rjs_index_d, &
         gpu_stream)

   end subroutine setup_gpu_xrd_forces
#endif

#ifdef _GPU
   subroutine free_gpu_xrd_forces(gpu_exp, gpu_host, n_dim_idx, gpu_stream)
      implicit none
      integer, intent(in) :: n_dim_idx
      type(gpu_storage_type), intent(inout) :: gpu_exp
      type(gpu_host_batch_storage_type), intent(inout), target :: gpu_host
      integer(c_size_t) :: st_rjs_index_d
      type(c_ptr) :: gpu_stream
      ! copy the xyz, j2 and k_index_d arrays

      call gpu_free_async(gpu_exp%k_index_d(n_dim_idx), gpu_stream)
      call gpu_free_async(gpu_exp%j2_index_d(n_dim_idx), gpu_stream)
      call gpu_free_async(gpu_exp%xyz_k_d(n_dim_idx), gpu_stream)
      call gpu_free_async(gpu_exp%pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

   end subroutine free_gpu_xrd_forces

   !--- GPU PAIR DISTRIBUTION FUNCTIONS ---!
#endif
#ifdef _GPU
   subroutine gpu_calculate_pair_distribution(n_dim_partial_out, params, x_pair_distribution&
        &, y_pair_distribution, y_pair_distribution_temp,&
        & pair_distribution_partial, pair_distribution_partial_temp, &
        & n_species, species_types, n_atoms_of_species, n_sites, a_box, b_box, c_box,&
        & indices, md_istep, mc_istep, i_beg, i_end, j_beg, j_end, ierr, rjs, xyz, &
        & neighbors_list, n_neigh, neighbor_species, species, rank,&
        & do_derivatives, pair_distribution_der, pair_distribution_partial_der, &
    & nk, pair_distribution_d, nk_d, k_index_d, j2_index_d, xyz_k_d, pair_distribution_partial_d, pair_distribution_partial_der_d, &
        & st_nk_d, st_k_index_d, st_j2_index_d, st_pair_distribution_partial_d, st_pair_distribution_partial_der_d,&
        & n_neigh_d, species_d, neighbor_species_d, neighbor_list_d, rjs_d, xyz_d, species_types_d, cublas_handle, gpu_stream, &
        gpu_host_storage, gpu_low_memory,&
        & pair_distribution_partial_temp_der, energies_pair_distribution, forces_pair_distribution, virial)
      implicit none
      type(input_parameters), intent(inout) :: params
      integer, intent(out) :: n_dim_partial_out
      real(dp), allocatable, intent(out), target :: x_pair_distribution(:)
      real(dp), allocatable, intent(out), target :: y_pair_distribution(:)
      real(dp), allocatable, intent(out), target :: pair_distribution_partial(:, :)
      real(dp), allocatable, intent(out), target :: n_atoms_of_species(:)
      real(dp), allocatable, intent(out), target :: pair_distribution_partial_temp(:, :)
      real(dp), allocatable, intent(out), target :: y_pair_distribution_temp(:)
      real(dp), allocatable, intent(out), target :: pair_distribution_der(:, :)
      real(dp), allocatable, intent(out), target :: pair_distribution_partial_der(:, :, :)
      real(dp), allocatable, intent(out), target :: pair_distribution_partial_temp_der(:, :, :)
      real(dp), allocatable, intent(out), target :: energies_pair_distribution(:)
      real(dp), allocatable, intent(out), target :: forces_pair_distribution(:, :)
      character*8, allocatable, intent(in) :: species_types(:)
      real(dp), intent(in), allocatable, target :: rjs(:)
      real(dp), intent(in), allocatable, target :: xyz(:, :)
      integer, intent(in), allocatable, target :: neighbors_list(:)
      integer, intent(in), allocatable, target :: n_neigh(:)
      integer, intent(in), allocatable, target :: neighbor_species(:)
      integer, intent(in), allocatable, target :: species(:)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      real(dp), intent(inout) :: virial(1:3, 1:3)
      real(dp) :: v_uc
      real(dp) :: f
!     -V dE/dV, the cell half of this observable's virial.
      real(dp) :: dedv
      real(dp) :: pdf_factor
      real(dp) :: der_factor
      real(dp) :: total_memory_usage = 0.d0
      integer, intent(in) :: n_species
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_beg
      integer, intent(in) :: i_end
      integer, intent(in) :: j_beg
      integer, intent(in) :: j_end
      integer, intent(in) :: indices(1:3)
      integer, intent(in) :: md_istep
      integer, intent(in) :: mc_istep
      integer, intent(in) :: rank
      integer, intent(inout) :: ierr
      real(dp), allocatable, target :: factors(:)
      real(dp), allocatable, target :: pair_distribution_der_temp(:)
      real(dp), allocatable, target :: dV(:)
      real(dp), allocatable, target :: pdf_gpu_check(:)
      real(dp), allocatable, target :: rjs_temp(:)
      real(dp), allocatable, target :: ders_temp(:, :)
      integer, allocatable, target :: ks_temp(:)
      integer, allocatable, target :: ksd_temp(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: i2
      integer :: n_dim_partial
      integer :: n_dim_idx
      logical, intent(in) :: do_derivatives
      real(dp), parameter :: pi = acos(-1.0)
      logical :: write_condition
      logical :: overwrite_condition
      character*1024 :: filename
      type(c_ptr) :: cublas_handle
      type(c_ptr) :: gpu_stream

      type(c_ptr) :: n_neigh_d
      type(c_ptr) :: species_d
      type(c_ptr) :: neighbor_species_d
      type(c_ptr) :: neighbor_list_d
      type(c_ptr) :: rjs_d
      type(c_ptr) :: xyz_d
      type(c_ptr) :: species_types_d
      type(c_ptr) :: x_d
      type(c_ptr) :: dV_d
      type(c_ptr) :: pair_distribution_d
      type(c_ptr), allocatable :: nk_d(:)
      type(c_ptr), allocatable :: nk_flags_d(:)
      type(c_ptr), allocatable :: nk_flags_sum_d(:)
      type(c_ptr), allocatable :: k_index_d(:)
      type(c_ptr), allocatable :: j2_index_d(:)
      type(c_ptr), allocatable :: rjs_index_d(:)
      type(c_ptr), allocatable :: xyz_k_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_d(:)
      type(c_ptr), allocatable :: pair_distribution_partial_der_d(:)
      integer(c_size_t), allocatable :: st_nk_d(:)
      integer(c_size_t), allocatable :: st_k_index_d(:)
      integer(c_size_t), allocatable :: st_j2_index_d(:)
      integer(c_size_t), allocatable :: st_rjs(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_d(:)
      integer(c_size_t), allocatable :: st_pair_distribution_partial_der_d(:)
      integer(c_size_t) :: st_nk_flags
      integer(c_size_t) :: st_nk_temp
      integer(c_size_t) :: st_x_d
      integer, allocatable :: nk(:)
      integer, allocatable :: k_index_single(:)
      integer, target :: nk_temp(1)

      type(gpu_host_storage_type), intent(inout), allocatable, target :: gpu_host_storage(:)
      logical, intent(in) :: gpu_low_memory
      integer(c_size_t) :: st_n_sites_int
      integer(c_size_t) :: st_n_atom_pairs_int
      integer(c_size_t) :: st_n_atom_pairs_double
      integer(c_size_t) :: st_species_types_d

!            st_species_types_d = n_species * c_int
!            call gpu_malloc_async(species_types_d,st_species_types_d,gpu_stream)
!            call cpy_htod(c_loc(species_types),species_types_d,st_species_types_d,gpu_stream)
! !           call gpu_device_sync()

      ! Seeing if my gpu helper mod actually helps

      ! Things that are allocated here:
      ! Always:
      !  > x_pair_distribution
      !  > y_pair_distribution
      ! if pair_distribution_partial == .true.
      !  > pair_distribution_partial( n_samples, n_spec * (n_spec + 1)/2 )
      !  if do_derivatives == .true.
      !    > pair_distribution_partial_der( n_samples, n_spec * (n_spec + 1)/2, j_beg : j_end )

      ! first allocate the necessary arrays for the
      ! calculation of the pair correlation function
      if (allocated(x_pair_distribution)) deallocate (x_pair_distribution)
      if (allocated(y_pair_distribution)) deallocate (y_pair_distribution)

      allocate (x_pair_distribution(1:params%pair_distribution_n_samples))
      allocate (y_pair_distribution(1:params%pair_distribution_n_samples))

      allocate (dV(1:params%pair_distribution_n_samples))

      if (params%n_exp > 0) then
         do i = 1, params%n_exp
            if (trim(params%exp_data(i)%label) == 'pair_distribution') then
               x_pair_distribution = params%exp_data(i)%x
            end if
         end do
      end if

      n_dim_partial_out = 0

      if (params%pair_distribution_partial) then
         n_dim_partial = n_species*(n_species + 1)/2
         n_dim_partial_out = n_dim_partial
         allocate (factors(1:n_dim_partial))

         n_dim_idx = 1
         outer: do i = 1, n_species
            do j = 1, n_species
               if (i > j) cycle

               if (i /= j) then
                  factors(n_dim_idx) = 2.d0
               else
                  factors(n_dim_idx) = 1.d0
               end if

               n_dim_idx = n_dim_idx + 1
               if (n_dim_idx > n_dim_partial) then
                  exit outer
               end if

            end do
         end do outer

         if (.not. allocated(pair_distribution_partial)) then   !deallocate(pair_distribution_partial)
            allocate (pair_distribution_partial(1:params%pair_distribution_n_samples,&
                 & 1:n_dim_partial))
         end if

         pair_distribution_partial = 0.d0

         if (params%do_forces .and. params%exp_forces) then
            allocate (pair_distribution_partial_der(1:params%pair_distribution_n_samples,&
              & 1:n_dim_partial, j_beg:j_end))
            pair_distribution_partial_der = 0.d0

            if (rank == 0 .and. md_istep == 0) write (*, '(A,1X,F7.4,1X,A)') "Gb/core: partial pdfder = ", dfloat(params&
                 &%pair_distribution_n_samples*n_dim_partial*j_end)&
                 & *8.d0/(dfloat(1024*1024*1024)), " Gb  |"
            if (rank == 0 .and. md_istep == 0) write (*, *) '                                       |'

         end if
      else
         if (params%do_forces .and. params%exp_forces) then
            allocate (pair_distribution_partial_der(1:params%pair_distribution_n_samples, 1:1, &
              &  j_beg:j_end))
            pair_distribution_partial_der = 0.d0
         end if

      end if

      if (allocated(n_atoms_of_species)) deallocate (n_atoms_of_species)
      allocate (n_atoms_of_species(1:n_species))

      do j = 1, n_species
         n_atoms_of_species(j) = 0.d0
         do i2 = 1, n_sites
            if (species(i2) == j) then
               n_atoms_of_species(j) = n_atoms_of_species(j) + 1.d0
            end if
         end do
      end do

#ifdef _MPIF90
      if (params%pair_distribution_partial) then
         allocate (pair_distribution_partial_temp(1:params%pair_distribution_n_samples, 1:n_dim_partial))

         pair_distribution_partial_temp = 0.0d0
      end if

      allocate (y_pair_distribution_temp(1:params%pair_distribution_n_samples))
      y_pair_distribution_temp = 0.d0

#endif
      v_uc = dot_product(cross_product(a_box,&
           & b_box), c_box)/(&
           & dfloat(indices(1)*indices(2)&
           &*indices(3)))

      !###---   Calculating the partial pair distribution functions   ---###!
!!    print *, " - Allocating pdf gpu pointers -"
      allocate (nk_d(1:n_dim_partial))
      allocate (nk_flags_d(1:n_dim_partial))
      allocate (nk_flags_sum_d(1:n_dim_partial))
      allocate (k_index_d(1:n_dim_partial))
      allocate (j2_index_d(1:n_dim_partial))
      allocate (rjs_index_d(1:n_dim_partial))
      allocate (xyz_k_d(1:n_dim_partial))
      allocate (pair_distribution_partial_d(1:n_dim_partial))
      allocate (pair_distribution_partial_der_d(1:n_dim_partial))
      allocate (st_nk_d(1:n_dim_partial))
      allocate (st_k_index_d(1:n_dim_partial))
      allocate (st_j2_index_d(1:n_dim_partial))
      allocate (st_rjs(1:n_dim_partial))
      allocate (st_pair_distribution_partial_d(1:n_dim_partial))
      allocate (st_pair_distribution_partial_der_d(1:n_dim_partial))
      allocate (nk(1:n_dim_partial))

      st_x_d = int(params%pair_distribution_n_samples, c_size_t)*c_double

      call setup_pdf_arrays(params%r_range_min, params%r_range_max,&
           & params%pair_distribution_rcut, params&
           &%pair_distribution_n_samples, x_pair_distribution, dV)

      call gpu_meminfo()
      call gpu_malloc_async(x_d, st_x_d, gpu_stream)
      call cpy_htod(c_loc(x_pair_distribution), x_d, st_x_d, gpu_stream)

      call gpu_malloc_async(dV_d, st_x_d, gpu_stream)
      call cpy_htod(c_loc(dV), dV_d, st_x_d, gpu_stream)

      if (params%pair_distribution_partial) then
         if (gpu_low_memory) then

            ! Allocate the host storage arrays

            allocate (gpu_host_storage(1:n_dim_partial))

            n_dim_idx = 1
            outera: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle ! We have already calculated the pair correlation function!

                  ! Note that with the calculation of the derivatives here,
                  ! this is without the -2 * delta_ik (r_j^alpha - r_i^alpha)
                  ! factor, which allows for some freeing of memory

                  ! Get all nk_flags
                  call gpu_meminfo()

                  call gpu_stream_sync(gpu_stream)

                  st_nk_temp = int(1, c_size_t)*c_int
                  call gpu_malloc_async(nk_d(n_dim_idx), st_nk_temp, gpu_stream)
                  st_nk_flags = int(j_end, c_size_t)*c_int
                  call gpu_malloc_async(nk_flags_d(n_dim_idx), st_nk_flags, gpu_stream)
                  call gpu_memset_async(nk_flags_d(n_dim_idx), 0, st_nk_flags, gpu_stream)
                  call gpu_malloc_async(nk_flags_sum_d(n_dim_idx), st_nk_flags, gpu_stream)

                  ! Here I need the neighbor_list_d, n_neigh_d, neighbor_species_d, species_d, rjs_d

                  st_n_sites_int = n_sites*sizeof(n_neigh(1))
                  call gpu_malloc_async(n_neigh_d, st_n_sites_int, gpu_stream)
                  call cpy_htod(c_loc(n_neigh), n_neigh_d, st_n_sites_int, gpu_stream)
                  call gpu_malloc_async(species_d, st_n_sites_int, gpu_stream)
                  call cpy_htod(c_loc(species), species_d, st_n_sites_int, gpu_stream)
                  st_n_atom_pairs_int = j_end*sizeof(neighbor_species(1))

                  call gpu_malloc_async(neighbor_species_d, st_n_atom_pairs_int, gpu_stream)
                  call cpy_htod(c_loc(neighbor_species), neighbor_species_d, st_n_atom_pairs_int, gpu_stream)
                  call gpu_malloc_async(neighbor_list_d, st_n_atom_pairs_int, gpu_stream)
                  call cpy_htod(c_loc(neighbors_list), neighbor_list_d, st_n_atom_pairs_int, gpu_stream)

                  st_n_atom_pairs_double = j_end*sizeof(rjs(1))
                  call gpu_malloc_async(rjs_d, st_n_atom_pairs_double, gpu_stream)
                  call cpy_htod(c_loc(rjs), rjs_d, st_n_atom_pairs_double, gpu_stream)

                  ! I don't actually need xyz_d here at all!
                  call gpu_get_pair_distribution_nk(i_beg, i_end, j_end, n_sites, neighbor_list_d, &
                       n_neigh_d, neighbor_species_d, species_d,&
                       & rjs_d, xyz_d, params%r_range_min, params%r_range_max, params%pair_distribution_rcut, 6.d0&
                       &*params%pair_distribution_kde_sigma,&
                       & nk_d(n_dim_idx), nk_flags_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), j, k, gpu_stream)

                  call gpu_free_async(n_neigh_d, gpu_stream)
                  call gpu_free_async(species_d, gpu_stream)
                  call gpu_free_async(neighbor_species_d, gpu_stream)
                  call gpu_free_async(neighbor_list_d, gpu_stream)
                  call gpu_free_async(rjs_d, gpu_stream)

                  call gpu_free_async(nk_flags_d(n_dim_idx), gpu_stream)

                  st_nk_temp = int(1, c_size_t)*c_int
                  call cpy_dtoh(nk_d(n_dim_idx), c_loc(nk_temp), st_nk_temp, gpu_stream)

                  nk(n_dim_idx) = nk_temp(1)

                  ! call estimate_device_memory_usage( n_sites, nk(n_dim_idx), params%pair_distribution_n_samples,&
                  !      params%structure_factor_n_samples, total_memory_usage )

                  ! Now we create temporary arrays for the k indices
                  ! -------------------- Setting k --------------------
                  st_k_index_d(n_dim_idx) = int(nk(n_dim_idx), c_size_t)*c_int
                  call gpu_malloc_async(k_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)
                  call gpu_memset_async(k_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)

                  call gpu_set_pair_distribution_k_index_only(j_end, k_index_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream)

                  allocate (gpu_host_storage(n_dim_idx)%k_index_h(1:nk(n_dim_idx)))
              call cpy_dtoh(k_index_d(n_dim_idx), c_loc(gpu_host_storage(n_dim_idx)%k_index_h), st_k_index_d(n_dim_idx), gpu_stream)
!                 The plain hipFree here synchronised the whole device, which is
!                 what made the host buffer above valid. hipFreeAsync does not,
!                 so the sync the author left commented out has to be real, and
!                 has to come before the free.
                  call gpu_stream_sync(gpu_stream)
                  call gpu_free_async(k_index_d(n_dim_idx), gpu_stream)

                  ! -------------------- Setting j2 --------------------
                  call gpu_malloc_async(j2_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)
                  call gpu_memset_async(j2_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)

                  call gpu_malloc_async(neighbor_list_d, st_n_atom_pairs_int, gpu_stream)
                  call cpy_htod(c_loc(neighbors_list), neighbor_list_d, st_n_atom_pairs_int, gpu_stream)

                call gpu_set_pair_distribution_j2_only(j_end, n_sites, neighbor_list_d, j2_index_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream )

                  call gpu_free_async(neighbor_list_d, gpu_stream)
                  allocate (gpu_host_storage(n_dim_idx)%j2_index_h(1:nk(n_dim_idx)))

                  call cpy_dtoh(j2_index_d(n_dim_idx), c_loc(gpu_host_storage(n_dim_idx)%j2_index_h), st_k_index_d(n_dim_idx), &
                                gpu_stream)
!                 The plain hipFree here synchronised the whole device, which is
!                 what made the host buffer above valid. hipFreeAsync does not,
!                 so the sync the author left commented out has to be real, and
!                 has to come before the free.
                  call gpu_stream_sync(gpu_stream)
                  call gpu_free_async(j2_index_d(n_dim_idx), gpu_stream)

                  ! -------------------- Setting xyz --------------------
                  call gpu_malloc_async(xyz_d, 3*st_n_atom_pairs_double, gpu_stream)
                  call cpy_htod(c_loc(xyz), xyz_d, 3*st_n_atom_pairs_double, gpu_stream)

                  st_rjs(n_dim_idx) = int(nk(n_dim_idx), c_size_t)*c_double
                  call gpu_malloc_async(xyz_k_d(n_dim_idx), 3*st_rjs(n_dim_idx), gpu_stream)
                  call gpu_memset_async(xyz_k_d(n_dim_idx), 0, 3*st_rjs(n_dim_idx), gpu_stream)

                  call gpu_set_pair_distribution_xyz_only(j_end, xyz_d, xyz_k_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream)

                  allocate (gpu_host_storage(n_dim_idx)%xyz_k_h(1:3, 1:nk(n_dim_idx)))
                  call cpy_dtoh(xyz_k_d(n_dim_idx), c_loc(gpu_host_storage(n_dim_idx)%xyz_k_h), 3*st_rjs(n_dim_idx), gpu_stream)
                  call gpu_free_async(xyz_d, gpu_stream)
                  call gpu_stream_sync(gpu_stream)

                  ! -------------------- Setting rjs --------------------
                  st_n_atom_pairs_double = j_end*sizeof(rjs(1))
                  call gpu_malloc_async(rjs_d, st_n_atom_pairs_double, gpu_stream)
                  call cpy_htod(c_loc(rjs), rjs_d, st_n_atom_pairs_double, gpu_stream)

                  st_rjs(n_dim_idx) = int(nk(n_dim_idx), c_size_t)*c_double
                  call gpu_malloc_async(rjs_index_d(n_dim_idx), st_rjs(n_dim_idx), gpu_stream)
                  call gpu_memset_async(rjs_index_d(n_dim_idx), 0, st_rjs(n_dim_idx), gpu_stream)

                call gpu_set_pair_distribution_rjs_only(j_end, rjs_d, rjs_index_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), gpu_stream)

                  call gpu_free_async(rjs_d, gpu_stream)

                  ! No need for device storage here

                  call gpu_meminfo()
                  call gpu_free_async(nk_flags_sum_d(n_dim_idx), gpu_stream)
                  call gpu_stream_sync(gpu_stream)

                  st_pair_distribution_partial_d(n_dim_idx) = int(params%pair_distribution_n_samples, c_size_t)*c_double
                call gpu_malloc_async(pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), gpu_stream)
             call gpu_memset_async(pair_distribution_partial_d(n_dim_idx), 0, st_pair_distribution_partial_d(n_dim_idx), gpu_stream)

                  pdf_factor = ((params%r_range_max - params%r_range_min)/ &
                                dfloat(params%pair_distribution_n_samples))/(sqrt(2.d0*pi)*params%pair_distribution_kde_sigma)

                  der_factor = v_uc/n_atoms_of_species(j)/ &
                       & n_atoms_of_species(k)/factors(n_dim_idx)

                  call gpu_stream_sync(gpu_stream)

                  call gpu_get_pair_distribution_only( &
                     pair_distribution_partial_d(n_dim_idx), &
                     nk(n_dim_idx), &
                     params%pair_distribution_n_samples, &
                     params%pair_distribution_kde_sigma, &
                     x_d, dV_d, &
                     rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)

                  call gpu_copy_pdf(params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                                    pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                                    gpu_stream)

                  call gpu_free_async(pair_distribution_partial_d(n_dim_idx), gpu_stream)

                  call gpu_stream_sync(gpu_stream)

            st_pair_distribution_partial_der_d(n_dim_idx) = int(params%pair_distribution_n_samples, c_size_t)*nk(n_dim_idx)*c_double
        call gpu_malloc_async(pair_distribution_partial_der_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)
     call gpu_memset_async(pair_distribution_partial_der_d(n_dim_idx), 0, st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

                  call gpu_get_pair_distribution_der_only( &
                     pair_distribution_partial_der_d(n_dim_idx), &
                     nk(n_dim_idx), &
                     params%pair_distribution_n_samples, &
                     params%pair_distribution_kde_sigma, &
                     x_d, dV_d, &
                     rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)

                  ! call gpu_get_pair_distribution_and_ders(&
                  !      pair_distribution_partial_d(n_dim_idx),&
                  !      pair_distribution_partial_der_d(n_dim_idx),&
                  !      nk(n_dim_idx), &
                  !      params%pair_distribution_n_samples, &
                  !      params%pair_distribution_kde_sigma, &
                  !      x_d, dV_d,&
                  !      rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)

                  ! We can deallocate the rjs as we don't need them any more
                  call gpu_free_async(rjs_index_d(n_dim_idx), gpu_stream)

                  call gpu_stream_sync(gpu_stream)

       allocate (gpu_host_storage(n_dim_idx)%pair_distribution_partial_der_h(1:params%pair_distribution_n_samples, 1:nk(n_dim_idx)))
                call cpy_dtoh( pair_distribution_partial_der_d(n_dim_idx), c_loc( gpu_host_storage( n_dim_idx ) % pair_distribution_partial_der_h ), st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream )
                  call gpu_free_async(pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

                  ! call gpu_copy_pdf_der( params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                  !      pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                  !      gpu_stream )

                  n_dim_idx = n_dim_idx + 1

                  call gpu_meminfo()

                  if (n_dim_idx > n_dim_partial) then
                     exit outera
                  end if

               end do
            end do outera

         else
            n_dim_idx = 1
            call gpu_stream_sync(gpu_stream)
            outer1: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle ! We have already calculated the pair correlation function!

                  ! Note that with the calculation of the derivatives here,
                  ! this is without the -2 * delta_ik (r_j^alpha - r_i^alpha)
                  ! factor, which allows for some freeing of memory

                  ! Get all nk_flags
                  call gpu_meminfo()

                  st_nk_temp = int(1, c_size_t)*c_int
                  call gpu_malloc_async(nk_d(n_dim_idx), st_nk_temp, gpu_stream)
                  st_nk_flags = int(j_end, c_size_t)*c_int
                  call gpu_malloc_async(nk_flags_d(n_dim_idx), st_nk_flags, gpu_stream)
                  call gpu_memset_async(nk_flags_d(n_dim_idx), 0, st_nk_flags, gpu_stream)
                  call gpu_malloc_async(nk_flags_sum_d(n_dim_idx), st_nk_flags, gpu_stream)

                  call gpu_get_pair_distribution_nk(i_beg, i_end, j_end, n_sites, neighbor_list_d, &
                       n_neigh_d, neighbor_species_d, species_d,&
                       & rjs_d, xyz_d, params%r_range_min, params%r_range_max, params%pair_distribution_rcut, 6.d0&
                       &*params%pair_distribution_kde_sigma,&
                       & nk_d(n_dim_idx), nk_flags_d(n_dim_idx), nk_flags_sum_d(n_dim_idx), j, k, gpu_stream)

                  call gpu_free_async(nk_flags_d(n_dim_idx), gpu_stream)

                  ! Now copy the value of nk from the gpu
!                print *, "out of pdf nk kernel "
                  st_nk_temp = int(1, c_size_t)*c_int
                  call cpy_dtoh(nk_d(n_dim_idx), c_loc(nk_temp), st_nk_temp, gpu_stream)

                  nk(n_dim_idx) = nk_temp(1)

                  call estimate_device_memory_usage(n_sites, 0, nk(n_dim_idx), params%pair_distribution_n_samples, &
                                                    params%structure_factor_n_samples, total_memory_usage, .false.)

                  ! Now we create temporary arrays for the k indices

                  st_k_index_d(n_dim_idx) = int(nk(n_dim_idx), c_size_t)*c_int
                  call gpu_malloc_async(k_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)
                  call gpu_memset_async(k_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)

                  call gpu_malloc_async(j2_index_d(n_dim_idx), st_k_index_d(n_dim_idx), gpu_stream)
                  call gpu_memset_async(j2_index_d(n_dim_idx), 0, st_k_index_d(n_dim_idx), gpu_stream)

                  st_rjs(n_dim_idx) = int(nk(n_dim_idx), c_size_t)*c_double
                  call gpu_malloc_async(rjs_index_d(n_dim_idx), st_rjs(n_dim_idx), gpu_stream)
                  call gpu_memset_async(rjs_index_d(n_dim_idx), 0, st_rjs(n_dim_idx), gpu_stream)

                  call gpu_malloc_async(xyz_k_d(n_dim_idx), 3*st_rjs(n_dim_idx), gpu_stream)
                  call gpu_memset_async(xyz_k_d(n_dim_idx), 0, 3*st_rjs(n_dim_idx), gpu_stream)

                  call gpu_meminfo()
                  call gpu_set_pair_distribution_k_index(i_beg, i_end, j_end, n_sites, neighbor_list_d,&
                       & rjs_d, xyz_d, k_index_d(n_dim_idx), j2_index_d(n_dim_idx),&
                       & rjs_index_d(n_dim_idx), xyz_k_d(n_dim_idx), nk_flags_d(n_dim_idx), nk_flags_sum_d(n_dim_idx),&
                       & gpu_stream)

                  call gpu_free_async(nk_flags_sum_d(n_dim_idx), gpu_stream)

                  !--- CALCULATING THE PAIR DISTRIBUTION FUNCTION ---!
                  !             print *, " Allocating pdf arrays   "
                  st_pair_distribution_partial_d(n_dim_idx) = int(params%pair_distribution_n_samples, c_size_t)*c_double
                call gpu_malloc_async(pair_distribution_partial_d(n_dim_idx), st_pair_distribution_partial_d(n_dim_idx), gpu_stream)
             call gpu_memset_async(pair_distribution_partial_d(n_dim_idx), 0, st_pair_distribution_partial_d(n_dim_idx), gpu_stream)

            st_pair_distribution_partial_der_d(n_dim_idx) = int(params%pair_distribution_n_samples, c_size_t)*nk(n_dim_idx)*c_double
        call gpu_malloc_async(pair_distribution_partial_der_d(n_dim_idx), st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)
     call gpu_memset_async(pair_distribution_partial_der_d(n_dim_idx), 0, st_pair_distribution_partial_der_d(n_dim_idx), gpu_stream)

                  pdf_factor = ((params%r_range_max - params%r_range_min)/ &
                                dfloat(params%pair_distribution_n_samples))/(sqrt(2.d0*pi)*params%pair_distribution_kde_sigma)

                  der_factor = v_uc/n_atoms_of_species(j)/ &
                       & n_atoms_of_species(k)/factors(n_dim_idx)

                  call gpu_get_pair_distribution_and_ders( &
                     pair_distribution_partial_d(n_dim_idx), &
                     pair_distribution_partial_der_d(n_dim_idx), &
                     nk(n_dim_idx), &
                     params%pair_distribution_n_samples, &
                     params%pair_distribution_kde_sigma, &
                     x_d, dV_d, &
                     rjs_index_d(n_dim_idx), pdf_factor, der_factor, gpu_stream)

                  ! call get_pair_distribution( n_sites,&
                  !      & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end)&
                  !      &, neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
                  !      & xyz(1:3,j_beg:j_end), params %r_range_min, params&
                  !      &%r_range_max, params%pair_distribution_n_samples,&
                  !      & x_pair_distribution,&
                  !      & pair_distribution_partial(1:params&
                  !      &%pair_distribution_n_samples, n_dim_idx), params &
                  !      &%pair_distribution_rcut, .false., params&
                  !      &%pair_distribution_partial, j, k, params&
                  !      &%pair_distribution_kde_sigma, dfloat(n_sites)/v_uc&
                  !      &,  params%exp_forces,&
                  !      & pair_distribution_partial_der, n_dim_idx, j_beg,&
                  !      & j_end )

                  ! -------- PDF CHECK ----------!
                  !--- CHECKING THAT IT WORKS ---!
                  ! print *, "checking pdf"
                  ! allocate(pdf_gpu_check(1:params%pair_distribution_n_samples))

                  !-------- PDF DER CHECK -------!

                  ! will need to do a separate function to copy the pair distribution
                  !             call gpu_device_sync()
                  call gpu_copy_pdf(params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                                    pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                                    gpu_stream)

                  ! call gpu_copy_pdf_der( params%pair_distribution_n_samples, pair_distribution_partial_d(n_dim_idx), &
                  !      pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), &
                  !      gpu_stream )

                  ! We can deallocate the rjs as we don't need them any more
                  call gpu_free_async(rjs_index_d(n_dim_idx), gpu_stream)

                  n_dim_idx = n_dim_idx + 1

                  call gpu_meminfo()

                  if (n_dim_idx > n_dim_partial) then
                     exit outer1
                  end if

               end do
            end do outer1
         end if
         call gpu_free_async(x_d, gpu_stream)
         call gpu_free_async(dV_d, gpu_stream)

         deallocate (rjs_index_d, st_rjs)
         deallocate (nk_flags_d)
         deallocate (nk_flags_sum_d)
         call gpu_meminfo()

      else
         call get_pair_distribution(n_sites, &
              & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
              & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3, j_beg:j_end),&
              & params%r_range_min, params%r_range_max, params &
              &%pair_distribution_n_samples, x_pair_distribution,&
              & y_pair_distribution, params &
              &%pair_distribution_rcut, .false., .false., 1, 1,&
              & params%pair_distribution_kde_sigma, dfloat(n_sites)&
              &/v_uc, params%do_forces .and. params%exp_forces, pair_distribution_partial_der, 1, &
              & j_beg, j_end)
      end if

      deallocate (dV)

      ! --- MPI communication is here  ---

      if (params%pair_distribution_partial) then
#ifdef _MPIF90
         call mpi_reduce(pair_distribution_partial,&
              & pair_distribution_partial_temp, params&
              &%pair_distribution_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
              & 0, MPI_COMM_WORLD, ierr)

         ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
         ! Note, we have only so far divided by 4 pi r^2 dr
         ! Therefore, we must scale by the density

         pair_distribution_partial = pair_distribution_partial_temp
         deallocate (pair_distribution_partial_temp)

         call mpi_bcast(pair_distribution_partial, params&
              &%pair_distribution_n_samples*n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
              & MPI_COMM_WORLD, ierr)

         ! Now, we have the derivatives of the partial pair distribution
         ! function with respect to the atom pairs in that rank
         !
         ! We can keep them in the rank and calculate forces

         ! call mpi_reduce(pair_distribution_partial_der,&
         !      & pair_distribution_partial_der_temp, params&
         !      &%pair_distribution_n_samples * n_species *&
         !      & n_species * 3 * n_pairs_tot, MPI_DOUBLE_PRECISION, MPI_SUM,&
         !      & 0, MPI_COMM_WORLD, ierr)

         ! ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
         ! ! Note, we have only so far divided by 4 pi r^2 dr
         ! ! Therefore, we must scale by the density

#endif

         if (params%valid_pdf) then
            allocate (energies_pair_distribution(1:n_sites))
            energies_pair_distribution = 0.d0

            if (params%do_forces .and. params%exp_forces) then
               allocate (forces_pair_distribution(1:3, 1:n_sites))
               forces_pair_distribution = 0.d0
            end if

         end if

         !###---   Accumulate the PDF   ---###!

         y_pair_distribution = 0.d0
         n_dim_idx = 1
         outer2: do j = 1, n_species
            do k = 1, n_species

               if (j > k) cycle

               pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) =&
                    & pair_distribution_partial(1:params&
                    &%pair_distribution_n_samples, n_dim_idx)*v_uc &
                    &  /n_atoms_of_species(j)/n_atoms_of_species(k)/factors(n_dim_idx) !real(n_sites)

               if (params%do_forces .and. params%exp_forces) then
                  pair_distribution_partial_der(1:params&
                       &%pair_distribution_n_samples, n_dim_idx, &
                       & j_beg:j_end) = pair_distribution_partial_der(1:params&
                       & %pair_distribution_n_samples, n_dim_idx, &
                       & j_beg:j_end)*v_uc/n_atoms_of_species(j)/ &
                       & n_atoms_of_species(k)/factors(n_dim_idx)!real(n_sites)

               end if

               y_pair_distribution(1:params%pair_distribution_n_samples) = &
                    & y_pair_distribution(1:params%pair_distribution_n_samples) +  &
                    &  factors(n_dim_idx)*(n_atoms_of_species(j)*n_atoms_of_species(k))* &
                    & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) &
                    &  /dfloat(n_sites)/dfloat(n_sites)

               n_dim_idx = n_dim_idx + 1

               if (n_dim_idx > n_dim_partial) then
                  exit outer2
               end if

            end do
         end do outer2

         ! --- Preprocess the pair distribution according to the output --- !
         if (trim(params%pair_distribution_output) == "D(r)") then
            y_pair_distribution = 4.d0*pi*(dfloat(n_sites)/v_uc)*x_pair_distribution*(y_pair_distribution - 1.d0)
         end if

         !###---   Calculate the forces   ---###!

         if (params%valid_pdf .and. allocated(params%exp_energy_scales)) then

            call get_energy_scale(params%do_md, params%do_mc,&
                 & md_istep, params%md_nsteps, mc_istep, params&
                 &%mc_nsteps, params &
                 &%exp_energy_scales_initial(params%pdf_idx), params &
                 &%exp_energy_scales_final(params%pdf_idx), params &
                 &%exp_energy_scales(params%pdf_idx))

            call get_exp_energies(params%exp_energy_scales(params&
                 &%pdf_idx), params%exp_data(params%pdf_idx)%y&
                 &, y_pair_distribution,&
                 & params%pair_distribution_n_samples, n_sites,&
                 & energies_pair_distribution(i_beg:i_end), params%exp_data(params%pdf_idx)%w)

            if (params%do_forces .and. params%exp_forces) then

               n_dim_idx = 1
               outerforces: do j = 1, n_species
                  do k = 1, n_species

                     if (j > k) cycle

                     call get_pair_distribution_forces(n_sites, params%exp_energy_scales(params%pdf_idx),&
                          & params%exp_data(params%pdf_idx)%x, params%exp_data(params%pdf_idx)%y,&
                          & forces_pair_distribution, virial,&
                          & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                          & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3, j_beg:j_end), params&
                          &%r_range_min, params%r_range_max, params&
                          &%pair_distribution_n_samples,&
                          & y_pair_distribution(1:params&
                          &%pair_distribution_n_samples), params%pair_distribution_rcut&
                          &, j, k, pair_distribution_partial_der(1:params &
                          &%pair_distribution_n_samples, n_dim_idx,&
                          & j_beg:j_end), params%pair_distribution_partial,&
                          & params%pair_distribution_kde_sigma,&
                          & ((n_atoms_of_species(j)*&
                          & n_atoms_of_species(k))/dfloat(n_sites)/&
                          & dfloat(n_sites)), (dfloat(n_sites)/v_uc), params%pair_distribution_output)

                     n_dim_idx = n_dim_idx + 1

                     if (n_dim_idx > n_dim_partial) then
                        exit outerforces
                     end if

                  end do
               end do outerforces

!              ---   The cell half of the virial   --- !
!
!              get_pair_distribution_forces differentiates the interatomic
!              distances; the pattern also depends on the cell directly,
!              because each partial is normalised by the number density, and a
!              homogeneous strain changes the volume as well as the distances.
!              That part is a property of the whole pattern rather than of one
!              (a,b) channel, so it is added here, once, and only by rank 0 --
!              the virial is summed over ranks afterwards.
!
!              g(r) is proportional to V outright, so V dy/dV is y itself.
!              D(r) = 4 pi rho r ( g(r) - 1 ) has the volume cancel in the g
!              term, leaving only the subtracted background, which goes as 1/V.
               if (rank == 0) then
                  if (trim(params%pair_distribution_output) == "D(r)") then
                     dedv = -params%exp_energy_scales(params%pdf_idx)*&
                          & sum((y_pair_distribution(1:params%pair_distribution_n_samples) - &
                          &      params%exp_data(params%pdf_idx)%y)*&
                          &     4.d0*pi*(dfloat(n_sites)/v_uc)*&
                          &     x_pair_distribution(1:params%pair_distribution_n_samples))
                  else
                     dedv = -params%exp_energy_scales(params%pdf_idx)*&
                          & sum((y_pair_distribution(1:params%pair_distribution_n_samples) - &
                          &      params%exp_data(params%pdf_idx)%y)*&
                          &     y_pair_distribution(1:params%pair_distribution_n_samples))
                  end if

                  do i = 1, 3
                     virial(i, i) = virial(i, i) + dedv
                  end do
               end if
            end if

            ! open(unit=1234, file="grad", status="unknown")
            ! do i = 100, 110
            !    write(1234,  '(A,1X,I8,1X,F20.8)'), "dg_dr_0^1 ", i, pair_distribution_der_temp( i )
            ! end do
            ! close(unit=1234)

         end if

         !###---   If not doing partial pair distribution functions   ---###!

      else
#ifdef _MPIF90
         call mpi_reduce(y_pair_distribution,&
              & y_pair_distribution_temp, params&
              &%pair_distribution_n_samples,&
              & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
              & MPI_COMM_WORLD, ierr)

         y_pair_distribution = y_pair_distribution_temp
         deallocate (y_pair_distribution_temp)

         call mpi_bcast(y_pair_distribution, params&
              &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
              & MPI_COMM_WORLD, ierr)

         !    call mpi_reduce(forces_pair_distribution,&
         !         & forces_pair_distribution_temp, 3 * n_sites,&
         !         & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         !         & MPI_COMM_WORLD, ierr)

         !    call mpi_bcast(forces_pair_distribution, 3*n_sites, MPI_DOUBLE_PRECISION, 0,&
         !         & MPI_COMM_WORLD, ierr)

#endif
         y_pair_distribution = y_pair_distribution* &
              & dot_product(cross_product(a_box, b_box),&
              & c_box)/(dfloat(indices(1)*indices(2)&
              &*indices(3)))/dfloat(n_sites)/dfloat(n_sites)

      end if

      if (params%pair_distribution_partial .and. allocated(factors)) deallocate (factors)

      ! Write out the partial pair distribution functions
      call get_write_condition(params%do_mc, params%do_md&
           &, mc_istep, md_istep, params%write_xyz,&
           & write_condition)

      if (rank == 0 .and. params%write_pair_distribution .and. write_condition) then
         ! call write_partial_exp(params%do_mc, params%do_md, mc_istep, md_istep,&
         !      & params%write_xyz, params%pair_distribution_partial,&
         !      & n_species, params%pair_distribution_n_samples,&
         !      & n_dim_partial , x_pair_distribution(1:params&
         !      &%pair_distribution_n_samples), y_pair_distribution(1:params&
         !      &%pair_distribution_n_samples), pair_distribution_partial(1:params &
         !      &%pair_distribution_n_samples, 1:n_dim_partial),&
         !      & species_types , 'pair_distribution')

         call get_overwrite_condition(params%do_mc, params%do_md,&
              & mc_istep, md_istep, params%write_xyz,&
              & overwrite_condition)

         if (params%pair_distribution_partial) then
            n_dim_idx = 1
            outer3: do j = 1, n_species
               do k = 1, n_species

                  if (j > k) cycle

                  write (filename, '(A)')&
                       & 'pair_distribution_'//trim(params&
                       &%species_types(j))//'_'//trim(params&
                       &%species_types(k))//&
                       & "_prediction.dat"
                  call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
                       & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx),&
                       & overwrite_condition, filename, 'pair_distribution')

                  n_dim_idx = n_dim_idx + 1
                  if (n_dim_idx > n_dim_partial) then
                     exit outer3
                  end if

               end do
            end do outer3
         end if

         write (filename, '(A)')&
              & "pair_distribution_total.dat"
         call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
              &y_pair_distribution(1:params%pair_distribution_n_samples),&
              & overwrite_condition, filename, "pair_distribution  output: "//trim(params&
              &%pair_distribution_output))

      end if

   end subroutine gpu_calculate_pair_distribution
#endif

end module exp_interface
