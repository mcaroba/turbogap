! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, gap.f90, is copyright (c) 2019-2022, Miguel A. Caro
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

   use kinds

   use splines
#ifdef _GPU
   use F_B_C
   use iso_c_binding
!  The default stream and cuBLAS handle come from gpu_context, not from the
!  argument list. The driver passes the GLOBAL ones here -- the SOAP loop is
!  not OpenMP-parallel, unlike the exp and estat batch loops -- so the names
!  are unchanged and the bodies below are untouched. Same trick as
!  gap_backend_gpu's ten buffers.
   use gpu_context, only: cublas_handle, gpu_stream
   use mpi

!  The pair->site index for the SOAP path. It is BUILT here and consumed by
!  local_properties, and until now it was allocated here and freed in
!  gap_interface -- split ownership across two modules. Now it is allocated and
!  released in one place, and reached from outside through soap_l_index_d()
!  rather than by being passed back out of the argument list.
!
!  It is private because n_neigh_d, soap_d and l_index_d are each declared in
!  several modules here (7, 4 and 4 respectively) and gap.f90 has no default
!  private, so publishing module state under those names collides on use.
   type(c_ptr) :: l_index_d

!  The sparse-set arrays for the descriptor currently being evaluated. They
!  were uploaded and released by the DRIVER, two levels above the only code
!  that reads them, and passed down through get_gap_soap purely to reach here.
!  soap_backend_begin/end take over that lifetime so neither the driver nor
!  get_gap_soap names a device pointer for them.
   type(c_ptr) :: alphas_d
   type(c_ptr) :: Qs_d

   private :: cublas_handle, gpu_stream, l_index_d, alphas_d, Qs_d

contains

   subroutine soap_backend_begin(hypers)
!     Upload the sparse set for one descriptor. Pairs with soap_backend_end.
      use types, only: soap_turbo
      implicit none
      type(soap_turbo), intent(in), target :: hypers
      integer(c_size_t) :: st

      st = hypers%n_sparse*sizeof(hypers%nf(1))
      call gpu_malloc_async(alphas_d, st, gpu_stream)
      call cpy_htod(c_loc(hypers%alphas), alphas_d, st, gpu_stream)

      st = hypers%n_sparse*hypers%dim*sizeof(hypers%nf(1))
      call gpu_malloc_async(Qs_d, st, gpu_stream)
      call cpy_htod(c_loc(hypers%Qs), Qs_d, st, gpu_stream)
   end subroutine soap_backend_begin

   subroutine soap_backend_end()
!     Release what soap_backend_begin uploaded.
      implicit none
      call gpu_free_async(alphas_d, gpu_stream)
      call gpu_free_async(Qs_d, gpu_stream)
   end subroutine soap_backend_end

   function soap_l_index_d() result(ptr)
!     The pair->site index built by get_soap_energy_and_forces, for the
!     local-properties path that consumes it. A function rather than public
!     module state, so the name cannot collide on use.
      implicit none
      type(c_ptr) :: ptr
      ptr = l_index_d
   end function soap_l_index_d

   subroutine soap_backend_release()
!     Release what the SOAP path allocated. Pairs with the allocation in
!     get_soap_energy_and_forces; gap_interface used to do this by hand.
      implicit none
      call gpu_free_async(l_index_d, gpu_stream)
   end subroutine soap_backend_release

   subroutine get_soap_energy_and_forces(n_sparse, soap, soap_der, delta, zeta0, e0, &
        n_neigh, neighbors_list, xyz, do_forces, do_timing, &
        energies, forces, virial, solo_time_soap, soap_d, &
        soap_der_d, n_neigh_d,&
        & n_pairs)
      !   soap(1:n_soap, 1:n_sites)

      !use mpi
      implicit none

      integer(c_int), intent(in) :: n_sparse
      !   real(c_double), intent(in),target :: soap(:,:), soap_der(:,:,:), alphas(:), delta, Qs(:,:), e0, zeta0, xyz(:,:)
      real(c_double), intent(in), target :: soap(:, :)
      real(c_double), intent(in), target :: soap_der(:, :, :)
      real(c_double), intent(in), target :: delta
      real(c_double), intent(in), target :: e0
      real(c_double), intent(in), target :: zeta0
      real(c_double), intent(in), target :: xyz(:, :)
      real(c_double), intent(out):: energies(:)
      real(c_double), intent(out):: forces(:, :)
      real(c_double), intent(out):: virial(1:3, 1:3)
      real(c_double), allocatable, target :: tmp_energies(:)
      real(c_double), allocatable, target :: tmp_forces(:, :)
      real(c_double), allocatable, target :: tmp_virial(:, :)
      integer(c_int), intent(in), target :: n_neigh(:)
      integer(c_int), intent(in), target :: neighbors_list(:)
      logical, intent(in) :: do_forces
      logical, intent(in) :: do_timing
      real(c_double), allocatable, target :: kernels(:, :)
      real(c_double), allocatable, target :: kernels_der(:, :)
      real(c_double), allocatable, target :: Qss(:, :)
      real(c_double), allocatable, target :: Qs_copy(:, :)
      real(c_double), allocatable, target :: this_Qss(:)
      real(c_double), allocatable, target :: kernels_copy(:, :)
      real(c_double), allocatable, target :: this_force_h(:, :)
      real(c_double) :: time1
      real(c_double) :: time2
      real(c_double) :: time3
      real(c_double) :: energies_time
      real(c_double) :: forces_time
      real(c_double) :: this_force(1:3)
      real(c_double) :: zeta
      real(c_double) :: cdelta_ene
      real(c_double) :: cdelta_force
      real(c_double) :: mzetam
      integer(c_int) :: n_sites
      integer(c_int) :: n_soap
      integer(c_int) :: i
      integer(c_int) :: j
      integer(c_int) :: k
      integer(c_int) :: l
      integer(c_int) :: j2
      integer(c_int) :: zeta_int
      integer(c_int) :: n_sites0
      integer(c_int) :: k1
      integer(c_int) :: k2
      logical :: is_zeta_int = .false.
      !   type(c_ptr) :: kernels_copy_d, kernels_d, Qs_d, energies_d, alphas_d
      type(c_ptr) :: kernels_copy_d
      type(c_ptr) :: kernels_d
      type(c_ptr) :: energies_d
      type(c_ptr) :: kernels_der_d
      type(c_ptr) :: Qss_d
      type(c_ptr) :: Qs_copy_d !
      type(c_ptr) :: this_Qss_d
      integer(c_int) :: size_kernels
      integer(c_int) :: size_soap
      integer(c_int) :: size_Qs
      integer(c_int) :: size_alphas
      integer(c_int) :: size_energies
      integer(c_int) :: maxnn
      integer(c_int) :: size_nnlist
      integer(c_int) :: size_xyz
      integer(c_int) :: n1xyz
      integer(c_int) :: n2xyz
      integer(c_int) :: n1forces
      integer(c_int) :: n2forces
      integer(c_int) :: n1virial
      integer(c_int) :: n2virial
      integer(c_int) :: size_forces
      integer(c_int) :: size_virial
      integer(c_int) :: size_soap_der
      integer(c_int) :: n1soap_der
      integer(c_int) :: n2soap_der
      integer(c_int) :: n3soap_der
      integer(c_int) :: rank
      integer(c_int) :: ierr
      integer(c_int), intent(inout) :: n_pairs
      integer(c_int), allocatable, target :: neighbors_beg(:)
      integer(c_int), allocatable, target :: neighbors_end(:)
      type(c_ptr) :: virial_d
      type(c_ptr) :: n_neigh_d
      type(c_ptr) :: this_force_d
      type(c_ptr) :: j2_index_d
      type(c_ptr) :: neighbors_beg_d
      type(c_ptr) :: neighbors_end_d
      type(c_ptr) :: xyz_d
      type(c_ptr) :: neighbors_list_d
      type(c_ptr) :: forces_d
      real(dp), intent(inout) :: solo_time_soap
      integer(c_int), allocatable, target :: j2_index(:)
      integer(c_int), allocatable, target :: l_index(:)
      type(c_ptr), intent(inout) :: soap_der_d
      type(c_ptr), intent(inout) :: soap_d
      real(dp) :: ttt(2)
      integer(c_size_t) :: st_alphas
      integer(c_size_t) :: st_Qs
      integer(c_size_t) :: st_kernels
      integer(c_size_t) :: st_energies
      integer(c_size_t) :: st_soap
      integer(c_size_t) :: st_n_pairs
      integer(c_size_t) :: st_xyz
      integer(c_size_t) :: st_nnlist
      integer(c_size_t) :: st_forces
      integer(c_size_t) :: st_virial
      integer(c_size_t) :: st_neigh
      integer(c_size_t) :: st_neigh_beg
      integer(c_size_t) :: st_neigh_end

#ifdef _MPIF90
      ttt(1) = MPI_Wtime()
#else
      call cpu_time(ttt(1))
#endif

      cdelta_ene = delta*delta
      if (dabs(zeta0 - dfloat(int(zeta0))) < 1.d-5) then
         is_zeta_int = .true.
         zeta_int = int(zeta0)
         zeta = dfloat(zeta_int)
      else
         zeta = zeta0
      end if

      !   Energies
      if (do_timing) then
         call cpu_time(time1)
         time3 = time1
      end if

      n_soap = size(soap, 1)
      n_sites = size(soap, 2)
      n_sites0 = size(forces, 2)

      allocate (kernels(1:n_sites, 1:n_sparse))
      kernels = 0.d0
      allocate (kernels_copy(1:n_sites, 1:n_sparse))

      size_kernels = n_sites*n_sparse
      size_soap = n_soap*n_sites
      size_Qs = n_soap*n_sparse

      size_energies = n_sites

      st_kernels = size_kernels*(sizeof(kernels(1, 1)))
      st_Qs = size_Qs*(sizeof(e0))
      st_energies = size_energies*(sizeof(energies(1)))

      call gpu_malloc_async(kernels_d, st_kernels, gpu_stream)
      call gpu_malloc_async(kernels_copy_d, st_kernels, gpu_stream)
      call gpu_malloc_async(energies_d, st_energies, gpu_stream)

      call gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites)
      call gpu_kernels_pow(kernels_d, kernels_copy_d, zeta, size_kernels, gpu_stream)
      call gpu_blas_mvmul_n(cublas_handle, kernels_copy_d, alphas_d, energies_d, n_sites, n_sparse)

      call gpu_axpe(energies_d, cdelta_ene, e0, size_energies, gpu_stream)

      allocate (tmp_energies(1:n_sites))
      energies = 0.d0

      if (do_forces) then

         st_soap = size_soap*sizeof(Qss(1, 1))

         call gpu_malloc_async(kernels_der_d, st_kernels, gpu_stream)
         call gpu_malloc_async(Qss_d, st_soap, gpu_stream)
         call gpu_malloc_async(Qs_copy_d, st_Qs, gpu_stream)

         call cpy_dtod(Qs_d, Qs_copy_d, st_Qs, gpu_stream)

         mzetam = zeta - 1
         call gpu_kernels_pow(kernels_d, kernels_der_d, mzetam, size_kernels, gpu_stream)

         if (n_sites < n_soap) then
            call gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse, gpu_stream)
         else
            call gpu_matvect(Qs_copy_d, alphas_d, n_soap, n_sparse, gpu_stream)
         end if

         cdelta_force = -zeta*delta**2
         call gpu_blas_mmul_n_t(cublas_handle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, &
                                n_soap, n_sites, cdelta_force)

         allocate (neighbors_beg(1:n_sites))
         allocate (neighbors_end(1:n_sites))
         l = 0
         maxnn = 0
         do i = 1, n_sites
            neighbors_beg(i) = l + 1
            do j = 1, n_neigh(i)
               l = l + 1
            end do
            neighbors_end(i) = l
            if (n_neigh(i) > maxnn) then
               maxnn = n_neigh(i)
            end if
         end do
         n_pairs = l

         allocate (l_index(1:n_pairs))

         l = 0
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               l = l + 1
               l_index(l) = i
            end do
         end do
         st_n_pairs = n_pairs*sizeof(l_index(1))
         call gpu_malloc_async(l_index_d, st_n_pairs, gpu_stream)
         call cpy_htod(c_loc(l_index), l_index_d, st_n_pairs, gpu_stream)

         allocate (j2_index(1:n_pairs))
         l = 0
         do i = 1, n_sites
            do j = 1, n_neigh(i)
               l = l + 1
               j2 = mod(neighbors_list(l) - 1, n_sites0) + 1
               j2_index(l) = j2
            end do
         end do
         st_n_pairs = n_pairs*sizeof(j2_index(1))
         call gpu_malloc_async(j2_index_d, st_n_pairs, gpu_stream)
         call cpy_htod(c_loc(j2_index), j2_index_d, st_n_pairs, gpu_stream)

         virial = 0.d0
         forces = 0.d0

         n1xyz = size(xyz, 1)
         n2xyz = size(xyz, 2)
         size_xyz = n1xyz*n2xyz

         size_nnlist = size(neighbors_list, 1)

         n1forces = size(forces, 1)
         n2forces = size(forces, 2)
         size_forces = n1forces*n2forces

         allocate (tmp_forces(1:n1forces, 1:n2forces))
         allocate (tmp_virial(1:3, 1:3))

         n1virial = size(virial, 1)
         n2virial = size(virial, 2)
         size_virial = n1virial*n2virial

         n1soap_der = size(soap_der, 1)
         n2soap_der = size(soap_der, 2)
         n3soap_der = size(soap_der, 3)
         size_soap_der = n1soap_der*n2soap_der*n3soap_der

         st_xyz = size_xyz*sizeof(xyz(1, 1))
         call gpu_malloc_async(xyz_d, st_xyz, gpu_stream)
         call cpy_htod(c_loc(xyz), xyz_d, st_xyz, gpu_stream)

         st_nnlist = size_nnlist*sizeof(neighbors_list(1))
         call gpu_malloc_async(neighbors_list_d, st_nnlist, gpu_stream)
         call cpy_htod(c_loc(neighbors_list), neighbors_list_d, st_nnlist, gpu_stream)

         st_forces = size_forces*sizeof(forces(1, 1))
         call gpu_malloc_async(forces_d, st_forces, gpu_stream)

         st_virial = size_virial*sizeof(virial(1, 1))
         call gpu_malloc_async(virial_d, st_virial, gpu_stream)

         st_neigh = n_sites*sizeof(n_neigh(1))

         st_neigh_end = n_sites*sizeof(neighbors_end(1))
         call gpu_malloc_async(neighbors_end_d, st_neigh_end, gpu_stream)
         call cpy_htod(c_loc(neighbors_end), neighbors_end_d, st_neigh_end, gpu_stream)

         st_neigh_beg = n_sites*sizeof(neighbors_beg(1))
         call gpu_malloc_async(neighbors_beg_d, st_neigh_beg, gpu_stream)
         call cpy_htod(c_loc(neighbors_beg), neighbors_beg_d, st_neigh_beg, gpu_stream)

         call gpu_final_soap_forces_virial(n_sites, &
                                           Qss_d, n_soap, l_index_d, j2_index_d, &
                                           soap_der_d, &
                                           xyz_d, virial_d, &
                                           n_sites0, &
                                           forces_d, &
                                           n_pairs, gpu_stream)

         call cpy_dtoh(forces_d, c_loc(tmp_forces), st_forces, gpu_stream)
         call cpy_dtoh(virial_d, c_loc(tmp_virial), st_virial, gpu_stream)

      end if

      call cpy_dtoh(energies_d, c_loc(tmp_energies), st_energies, gpu_stream)

      if (do_timing) then
         call cpu_time(time2)
         forces_time = time2 - time1
      end if

      ! Now moving deallocation outside of routine so that local properties can use them.

      call gpu_free_async(neighbors_list_d, gpu_stream)
      call gpu_free_async(neighbors_end_d, gpu_stream)
      call gpu_free_async(neighbors_beg_d, gpu_stream)
      call gpu_free_async(kernels_der_d, gpu_stream)
      call gpu_free_async(Qss_d, gpu_stream)
      call gpu_free_async(Qs_copy_d, gpu_stream)
      call gpu_free_async(xyz_d, gpu_stream)
      call gpu_free_async(kernels_d, gpu_stream)
      call gpu_free_async(kernels_copy_d, gpu_stream)

      if (do_forces) then
         call gpu_free_async(j2_index_d, gpu_stream)
         call gpu_free_async(forces_d, gpu_stream)
         call gpu_free_async(virial_d, gpu_stream)
      end if

      call gpu_free_async(energies_d, gpu_stream)

      energies = tmp_energies

      if (do_forces) then
         forces = tmp_forces
         virial = tmp_virial
      end if

      deallocate (kernels, kernels_copy, tmp_energies)

      if (allocated(j2_index)) deallocate (j2_index)
      if (allocated(l_index)) deallocate (l_index)

      if (do_forces) then
         deallocate (tmp_forces, tmp_virial)
         deallocate (neighbors_beg, neighbors_end)
      end if

#ifdef _MPIF90
      ttt(2) = MPI_Wtime()
#else
      call cpu_time(ttt(2))
#endif

      ttt(2) = MPI_Wtime()
      solo_time_soap = solo_time_soap + ttt(2) - ttt(1)

   end subroutine get_soap_energy_and_forces

   subroutine get_soap_dipole(n_sparse, soap, delta, zeta0, n_neigh, dipoles, energies, &
                              soap_d, soap_der_d)
      !   Local dipoles from a dipole GAP, on the device.
      !
      !   The model is a GAP whose local "energy" E_i is a fictitious scalar,
      !   fitted so that its derivative with respect to the central atom's OWN
      !   position is the local dipole:
      !
      !       mu_i = dE_i/dr_i   and   mu = sum_i mu_i
      !
      !   Everything up to Qss is identical to get_soap_energy_and_forces; the
      !   two differences are the sign (mu is a gradient, not a force, so the
      !   prefactor is +zeta*delta**2 rather than -zeta*delta**2) and the fact
      !   that only the SELF pair of each site is contracted rather than every
      !   pair scattered onto neighbours. soap_turbo builds that self pair as
      !   -sum_{j/=1} of the neighbour terms, which is exactly d(soap_i)/d(r_i),
      !   and build_neighbors_list puts it first for every site -- so the pair
      !   index wanted for site i is neighbors_beg(i).
      !
      !   E_i is returned too, because the kernel matrix it needs has already
      !   been formed here. It is a fitting artefact with no physical meaning:
      !   it is reported as energy_dipole and must never be added to the total
      !   energy, and its gradient must never be added to the forces.

      implicit none

      integer(c_int), intent(in) :: n_sparse
      real(c_double), intent(in), target :: soap(:, :)
      real(c_double), intent(in), target :: delta
      real(c_double), intent(in), target :: zeta0
      integer(c_int), intent(in), target :: n_neigh(:)
      real(c_double), intent(out), target :: dipoles(:, :)
      real(c_double), intent(out), target :: energies(:)
      type(c_ptr), intent(inout) :: soap_d
      type(c_ptr), intent(inout) :: soap_der_d

      real(c_double), allocatable, target :: kernels(:, :)
      real(c_double) :: zeta
      real(c_double) :: cdelta_ene
      real(c_double) :: cdelta_dipole
      real(c_double) :: mzetam
      integer(c_int) :: n_sites
      integer(c_int) :: n_soap
      integer(c_int) :: i
      integer(c_int) :: j
      integer(c_int) :: l
      integer(c_int) :: zeta_int
      integer(c_int) :: size_kernels
      integer(c_int) :: size_soap
      integer(c_int) :: size_Qs
      integer(c_int) :: size_energies
      logical :: is_zeta_int = .false.
      type(c_ptr) :: kernels_d
      type(c_ptr) :: kernels_copy_d
      type(c_ptr) :: kernels_der_d
      type(c_ptr) :: energies_d
      type(c_ptr) :: Qss_d
      type(c_ptr) :: Qs_copy_d
      type(c_ptr) :: dipoles_d
      type(c_ptr) :: beg_index_d
      integer(c_int), allocatable, target :: beg_index(:)
      integer(c_size_t) :: st_kernels
      integer(c_size_t) :: st_Qs
      integer(c_size_t) :: st_energies
      integer(c_size_t) :: st_soap
      integer(c_size_t) :: st_dipoles
      integer(c_size_t) :: st_beg_index

      n_soap = size(soap, 1)
      n_sites = size(soap, 2)

      dipoles = 0.d0
      energies = 0.d0
      if (n_sites == 0) then
         return
      end if

      cdelta_ene = delta*delta
      if (dabs(zeta0 - dfloat(int(zeta0))) < 1.d-5) then
         is_zeta_int = .true.
         zeta_int = int(zeta0)
         zeta = dfloat(zeta_int)
      else
         zeta = zeta0
      end if

      allocate (kernels(1:n_sites, 1:n_sparse))

      size_kernels = n_sites*n_sparse
      size_soap = n_soap*n_sites
      size_Qs = n_soap*n_sparse
      size_energies = n_sites

      st_kernels = size_kernels*(sizeof(kernels(1, 1)))
      st_Qs = size_Qs*(sizeof(delta))
      st_energies = size_energies*(sizeof(energies(1)))
      st_soap = size_soap*(sizeof(delta))
      st_dipoles = 3*n_sites*(sizeof(delta))

      call gpu_malloc_async(kernels_d, st_kernels, gpu_stream)
      call gpu_malloc_async(kernels_copy_d, st_kernels, gpu_stream)
      call gpu_malloc_async(energies_d, st_energies, gpu_stream)

      !   The fictitious energy, formed from the same kernel matrix the dipole needs
      call gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites)
      call gpu_kernels_pow(kernels_d, kernels_copy_d, zeta, size_kernels, gpu_stream)
      call gpu_blas_mvmul_n(cublas_handle, kernels_copy_d, alphas_d, energies_d, n_sites, n_sparse)
      call gpu_axpe(energies_d, cdelta_ene, 0.d0, size_energies, gpu_stream)
      call cpy_dtoh(energies_d, c_loc(energies), st_energies, gpu_stream)

      !   Qss, exactly as for the forces but for the sign
      call gpu_malloc_async(kernels_der_d, st_kernels, gpu_stream)
      call gpu_malloc_async(Qss_d, st_soap, gpu_stream)
      call gpu_malloc_async(Qs_copy_d, st_Qs, gpu_stream)
      call cpy_dtod(Qs_d, Qs_copy_d, st_Qs, gpu_stream)

      mzetam = zeta - 1
      call gpu_kernels_pow(kernels_d, kernels_der_d, mzetam, size_kernels, gpu_stream)

      if (n_sites < n_soap) then
         call gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse, gpu_stream)
      else
         call gpu_matvect(Qs_copy_d, alphas_d, n_soap, n_sparse, gpu_stream)
      end if

      !   Plus, not minus: a dipole is a gradient, where a force is its negative
      cdelta_dipole = zeta*delta**2
      call gpu_blas_mmul_n_t(cublas_handle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, &
                             n_soap, n_sites, cdelta_dipole)

      !   Site -> its own self pair. This is neighbors_beg, and it is the only
      !   pair index the dipole needs, so the kernel runs over n_sites blocks
      !   rather than n_pairs.
      allocate (beg_index(1:n_sites))
      l = 0
      do i = 1, n_sites
         beg_index(i) = l + 1
         do j = 1, n_neigh(i)
            l = l + 1
         end do
      end do
      st_beg_index = n_sites*sizeof(beg_index(1))
      call gpu_malloc_async(beg_index_d, st_beg_index, gpu_stream)
      call cpy_htod(c_loc(beg_index), beg_index_d, st_beg_index, gpu_stream)

      call gpu_malloc_async(dipoles_d, st_dipoles, gpu_stream)
      call gpu_soap_dipole(n_sites, Qss_d, n_soap, beg_index_d, soap_der_d, dipoles_d, gpu_stream)
      call cpy_dtoh(dipoles_d, c_loc(dipoles), st_dipoles, gpu_stream)

      call gpu_free_async(dipoles_d, gpu_stream)
      call gpu_free_async(beg_index_d, gpu_stream)
      call gpu_free_async(kernels_der_d, gpu_stream)
      call gpu_free_async(Qss_d, gpu_stream)
      call gpu_free_async(Qs_copy_d, gpu_stream)
      call gpu_free_async(kernels_d, gpu_stream)
      call gpu_free_async(kernels_copy_d, gpu_stream)
      call gpu_free_async(energies_d, gpu_stream)

      deallocate (kernels, beg_index)

   end subroutine get_soap_dipole

#else

!  The sparse set for the descriptor currently being evaluated. It is POINTED
!  at, never copied: Qs feeds a dgemm and is n_sparse x dim, so a copy per
!  descriptor per step would be real time. soap_backend_begin associates them
!  and soap_backend_end nullifies.
!
!  They keep the names the dummy arguments had, so the body of
!  get_soap_energy_and_forces is untouched -- the same trick gap_backend_gpu
!  uses for its ten buffers. Private, because these are common names and
!  gap.f90 has no default private, so publishing them collides on use.
   real(dp), pointer :: alphas(:) => null()
   real(dp), pointer :: Qs(:, :) => null()
   private :: alphas, Qs

contains

   subroutine soap_backend_begin(hypers)
!     Point at one descriptor's sparse set. The GPU branch's version of this
!     uploads it to the device instead; the driver calls the same name on both.
      use types, only: soap_turbo
      implicit none
      type(soap_turbo), intent(in), target :: hypers

      alphas => hypers%alphas
      Qs => hypers%Qs
   end subroutine soap_backend_begin

   subroutine soap_backend_end()
      implicit none
      nullify (alphas)
      nullify (Qs)
   end subroutine soap_backend_end

   subroutine get_soap_energy_and_forces(soap, soap_der, delta, zeta0, e0, &
                                         n_neigh, neighbors_list, xyz, do_forces, do_timing, &
                                         energies, forces, virial)
!   soap(1:n_soap, 1:n_sites)

      implicit none

      real(dp), intent(in) :: soap(:, :)
      real(dp), intent(in) :: soap_der(:, :, :)
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: e0
      real(dp), intent(in) :: zeta0
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(out) :: energies(:)
      real(dp), intent(out) :: forces(:, :)
      real(dp), intent(out) :: virial(1:3, 1:3)
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: neighbors_list(:)
      logical, intent(in) :: do_forces
      logical, intent(in) :: do_timing
      real(dp), allocatable :: kernels(:, :)
      real(dp), allocatable :: kernels_der(:, :)
      real(dp), allocatable :: Qss(:, :)
      real(dp), allocatable :: Qs_copy(:, :)
      real(dp), allocatable :: this_Qss(:)
      real(dp), allocatable :: kernels_copy(:, :)
      real(dp) :: time1
      real(dp) :: time2
      real(dp) :: time3
      real(dp) :: energies_time
      real(dp) :: forces_time
      real(dp) :: zeta
      real(dp) :: this_force(1:3)
      integer :: n_sites
      integer :: n_sparse
      integer :: n_soap
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: j2
      integer :: zeta_int
      integer :: n_sites0
      integer :: k1
      integer :: k2
      logical :: is_zeta_int = .false.
!    integer, allocatable :: neighbors_beg(:), neighbors_end(:)

      if (dabs(zeta0 - dfloat(int(zeta0))) < 1.d-5) then
         is_zeta_int = .true.
         zeta_int = int(zeta0)
         zeta = dfloat(zeta_int)
      else
         zeta = zeta0
      end if

!   Energies
      if (do_timing) then
         call cpu_time(time1)
         time3 = time1
      end if
      n_sparse = size(alphas)
      n_soap = size(soap, 1)
      n_sites = size(soap, 2)
      n_sites0 = size(forces, 2)

      allocate (kernels(1:n_sites, 1:n_sparse))
      kernels = 0.d0
      allocate (kernels_copy(1:n_sites, 1:n_sparse))
      if (n_sites > 0) then
         call dgemm("t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
                    kernels, n_sites)
      end if
!   We copy the kernels because it makes the matmul() operation (WHICH SHOULD BY THE WAY BE WRITTEN
!   USING LAPACK ROUTINES) a lot faster
      if (is_zeta_int) then
         kernels_copy = kernels**zeta_int
         energies = matmul(kernels_copy, alphas)
      else
         kernels_copy = kernels**zeta
         energies = matmul(kernels_copy, alphas)
      end if
      energies = delta**2*energies + e0
      if (do_timing) then
         call cpu_time(time2)
         energies_time = time2 - time1
      end if

!   Forces
      if (do_forces) then
         if (do_timing) then
            call cpu_time(time1)
         end if
         allocate (kernels_der(1:n_sites, 1:n_sparse))
         allocate (Qss(1:n_sites, 1:n_soap))
         Qss = 0.d0
         allocate (Qs_copy(1:n_soap, 1:n_sparse))
         allocate (this_Qss(1:n_soap))

         Qs_copy = Qs

         if (is_zeta_int) then
            kernels_der = kernels**(zeta_int - 1)
         else
            kernels_der = kernels**(zeta - 1.d0)
         end if
         if (n_sites < n_soap) then
            do i = 1, n_sites
               kernels_der(i, :) = kernels_der(i, :)*alphas(:)
            end do
         else
            do i = 1, n_soap
               Qs_copy(i, :) = Qs(i, :)*alphas(:)
            end do
         end if

         if (n_sites > 0) then
            call dgemm("n", "t", n_sites, n_soap, n_sparse, -zeta*delta**2, kernels_der, n_sites, &
                       Qs_copy, n_soap, 0.d0, Qss, n_sites)
         end if

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
            this_Qss = Qss(i, 1:n_soap)
            do j = 1, n_neigh(i)
               l = l + 1
               j2 = mod(neighbors_list(l) - 1, n_sites0) + 1
               do k = 1, 3
                  this_force(k) = dot_product(this_Qss, soap_der(k, :, l))
                  forces(k, j2) = forces(k, j2) + this_force(k)
               end do
!         This is a many body potential, so there's no factor of 1/2 here
!          virial = virial + dot_product( this_force(1:3), xyz(1:3,l) )
               do k1 = 1, 3
                  do k2 = 1, 3
                     virial(k1, k2) = virial(k1, k2) + 0.5d0*(this_force(k1)*xyz(k2, l) + this_force(k2)*xyz(k1, l))
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

         if (do_timing) then
            call cpu_time(time2)
            forces_time = time2 - time1
         end if
      end if

!   Wrap it up
      deallocate (kernels, kernels_copy)
      if (do_forces) then
         deallocate (kernels_der, Qs_copy, Qss, this_Qss)
      end if

      if (do_timing) then
         call cpu_time(time2)
         write (*, *) '                                       |'
         write (*, *) 'Prediction timings (SOAP):             |'
         write (*, *) '                                       |'
         write (*, '(A, F7.3, A)') '  *) Energy prediction: ', energies_time, ' seconds |'
         if (do_forces) then
            write (*, '(A, F7.3, A)') '  *) Forces prediction: ', forces_time, ' seconds |'
         end if
         write (*, '(A, F8.3, A)') '  *) Total prediction: ', time2 - time3, ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if

   end subroutine

   subroutine get_soap_dipole(soap, soap_cart_der, delta, zeta0, n_neigh, do_timing, dipoles, energies)
!   Local dipoles from a dipole GAP.
!
!   The model is a GAP whose local "energy" E_i is a fictitious scalar, fitted so
!   that its derivative with respect to the central atom's OWN position is the
!   local dipole:
!
!       mu_i = dE_i/dr_i   and   mu = sum_i mu_i
!
!   This is the same quantity get_soap_energy_and_forces already contracts, with
!   two differences: the sign (mu is a gradient, not a force, so the prefactor is
!   +zeta*delta**2 rather than -zeta*delta**2), and the fact that only the SELF
!   pair of each site is wanted rather than the scatter over all its neighbors.
!   soap_turbo builds that self pair as -sum_{j/=1} of the neighbor terms, which
!   is exactly d(soap_i)/d(r_i), and the neighbor list puts it first for every
!   site (see build_neighbors_list).
!
!   E_i itself is returned too, because the kernel matrix it needs has already
!   been formed here. It is a fitting artefact with no physical meaning: it is
!   reported as energy_dipole and must never be added to the total energy, and
!   its gradient must never be added to the forces.
!
!   soap(1:n_soap, 1:n_sites), dipoles(1:3, 1:n_sites)

      implicit none

      real(dp), intent(in) :: soap(:, :)
      real(dp), intent(in) :: soap_cart_der(:, :, :)
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: zeta0
      integer, intent(in) :: n_neigh(:)
      logical, intent(in) :: do_timing
      real(dp), intent(out) :: dipoles(:, :)
      real(dp), intent(out) :: energies(:)
      real(dp), allocatable :: kernels(:, :)
      real(dp), allocatable :: kernels_copy(:, :)
      real(dp), allocatable :: kernels_der(:, :)
      real(dp), allocatable :: Qss(:, :)
      real(dp), allocatable :: Qs_copy(:, :)
      real(dp) :: time1
      real(dp) :: time2
      real(dp) :: zeta
      integer :: n_sites
      integer :: n_sparse
      integer :: n_soap
      integer :: i
      integer :: k
      integer :: l
      integer :: zeta_int
      logical :: is_zeta_int = .false.

      if (do_timing) then
         call cpu_time(time1)
      end if

      if (dabs(zeta0 - dfloat(int(zeta0))) < 1.d-5) then
         is_zeta_int = .true.
         zeta_int = int(zeta0)
         zeta = dfloat(zeta_int)
      else
         zeta = zeta0
      end if

      n_sparse = size(alphas)
      n_soap = size(soap, 1)
      n_sites = size(soap, 2)

      dipoles = 0.d0
      energies = 0.d0
      if (n_sites == 0) then
         return
      end if

      allocate (kernels(1:n_sites, 1:n_sparse))
      kernels = 0.d0
      call dgemm("t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
                 kernels, n_sites)

!   The fictitious energy, taken before the alphas are folded into kernels_der
!   below. No e0: a dipole model has no reference energy to speak of.
      allocate (kernels_copy(1:n_sites, 1:n_sparse))
      if (is_zeta_int) then
         kernels_copy = kernels**zeta_int
      else
         kernels_copy = kernels**zeta
      end if
      energies = delta**2*matmul(kernels_copy, alphas)
      deallocate (kernels_copy)

      allocate (kernels_der(1:n_sites, 1:n_sparse))
      allocate (Qss(1:n_sites, 1:n_soap))
      Qss = 0.d0
      allocate (Qs_copy(1:n_soap, 1:n_sparse))
      Qs_copy = Qs

      if (is_zeta_int) then
         kernels_der = kernels**(zeta_int - 1)
      else
         kernels_der = kernels**(zeta - 1.d0)
      end if
!   Fold the alphas into whichever of the two factors is cheaper, exactly as
!   get_soap_energy_and_forces does
      if (n_sites < n_soap) then
         do i = 1, n_sites
            kernels_der(i, :) = kernels_der(i, :)*alphas(:)
         end do
      else
         do i = 1, n_soap
            Qs_copy(i, :) = Qs(i, :)*alphas(:)
         end do
      end if

      call dgemm("n", "t", n_sites, n_soap, n_sparse, zeta*delta**2, kernels_der, n_sites, &
                 Qs_copy, n_soap, 0.d0, Qss, n_sites)

      l = 0
      do i = 1, n_sites
!       j = 1 is the central atom itself, so soap_cart_der(:,:,l) here is
!       d(soap_i)/d(r_i). The remaining n_neigh(i)-1 pairs are skipped.
         l = l + 1
         do k = 1, 3
            dipoles(k, i) = dot_product(Qss(i, 1:n_soap), soap_cart_der(k, :, l))
         end do
         l = l + n_neigh(i) - 1
      end do

      deallocate (kernels, kernels_der, Qs_copy, Qss)

      if (do_timing) then
         call cpu_time(time2)
         write (*, *) '                                       |'
         write (*, '(A, F7.3, A)') '  *) Dipole prediction: ', time2 - time1, ' seconds |'
         write (*, *) '.......................................|'
      end if

   end subroutine

   subroutine get_soap_dipole_weights(soap, soap_cart_der, delta, zeta0, n_neigh, w, V)
!   The model-side half of a dipole gradient.
!
!   mu_i = w_i . dq/dr_i with w_i = dE_i/dq, so
!
!     d mu_ia / d r_jb = sum_de (d2E_i/dq_d dq_e)(dq_d/dr_ia)(dq_e/dr_jb)   [A]
!                      + sum_d (dE_i/dq_d)(d2 q_d/dr_ia dr_jb)              [B]
!
!   This routine returns what both terms need. w is term B's weight vector, to
!   be handed to get_soap_central_hessian as its vecs argument. V is term A,
!   pre-contracted.
!
!   Term A is the one that has to be written carefully. For the dot_product
!   covariance,
!
!     d2E_i/dq_d dq_e = zeta(zeta-1) delta^2 sum_s alpha_s k_s^(zeta-2) q_sd q_se
!
!   so term A is sum_s beta_s (q_s . dq/dr_ia)(q_s . dq/dr_jb). Evaluated as
!   written that needs the sparse-space projection of EVERY pair's derivative,
!   which is n_sparse times the cost of the force contraction and is what makes
!   a naive dipole gradient unaffordable. But the first factor depends only on
!   the site, so folding the sparse sum into
!
!     V_a(i) = sum_s beta_s (q_s . dq/dr_ia) q_s
!
!   leaves term A = V_a . dq/dr_jb: three descriptor-space vectors per site,
!   built once, and 9*n_soap per pair afterwards -- the same order as the
!   forces, and linear in the number of pairs.
!
!   zeta = 1 has no term A at all (the kernel is linear in q); it returns V = 0
!   rather than going through k^(-1).
!
!   soap(1:n_soap, 1:n_sites), w(1:n_soap, 1:n_sites), V(1:n_soap, 1:3, 1:n_sites)

      implicit none

      real(dp), intent(in) :: soap(:, :)
      real(dp), intent(in) :: soap_cart_der(:, :, :)
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: zeta0
      integer, intent(in) :: n_neigh(:)
      real(dp), intent(out) :: w(:, :)
      real(dp), intent(out) :: V(:, :, :)
      real(dp), allocatable :: kernels(:, :)
      real(dp) :: zeta
      real(dp) :: beta
      real(dp) :: gsa(1:3)
      integer :: n_sites
      integer :: n_sparse
      integer :: n_soap
      integer :: i
      integer :: s
      integer :: k
      integer :: a
      integer :: zeta_int
      logical :: is_zeta_int

      is_zeta_int = .false.
      if (dabs(zeta0 - dfloat(int(zeta0))) < 1.d-5) then
         is_zeta_int = .true.
         zeta_int = int(zeta0)
         zeta = dfloat(zeta_int)
      else
         zeta = zeta0
      end if

      n_sparse = size(alphas)
      n_soap = size(soap, 1)
      n_sites = size(soap, 2)

      w = 0.d0
      V = 0.d0
      if (n_sites == 0) then
         return
      end if

      allocate (kernels(1:n_sites, 1:n_sparse))
      kernels = 0.d0
      call dgemm("t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
                 kernels, n_sites)

      do i = 1, n_sites
         do s = 1, n_sparse
            if (is_zeta_int) then
               beta = zeta*delta**2*alphas(s)*kernels(i, s)**(zeta_int - 1)
            else
               beta = zeta*delta**2*alphas(s)*kernels(i, s)**(zeta - 1.d0)
            end if
            w(1:n_soap, i) = w(1:n_soap, i) + beta*Qs(1:n_soap, s)
         end do
      end do

      if (dabs(zeta - 1.d0) < 1.d-10) then
         deallocate (kernels)
         return
      end if

!   V_a(i). j = 1 is the self pair, so soap_cart_der(:,:,k) there is dq/dr_i.
      k = 0
      do i = 1, n_sites
         k = k + 1
         do s = 1, n_sparse
            if (is_zeta_int) then
               beta = zeta*(zeta - 1.d0)*delta**2*alphas(s)*kernels(i, s)**(zeta_int - 2)
            else
               beta = zeta*(zeta - 1.d0)*delta**2*alphas(s)*kernels(i, s)**(zeta - 2.d0)
            end if
            do a = 1, 3
               gsa(a) = dot_product(Qs(1:n_soap, s), soap_cart_der(a, 1:n_soap, k))
            end do
            do a = 1, 3
               V(1:n_soap, a, i) = V(1:n_soap, a, i) + (beta*gsa(a))*Qs(1:n_soap, s)
            end do
         end do
         k = k + n_neigh(i) - 1
      end do

      deallocate (kernels)

   end subroutine

   subroutine assemble_soap_dipole_gradient(V, soap_cart_der, hess, n_neigh, dipole_der)
!   The two halves put together.
!
!     dipole_der(a, b, k2) = d mu_ia / d r_jb   for pair k2 = (site i, neighbour j)
!
!   with V from get_soap_dipole_weights and hess from get_soap_central_hessian
!   called with vecs = w from the same routine. For k2 the self pair the block
!   returned is d mu_ia / d r_ib, the central atom's own.
!
!   The caller's sequence is:
!
!       soap_hessian_enabled = .true.
!       call get_soap(..., soap, soap_cart_der)
!       call get_soap_dipole_weights(soap, soap_cart_der, delta, zeta, n_neigh, w, V)
!       call get_soap_central_hessian(..., soap, 1, w, hess)
!       call assemble_soap_dipole_gradient(V, soap_cart_der, hess, n_neigh, dipole_der)
!
!   and it needs soap_radial_legacy_filter = .false., because the radial second
!   derivatives do.
!
!   A rigid translation cannot change a dipole, so the blocks of one site sum to
!   zero for each (a,b). That holds term by term here, which makes it a cheap
!   check on the caller's pair indexing rather than on the physics.

      implicit none

      real(dp), intent(in) :: V(:, :, :)
      real(dp), intent(in) :: soap_cart_der(:, :, :)
      real(dp), intent(in) :: hess(:, :, :, :)
      integer, intent(in) :: n_neigh(:)
      real(dp), intent(out) :: dipole_der(:, :, :)
      integer :: n_sites
      integer :: n_soap
      integer :: i
      integer :: j
      integer :: k
      integer :: a
      integer :: b

      n_soap = size(V, 1)
      n_sites = size(V, 3)

      k = 0
      do i = 1, n_sites
         do j = 1, n_neigh(i)
            k = k + 1
            do b = 1, 3
               do a = 1, 3
                  dipole_der(a, b, k) = dot_product(V(1:n_soap, a, i), soap_cart_der(b, 1:n_soap, k)) &
                                        + hess(a, b, 1, k)
               end do
            end do
         end do
      end do

   end subroutine

   subroutine accumulate_dmu_dr(V, soap_cart_der, hess, n_neigh, neighbors_list, n_sites0, dmu_dr)
!   d(total dipole)/d(atom position), accumulated per atom.
!
!     dmu_dr(a, b, j) = d mu_a / d r_jb,     mu = sum_i mu_i
!
!   This is the form a MAD IR force wants. The per-pair blocks are never
!   stored: (3,3,n_atom_pairs) is 52 MB for 7000 atoms and grows with the
!   system, whereas the scatter target is 9*n_atoms -- half a megabyte at the
!   same size, and it does not care how many neighbours each atom has.
!
!   The force then costs one contraction with a single 3-vector,
!   f_jb = -sum_a lambda_a dmu_dr(a,b,j), which is why lambda does not have to
!   be known before the descriptor pass. Contracting it in earlier would make
!   the Hessian pass about three times cheaper, but it would force either a
!   second get_soap pass per step or a one-step-lagged lambda, and neither is
!   worth it against half a megabyte.
!
!   Accumulates rather than assigns, so several dipole descriptors, or several
!   batches of sites, can add into the same array. Zero it before the first
!   call.
!
!   V and hess come from get_soap_dipole_weights and get_soap_central_hessian.

      implicit none

      real(dp), intent(in) :: V(:, :, :)
      real(dp), intent(in) :: soap_cart_der(:, :, :)
      real(dp), intent(in) :: hess(:, :, :, :)
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: neighbors_list(:)
      integer, intent(in) :: n_sites0
      real(dp), intent(inout) :: dmu_dr(:, :, :)
      integer :: n_sites
      integer :: n_soap
      integer :: i
      integer :: j
      integer :: k
      integer :: a
      integer :: b
      integer :: j2

      n_soap = size(V, 1)
      n_sites = size(V, 3)

      k = 0
      do i = 1, n_sites
         do j = 1, n_neigh(i)
            k = k + 1
!           the atom this pair points at, folded back into the primitive cell
            j2 = mod(neighbors_list(k) - 1, n_sites0) + 1
            do b = 1, 3
               do a = 1, 3
                  dmu_dr(a, b, j2) = dmu_dr(a, b, j2) &
                                     + dot_product(V(1:n_soap, a, i), soap_cart_der(b, 1:n_soap, k)) &
                                     + hess(a, b, 1, k)
               end do
            end do
         end do
      end do

   end subroutine

   subroutine accumulate_mu_weighted_force(wsite, V, soap_cart_der, hess, n_neigh, &
                                           neighbors_list, n_sites0, in_to_out_site, force)
!   The MAD IR force when the weight is already known: a per-site 3-vector
!   contracted against the local dipole gradient as the descriptor is walked.
!
!     force(b, j) += sum_i sum_a wsite(a, i) d mu_ia / d r_jb
!
!   This is what an EXTENDED-LAGRANGIAN bias wants (mad_ir_xl.f90, equation 3),
!   and it is a different shape of calculation from accumulate_dmu_dr's even
!   though it walks the same pairs. There, the weight lambda is a function of
!   the dipole of the current configuration, so it cannot exist until the
!   descriptor pass is over, and the pass has to leave behind the whole
!   (3,3,n_atoms) tensor for lambda to meet later. Here the weight is
!   w_i = sum_k g x_{k,i}, a dynamical variable carried in from the previous
!   step, so it is known BEFORE the pass and is contracted inside it. Nothing
!   of size 9*n_atoms is ever formed, and there is 3*n_atoms rather than
!   9*n_atoms to reduce across ranks.
!
!   It is also cheaper by about a factor of three in the term that dominates.
!   Contracting the weight into V first,
!
!     Vw(:, i) = sum_a wsite(a, i) V(:, a, i)
!
!   turns the three dot products of length n_soap that accumulate_dmu_dr does
!   per (pair, direction) into one. The Hessian term does not compress the same
!   way -- hess carries both indices -- but it is a 3x3 contraction rather than
!   an n_soap one.
!
!   THE WEIGHT IS INDEXED BY GLOBAL ATOM, THE DESCRIPTOR BY BATCH SITE. wsite
!   is (3, n_atoms) in the caller's numbering and in_to_out_site maps this
!   batch's local site index to it, exactly as the dipole accumulation above
!   does for local_dipoles0. Getting this wrong gives a force that is the right
!   size and attached to the wrong atoms.
!
!   Accumulates rather than assigns, so several dipole descriptors, or several
!   batches, add into the same array. force is the run's force array and is
!   already non-zero: this ADDS the bias to it.

      implicit none

      real(dp), intent(in) :: wsite(:, :)
      real(dp), intent(in) :: V(:, :, :)
      real(dp), intent(in) :: soap_cart_der(:, :, :)
      real(dp), intent(in) :: hess(:, :, :, :)
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: neighbors_list(:)
      integer, intent(in) :: n_sites0
      integer, intent(in) :: in_to_out_site(:)
      real(dp), intent(inout) :: force(:, :)
      real(dp), allocatable :: Vw(:, :)
      real(dp) :: wi(1:3), acc
      integer :: n_sites, n_soap, i, i2, j, k, a, b, j2

      n_soap = size(V, 1)
      n_sites = size(V, 3)

      allocate (Vw(1:n_soap, 1:n_sites))
      do i = 1, n_sites
         i2 = in_to_out_site(i)
         Vw(1:n_soap, i) = wsite(1, i2)*V(1:n_soap, 1, i) &
                           + wsite(2, i2)*V(1:n_soap, 2, i) &
                           + wsite(3, i2)*V(1:n_soap, 3, i)
      end do

      k = 0
      do i = 1, n_sites
         i2 = in_to_out_site(i)
         wi(1:3) = wsite(1:3, i2)
         do j = 1, n_neigh(i)
            k = k + 1
!           the atom this pair points at, folded back into the primitive cell
            j2 = mod(neighbors_list(k) - 1, n_sites0) + 1
            j2 = in_to_out_site(j2)
            do b = 1, 3
               acc = dot_product(Vw(1:n_soap, i), soap_cart_der(b, 1:n_soap, k))
               do a = 1, 3
                  acc = acc + wi(a)*hess(a, b, 1, k)
               end do
               force(b, j2) = force(b, j2) + acc
            end do
         end do
      end do

      deallocate (Vw)

   end subroutine

#endif
   subroutine get_2b_energy_and_forces(rjs, xyz, alphas, cutoff, rcut, buffer, delta, sigma, e0, Qs, &
                                       n_neigh, do_forces, do_timing, species, neighbor_species, &
                                       species1, species2, species_types, energies, forces, virial)

      implicit none

      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(in) :: alphas(:)
      real(dp), intent(in) :: cutoff(:)
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: sigma
      real(dp), intent(in) :: e0
      real(dp), intent(in) :: Qs(:)
      real(dp), intent(in) :: rcut
      real(dp), intent(in) :: buffer
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: species(:)
      integer, intent(in) :: neighbor_species(:)
      character*8, intent(in) :: species_types(:)
      character*8, intent(in) :: species1
      character*8, intent(in) :: species2
      logical, intent(in) :: do_forces
      logical, intent(in) :: do_timing

      real(dp), intent(out) :: energies(:)
      real(dp), intent(out) :: forces(:, :)
      real(dp), intent(out) :: virial(1:3, 1:3)

      real(dp) :: time1
      real(dp) :: time2
      real(dp) :: fcut
      real(dp) :: pi
      real(dp) :: dfcut
      real(dp) :: this_force(1:3)
      integer :: n_sparse
      integer :: i
      integer :: j
      integer :: k
      integer :: n_sites
      integer :: n_atom_pairs
      integer :: s
      integer :: sp1
      integer :: sp2
      integer :: n_sites0
      integer :: k1
      integer :: k2

      if (do_timing) then
         call cpu_time(time1)
      end if

      pi = dacos(-1.d0)

      n_sparse = size(alphas)
      n_sites = size(n_neigh)
      n_atom_pairs = size(rjs)
      n_sites0 = size(forces, 2)

!   Map species to index
      do i = 1, size(species_types)
         if (species1 == species_types(i)) then
            sp1 = i
            exit
         end if
      end do
      do i = 1, size(species_types)
         if (species2 == species_types(i)) then
            sp2 = i
            exit
         end if
      end do

!   Energy calculation
      energies = e0
      k = 0
      do i = 1, n_sites
         if (species(i) /= sp1 .and. species(i) /= sp2) then
            k = k + n_neigh(i)
            cycle
         end if
         k = k + 1
         do j = 2, n_neigh(i)
            k = k + 1
            if ((species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
                (species(i) == sp2 .and. neighbor_species(k) == sp1)) then
               continue
            else
               cycle
            end if
            if (rjs(k) < rcut) then
               if (rjs(k) < rcut - buffer) then
                  fcut = 1.d0
               else
                  fcut = (dcos(pi*(rjs(k) - rcut + buffer)/buffer) + 1.d0)/2.d0
               end if
               do s = 1, n_sparse
                  energies(i) = energies(i) + delta**2*alphas(s)*cutoff(s)*fcut* &
                                dexp(-0.5d0*(rjs(k) - Qs(s))**2/sigma**2)
               end do
            end if
         end do
      end do

!   Force calculation
      if (do_forces) then
         forces = 0.d0
         virial = 0.d0
         k = 0
         do i = 1, n_sites
            if (species(i) /= sp1 .and. species(i) /= sp2) then
               k = k + n_neigh(i)
               cycle
            end if
            k = k + 1
            do j = 2, n_neigh(i)
               k = k + 1
               if ((species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
                   (species(i) == sp2 .and. neighbor_species(k) == sp1)) then
                  continue
               else
                  cycle
               end if
               if (rjs(k) < rcut) then
                  if (rjs(k) < rcut - buffer) then
                     fcut = 1.d0
                     dfcut = 0.d0
                  else
                     fcut = (dcos(pi*(rjs(k) - rcut + buffer)/buffer) + 1.d0)/2.d0
                     dfcut = pi/2.d0/buffer*dsin(pi*(rjs(k) - rcut + buffer)/buffer)
                  end if
                  do s = 1, n_sparse
                     this_force(1:3) = -2.d0*delta**2*alphas(s)*cutoff(s)* &
                                       dexp(-0.5d0*(rjs(k) - Qs(s))**2/sigma**2)* &
                                       xyz(1:3, k)/rjs(k)*((rjs(k) - Qs(s))/sigma**2*fcut + dfcut)
                     forces(1:3, i) = forces(1:3, i) + this_force(1:3)
                     do k1 = 1, 3
                        do k2 = 1, 3
                           virial(k1, k2) = virial(k1, k2) - 0.5d0*(this_force(k1)*xyz(k2, k) + this_force(k2)*xyz(k1, k))
                        end do
                     end do
                  end do
               end if
            end do
         end do
!     Half contribution to the virial for pair potentials
         virial = 0.5d0*virial
      end if

      if (do_timing) then
         call cpu_time(time2)
         write (*, *) '                                       |'
         write (*, *) 'Prediction timings (2b):               |'
         write (*, *) '                                       |'
         write (*, '(A, F8.3, A)') '  *) Total prediction: ', time2 - time1, ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if

   end subroutine

   subroutine get_core_pot_energy_and_forces(rjs, xyz, x, V, yp1, ypn, dVdx2, n_neigh, do_forces, do_timing, species, &
                                             neighbor_species, species1, species2, species_types, energies, forces, &
                                             virial)

      implicit none

      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: V(:)
      real(dp), intent(in) :: yp1
      real(dp), intent(in) :: ypn
      real(dp), intent(in) :: dVdx2(:)
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: species(:)
      integer, intent(in) :: neighbor_species(:)
      character*8, intent(in) :: species_types(:)
      character*8, intent(in) :: species1
      character*8, intent(in) :: species2
      logical, intent(in) :: do_forces
      logical, intent(in) :: do_timing

      real(dp), intent(out) :: energies(:)
      real(dp), intent(out) :: forces(:, :)
      real(dp), intent(out) :: virial(1:3, 1:3)

!   There are two ways of doing the core_pot interpolation; most efficient probably depends on
!   whether a small subset or big subset of the total number of atom pairs has a core potential
!   term associated to it. The current implementation is fast going over pairs, but slow computing
!   the splines. This will give the best performance for systems with many atom types
      real(dp) :: V_int(1:1)
      real(dp) :: dV_int(1:1)
      real(dp) :: this_force(1:3)
!    real(dp), allocatable :: V_int(:), dV_int(:)
      real(dp) :: time1
      real(dp) :: time2
      real(dp) :: rcut
      integer :: n_sparse
      integer :: i
      integer :: j
      integer :: k
      integer :: n_sites
      integer :: n_atom_pairs
      integer :: s
      integer :: sp1
      integer :: sp2
      integer :: n_sites0
      integer :: k1
      integer :: k2

      if (do_timing) then
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
         if (species1 == species_types(i)) then
            sp1 = i
            exit
         end if
      end do
      do i = 1, size(species_types)
         if (species2 == species_types(i)) then
            sp2 = i
            exit
         end if
      end do

!   Energy calculation
      energies = 0.d0
      k = 0
      do i = 1, n_sites
         if (species(i) /= sp1 .and. species(i) /= sp2) then
            k = k + n_neigh(i)
            cycle
         end if
         k = k + 1
         do j = 2, n_neigh(i)
            k = k + 1
            if ((species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
                (species(i) == sp2 .and. neighbor_species(k) == sp1)) then
               continue
            else
               cycle
            end if
            if (rjs(k) < rcut) then
               V_int = spline(x, V, dVdx2, yp1, ypn, rjs(k:k), rcut)
               energies(i) = energies(i) + 0.5d0*V_int(1)
            end if
         end do
      end do

!   Force calculation
      if (do_forces) then
         forces = 0.d0
         virial = 0.d0
         k = 0
         do i = 1, n_sites
            if (species(i) /= sp1 .and. species(i) /= sp2) then
               k = k + n_neigh(i)
               cycle
            end if
            k = k + 1
            do j = 2, n_neigh(i)
               k = k + 1
               if ((species(i) == sp1 .and. neighbor_species(k) == sp2) .or. &
                   (species(i) == sp2 .and. neighbor_species(k) == sp1)) then
                  continue
               else
                  cycle
               end if
               if (rjs(k) < rcut) then
                  dV_int = spline_der(x, V, dVdx2, yp1, ypn, rjs(k:k), rcut)
                  this_force(1:3) = dV_int(1)*xyz(1:3, k)/rjs(k)
                  forces(1:3, i) = forces(1:3, i) + this_force(1:3)
                  do k1 = 1, 3
                     do k2 = 1, 3
                        virial(k1, k2) = virial(k1, k2) - 0.5d0*(this_force(k1)*xyz(k2, k) + this_force(k2)*xyz(k1, k))
                     end do
                  end do
               end if
            end do
         end do
!     Half contribution to the virial for pair potentials
         virial = 0.5d0*virial
      end if

      if (do_timing) then
         call cpu_time(time2)
         write (*, *) '                                       |'
         write (*, *) 'Prediction timings (core_pot):         |'
         write (*, *) '                                       |'
         write (*, '(A, F8.3, A)') '  *) Total prediction: ', time2 - time1, ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if

   end subroutine

#ifdef _GPU
   subroutine setup_3b_gpu(kernel_type, species_center, species1, species2, species_types, c_name, sp0, sp1, sp2)
      implicit none
      character(kind=c_char, len=4) :: c_name
      character*3, intent(in) :: kernel_type
      character*8, intent(in) :: species_center
      character*8, intent(in) :: species1
      character*8, intent(in) :: species2
      character*8, intent(in) :: species_types(:)
      integer, intent(out) :: sp0
      integer, intent(out) :: sp1
      integer, intent(out) :: sp2
      integer :: i

      c_name = trim(kernel_type)//c_null_char
      !Map species to index
      do i = 1, size(species_types)
         if (species_center == species_types(i)) then
            sp0 = i
            exit
         end if
      end do
      do i = 1, size(species_types)
         if (species1 == species_types(i)) then
            sp1 = i
            exit
         end if
      end do
      do i = 1, size(species_types)
         if (species2 == species_types(i)) then
            sp2 = i
            exit
         end if
      end do

   end subroutine
#endif

   subroutine get_3b_energy_and_forces(rjs, xyz, alphas, cutoff, rcut, buffer, delta, sigma, e0, Qs, &
                                       n_neigh, neighbors_list, do_forces, do_timing, kernel_type, &
                                       species, neighbor_species, species_center, species1, species2, &
                                       species_types, energies, forces, virial)

      implicit none

      real(dp), intent(in) :: rjs(:)
      real(dp), intent(in) :: xyz(:, :)
      real(dp), intent(in) :: alphas(:)
      real(dp), intent(in) :: cutoff(:)
      real(dp), intent(in) :: rcut
      real(dp), intent(in) :: buffer
      real(dp), intent(in) :: delta
      real(dp), intent(in) :: sigma(:)
      real(dp), intent(in) :: e0
      real(dp), intent(in) :: Qs(:, :)
      integer, intent(in) :: n_neigh(:)
      integer, intent(in) :: neighbors_list(:)
      integer, intent(in) :: species(:)
      integer, intent(in) :: neighbor_species(:)
      logical, intent(in) :: do_timing
      logical, intent(in) :: do_forces
      character*3, intent(in) :: kernel_type
      character*8, intent(in) :: species_center
      character*8, intent(in) :: species1
      character*8, intent(in) :: species2
      character*8, intent(in) :: species_types(:)

      real(dp), intent(out) :: energies(:)
      real(dp), intent(out) :: forces(:, :)
      real(dp), intent(out) :: virial(1:3, 1:3)

      real(dp) :: time1
      real(dp) :: time2
      real(dp) :: fcut
      real(dp) :: pi
      real(dp) :: r12
      real(dp) :: r13
      real(dp) :: r23
      real(dp) :: xyz12(1:3)
      real(dp) :: xyz13(1:3)
      real(dp) :: xyz23(1:3)
      real(dp) :: q(1:3)
      real(dp) :: fcut12
      real(dp) :: fcut13
      real(dp) :: dfcut12(1:3)
      real(dp) :: dfcut13(1:3)
      real(dp) :: force1(1:3)
      real(dp) :: force2(1:3)
      real(dp) :: force3(1:3)
      real(dp) :: dfcut(1:3)
      real(dp) :: xyz12_red(1:3)
      real(dp) :: xyz13_red(1:3)
      real(dp) :: xyz23_red(1:3)
      real(dp) :: this_force(1:3)
      real(dp), allocatable :: r(:)
      real(dp), allocatable :: drdq(:, :)
      real(dp), allocatable :: kernel(:)
      real(dp), allocatable :: drdx1(:, :)
      real(dp), allocatable :: drdx2(:, :)
      real(dp), allocatable :: drdx3(:, :)
      real(dp), allocatable :: pref(:)
      real(dp), allocatable :: kernel_der(:)
      integer :: n_sparse
      integer :: i
      integer :: j
      integer :: k
      integer :: k2
      integer :: n_sites
      integer :: n_atom_pairs
      integer :: s
      integer :: j2
      integer :: i3
      integer :: j3
      integer :: k3
      integer :: l
      integer :: sp0
      integer :: sp1
      integer :: sp2
      integer :: n_sites0
      integer :: k1
      integer :: k4

      if (do_timing) then
         call cpu_time(time1)
      end if

      pi = dacos(-1.d0)

      n_sparse = size(alphas)
      n_sites = size(n_neigh)
      n_atom_pairs = size(rjs)
      n_sites0 = size(forces, 2)

!   Map species to index
      do i = 1, size(species_types)
         if (species_center == species_types(i)) then
            sp0 = i
            exit
         end if
      end do
      do i = 1, size(species_types)
         if (species1 == species_types(i)) then
            sp1 = i
            exit
         end if
      end do
      do i = 1, size(species_types)
         if (species2 == species_types(i)) then
            sp2 = i
            exit
         end if
      end do

      allocate (r(1:n_sparse))
      allocate (kernel(1:n_sparse))
      allocate (pref(1:n_sparse))
      if (do_forces) then
         allocate (drdq(1:n_sparse, 1:3))
         allocate (drdx1(1:n_sparse, 1:3))
         allocate (drdx2(1:n_sparse, 1:3))
         allocate (drdx3(1:n_sparse, 1:3))
         if (kernel_type == "pol") then
            allocate (kernel_der(1:n_sparse))
         end if
      end if

!   NOTE: these loops can be made between 2 and 3 times faster by summing over triplets instead of
!   over sites, however it may be a better strategy with parallelization in mind to sum over sites
!   (or maybe not...)

      pref(1:n_sparse) = 2.d0*delta**2*alphas(1:n_sparse)*cutoff(1:n_sparse)

!   Energy and force calculation
      energies = e0
      if (do_forces) then
         forces = 0.d0
         virial = 0.d0
      end if
      k = 0
      do i = 1, n_sites
         if (species(i) /= sp0) then
            k = k + n_neigh(i)
            cycle
         end if
         k = k + 1
!      i3 = neighbors_list(k)
         i3 = mod(neighbors_list(k) - 1, n_sites0) + 1
         do j = 2, n_neigh(i)
            k = k + 1
            r12 = rjs(k)
            if ((neighbor_species(k) /= sp1 .and. neighbor_species(k) /= sp2) .or. r12 > rcut) then
               cycle
            end if
!        j3 = neighbors_list(k)
            j3 = mod(neighbors_list(k) - 1, n_sites0) + 1
            xyz12 = xyz(1:3, k)
            xyz12_red = xyz12/r12
            k2 = k
            do j2 = j + 1, n_neigh(i)
               k2 = k2 + 1
               r13 = rjs(k2)
               if (((neighbor_species(k) == sp1 .and. neighbor_species(k2) == sp2) .or. &
                    (neighbor_species(k) == sp2 .and. neighbor_species(k2) == sp1)) .and. r13 < rcut) then
                  continue
               else
                  cycle
               end if
!          k3 = neighbors_list(k2)
               k3 = mod(neighbors_list(k2) - 1, n_sites0) + 1
               xyz13 = xyz(1:3, k2)
               xyz13_red = xyz13/r13
               xyz23 = xyz13 - xyz12
               r23 = dsqrt(sum(xyz23**2))
               xyz23_red = xyz23/r23
!         It would be nice that if r23 < rcut, we only evaluate the expression if k3 > i3 and j3 > i3,
!         however that gets messy because the cutoff functions are not the same for each of the 3
!         evaluations of the triplet, i.e., atom 1 sees two cutoff functions, atom 2 sees two *different*
!         cutoff functions and 3 sees another 2 differetn cutoff functions. I actually tried and it
!         reduced the calculation time by about 40% for the large random C system with a 3 Angstrom cutoff
               if (r12 < rcut .and. r13 < rcut) then
                  if (r12 < rcut - buffer) then
                     fcut12 = 1.d0
                     dfcut12(1:3) = 0.d0
                  else
                     fcut12 = (dcos(pi*(r12 - rcut + buffer)/buffer) + 1.d0)/2.d0
                     dfcut12(1:3) = dsin(pi*(r12 - rcut + buffer)/buffer)/2.d0*pi/buffer*xyz12_red(1:3)
                  end if
                  if (r13 < rcut - buffer) then
                     fcut13 = 1.d0
                     dfcut13(1:3) = 0.d0
                  else
                     fcut13 = (dcos(pi*(r13 - rcut + buffer)/buffer) + 1.d0)/2.d0
                     dfcut13(1:3) = dsin(pi*(r13 - rcut + buffer)/buffer)/2.d0*pi/buffer*xyz13_red(1:3)
                  end if
                  fcut = fcut12*fcut13
                  dfcut(1:3) = fcut12*dfcut13(1:3) + fcut13*dfcut12(1:3)
!           This builds the actual descriptor
                  q = [r12 + r13, (r12 - r13)**2, r23]
!           This gets the Euclidean distance between q and all the qs
                  r(1:n_sparse) = (q(1) - Qs(1:n_sparse, 1))**2/sigma(1)**2
                  r(1:n_sparse) = r(1:n_sparse) + (q(2) - Qs(1:n_sparse, 2))**2/sigma(2)**2
                  r(1:n_sparse) = r(1:n_sparse) + (q(3) - Qs(1:n_sparse, 3))**2/sigma(3)**2
                  r(1:n_sparse) = dsqrt(r(1:n_sparse))
!           This gets the derivatives of the distances wrt the descriptors (the 1/r factor is not included)
                  if (do_forces) then
                     drdq(1:n_sparse, 1) = (q(1) - Qs(1:n_sparse, 1))/sigma(1)**2
                     drdq(1:n_sparse, 2) = (q(2) - Qs(1:n_sparse, 2))/sigma(2)**2
                     drdq(1:n_sparse, 3) = (q(3) - Qs(1:n_sparse, 3))/sigma(3)**2
                  end if
!           Evaluate the kernels
                  if (kernel_type == "exp") then
!             This kernel already contains the prefactor
                     kernel(1:n_sparse) = pref(1:n_sparse)*dexp(-0.5d0*r(1:n_sparse)**2)
                     energies(i) = energies(i) + fcut*sum(kernel(1:n_sparse))
                     if (do_forces) then
                        force1 = 0.d0
                        force2 = 0.d0
                        force3 = 0.d0
                        do l = 1, 3
!                 For atom 1
                           drdx1(1:n_sparse, l) = drdq(1:n_sparse, 1)*(xyz12_red(l) + xyz13_red(l)) &
                                                  + drdq(1:n_sparse, 2)*2.d0*(r12 - r13)*(xyz12_red(l) - xyz13_red(l))
!                 For atom 2
                           drdx2(1:n_sparse, l) = -drdq(1:n_sparse, 1)*xyz12_red(l) &
                                                  - drdq(1:n_sparse, 2)*2.d0*(r12 - r13)*xyz12_red(l) &
                                                  + drdq(1:n_sparse, 3)*xyz23_red(l)
!                 For atom 3
                           drdx3(1:n_sparse, l) = -drdq(1:n_sparse, 1)*xyz13_red(l) &
                                                  - drdq(1:n_sparse, 2)*2.d0*(r13 - r12)*xyz13_red(l) &
                                                  - drdq(1:n_sparse, 3)*xyz23_red(l)
                           force1(l) = -sum(kernel(1:n_sparse)*(dfcut(l) + fcut*drdx1(1:n_sparse, l)))
                           force2(l) = -sum(kernel(1:n_sparse)*(-fcut13*dfcut12(l) + fcut*drdx2(1:n_sparse, l)))
                           force3(l) = -sum(kernel(1:n_sparse)*(-fcut12*dfcut13(l) + fcut*drdx3(1:n_sparse, l)))
                        end do
                        forces(1:3, i3) = forces(1:3, i3) + force1(1:3)
                        forces(1:3, j3) = forces(1:3, j3) + force2(1:3)
                        forces(1:3, k3) = forces(1:3, k3) + force3(1:3)
!               force1 acting on i3 does not contribute to the virial
!                virial = virial + dot_product( force2(1:3), xyz12(1:3) )
!                virial = virial + dot_product( force3(1:3), xyz13(1:3) )
                        do k1 = 1, 3
                           do k4 = 1, 3
                              virial(k1, k4) = virial(k1, k4) + 0.5d0*(force2(k1)*xyz12(k4) + force2(k4)*xyz12(k1))
                              virial(k1, k4) = virial(k1, k4) + 0.5d0*(force3(k1)*xyz13(k4) + force3(k4)*xyz13(k1))
                           end do
                        end do
                     end if
                  else if (kernel_type == "pol") then
!             This kernel already contains the prefactor
                     kernel(1:n_sparse) = pref(1:n_sparse)*cov_pp(r(1:n_sparse), 3, 1)
                     energies(i) = energies(i) + fcut*sum(kernel(1:n_sparse))
                     if (do_forces) then
                        force1 = 0.d0
                        force2 = 0.d0
                        force3 = 0.d0
!               This derivative contains the 1/r term
                        kernel_der(1:n_sparse) = pref(1:n_sparse)*cov_pp_der(r(1:n_sparse), 3, 1)
                        do l = 1, 3
!                 For atom 1
                           drdx1(1:n_sparse, l) = drdq(1:n_sparse, 1)*(xyz12_red(l) + xyz13_red(l)) &
                                                  + drdq(1:n_sparse, 2)*2.d0*(r12 - r13)*(xyz12_red(l) - xyz13_red(l))
!                 For atom 2
                           drdx2(1:n_sparse, l) = -drdq(1:n_sparse, 1)*xyz12_red(l) &
                                                  - drdq(1:n_sparse, 2)*2.d0*(r12 - r13)*xyz12_red(l) &
                                                  + drdq(1:n_sparse, 3)*xyz23_red(l)
!                 For atom 3
                           drdx3(1:n_sparse, l) = -drdq(1:n_sparse, 1)*xyz13_red(l) &
                                                  - drdq(1:n_sparse, 2)*2.d0*(r13 - r12)*xyz13_red(l) &
                                                  - drdq(1:n_sparse, 3)*xyz23_red(l)
                           force1(l) = -sum(kernel(1:n_sparse)*dfcut(l) &
                                            - kernel_der(1:n_sparse)*fcut*drdx1(1:n_sparse, l))
                           force2(l) = -sum(-kernel(1:n_sparse)*fcut13*dfcut12(l) &
                                            - kernel_der(1:n_sparse)*fcut*drdx2(1:n_sparse, l))
                           force3(l) = -sum(-kernel(1:n_sparse)*fcut12*dfcut13(l) &
                                            - kernel_der(1:n_sparse)*fcut*drdx3(1:n_sparse, l))
                        end do
                        forces(1:3, i3) = forces(1:3, i3) + force1(1:3)
                        forces(1:3, j3) = forces(1:3, j3) + force2(1:3)
                        forces(1:3, k3) = forces(1:3, k3) + force3(1:3)
!               force1 acting on i3 does not contribute to the virial
!                virial = virial + dot_product( force2(1:3), xyz12(1:3) )
!                virial = virial + dot_product( force3(1:3), xyz13(1:3) )
                        do k1 = 1, 3
                           do k4 = 1, 3
                              virial(k1, k4) = virial(k1, k4) + 0.5d0*(force2(k1)*xyz12(k4) + force2(k4)*xyz12(k1))
                              virial(k1, k4) = virial(k1, k4) + 0.5d0*(force3(k1)*xyz13(k4) + force3(k4)*xyz13(k1))
                           end do
                        end do
                     end if
                  end if
               end if
            end do
         end do
      end do

      deallocate (kernel, pref, r)
      if (do_forces) then
         deallocate (drdq, drdx1, drdx2, drdx3)
         if (kernel_type == "pol") then
            deallocate (kernel_der)
         end if
      end if

      if (do_timing) then
         call cpu_time(time2)
         write (*, *) '                                       |'
         write (*, *) 'Prediction timings (3b):               |'
         write (*, *) '                                       |'
         write (*, '(A, F8.3, A)') '  *) Total prediction: ', time2 - time1, ' seconds |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if

   end subroutine

!
! This is two functions
   function cov_pp(r, d, q) result(cov)

      implicit none

      real(dp), intent(in) :: r(:)
      integer, intent(in) :: d
      integer, intent(in) :: q
      real(dp), dimension(1:size(r)) :: cov

      real(dp) :: j
      integer :: j_int
      integer :: i

      j_int = d/2 + q + 1
      j = dfloat(j_int)

      if (d == 3 .and. q == 1) then
         do i = 1, size(r)
            if (r(i) >= 1.d0) then
               cov(i) = 0.d0
            else if (r(i) <= 0.d0) then
               cov(i) = 1.d0
            else
               cov(i) = (1.d0 - r(i))**(j_int + 1)*((j + 1.d0)*r(i) + 1.d0)
            end if
         end do
      else
         write (*, *) "ERROR: This combination of input parameters is not currently supported for cov_pp!"
         stop
      end if

   end function
   function cov_pp_der(r, d, q) result(cov_der)

      implicit none

      real(dp), intent(in) :: r(:)
      integer, intent(in) :: d
      integer, intent(in) :: q
      real(dp), dimension(1:size(r)) :: cov_der

      real(dp) :: j
      integer :: j_int
      integer :: i

      j_int = d/2 + q + 1
      j = dfloat(j_int)

      if (d == 3 .and. q == 1) then
         do i = 1, size(r)
            if (r(i) >= 1.d0 .or. r(i) <= 0.d0) then
               cov_der(i) = 0.d0
            else
!         This expression contains the 1/r factor
               cov_der(i) = (-(j + 1.d0)*(1.d0 - r(i))**j_int*((j + 1.d0)*r(i) + 1.d0) + &
                             (1.d0 - r(i))**(j_int + 1)*(j + 1.d0))/r(i)
            end if
         end do
      else
         write (*, *) "ERROR: This combination of input parameters is not currently supported for cov_pp!"
         stop
      end if

   end function

end module gap
