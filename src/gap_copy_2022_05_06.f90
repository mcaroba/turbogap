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

! module load gcc/10.3.0 openblas openmpi cuda
! rm cuda_wrappers.o ; nvcc -arch=sm_70 -c src/cuda_wrappers.cu ; make clean; make -B
!  srun  --time=00:10:00 --partition=gputest --account=project_2000634 --nodes=1 --ntasks-per-node=4  --cpus-per-task=1 --gres=gpu:a100:4  ../turbogap_dev/bin/turbogap md
! gnupot
!plot "< paste trajectory_out.xyz_64k_oiriginal trajectory_out.xyz | awk 'NR>2{print $5,$13}'", x
MODULE F_B_C
    INTERFACE
      subroutine gpu_malloc_double(a_d,n) bind(C,name="cuda_malloc_double")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d
        integer(c_int),value :: n
      end subroutine

      subroutine gpu_malloc_int(a_d,n) bind(C,name="cuda_malloc_int")
        use iso_c_binding
        implicit none
          type(c_ptr) :: a_d
          integer(c_int),value :: n
      end subroutine

      subroutine gpu_free(a_d) bind(C,name="cuda_free")
        use iso_c_binding
        implicit none
        type(c_ptr) :: a_d
      end subroutine

      subroutine cpy_double_htod(a,a_d,n) bind(C,name="cuda_cpy_double_htod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine cpy_int_htod(a,a_d,n) bind(C,name="cuda_cpy_int_htod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine cpy_double_dtod(b_d,a_d,n) bind(C,name="cuda_cpy_double_dtod")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,b_d
        integer(c_int),value :: n
      end subroutine


      subroutine cpy_double_dtoh(a_d,a,n) bind(C,name="cuda_cpy_double_dtoh")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d,a
        integer(c_int),value :: n
      end subroutine

      subroutine gpu_vector_fill_curand(a_d,n,c) bind(C,name="GPU_fill_rand")
        use iso_c_binding
        implicit none
        type(c_ptr),value :: a_d
        integer(c_int),value :: n,c
      end subroutine

      subroutine create_cublas_handle(cubhandle)bind(C,name="create_cublas_handle")
        use iso_c_binding
        implicit none
        type(c_ptr) :: cubhandle
      end subroutine

      subroutine destroy_cublas_handle(cubhandle)bind(C,name="destroy_cublas_handle")
        use iso_c_binding
        implicit none
        type(c_ptr) :: cubhandle
      end subroutine

      subroutine gpu_blas_mmul_t_n(cubhandle, Qs_d, soap_d, kernels_d, &
                         n_sparse, n_soap, n_sites)  bind(C,name="gpu_blas_mmul_t_n")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cubhandle,kernels_d,Qs_d,soap_d
        integer(c_int),value :: n_sparse, n_soap, n_sites
      end subroutine

      subroutine gpu_blas_mmul_n_t(cubhandle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, &
                                  n_soap, n_sites, cdelta) bind(C,name="gpu_blas_mmul_n_t")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cubhandle,kernels_der_d, Qs_copy_d, Qss_d
        integer(c_int),value :: n_sparse,n_soap, n_sites
        real(c_double), value :: cdelta
      end subroutine

      subroutine gpu_blas_mvmul_n(cubhandle, A, B, C, nAx, nAy) bind(C,name="gpu_blas_mvmul_n")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: cubhandle,A,B,C
        integer(c_int),value :: nAx,nAy
      end subroutine

      subroutine gpu_kernels_pow( A,B,  zeta, N) bind(C,name="gpu_kernels_pow")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: A,B
        integer(c_int),value :: N
        real(c_double), value :: zeta
      end subroutine


      subroutine gpu_axpe( A,  dccc,e0, N) bind(C,name="gpu_axpc")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: A
        integer(c_int),value :: N
        real(c_double), value :: dccc,e0
      end subroutine

      subroutine one_for_all(soap, kernels, kernels_copy, Qs, energies, &
            delta, zeta, e0, &
            n_sites, n_soap, n_sparse, &
            size_kernels, size_soap, size_Qs, size_alphas, size_energies) bind(C, name="wrappers_all")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: soap, kernels,kernels_copy, Qs, energies
        real(c_double), value :: delta, zeta, e0
        integer(c_int), value :: n_sites, n_soap, n_sparse
        integer(c_int), value :: size_kernels, size_soap, size_Qs, size_alphas, size_energies
      end subroutine

      subroutine gpu_set_device(my_rank) bind(C, name="cuda_set_device")
        use iso_c_binding
        integer(c_int), value :: my_rank
      end subroutine

      subroutine gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse) bind(C,name="cuda_matvect_kernels")
        use iso_c_binding
        type(c_ptr), value :: kernels_der_d, alphas_d
        integer(c_int), value :: n_sites,n_sparse
      end subroutine

      subroutine gpu_final_soap_forces_virial(n_neigh_d,n_sites,maxnn, &
                                              Qss_d,n_soap,neighbors_beg_d, &
                                              soap_der_d, &
                                              xyz_d, virial_d, &
                                              neighbors_list_d, n_sites0, forces_d) &
                  bind(C,name="gpu_final_soap_forces_virial")
        use iso_c_binding
        type(c_ptr), value :: Qss_d, n_neigh_d,neighbors_beg_d, neighbors_list_d
        type(c_ptr), value :: forces_d, soap_der_d, xyz_d, virial_d
        integer(c_int), value :: n_sites,n_soap,maxnn, n_sites0
      end subroutine

      subroutine gpu_soap_energies_forces_virial(n_neigh_d,n_sites,maxnn, &
                                              Qss_d,n_soap,neighbors_beg_d, &
                                              soap_der_d, &
                                              xyz_d, virial_d, &
                                              neighbors_list_d, n_sites0, forces_d, &
                                              cuhandle, kernels_der_d, Qs_copy_d, &
                                              n_sparse, cdelta, &
                                              alphas_d, &
                                              kernels_d, mzetam, size_kernels, &
                                              do_forces ) &
                  bind(C,name="gpu_soap_energies_forces_virial")
        use iso_c_binding
        type(c_ptr), value :: Qss_d, n_neigh_d,neighbors_beg_d, neighbors_list_d
        type(c_ptr), value :: forces_d, soap_der_d, xyz_d, virial_d
        integer(c_int), value :: n_sites,n_soap,maxnn, n_sites0
        type(c_ptr), value :: cuhandle,kernels_der_d, Qs_copy_d, alphas_d, kernels_d
        integer(c_int),value :: n_sparse, size_kernels
        logical, value:: do_forces
        real(c_double), value :: cdelta, mzetam
      end subroutine

    END INTERFACE
  END MODULE F_B_C

module gap

  use splines
  use F_B_C
  use iso_c_binding
  use mpi

  contains

  subroutine get_soap_energy_and_forces(soap, soap_der, alphas, delta, zeta0, e0, Qs, &
                                        n_neigh, neighbors_list, xyz, do_forces, do_timing, &
                                        energies, forces, virial, solo_time_soap)
!   **********************************************
!   soap(1:n_soap, 1:n_sites)

    !use mpi
    implicit none

    real(c_double), intent(in),target :: soap(:,:), soap_der(:,:,:), alphas(:), delta, Qs(:,:), e0, zeta0, xyz(:,:)
    real(c_double), intent(out), target :: energies(:), forces(:,:), virial(1:3,1:3)
    integer(c_int), intent(in), target :: n_neigh(:), neighbors_list(:)
    logical, intent(in) :: do_forces, do_timing
    real(c_double), allocatable,target :: kernels(:,:), kernels_der(:,:), &
                            Qss(:,:), Qs_copy(:,:), this_Qss(:), &
                           kernels_copy(:,:)
    real(c_double) :: time1, time2, time3, energies_time, forces_time, this_force(1:3)
    real(c_double) ::  zeta, cdelta, mzetam
    integer(c_int) :: n_sites, n_sparse, n_soap, i, j, k, l, j2, zeta_int, n_sites0, k1, k2
    logical :: is_zeta_int = .false.
    type(c_ptr) :: cublas_handle, kernels_copy_d, kernels_d, soap_d, Qs_d, energies_d, alphas_d
    type(c_ptr) :: kernels_der_d, Qss_d, Qs_copy_d, this_Qss_d
    integer(c_int) :: size_kernels, size_soap, size_Qs, size_alphas, size_energies,maxnn
    integer(c_int) :: size_nnlist, size_xyz, n1xyz, n2xyz, n1forces,n2forces, n1virial,n2virial
    integer(c_int) :: size_forces, size_virial, size_soap_der,n1soap_der,n2soap_der,n3soap_der
    integer(c_int) :: rank, ierr
    integer(c_int), allocatable, target :: neighbors_beg(:), neighbors_end(:)
    type(c_ptr) :: soap_der_d, virial_d, n_neigh_d
    type(c_ptr) :: neighbors_beg_d, neighbors_end_d, xyz_d,  neighbors_list_d, forces_d
    real*8 :: solo_time_soap(1:3)

    call cpu_time(solo_time_soap(1))
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call gpu_set_device(mod(rank,4)) ! Every node has 4 GPUs. Even if there are more than 1 nodes used. This will assing th ranks to GPU in a roundbin fashion
    cdelta=delta*delta
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

    !if(rank==0) then
    !write(*,*)
    !write(*,*) "n_sparse= ", n_sparse
    !write(*,*) "n_soap= ", n_soap
    !write(*,*) "n_sites= ", n_sites
    !write(*,*) "n_sites0= ", n_sites0
    !write(*,*)
    !endif
    !stop
    allocate( kernels(1:n_sites, 1:n_sparse) )
    kernels = 0.d0
    allocate( kernels_copy(1:n_sites, 1:n_sparse) )

    size_kernels=n_sites*n_sparse
    size_soap=n_soap*n_sites
    size_Qs = n_soap*n_sparse
    size_alphas=n_sparse
    size_energies=n_sites
    !write(*,*)
    !write(*,*)  "n_sites", n_sites, "n_soap", n_soap, "n_sparse", n_sparse
    !write(*,*) size_kernels,size_soap,size_Qs,size_alphas,size_energies
    call create_cublas_handle(cublas_handle)
    call gpu_malloc_double(kernels_d,size_kernels)
    call gpu_malloc_double(kernels_copy_d,size_kernels)
    call gpu_malloc_double(soap_d,size_soap)
    call gpu_malloc_double(Qs_d,size_Qs)

    call cpy_double_htod(c_loc(soap),soap_d,size_soap)
    call cpy_double_htod(c_loc(Qs),   Qs_d ,size_Qs)
        !  gpu_blas_mmul_t_n (cubhandle, Qs_d, soap_d, kernels_d,  n_sparse, n_soap, n_sites)
    call gpu_blas_mmul_t_n(cublas_handle, Qs_d, soap_d, kernels_d, n_sparse, n_soap, n_sites) ! Here the result is in variables kernels_d
    !call cpy_double_dtoh(kernels_d,c_loc(kernels),size_kernels)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !stop
    ! soap(n_soap,n_sites) nAx=n_sites, nAy=n_soap
    ! Qs(1:n_soap, 1:n_sparse) nBx=nsoap, nBy=n_sparse
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !call dgemm( "t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
    !          kernels, n_sites)
    !call one_for_all(c_loc(soap), c_loc(kernels),c_loc(kernels_copy), c_loc(Qs), c_loc(energies), &
    !  delta, zeta, e0, &
    !  n_sites, n_soap,n_sparse, &
    !  size_kernels, size_soap, size_Qs, size_alphas, size_energies)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   We copy the kernels because it makes the matmul() operation (WHICH SHOULD BY THE WAY BE WRITTEN
!   USING LAPACK ROUTINES) a lot faster
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if( is_zeta_int )then
!      kernels_copy = kernels**zeta_int
!    else
!      kernels_copy = kernels**zeta
!    end if

    call gpu_kernels_pow( kernels_d, kernels_copy_d,zeta, size_kernels)
    !call cpy_double_dtoh(kernels_copy_d,c_loc(kernels_copy),size_kernels)


    !write(*,*) zeta, zeta_int, is_zeta_int
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kernels_copy(1:n_sites, 1:n_sparse)
    ! alphas(n_sparse)
    ! energies (n_sites)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gpu_malloc_double(alphas_d,size_alphas)
    call gpu_malloc_double(energies_d,size_energies)
    call cpy_double_htod(c_loc(alphas),alphas_d,size_alphas)

    call gpu_blas_mvmul_n(cublas_handle, kernels_copy_d, alphas_d, energies_d, n_sites, n_sparse)
!!    call gpu_axpe( energies_d,cdelta,e0, size_energies)
!!    call cpy_double_dtoh(energies_d,c_loc(energies),size_energies)
    !energies = matmul(kernels_copy, alphas)
    !energies = delta**2 * energies + e0

    if( do_timing )then
      call cpu_time(time2)
      energies_time = time2 - time1
    end if



!   Forces
!    if( do_forces )then
      if( do_timing )then
        call cpu_time(time1)
      end if
!!      allocate( kernels_der(1:n_sites, 1:n_sparse) )
!!      allocate( Qss(1:n_sites, 1:n_soap) )
!!      Qss = 0.d0
!!      allocate( Qs_copy(1:n_soap, 1:n_sparse) )
!!      allocate(this_Qss(1:n_soap))
!!      Qs_copy = Qs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    size_kernels=n_sites*n_sparse
!    size_soap=n_soap*n_sites
!    size_Qs = n_soap*n_sparse
!    size_alphas=n_sparse
!    size_energies=n_sites
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call gpu_malloc_double(kernels_der_d,size_kernels)
      call gpu_malloc_double(Qss_d,size_soap)
      call gpu_malloc_double(Qs_copy_d,size_Qs)
      call gpu_malloc_double(this_Qss_d,n_soap)

      call cpy_double_dtod(Qs_d,   Qs_copy_d ,size_Qs)

      !if( is_zeta_int )then
      !  kernels_der = kernels**(zeta_int-1)
      !else
      !  kernels_der = kernels**(zeta-1.d0)
      !end if
      mzetam=zeta-1
!!      call gpu_kernels_pow( kernels_d, kernels_der_d,mzetam, size_kernels)
      !call cpy_double_dtoh(kernels_der_d,c_loc(kernels_der),size_kernels)

!!      if( n_sites < n_soap )then
!!        call gpu_matvect(kernels_der_d, alphas_d, n_sites, n_sparse)
!!        !call cpy_double_dtoh(kernels_der_d,c_loc(kernels_der),size_kernels)
!!        !allocate( kernels_der(1:n_sites, 1:n_sparse) )
!!        !do i = 1, n_sites
!!        !  kernels_der(i,:) = kernels_der(i,:)*alphas(:)
!!        !end do
!!      else
!!
!!        call gpu_matvect(Qs_copy_d, alphas_d, n_soap, n_sparse)
!!        !call cpy_double_dtoh(Qs_copy_d,c_loc(Qs_copy),size_Qs)
!!        !allocate( Qs_copy(1:n_soap, 1:n_sparse) )!!
!!
!!        !do i = 1, n_soap
!!        !  Qs_copy(i,:) = Qs(i,:)*alphas(:)
!!        !end do
!!      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       cdelta=-zeta*delta**2
      !call dgemm("n", "t", n_sites, n_soap, n_sparse, cdelta, kernels_der, n_sites, &
      !           Qs_copy, n_soap, 0.d0, Qss, n_sites)
      !call gpu_blas_mmul_n_t(cublas_handle, kernels_der_d, Qs_copy_d, Qss_d, n_sparse, n_soap, n_sites, cdelta)
      !call cpy_double_dtoh( Qss_d,c_loc(Qss),size_soap)

! EXPERIMENTAL CODE
      call cpu_time(time1)
      allocate( neighbors_beg(1:n_sites) )
      allocate( neighbors_end(1:n_sites) )
      l = 0
      maxnn=0
      do i = 1, n_sites
        neighbors_beg(i) = l + 1
        do j = 1, n_neigh(i)
          l = l + 1
        end do
        neighbors_end(i) = l
        if(n_neigh(i)>maxnn) then
          maxnn=n_neigh(i)
        end if
      end do
! END EXPERIMENTAL CODE
    !write(*,*) size(n_neigh), size(n_neigh,1)
    !stop
     virial = 0.d0
     forces = 0.d0

     n1xyz=size(xyz,1)
     n2xyz=size(xyz,2)
     size_xyz=n1xyz*n2xyz

     call gpu_malloc_double(xyz_d,size_xyz)

     size_nnlist=size(neighbors_list,1)

    n1forces = size(forces, 1)
    n2forces = size(forces, 2)
    size_forces=n1forces*n2forces

    n1virial = size(virial, 1)
    n2virial = size(virial, 2)
    size_virial=n1virial*n2virial

    n1soap_der = size(soap_der, 1)
    n2soap_der = size(soap_der, 2)
    n3soap_der = size(soap_der, 3)
    size_soap_der=n1soap_der*n2soap_der*n3soap_der
    call gpu_malloc_int(neighbors_list_d,size_nnlist)
    call gpu_malloc_double(forces_d,size_forces)
    call gpu_malloc_double(virial_d,size_virial)
    call gpu_malloc_int(n_neigh_d,n_sites)
    call gpu_malloc_int(neighbors_end_d,n_sites)
    call gpu_malloc_int(neighbors_beg_d,n_sites)
    call gpu_malloc_double(soap_der_d,size_soap_der)
    call cpy_double_htod(c_loc( xyz), xyz_d,size_xyz)
    call cpy_int_htod(c_loc( n_neigh), n_neigh_d,n_sites)
    call cpy_int_htod(c_loc( neighbors_beg), neighbors_beg_d,n_sites)
    call cpy_int_htod(c_loc( neighbors_end), neighbors_end_d,n_sites)
    call cpy_int_htod(c_loc(neighbors_list), neighbors_list_d, size_nnlist)
    call cpy_double_htod(c_loc(soap_der), soap_der_d, size_soap_der)
!    write(*,*) "n_sites ", n_sites, "n_sites0 ", n_sites0, "maxnn ", maxnn, "size_Qs ", size_Qs, &
!     "size_soap_der ", size_soap_der, "n3soap_der ",n3soap_der, "n2soap_der ",n2soap_der, "n1soap_der ",n1soap_der, &
!     "size_nnlist ", size_nnlist, "n_soap" , n_soap, "size_forces ", size_forces, "size_soap", size_soap
!    stop

     call gpu_soap_energies_forces_virial(n_neigh_d,n_sites, maxnn, &
                                      Qss_d,n_soap, neighbors_beg_d, &
                                      soap_der_d,  &
                                      xyz_d, virial_d, &
                                      neighbors_list_d, n_sites0, forces_d, &
                                      cublas_handle, kernels_der_d, Qs_copy_d, &
                                      n_sparse, cdelta,&
                                      alphas_d, &
                                      kernels_d, mzetam, size_kernels, &
                                      do_forces)

      call cpy_double_dtoh(forces_d,c_loc(forces), size_forces)
      call cpy_double_dtoh(virial_d,c_loc(virial), size_virial)


!!$OMP parallel do private(i,j,l,j2,this_Qss)



 !!!!!      l = 0
 !      do i = 1, n_sites
 !        this_Qss = Qss(i,1:n_soap)
 !!!!!        do j = 1, n_neigh(i)
 !!!!!          l = l + 1
 !!!!!         do l = neighbors_beg(i), neighbors_end(i) !!
 !          do l= neighbors_beg(i), neighbors_end(i) !,neighbors_beg(i)-1+n_neigh(i)
!
 !          j2 = mod(neighbors_list(l)-1, n_sites0) + 1
 !          do k = 1, 3
 !              this_force(k) = dot_product(this_Qss, soap_der(k,:,l))
 !            forces(k, j2) = forces(k, j2) + this_force(k)
 !          end do
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!         This is a many body potential, so there's no factor of 1/2 here
 !!!          virial = virial + dot_product( this_force(1:3), xyz(1:3,l) )
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 !          do k1 = 1, 3
 !            do k2 =1, 3
 !              virial(k1, k2) = virial(k1, k2) + 0.5d0 * (this_force(k1)*xyz(k2,l) + this_force(k2)*xyz(k1,l))
 !            end do
 !          end do
 !        end do
 !      end do
 !!!$OMP end parallel do

      if( do_timing )then
        call cpu_time(time2)
        forces_time = time2 - time1
      end if
!    end if


! EXPERIMENTAL CODE
      deallocate( neighbors_beg, neighbors_end )
!      call cpu_time(time2)
!      write(*,*) time2-time1
! END EXPERIMENTAL CODE


!   Wrap it up
    deallocate(kernels, kernels_copy)
!    if( do_forces )then
!      deallocate(kernels_der, Qs_copy, Qss, this_Qss)

      call gpu_free(neighbors_list_d)
      call gpu_free(neighbors_end_d)
      call gpu_free(neighbors_beg_d)
      call gpu_free(kernels_der_d)
      call gpu_free(Qss_d)
      call gpu_free(Qs_copy_d)
      call gpu_free(forces_d)
      call gpu_free(neighbors_list_d)
      call gpu_free(virial_d)
!    end if


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

    call cpy_double_dtoh(energies_d,c_loc(energies),size_energies)
  call destroy_cublas_handle(cublas_handle)
  call gpu_free(kernels_d)
  call gpu_free(kernels_d)
  call gpu_free(kernels_copy_d)
  call gpu_free(soap_d)
  call gpu_free(Qs_d)
  call gpu_free(energies_d)
  call gpu_free(alphas_d)
  call gpu_free(this_Qss_d)
  call cpu_time(solo_time_soap(2))
  solo_time_soap(3)=solo_time_soap(3)+solo_time_soap(2)-solo_time_soap(1)
  !stop
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
