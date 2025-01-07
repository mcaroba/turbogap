! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, get_linearized_sparse_model.f90, is copyright (c) 2019-2024, Miguel A. Caro and
! HND X   Patricia Hernandez-Leon
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

! This program checks whether kernel linearization is faster than usual kernel
! computation, given a trained GAP with soap_turbo descriptors.
! If favorable, the program also generates any file related to the linearized
! model, needed for TurboGAP calculations with kernel linearization support.

program get_linear_model

  use kernel_linearization

  implicit none

  integer :: n_soap, n_linear, n_sparse, n_soap_turbo = 0, zeta, iostatus, i, j, n
  integer :: n_sites = 1000 ! n. atoms used for timings tests
  character*1024 :: file_gap, filename, temp
  character*64 :: keyword
  character*128 :: cjunk
  character*1024, allocatable :: file_alphas(:), file_desc(:), file_gap_updated(:)
  real*8, allocatable :: zetas(:), alphas(:), alphas_linear(:)
  real*8, allocatable :: Qs(:,:), Ms_linear(:,:)
  logical :: do_check = .true., do_linear = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters that should be provided by the user
! path to the GAP file
  call get_command_argument(1, file_gap)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Extract the info about zeta, alphas file and descriptor files
  call read_gap_hypers(file_gap, n_soap_turbo, zetas, file_alphas, file_desc)

  do n = 1, n_soap_turbo
    write(*,*) n
!   Check that zeta is integer, otherwise continue
    if( dabs(zetas(n) - dfloat(int(zetas(n)))) < 1.d-5 )then
      zeta = int(zetas(n))
    else
      continue
    end if
!   Get alphas, sparse descriptors and their size
!   Here size of Qs is (1:n_soap, 1:n_sparse)
    call read_alphas_and_descriptors(file_desc(n), file_alphas(n), n_sparse, n_soap, alphas, Qs)
    write(*,*) file_desc(n)
!   Compute the size of the linearized model
    n_linear = 1
    do i = 1, zeta
      n_linear = n_linear*(n_soap + i - 1)/i
    end do
    write(*,*) "n_sites = ", n_sites
    write(*,*) "zeta = ", zeta
    write(*,*) "n_sparse = ", size(Qs, 1)
    write(*,*) "n_soap = ", size(Qs, 2)
    write(*,*) "n_linear = ", n_linear

!   Check whether kernel linearization is a better option
    if( do_check )then
      do_linear = check_kernel_vs_linear(zeta, n_sites, n_sparse, n_soap, n_linear, .true.)
    else
      do_linear = .true.
    end if

!   Compute the linearized sparse contributions if favourable
    if( do_linear )then
!     ENERGIES
      alphas_linear = get_linearized_sparse(zeta, alphas, transpose(Qs))
!     save them
      filename = trim(file_alphas(n))
      filename = filename(1:29) // "_linear.dat"
      open(unit=10, file=filename, status="unknown")
      do i = 1, n_linear
        write(10,*) alphas_linear(i), 1.
      end do
      close(10)

!     FORCES
      allocate( Ms_linear(1:n_soap,1:n_linear) )
      Ms_linear = 1.d0
!      Ms_linear = get_linearized_sparse_forces(zeta, alphas, Qs)
!     add the file for forces contributions
      filename = trim(file_desc(n)) // "_linear"
      open(unit=10, file=filename, status="unknown")
      do i = 1, n_linear
        do j = 1, n_soap
          write(10,*) Ms_linear(j, i)
        end do
      end do
      close(10)

    end if

  end do


contains
!**************************************************************************
!
! This subroutine is a much shorter version of src/read_files.f90 subroutine
! with same name, modified to only read the information of soap_turbo descriptors.
! It only looks for the zetas as well as filenames of alphas, descriptors
!
  subroutine read_gap_hypers(file_gap, n_soap_turbo, zeta, file_alphas, file_desc)

    implicit none

!   Input variables
    character(len=*), intent(in) :: file_gap
!   Output variables
    real*8, allocatable, intent(out) :: zeta(:)
    character(len=*), allocatable, intent(out) :: file_alphas(:), file_desc(:)
    integer, intent(out) :: n_soap_turbo
!   Internal variables
    integer :: iostatus
    character*64 :: keyword, cjunk
    character*1 :: keyword_first

    open(unit=10, file=file_gap, status="old", iostat=iostatus)
!   Look for the number of instances of each GAP
    n_soap_turbo = 0
    iostatus = 0
    do while(iostatus==0)
      read(10, *, iostat=iostatus) keyword
      keyword = trim(keyword)
      if(iostatus/=0)then
        exit
      end if
      keyword_first = keyword(1:1)
      if(keyword_first=='#' .or. keyword_first=='!')then
        continue
      else if(keyword=='gap_beg')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, keyword
        if( keyword == "soap_turbo" )then
          n_soap_turbo = n_soap_turbo + 1
        end if
      end if
    end do
!   Allocate the variables
    if( n_soap_turbo > 0 )then
      allocate( zeta(1:n_soap_turbo) )
      allocate( file_alphas(1:n_soap_turbo) )
      allocate( file_desc(1:n_soap_turbo) )
    end if

!   Now record the hypers
    rewind(10)
    n_soap_turbo = 0
    iostatus = 0
    do while(iostatus==0)
      read(10, *, iostat=iostatus) keyword
      keyword = trim(keyword)
      if(iostatus/=0)then
        exit
      end if
      keyword_first = keyword(1:1)
      if(keyword_first=='#' .or. keyword_first=='!')then
        continue
      else if(keyword=='gap_beg')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, keyword
!       soap_turbo definitions here
        if( keyword == "soap_turbo" )then
          n_soap_turbo = n_soap_turbo + 1
          iostatus = 0
          do while( keyword /= "gap_end" .and. iostatus == 0 )
            read(10, *, iostat=iostatus) keyword
            if( keyword == "zeta" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, zeta(n_soap_turbo)
            else if( keyword == "alphas_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, file_alphas(n_soap_turbo)
            else if( keyword == "desc_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, file_desc(n_soap_turbo)
            end if
          end do
        end if
      end if
    end do
    close(10)

  end subroutine
!**************************************************************************

!**************************************************************************
!
! This subroutine is a shorter version of src/read_files.f90 subroutine with
! same name, which assumes the input is soap_turbo descriptors
! Note that the dimensions of Qs are changed, i.e., Qs(1:n_sparse, 1:n_soap)
!
  subroutine read_alphas_and_descriptors(file_desc, file_alphas, n_sparse, &
                                         n_soap, alphas, Qs)

    implicit none

!   Input variables
    character*1024, intent(in) :: file_desc, file_alphas
!   Output variables
    real*8, allocatable, intent(out) :: alphas(:), Qs(:,:)
    integer, intent(out) :: n_sparse, n_soap
!   Internal variables
    integer :: i, j, iostatus, unit_number

!   Read alphas to figure out sparse set size
    open(newunit=unit_number, file=file_alphas, status="old")
    iostatus = 0
    n_sparse = -1
    do while(iostatus == 0)
      read(unit_number, *, iostat=iostatus)
      n_sparse = n_sparse + 1
    end do
    close(unit_number)

!   Get SOAP descriptor vectors size
    open(newunit=unit_number, file=file_desc, status="old")
    iostatus = 0
    i = -1
    do while(iostatus == 0)
      read(unit_number, *, iostat=iostatus)
      i = i + 1
    end do
    n_soap = i / n_sparse
    close(unit_number)

!   Read descriptor vectors in spare set
!   Read alphas SOAP
    allocate( alphas(1:n_sparse) )
    open(newunit=unit_number, file=file_alphas, status="old")
    do i = 1, n_sparse
      read(unit_number, *) alphas(i)
    end do
    close(unit_number)
!   Read sparse set descriptors
    allocate( Qs(1:n_soap, 1:n_sparse) )
    open(newunit=unit_number, file=file_desc, status="old")
    do i = 1, n_sparse
      do j = 1, n_soap
        read(unit_number, *) Qs(j, i)
      end do
    end do
    close(unit_number)

  end subroutine
!**************************************************************************


end program get_linear_model
