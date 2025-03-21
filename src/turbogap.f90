! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap.f90, is copyright (c) 2019-2023, Miguel A. Caro and
! HND X   Tigany Zarrouk
! HND X   Uttiyoarnab saha
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

program turbogap

  use turbogap_wrap, only: turbogap_routine
#ifdef _MPIF90
  use mpi
#endif
  implicit none
  integer :: ierr, mpierr
  character*16 :: mode="none"

#ifdef _MPIF90
  ! initialize mpi
  call mpi_init(ierr)
#endif

  !**************************************************************************
  ! Read the mode. It should be "soap", "predict" or "md"
  !
  call get_command_argument(1,mode)
  if( mode == "" .or. mode == "none" )then
     write(*,*) "ERROR: you need to run 'turbogap md' or 'turbogap predict'"
     stop
     ! THIS SHOULD BE FIXED, IN CASE THE USER JUST WANT TO OUTPUT THE SOAP DESCRIPTORS
     mode = "soap"
  end if
  !**************************************************************************


  ! execute turbogap
  call turbogap_routine( MPI_COMM_WORLD, mode, "input", "trajectory_out.xyz", ierr )

  ! check if there is an error
  if( ierr/= 0 ) then
     ! handle the error in some way, if needed
     stop
  end if

#ifdef _MPIF90
  ! succesful termination
  call mpi_finalize(ierr)
#endif

end program turbogap
