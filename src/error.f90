! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, error.f90, is copyright (c) 2019-2025, Miguel A. Caro and
! HND X   Tigany Zarrouk
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

module error
   use kinds, only: dp
#ifdef _MPIF90
   use mpi
#endif

contains

!  Exit non-zero, and take the whole job down under MPI.
!
!  The original exited via a bare `stop`, which is status 0: a fatal input error
!  reported itself as success, so nothing checking a return code -- a script, a
!  queue system, the regression suite -- could tell the run had failed. That is
!  the same trap as the debugging exit(0) recorded in the GPU branch's handoff
!  section 6b, which produced no trajectory and still looked like a clean run.
!
!  MPI_abort rather than MPI_finalize: finalize is the *orderly* shutdown and
!  expects every rank to call it. One rank hitting a bad keyword and finalizing
!  alone leaves the others waiting in the next collective.
   subroutine turbogap_abort()
      integer :: ierror

!     Flush before aborting. MPI_abort kills the process outright, and Fortran
!     unit 6 is buffered, so the diagnostic every caller writes just before
!     getting here -- which keyword was wrong, and what the valid values are --
!     was discarded along with the rest of the buffer. What the user saw was
!     "with errorcode 1" and nothing else.
      flush (6)
#ifdef _MPIF90
      call MPI_abort(MPI_COMM_WORLD, 1, ierror)
#endif
      stop 1
   end subroutine turbogap_abort

end module error
