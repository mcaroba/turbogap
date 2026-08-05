! Wall-clock timing.
!
! The GPU branch replaced every cpu_time call in the driver with a get_time
! wrapper, because under MPI cpu_time measures per-rank processor time, which
! is not what you want to report for a run whose cost is dominated by device
! work and communication. Adopting the same primitive here removes ~76
! gratuitous differences between the two turbogap.f90 files and makes the two
! branches' timing numbers directly comparable.
!
! Timings are printed, never written to any file compared by
! tests/regression/run.sh, so switching the primitive cannot move a result.

module timing

  use kinds

#ifdef _MPIF90
  use mpi
#endif

  implicit none

contains

  !**************************************************************************
  subroutine get_time(time)
    implicit none
 real(dp), intent(out) :: time

#ifdef _MPIF90
    time = MPI_Wtime()
#else
    call cpu_time(time)
#endif
  end subroutine get_time
  !**************************************************************************

end module timing
