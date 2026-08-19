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
!
! -- times_t ----------------------------------------------------------------
!
! Every bucket is a three-element accumulator: (1) the start stamp, (2) the
! end stamp, (3) the running total. time_start and time_end are the only
! things that touch (1) and (2); everything that reports reads (3). Writing
! the accumulation by hand at each site is how a bucket silently loses its
! `+`, which is the same "one statement, N copies that must agree" shape that
! produced the ts+mbd, MPI-reduce and electrostatics defects.
!
! The type carries the union of both branches' buckets, so this file is
! byte-identical on the CPU and GPU trees and the extracted modules
! (turbogap_setup, turbogap_exp, turbogap_md, gap_backend_*) take one
! argument instead of thirteen. A bucket the local branch never starts stays
! at zero and costs three doubles.
!
! -- the accounting ---------------------------------------------------------
!
! Buckets nest. sum_times adds only the PARENT buckets, and Miscellaneous is
! total minus that sum. A child bucket added to the sum is subtracted twice
! and Miscellaneous goes negative, which is exactly the defect this replaces:
! time_gap was accumulated once per SOAP descriptor from a stamp taken before
! the loop, so an N-descriptor potential counted the same interval N times,
! and it enclosed an mpi_reduce that was also charged to time_mpi_ef.
!
! The rule when a new bucket's nesting is not obvious: do NOT add it to
! sum_times. An over-large Miscellaneous is a correct-signed statement of
! ignorance; a negative one is a lie, and it is the reason the printed totals
! were unusable as refactor instrumentation.

module timing

   use kinds

#ifdef _MPIF90
   use mpi
#endif

   implicit none

   !**************************************************************************
   type times_t
      !! Wall-clock accumulators. (1) start, (2) end, (3) running total.
      !!
      !! To add a bucket: declare it here, decide whether it is a parent or a
      !! child, and add it to sum_times ONLY if it is a parent (see above).
      !! Then print it in print_times.

      !----- the reference -----------------------------------------------
      real(dp) :: total(3) = 0.0_dp

      !----- parents: disjoint, and summed by sum_times -------------------
      !! setup is everything before the first pass of the main loop: the input
      !! file, the potential files, and the allocation and broadcast that
      !! follow them. read_input measures part of that and is reported under
      !! it, never summed -- see below.
      real(dp) :: setup(3) = 0.0_dp
      real(dp) :: read_xyz(3) = 0.0_dp
      real(dp) :: neigh(3) = 0.0_dp
      real(dp) :: gap(3) = 0.0_dp
      real(dp) :: vdw(3) = 0.0_dp
      real(dp) :: estat(3) = 0.0_dp
      real(dp) :: pdf(3) = 0.0_dp
      real(dp) :: sf(3) = 0.0_dp
      real(dp) :: xrd(3) = 0.0_dp
      real(dp) :: nd(3) = 0.0_dp
      real(dp) :: xps(3) = 0.0_dp
      !! The MAD IR bias: everything between the dipole arriving and the
      !! biased forces leaving. The MPI all-reduce that precedes it is charged
      !! to `mpi`, so this bucket starts after it and the two stay disjoint.
      real(dp) :: ir(3) = 0.0_dp
      real(dp) :: md(3) = 0.0_dp
      real(dp) :: mc(3) = 0.0_dp
      real(dp) :: mpi(3) = 0.0_dp
      real(dp) :: mpi_positions(3) = 0.0_dp
      real(dp) :: mpi_ef(3) = 0.0_dp

      !----- children of setup: reported, never summed ---------------------
      !! read_input is the input and potential files. mpi_setup is the
      !! broadcast of what they produced, which happens inside read_input and
      !! is why it cannot share the `mpi` bucket with the per-step ones: one
      !! bucket used both nested and standalone is counted twice by sum_times
      !! however carefully the list is maintained.
      real(dp) :: read_input(3) = 0.0_dp
      real(dp) :: mpi_setup(3) = 0.0_dp

      !----- children of gap: reported, never summed ----------------------
      real(dp) :: soap(3) = 0.0_dp
      real(dp) :: gap_2b(3) = 0.0_dp
      real(dp) :: gap_3b(3) = 0.0_dp
      real(dp) :: gap_core_pot(3) = 0.0_dp
      real(dp) :: get_soap(3) = 0.0_dp        ! GPU branch
      real(dp) :: local_prop(3) = 0.0_dp      ! GPU branch
      real(dp) :: soap_solo(3) = 0.0_dp       ! GPU branch, "lolo__soap"
      real(dp) :: soap_lin(3) = 0.0_dp        ! GPU branch, "lin__turbo"

      !----- children of ir: reported, never summed ------------------------
      !! predict is the autocorrelation, the cosine transform, the loss and
      !! lambda -- mad_ir_evaluate. forces is the 9*n_atoms contraction of
      !! lambda with dmu/dr. io is the restart buffer and the spectrum files.
      !!
      !! What none of them measure is the descriptor second derivatives that
      !! produce dmu/dr in the first place: those are built inside get_soap
      !! and are charged to `soap`, and they are paid on every collected step
      !! whether or not the ensemble is full yet. So ir_predict rising from
      !! zero is the cost of the spectrum, not the cost of the bias.
      real(dp) :: ir_predict(3) = 0.0_dp
      real(dp) :: ir_forces(3) = 0.0_dp
      real(dp) :: ir_io(3) = 0.0_dp

      !----- GPU-only, nesting not established: never summed --------------
      real(dp) :: batch_alloc(3) = 0.0_dp
      real(dp) :: batch_pdf(3) = 0.0_dp
      real(dp) :: batch_xrd(3) = 0.0_dp
      real(dp) :: create_streams(3) = 0.0_dp
      real(dp) :: exp_batched(3) = 0.0_dp
   end type times_t
   !**************************************************************************

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

   !**************************************************************************
   subroutine time_start(bucket)
      !! Stamp the start of an interval. Pairs with time_end on the same bucket.
      implicit none
      real(dp), intent(inout) :: bucket(3)

      call get_time(bucket(1))
   end subroutine time_start
   !**************************************************************************

   !**************************************************************************
   subroutine time_end(bucket)
      !! Close the interval opened by time_start and add it to the total.
      implicit none
      real(dp), intent(inout) :: bucket(3)

      call get_time(bucket(2))
      bucket(3) = bucket(3) + bucket(2) - bucket(1)
   end subroutine time_end
   !**************************************************************************

   !**************************************************************************
   pure function sum_times(time) result(total)
      !! The sum of the PARENT buckets only. Adding a child here makes
      !! Miscellaneous negative; see the header.
      implicit none
      type(times_t), intent(in) :: time
      real(dp) :: total

      total = time%setup(3) + time%read_xyz(3) + time%neigh(3) &
              + time%gap(3) + time%vdw(3) + time%estat(3) &
              + time%pdf(3) + time%sf(3) + time%xrd(3) + time%nd(3) + time%xps(3) &
              + time%ir(3) &
              + time%md(3) + time%mc(3) &
              + time%mpi(3) + time%mpi_positions(3) + time%mpi_ef(3)
   end function sum_times
   !**************************************************************************

end module timing
