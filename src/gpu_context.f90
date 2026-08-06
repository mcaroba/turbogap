! Device context: the CPU implementation.
!
! One module name, two implementation files, chosen in the Makefile -- the same
! arrangement as gap_backend_cpu / gap_backend_gpu. The GPU branch's
! src/gpu_context.f90 owns the default stream, the cuBLAS handles, the
! per-OpenMP-task stream arrays and the device-side batch storage, and brings
! them up and down in these two procedures. Here there is nothing to bring up,
! so both bodies are empty.
!
! The point is not the eight lines it saves. It is that ~50 lines of device
! set-up and teardown stop being a difference between the two turbogap.f90
! files: the driver calls gpu_context_init and gpu_context_finalize on both
! branches and never learns which implementation it got.
!
! n_omp is an argument rather than module state on either side. The CPU has a
! correct answer for it -- one stream's worth of work -- and on the GPU branch
! keeping it out of the module preserves the exclusion documented at the top of
! that file: the arrays indexed by an OpenMP-private index are shared context,
! the index itself is not.

module gpu_context

   use types, only: input_parameters

   implicit none

   public

contains

   !**************************************************************************
   subroutine gpu_context_init(params, rank, n_omp)
      implicit none
      type(input_parameters), intent(in) :: params
      integer, intent(in) :: rank
      integer, intent(out) :: n_omp

      n_omp = 1

   end subroutine gpu_context_init
   !**************************************************************************

   !**************************************************************************
   subroutine gpu_context_finalize(params, n_omp)
      implicit none
      type(input_parameters), intent(in) :: params
      integer, intent(in) :: n_omp

   end subroutine gpu_context_finalize
   !**************************************************************************

end module gpu_context
