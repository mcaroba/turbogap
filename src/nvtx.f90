! HND XX
! HND XX This file is part of the TurboGAP code
! HND XX
! HND XX nvtx.f90, host-side range markers for the NVIDIA profilers
! HND XX
module nvtx

! Why this exists
! ---------------
! nsys shows two rows that are hard to join by eye: what the CPU is doing, and
! what the GPU is doing. TurboGAP's phases (neighbour build, SOAP, 2b, 3b, vdW,
! electrostatics, the batched pdf/xrd loops) are host-side control flow, so
! without markers a timeline is a wall of kernel names with no statement of
! which phase launched them. An NVTX range makes the phase a first-class row,
! and -- the reason this was written -- it makes GPU IDLE inside a phase
! visible, which is the measurement the stream work depends on.
!
! It is deliberately tied to the timing buckets rather than being a second,
! independent set of markers. times_t already names every phase anyone cares
! about, and every phase already has exactly one time_start/time_end pair. A
! separate set of push/pop calls would be a second thing to keep in sync with
! the first, and the two would disagree the first time a bucket moved.
!
! Cost when off
! -------------
! Without -D_NVTX the two routines have empty bodies and gfortran inlines them
! away. Without the profiler attached, libnvToolsExt's push/pop are a load and
! a predictable branch: measured below the noise of the wall-clock call that
! time_start already makes. So the markers are compiled in for PROFILE=1 builds
! only, and cost nothing in either case.
!
! Nesting
! -------
! NVTX ranges are a per-thread stack: every push must be popped on the same
! thread, in order. time_start/time_end take the label as an OPTIONAL argument
! precisely so this cannot break -- a call site that labels its start must
! label its end, and one that labels neither pushes nothing. A partially
! labelled pair would unbalance the stack, so do not add one.

   use iso_c_binding

   implicit none

   private
   public :: nvtx_push, nvtx_pop

#ifdef _NVTX
!  libnvToolsExt exports these as ordinary C symbols (the v2 ABI), so they bind
!  directly and no C shim is needed. -lnvToolsExt is already on the link line of
!  every CUDA architecture makefile.
   interface
      subroutine nvtxRangePushA(name) bind(C, name="nvtxRangePushA")
         import :: c_char
         character(kind=c_char), dimension(*), intent(in) :: name
      end subroutine nvtxRangePushA

      subroutine nvtxRangePop() bind(C, name="nvtxRangePop")
      end subroutine nvtxRangePop
   end interface
#endif

contains

!**************************************************************************
   subroutine nvtx_push(label)
      !! Open a named range on the calling thread's NVTX stack.
      implicit none
      character(len=*), intent(in) :: label

#ifdef _NVTX
      call nvtxRangePushA(trim(label)//c_null_char)
#endif
   end subroutine nvtx_push
!**************************************************************************

!**************************************************************************
   subroutine nvtx_pop()
      !! Close the innermost range on the calling thread's NVTX stack.
      implicit none

#ifdef _NVTX
      call nvtxRangePop()
#endif
   end subroutine nvtx_pop
!**************************************************************************

end module nvtx
