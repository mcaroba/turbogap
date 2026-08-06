! HND XX
! HND XX This file is part of the TurboGAP code
! HND XX
! HND XX gpu_context.f90, the shared GPU device context
! HND XX
module gpu_context

! The device state that the electrostatics, SOAP and experimental-observable
! paths all touch. It used to be declared in turbogap.f90 and reached by host
! association, which is why no one of those blocks could be lifted out of the
! driver: each one names the stream and the batch storage, and so does the code
! around it.
!
! Holding it here rather than passing it means an extraction stays a verbatim
! move -- the lifted code keeps the same names, and the driver's existing uses
! are untouched. That is the same reason gap_backend_gpu's ten buffers became
! module state.
!
! Deliberately NOT here:
!   * omp_task, n_omp -- omp_task is intended to be OpenMP-private (see the
!     '!$omp parallel private(omp_task)' and '!$OMP PRIVATE(i, omp_task, ...)'
!     directives in turbogap.f90). Making a thread-private index into shared
!     module state would be a correctness bug the moment -fopenmp is enabled.
!     The arrays indexed BY it (gpu_streams, cublas_handles) are shared and do
!     belong here; the index does not.
!   * i_beg_list / i_end_list / j_beg_list / j_end_list -- the batch
!     decomposition. These are recomputed per snapshot and handed to
!     neighbors.f90 as intent(out); they read as data flowing through the
!     driver, not as context, so they are passed as arguments.
!   * params%gpu_batched, params%gpu_n_batches, params%gpu_max_batch_size --
!     input parameters, already carried by params.

   use kinds
   use iso_c_binding
   use types, only: gpu_storage_type, gpu_neigh_storage_type, &
                    gpu_host_batch_storage_type, input_parameters
   use F_B_C
!$ use omp_lib

   implicit none

   public

! The default stream and cuBLAS handle, created once at start-up.
   type(c_ptr) :: cublas_handle
   type(c_ptr) :: gpu_stream

! Per-OpenMP-task streams and handles. Shared arrays, indexed by the
! thread-private omp_task.
   type(c_ptr), allocatable :: cublas_handles(:)
   type(c_ptr), allocatable :: gpu_streams(:)

! Device-side storage for the batched experimental-observable paths.
! gpu_batch_storage is indexed
!   gpu_batch_storage(i_batch) % host(i_n_dim_idx) % pair_distribution_h(1:n_samples)
   type(gpu_storage_type), allocatable :: gpu_exp(:)
   type(gpu_neigh_storage_type), allocatable :: gpu_neigh(:)
   type(gpu_host_batch_storage_type), allocatable :: gpu_batch_storage(:)

! Running total of device memory handed out, for the end-of-run report.
   real(dp) :: gpu_memory_usage = 0.d0

contains

!**************************************************************************
!  Bring the device up, and take it down.
!
!  These two exist on both branches under the same names and the same argument
!  lists; the CPU branch's src/gpu_context.f90 implements them as empty stubs
!  and the Makefile compiles one or the other, exactly as it does for
!  gap_backend_cpu / gap_backend_gpu. The driver calls the names and never
!  learns which implementation it got, so ~50 lines of device set-up and
!  teardown stop being a difference between the two turbogap.f90 files.
!
!  n_omp is an OUT argument rather than another module variable, deliberately:
!  the CPU stub has a correct answer for it (one stream's worth of work), and
!  keeping it out of the module keeps the exclusion documented at the top of
!  this file intact.
   subroutine gpu_context_init(params, rank, n_omp)
      type(input_parameters), intent(in) :: params
      integer, intent(in) :: rank
      integer, intent(out) :: n_omp

      integer :: omp_task, n_omp_temp

      n_omp = 1

!     One visible device per rank; the slurm submission script arranges that.
      call gpu_set_device(rank)

      call create_cublas_handle(cublas_handle, gpu_stream)

      if (params%gpu_batched) then
!        Do not ask for more threads than there are batches to give them.
#ifdef _OPENMP
         n_omp_temp = omp_get_max_threads()

         !$omp parallel DEFAULT(SHARED)
         n_omp_temp = omp_get_num_threads()
         if (n_omp_temp > params%gpu_n_batches) n_omp_temp = params%gpu_n_batches
         call omp_set_num_threads(n_omp_temp)
         !$omp end parallel

         n_omp = omp_get_max_threads()
#endif
         allocate (cublas_handles(1:n_omp))
         allocate (gpu_streams(1:n_omp))
         do omp_task = 1, n_omp
            call create_cublas_handle(cublas_handles(omp_task), gpu_streams(omp_task))
         end do

         if (rank == 0) write (*, '(A,I0,A)') &
            ' <<<< OPENMP >>>> created ', n_omp, ' streams for the batched gpu calculation'
      end if

   end subroutine gpu_context_init
!**************************************************************************

!**************************************************************************
   subroutine gpu_context_finalize(params, n_omp)
      type(input_parameters), intent(in) :: params
      integer, intent(in) :: n_omp

      integer :: i

      if (params%gpu_batched) then
         do i = 1, n_omp
            call destroy_cublas_handle(cublas_handles(i), gpu_streams(i))
         end do
      end if

      call destroy_cublas_handle(cublas_handle, gpu_stream)
      call gpu_device_reset()

   end subroutine gpu_context_finalize
!**************************************************************************

end module gpu_context
