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
                   gpu_host_batch_storage_type

  implicit none

  public

! The default stream and cuBLAS handle, created once at start-up.
  type(c_ptr) :: cublas_handle, gpu_stream

! Per-OpenMP-task streams and handles. Shared arrays, indexed by the
! thread-private omp_task.
  type(c_ptr), allocatable :: cublas_handles(:), gpu_streams(:)

! Device-side storage for the batched experimental-observable paths.
! gpu_batch_storage is indexed
!   gpu_batch_storage(i_batch) % host(i_n_dim_idx) % pair_distribution_h(1:n_samples)
  type( gpu_storage_type ),            allocatable :: gpu_exp(:)
  type( gpu_neigh_storage_type ),      allocatable :: gpu_neigh(:)
  type( gpu_host_batch_storage_type ), allocatable :: gpu_batch_storage(:)

! Running total of device memory handed out, for the end-of-run report.
  real(dp) :: gpu_memory_usage = 0.d0

end module gpu_context
