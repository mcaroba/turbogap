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

! The device memory budget this rank may use, in bytes. Set by
! gpu_memory_budget_init from the device's own free memory; zero means it was
! never established and every caller falls back to the input's fixed keywords.
   integer(c_size_t) :: gpu_mem_budget = 0_c_size_t

contains

!**************************************************************************
!  Establish the device memory budget, and report it.
!
!  Called once, after gpu_set_device and before anything allocates. It has to be
!  after the device is chosen -- hipMemGetInfo on an uninitialised context
!  returns the wrong card -- and before the first allocation, or the "free"
!  figure already has TurboGAP's own buffers subtracted from it.
!
!  Why a FRACTION of free rather than all of it. Three things are not visible
!  here and all of them consume device memory: cuBLAS/cutlass workspaces, the
!  driver's own context growth, and fragmentation inside the stream-ordered
!  pool, which can refuse a large contiguous request while reporting plenty
!  free. Sizing batches to 100% of free reliably produces a run that dies near
!  the end, which is the worst time.
!
!  n_ranks_on_device divides the budget. One rank per GPU is the arrangement the
!  slurm scripts make, but nothing enforces it, and two ranks each sizing their
!  batches to the whole card is a guaranteed failure that looks like a code bug.
   subroutine gpu_memory_budget_init(params, rank, ntasks)
      type(input_parameters), intent(inout) :: params
      integer, intent(in) :: rank
      integer, intent(in) :: ntasks

      integer(c_size_t) :: dev_free, dev_total, live, peak
      real(dp) :: budget_gb
      real(dp) :: free_gb
      integer :: n_devices
      integer :: n_ranks_on_device

      call gpu_mem_stats(dev_free, dev_total, live, peak)
      free_gb = real(dev_free, dp)/1024.d0**3

!     cuda_set_device assigns devices round-robin (rank mod n_devices), so this
!     many ranks are looking at the same free memory and would each budget all
!     of it. One rank per GPU is the usual slurm arrangement and makes this 1.
      n_devices = max(1, int(gpu_device_count()))
      n_ranks_on_device = max(1, (ntasks + n_devices - 1)/n_devices)

      if (params%gpu_mem_fraction <= 0.d0) then
!        Opted out: the input's own max_Gbytes_per_process and gpu_n_batches
!        stand exactly as written. Still report, because knowing the card has
!        5.7 GB while max_Gbytes_per_process says 1.0 is most of the diagnosis.
         gpu_mem_budget = 0_c_size_t
         if (rank == 0) then
            write (*, '(A,F8.3,A,F8.3,A)') ' <<<< GPU MEM >>>> device has ', free_gb, &
               ' GB free of ', real(dev_total, dp)/1024.d0**3, ' GB'
            write (*, '(A,F8.3,A)') ' <<<< GPU MEM >>>> budgeting is OFF; max_Gbytes_per_process = ', &
               params%max_Gbytes_per_process, ' as given'
            write (*, '(A)') ' <<<< GPU MEM >>>> set gpu_mem_fraction (e.g. 0.8) to size it from the device'
         end if
      else
         budget_gb = params%gpu_mem_fraction*free_gb/dfloat(max(1, n_ranks_on_device))
         gpu_mem_budget = int(budget_gb*1024.d0**3, c_size_t)

!        max_Gbytes_per_process is the SOAP loop's divisor: get_number_of_atom_pairs
!        splits the descriptor loop until an estimated batch fits inside it. It
!        defaults to 1.0 GB, a number chosen with no reference to any device, so
!        on this card it asks for 6x more batches than necessary and on a 80 GB
!        card 80x. Overwriting it here is the whole point of the keyword.
         params%max_Gbytes_per_process = budget_gb

         if (rank == 0) then
            write (*, '(A,F8.3,A,F8.3,A)') ' <<<< GPU MEM >>>> device has ', free_gb, &
               ' GB free of ', real(dev_total, dp)/1024.d0**3, ' GB'
            write (*, '(A,F6.3,A,I0,A,F8.3,A)') ' <<<< GPU MEM >>>> budget = ', params%gpu_mem_fraction, &
               ' x free / ', max(1, n_ranks_on_device), ' rank(s) = ', budget_gb, ' GB per rank'
            write (*, '(A,F8.3)') ' <<<< GPU MEM >>>> max_Gbytes_per_process set from the device to ', budget_gb
         end if
      end if

!     So an out-of-memory abort can print the keyword and its value rather than
!     just "out of memory".
      call gpu_mem_set_hint(params%max_Gbytes_per_process, params%gpu_n_batches)

   end subroutine gpu_memory_budget_init
!**************************************************************************

!**************************************************************************
!  How many batches a loop needs, given what one unbatched pass would cost.
!
!  need_gb comes from the caller because only the caller knows what its kernels
!  allocate. The tree already contains exactly this kind of estimate --
!  estimate_max_exp_forces_device_memory_usage enumerates the pdf/xrd buffers
!  one by one, and get_number_of_atom_pairs has the same thing for SOAP folded
!  into a 150-byte-per-pair constant -- so this is the consumer for them rather
!  than a competing model.
!
!  Two properties matter more than the arithmetic:
!
!  * It returns n_batches_in unchanged when budgeting is off. gpu_mem_fraction
!    defaults to zero, so nothing about an existing input changes until it is
!    set.
!
!  * It never returns FEWER batches than asked for. Raising the count is a
!    safety measure that cannot make a run fail; lowering it would override a
!    setting the user may have arrived at precisely because this estimate is
!    wrong for their system. The user's number is a floor, not a suggestion --
!    which is the answer to "it is likely that a manual setting is necessary".
   integer function gpu_batches_for_gb(need_gb, n_batches_in, n_sites, rank, what)
      real(dp), intent(in) :: need_gb
      integer, intent(in) :: n_batches_in
      integer, intent(in) :: n_sites
      integer, intent(in) :: rank
      character(len=*), intent(in) :: what

      real(dp) :: budget_gb

      gpu_batches_for_gb = n_batches_in

      if (gpu_mem_budget <= 0_c_size_t) return
      if (need_gb <= 0.d0) return

      budget_gb = real(gpu_mem_budget, dp)/1024.d0**3

      gpu_batches_for_gb = max(1, ceiling(need_gb/budget_gb))
!     A batch cannot be smaller than one site, so this is the hard ceiling. If
!     it is reached the estimate says one site's worth still does not fit, and
!     the run will fail however it is batched -- worth saying so rather than
!     letting it look like the batching handled it.
      if (gpu_batches_for_gb >= n_sites) then
         gpu_batches_for_gb = max(1, n_sites)
         if (rank == 0) then
            write (*, '(A,A,A)') ' <<<< GPU MEM >>>> WARNING: ', trim(what), &
               ' needs more memory than the device has even at one site per batch.'
            write (*, '(A)') ' <<<< GPU MEM >>>> This run is expected to fail. Reduce the cutoff, the'
            write (*, '(A)') ' <<<< GPU MEM >>>> number of samples, or the system size, or use more ranks.'
         end if
      end if

      if (gpu_batches_for_gb < n_batches_in) gpu_batches_for_gb = n_batches_in

      if (rank == 0 .and. gpu_batches_for_gb /= n_batches_in) then
         write (*, '(A,A,A,F12.3,A,F8.3,A,I0,A,I0)') ' <<<< GPU MEM >>>> ', trim(what), &
            ' needs ~', need_gb, ' GB, budget ', budget_gb, &
            ' GB/rank -> batches ', n_batches_in, ' -> ', gpu_batches_for_gb
      end if

   end function gpu_batches_for_gb
!**************************************************************************

!**************************************************************************
!  Print the ledger next to the device's own figures.
   subroutine gpu_memory_report(label)
      character(len=*), intent(in) :: label

      call gpu_mem_report(trim(label)//c_null_char)
   end subroutine gpu_memory_report
!**************************************************************************

!**************************************************************************
!  Which stream the calling thread owns.
!
!  The only correct answer is the thread's own index. Two threads running at
!  the same time must never be handed the same stream: the work they enqueue
!  would be serialised against each other (which defeats the point) and, worse,
!  a buffer one thread allocates stream-ordered could be freed by the other's
!  hipFreeAsync on the same stream while the first is still using it.
!
!  Both batched loops got this wrong in different ways, which is why it is a
!  function now rather than an expression repeated at each site:
!
!    * the pdf loop computed mod(i-1, n_omp)+1 and then assigned 1 over the
!      top unconditionally, so every batch used stream 1 whatever the build;
!    * the xrd loop kept mod(i-1, n_omp)+1, which is a property of the
!      ITERATION, not of the thread. Under any schedule that does not hand out
!      iterations strictly round-robin -- and static, the default, does not --
!      two concurrent threads can draw the same index.
!
!  Outside a parallel region omp_get_thread_num() is 0, so this returns 1 and
!  the serial path uses gpu_streams(1), exactly as before.
   integer function gpu_omp_task()
      implicit none

      gpu_omp_task = 1
!$    gpu_omp_task = omp_get_thread_num() + 1
   end function gpu_omp_task
!**************************************************************************

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
