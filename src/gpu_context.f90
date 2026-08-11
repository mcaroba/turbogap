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
!
! gpu_memory_budget_init is the third procedure the driver calls by the same
! name on both branches. It answers the same question from a different place:
! how much memory may one rank give the SOAP descriptor loop. The GPU branch
! asks the device; here we ask the node.

module gpu_context

   use kinds
   use types, only: input_parameters
#ifdef _MPIF90
   use mpi
#endif

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

   !**************************************************************************
   !  Establish this rank's host memory budget for the SOAP descriptor loop,
   !  and report it.
   !
   !  max_Gbytes_per_process is get_number_of_atom_pairs's divisor: it splits
   !  the descriptor loop until an estimated batch fits inside it. Its default
   !  of 1.0 GB was chosen with no machine in mind. On a 256 GB node that asks
   !  for roughly sixty times more batches than necessary, and every extra batch
   !  repeats the whole per-batch setup -- assign_species_multiplicity, the
   !  neighbour-subset gather, the allocation and the scatter back -- for the
   !  same total work.
   !
   !  So this reads the node's memory and sets the keyword from it. The GPU
   !  branch does exactly this from hipMemGetInfo under the name
   !  gpu_memory_budget_init; the driver's call is identical on both branches
   !  and neither driver learns which one it got.
   !
   !  Three decisions worth stating, because each one is a place this could be
   !  wrong in a way that looks like a code bug:
   !
   !  * It reads the CGROUP limit before /proc/meminfo. Under slurm, --mem puts
   !    the job in a cgroup whose limit is a fraction of the node, and
   !    MemTotal reports the whole node regardless. Budgeting from MemTotal on a
   !    shared node is a reliable way to be OOM-killed while the code believes
   !    it left 75% of memory spare.
   !
   !  * It uses the TOTAL, not what is free right now. Free memory moves between
   !    runs and between MD steps, and the batch count changes the order the
   !    per-batch force contributions are summed -- so a budget read from live
   !    memory makes two runs of the same input differ in the last bits. A
   !    ceiling that is a property of the machine and the rank count gives the
   !    same split every time.
   !
   !  * The fraction is well under 1. The SOAP batch is not the only thing
   !    resident: the global neighbour arrays (rjs, thetas, phis, xyz and the
   !    list itself, ~52 bytes per pair over ALL pairs, not just this batch),
   !    the positions, forces and velocities, and the descriptor outputs are all
   !    live at the same time and none of them is counted by the estimate this
   !    keyword bounds.
   subroutine gpu_memory_budget_init(params, rank, ntasks)
      implicit none
      type(input_parameters), intent(inout) :: params
      integer, intent(in) :: rank
      integer, intent(in) :: ntasks

      integer(8) :: mem_bytes
      integer :: n_ranks_on_node
      real(dp) :: total_gb, budget_gb
      character(len=32) :: source
      character(len=64) :: line

      call host_memory_limit(mem_bytes, source)
      n_ranks_on_node = ranks_sharing_node(ntasks)

      call banner(rank, '')
      call banner(rank, 'Host memory budget:')

      !  Nothing to go on. Leave the keyword exactly as it was -- a guess here
      !  would be worse than the documented default, because it would look
      !  like a measurement.
      if (mem_bytes <= 0_8) then
         call banner(rank, 'could not read the node memory')
         write (line, '(A,F10.3,A)') 'max_Gbytes_per_process = ', params%max_Gbytes_per_process, ' GB'
         call banner(rank, line)
         call banner(rank, '')
         return
      end if

      total_gb = real(mem_bytes, dp)/1024.d0**3
      write (line, '(A,F10.1,A,A,A)') 'node has ', total_gb, ' GB (', trim(source), ')'
      call banner(rank, line)

      !  An explicit value in the input wins. Someone who set this keyword did
      !  so because the estimate was wrong for their system, and silently
      !  overriding them would be worse than never estimating at all.
      if (params%max_Gbytes_set) then
         write (line, '(A,F10.3,A)') 'max_Gbytes_per_process = ', params%max_Gbytes_per_process, ' GB'
         call banner(rank, line)
         call banner(rank, 'as given in the input, not sized here')
         call banner(rank, '')
         return
      end if

      if (params%mem_fraction <= 0.d0) then
         call banner(rank, 'sizing is OFF (mem_fraction <= 0)')
         write (line, '(A,F10.3,A)') 'max_Gbytes_per_process = ', params%max_Gbytes_per_process, ' GB'
         call banner(rank, line)
         call banner(rank, '')
         return
      end if

      budget_gb = params%mem_fraction*total_gb/dfloat(max(1, n_ranks_on_node))
      params%max_Gbytes_per_process = budget_gb

      write (line, '(A,F6.3,A,I0,A)') 'budget = ', params%mem_fraction, ' x total / ', &
         max(1, n_ranks_on_node), ' rank(s)'
      call banner(rank, line)
      write (line, '(A,F10.3,A)') 'max_Gbytes_per_process = ', budget_gb, ' GB'
      call banner(rank, line)
      call banner(rank, '')

   end subroutine gpu_memory_budget_init
   !**************************************************************************

   !**************************************************************************
   !  One line of the start-up box: 39 columns of text, then the wall.
   !
   !  The surrounding printouts draw that wall by hand-counting spaces into the
   !  format string, which is why they go crooked whenever a number gains a
   !  digit. Padding a buffer instead keeps it straight for any value.
   subroutine banner(rank, text)
      implicit none
      integer, intent(in) :: rank
      character(len=*), intent(in) :: text

      character(len=39) :: padded

      if (rank /= 0) return
      padded = adjustl(text)
      write (*, '(1X,A39,A)') padded, '|'

   end subroutine banner
   !**************************************************************************

   !**************************************************************************
   !  How many ranks share this node's memory.
   !
   !  Every rank on a node budgets from the same total, so without this each of
   !  them would claim all of it. MPI answers this exactly --
   !  MPI_COMM_TYPE_SHARED splits the world into groups that can share memory --
   !  which is better than dividing by ntasks: ntasks is right only when the job
   !  runs on one node.
   integer function ranks_sharing_node(ntasks)
      implicit none
      integer, intent(in) :: ntasks
#ifdef _MPIF90
      integer :: shared_comm, ierr, n
#endif

      ranks_sharing_node = 1

#ifdef _MPIF90
      call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                               MPI_INFO_NULL, shared_comm, ierr)
      if (ierr == MPI_SUCCESS) then
         call mpi_comm_size(shared_comm, n, ierr)
         if (ierr == MPI_SUCCESS .and. n > 0) ranks_sharing_node = n
         call mpi_comm_free(shared_comm, ierr)
      else
         !  A split that fails leaves us knowing nothing about the layout.
         !  ntasks is the pessimistic answer -- it is right on one node and
         !  merely conservative on several -- and being conservative here costs
         !  batches, while being wrong costs the run.
         ranks_sharing_node = max(1, ntasks)
      end if
#endif

   end function ranks_sharing_node
   !**************************************************************************

   !**************************************************************************
   !  The memory ceiling this process must live under, in bytes, and where the
   !  number came from. Zero means we could not find out.
   !
   !  cgroup first, /proc/meminfo second: see the note on gpu_memory_budget_init.
   subroutine host_memory_limit(bytes, source)
      implicit none
      integer(8), intent(out) :: bytes
      character(len=*), intent(out) :: source

      bytes = 0_8
      source = 'unknown'

      call cgroup_memory_limit(bytes)
      if (bytes > 0_8) then
         source = 'cgroup limit'
         return
      end if

      call meminfo_memtotal(bytes)
      if (bytes > 0_8) source = 'MemTotal'

   end subroutine host_memory_limit
   !**************************************************************************

   !**************************************************************************
   !  MemTotal from /proc/meminfo, in bytes. Zero if the file is not there,
   !  which is the normal answer anywhere that is not Linux.
   subroutine meminfo_memtotal(bytes)
      implicit none
      integer(8), intent(out) :: bytes

      integer :: unit, iostatus
      character(len=256) :: line
      character(len=32) :: label, units
      integer(8) :: value

      bytes = 0_8

      open (newunit=unit, file='/proc/meminfo', status='old', action='read', iostat=iostatus)
      if (iostatus /= 0) return

      do
         read (unit, '(A)', iostat=iostatus) line
         if (iostatus /= 0) exit
         if (line(1:9) /= 'MemTotal:') cycle
         read (line, *, iostat=iostatus) label, value, units
         if (iostatus == 0) bytes = value*1024_8
         exit
      end do

      close (unit)

   end subroutine meminfo_memtotal
   !**************************************************************************

   !**************************************************************************
   !  This process's cgroup memory limit in bytes, or zero if there is none.
   !
   !  The limit that applies is not the one at the root of the hierarchy: slurm
   !  puts each job step in a nested cgroup, and only some levels carry a limit
   !  ('max' means unlimited at that level). So we read our own cgroup path from
   !  /proc/self/cgroup and walk UP it, taking the first numeric limit we find,
   !  which is the tightest one that applies to us.
   !
   !  cgroup v2 only. v1's memory.limit_in_bytes reports a number near 2^63 for
   !  "no limit" rather than a word, and telling that apart from a real limit
   !  means knowing the page size and the kernel's rounding; every machine this
   !  code targets has been on v2 for years, and a wrong limit is worse than no
   !  limit because the run dies believing it had room.
   subroutine cgroup_memory_limit(bytes)
      implicit none
      integer(8), intent(out) :: bytes

      integer :: unit, iostatus, i, last
      character(len=512) :: line, path, fname

      bytes = 0_8

      !  /proc/self/cgroup: the v2 entry is the line beginning "0::".
      path = ''
      open (newunit=unit, file='/proc/self/cgroup', status='old', action='read', iostat=iostatus)
      if (iostatus /= 0) return
      do
         read (unit, '(A)', iostat=iostatus) line
         if (iostatus /= 0) exit
         if (line(1:3) == '0::') then
            path = adjustl(line(4:))
            exit
         end if
      end do
      close (unit)

      if (len_trim(path) == 0) return

      !  Walk up: /a/b/c, then /a/b, then /a, then the root.
      do
         fname = '/sys/fs/cgroup'//trim(path)//'/memory.max'
         call read_integer_file(fname, bytes)
         if (bytes > 0_8) return

         if (len_trim(path) == 0) exit
         last = 0
         do i = 1, len_trim(path)
            if (path(i:i) == '/') last = i
         end do
         if (last <= 1) then
            path = ''
         else
            path = path(1:last - 1)
         end if
         if (len_trim(path) == 0) then
            !  One last look at the root cgroup, then stop.
            call read_integer_file('/sys/fs/cgroup/memory.max', bytes)
            exit
         end if
      end do

   end subroutine cgroup_memory_limit
   !**************************************************************************

   !**************************************************************************
   !  Read a file holding a single integer. Zero if it is absent, unreadable, or
   !  holds 'max' -- which is how cgroup v2 spells "no limit at this level".
   subroutine read_integer_file(fname, value)
      implicit none
      character(len=*), intent(in) :: fname
      integer(8), intent(out) :: value

      integer :: unit, iostatus
      character(len=64) :: buf

      value = 0_8

      open (newunit=unit, file=trim(fname), status='old', action='read', iostat=iostatus)
      if (iostatus /= 0) return
      read (unit, '(A)', iostat=iostatus) buf
      close (unit)
      if (iostatus /= 0) return

      read (buf, *, iostat=iostatus) value
      if (iostatus /= 0) value = 0_8

   end subroutine read_integer_file
   !**************************************************************************

end module gpu_context
