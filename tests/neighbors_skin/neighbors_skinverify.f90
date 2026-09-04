program neighbors_skinverify

!  Standalone driver for src/neighbors_skin.f90.
!
!  The question the module has to answer is whether a list built once, with a
!  padded cutoff, still holds every pair inside the real cutoff some steps
!  later. That is checkable without the list machinery: build both sets by
!  brute force over all N^2 pairs and ask whether one contains the other. The
!  reference is therefore not this module, not the neighbours builder, and not
!  a stored baseline -- it is the definition of a neighbour.
!
!  Usage: neighbors_skinverify [seed]

   use kinds, only: dp
   use neighbors_skin, only: skin_accumulate, skin_needs_rebuild, skin_two_largest

   implicit none

   integer :: seed_arg
   integer :: nfail
   character(len=32) :: arg

   nfail = 0
   seed_arg = 20260904
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      read (arg, *) seed_arg
   end if
   call seed_rng(seed_arg)

   write (*, '(A)') ''
   write (*, '(A)') '==== 1. the two largest displacements, against a full selection ===='
   call check_two_largest(nfail)

   write (*, '(A)') ''
   write (*, '(A)') '==== 2. the accumulator under wrapping ===='
   call check_accumulate_mic(nfail)

   write (*, '(A)') ''
   write (*, '(A)') '==== 3. safety: no pair inside rcut escapes the padded list ===='
   write (*, '(A)') '     orthorhombic cell'
   call check_safety(.false., nfail)
   write (*, '(A)') '     triclinic cell'
   call check_safety(.true., nfail)

   write (*, '(A)') ''
   write (*, '(A)') '==== 4. the bound itself, on a pair closed by hand ===='
   call check_bound(nfail)

   write (*, '(A)') ''
   if (nfail == 0) then
      write (*, '(A)') 'neighbors_skin: all checks passed'
   else
      write (*, '(A,I0,A)') 'neighbors_skin: ', nfail, ' check(s) FAILED'
      stop 1
   end if

contains

!**************************************************************************
   subroutine seed_rng(s)

      implicit none

      integer, intent(in) :: s
      integer, allocatable :: seed(:)
      integer :: n
      integer :: i

      call random_seed(size=n)
      allocate (seed(1:n))
      do i = 1, n
         seed(i) = s + 977*i
      end do
      call random_seed(put=seed)
      deallocate (seed)

   end subroutine seed_rng
!**************************************************************************

!**************************************************************************
!  skin_two_largest selects in the square and never sorts. Compare against a
!  selection that does sort, on fields that include ties and a single mover.
!
   subroutine check_two_largest(nfail)

      implicit none

      integer, intent(inout) :: nfail
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: mag(:)
      real(dp) :: d1
      real(dp) :: d2
      real(dp) :: r1
      real(dp) :: r2
      real(dp) :: worst
      integer :: trial
      integer :: n
      integer :: i
      integer :: j

      worst = 0.d0
      do trial = 1, 200
         n = 1 + modulo(trial*37, 64)
         allocate (d(1:3, 1:n))
         allocate (mag(1:n))
         call random_number(d)
         d = d - 0.5d0
!        A tie and a lone large mover, so the branch that keeps the runner-up
!        is actually taken.
         if (n >= 4) then
            d(1:3, 2) = d(1:3, 1)
            d(1:3, n) = 10.d0*d(1:3, n)
         end if
         do i = 1, n
            mag(i) = dsqrt(d(1, i)**2 + d(2, i)**2 + d(3, i)**2)
         end do
!        Reference: full descending sort.
         do i = 1, n - 1
            do j = i + 1, n
               if (mag(j) > mag(i)) call swap(mag(i), mag(j))
            end do
         end do
         r1 = mag(1)
         if (n >= 2) then
            r2 = mag(2)
         else
            r2 = 0.d0
         end if
         call skin_two_largest(d, d1, d2)
         worst = max(worst, abs(d1 - r1), abs(d2 - r2))
         deallocate (d, mag)
      end do
      write (*, '(A,ES12.4)') '     worst absolute difference from the sorted answer: ', worst
      if (worst > 1.d-14) then
         write (*, '(A)') '     FAIL'
         nfail = nfail + 1
      else
         write (*, '(A)') '     PASS'
      end if

   end subroutine check_two_largest
!**************************************************************************

!**************************************************************************
!  An atom walking steadily across a periodic boundary. The accumulated
!  displacement must be its true path length, not a box length. The naive
!  difference is reported too: if it is not O(L) the test has stopped
!  exercising the wrap and is no longer discriminating.
!
   subroutine check_accumulate_mic(nfail)

      implicit none

      integer, intent(inout) :: nfail
      real(dp) :: cell(1:3, 1:3)
      real(dp) :: cell_inv(1:3, 1:3)
      real(dp) :: pos(1:3, 1:1)
      real(dp) :: prev(1:3, 1:1)
      real(dp) :: diff(1:3, 1:1)
      real(dp) :: naive(1:3)
      real(dp) :: step(1:3)
      real(dp) :: expect
      real(dp) :: err
      integer :: k
      integer :: nsteps

      cell = 0.d0
      cell(1, 1) = 12.d0
      cell(2, 2) = 13.d0
      cell(3, 3) = 14.d0
      call inv3(cell, cell_inv)

      step = (/0.07d0, -0.03d0, 0.11d0/)
      nsteps = 400
      pos(1:3, 1) = (/11.9d0, 0.1d0, 13.8d0/)
      diff = 0.d0
      naive = 0.d0
      do k = 1, nsteps
         prev(1:3, 1) = pos(1:3, 1)
         pos(1:3, 1) = pos(1:3, 1) + step(1:3)
         call wrap(pos, cell, cell_inv)
         naive(1:3) = naive(1:3) + pos(1:3, 1) - prev(1:3, 1)
         call skin_accumulate(pos, prev, cell, diff)
      end do

      expect = real(nsteps, dp)*dsqrt(step(1)**2 + step(2)**2 + step(3)**2)
      err = abs(dsqrt(diff(1, 1)**2 + diff(2, 1)**2 + diff(3, 1)**2) - expect)
      write (*, '(A,F10.4)') '     true path length:                     ', expect
      write (*, '(A,F10.4)') '     accumulated with the minimum image:   ', &
         dsqrt(diff(1, 1)**2 + diff(2, 1)**2 + diff(3, 1)**2)
      write (*, '(A,F10.4)') '     accumulated without it (the old way): ', &
         dsqrt(naive(1)**2 + naive(2)**2 + naive(3)**2)
      if (err > 1.d-10) then
         write (*, '(A,ES12.4)') '     FAIL: error ', err
         nfail = nfail + 1
      else if (abs(dsqrt(naive(1)**2 + naive(2)**2 + naive(3)**2) - expect) < 1.d0) then
         write (*, '(A)') '     FAIL: the naive accumulator agrees, so no wrap was crossed'
         nfail = nfail + 1
      else
         write (*, '(A)') '     PASS'
      end if

   end subroutine check_accumulate_mic
!**************************************************************************

!**************************************************************************
!  The safety check. Build the padded pair set once by brute force, then walk
!  the atoms. For as long as skin_needs_rebuild says the list still stands,
!  every pair inside rcut must be in that padded set.
!
!  The criterion the module replaces is walked along the same trajectory and
!  its losses are reported for contrast. They are not asserted on here: whether
!  a random walk happens to produce one is chance, and check 4 settles the same
!  question by construction instead.
!
   subroutine check_safety(triclinic, nfail)

      implicit none

      logical, intent(in) :: triclinic
      integer, intent(inout) :: nfail
      real(dp), parameter :: rcut = 4.0d0
      real(dp), parameter :: buffer = 0.8d0
      integer, parameter :: n = 250
      integer, parameter :: ntrial = 20
      integer, parameter :: nstage = 200
      real(dp) :: cell(1:3, 1:3)
      real(dp) :: cell_inv(1:3, 1:3)
      real(dp) :: pos(1:3, 1:n)
      real(dp) :: prev(1:3, 1:n)
      real(dp) :: diff(1:3, 1:n)
      real(dp) :: diff_old(1:3, 1:n)
      real(dp) :: frac(1:3, 1:n)
      real(dp) :: r(1:3, 1:n)
      real(dp) :: delta
      real(dp) :: mean_stages
      logical, allocatable :: padded(:, :)
      integer :: trial
      integer :: stage
      integer :: escapes_new
      integer :: escapes_old
      integer :: stages_new
      integer :: total_stages_new
      logical :: old_alive
      logical :: new_alive

      allocate (padded(1:n, 1:n))
      escapes_new = 0
      escapes_old = 0
      total_stages_new = 0

      do trial = 1, ntrial
         cell = 0.d0
         cell(1, 1) = 16.d0
         cell(2, 2) = 17.d0
         cell(3, 3) = 18.d0
         if (triclinic) then
            cell(1, 2) = 1.8d0
            cell(1, 3) = -1.2d0
            cell(2, 3) = 2.2d0
         end if
         call inv3(cell, cell_inv)

         call random_number(frac)
         pos = matmul(cell, frac)
         call pair_set(pos, cell, cell_inv, rcut + buffer, padded)

!        A displacement scale that puts the rebuild a handful of stages away,
!        swept across trials so the criterion is approached from both sides.
         delta = 0.03d0 + 0.08d0*real(trial, dp)/real(ntrial, dp)
         diff = 0.d0
         diff_old = 0.d0
         new_alive = .true.
         old_alive = .true.
         stages_new = 0

         do stage = 1, nstage
            prev = pos
            call random_number(r)
            pos = pos + delta*(r - 0.5d0)
            call wrap(pos, cell, cell_inv)
            call skin_accumulate(pos, prev, cell, diff)
            call skin_accumulate(pos, prev, cell, diff_old)

            if (new_alive) then
               new_alive = .not. skin_needs_rebuild(diff, buffer)
               if (new_alive) then
                  stages_new = stages_new + 1
                  if (.not. contained(pos, cell, cell_inv, rcut, padded)) escapes_new = escapes_new + 1
               end if
            end if
            if (old_alive) then
               old_alive = .not. old_needs_rebuild(diff_old, buffer)
               if (old_alive) then
                  if (.not. contained(pos, cell, cell_inv, rcut, padded)) escapes_old = escapes_old + 1
               end if
            end if
            if (.not. new_alive .and. .not. old_alive) exit
         end do
         total_stages_new = total_stages_new + stages_new
      end do

      mean_stages = real(total_stages_new, dp)/real(ntrial, dp)
      write (*, '(A,F8.2)') '       mean steps survived before a rebuild: ', mean_stages
      write (*, '(A,I0)') '       pairs lost while this criterion said the list stood: ', escapes_new
      write (*, '(A,I0)') '       pairs lost under the criterion it replaces:          ', escapes_old

      if (escapes_new > 0) then
         write (*, '(A)') '       FAIL: the skin let a pair inside rcut escape the list'
         nfail = nfail + 1
      else if (mean_stages < 2.d0) then
         write (*, '(A)') '       FAIL: rebuilds every step or two -- safe, but the skin is doing no work'
         nfail = nfail + 1
      else
         write (*, '(A)') '       PASS'
      end if
      deallocate (padded)

   end subroutine check_safety
!**************************************************************************

!**************************************************************************
!  Two atoms placed a hair outside the padded cutoff, then closed symmetrically
!  one small increment at a time. The pair is not in the list. The criterion
!  has to call for a rebuild while the pair is still outside rcut, because once
!  it is inside and the list has not been rebuilt the pair is simply absent and
!  nothing downstream can notice.
!
!  Swept over the buffer, which is what makes the defect legible rather than
!  anecdotal: the replaced criterion fires at a closing distance of
!  sqrt(2*buffer) instead of buffer, so it is late for every buffer below 2 A
!  and only accidentally safe above it.
!
   subroutine check_bound(nfail)

      implicit none

      integer, intent(inout) :: nfail
      real(dp), parameter :: rcut = 4.0d0
      real(dp), parameter :: dx = 1.0d-4
      real(dp) :: cell(1:3, 1:3)
      real(dp) :: pos(1:3, 1:2)
      real(dp) :: prev(1:3, 1:2)
      real(dp) :: diff(1:3, 1:2)
      real(dp) :: diff_old(1:3, 1:2)
      real(dp) :: buffer
      real(dp) :: sep
      real(dp) :: sep_new
      real(dp) :: sep_old
      logical :: fired_new
      logical :: fired_old
      integer :: ib
      integer :: k
      integer :: bad_new
      integer :: bad_old

      cell = 0.d0
      cell(1, 1) = 60.d0
      cell(2, 2) = 60.d0
      cell(3, 3) = 60.d0
      bad_new = 0
      bad_old = 0

      write (*, '(A)') '     buffer   separation when a rebuild is called      lost?'
      write (*, '(A)') '        (A)        this module   the one replaced'
      do ib = 1, 15
         buffer = 0.2d0*real(ib, dp)
         pos = 0.d0
         pos(1, 1) = 30.d0 - 0.5d0*(rcut + buffer) - 1.d-6
         pos(1, 2) = 30.d0 + 0.5d0*(rcut + buffer) + 1.d-6
         diff = 0.d0
         diff_old = 0.d0
         fired_new = .false.
         fired_old = .false.
         sep_new = -1.d0
         sep_old = -1.d0
         do k = 1, 40000
            prev = pos
            pos(1, 1) = pos(1, 1) + dx
            pos(1, 2) = pos(1, 2) - dx
            call skin_accumulate(pos, prev, cell, diff)
            call skin_accumulate(pos, prev, cell, diff_old)
            sep = pos(1, 2) - pos(1, 1)
            if (.not. fired_new) then
               if (skin_needs_rebuild(diff, buffer)) then
                  fired_new = .true.
                  sep_new = sep
               end if
            end if
            if (.not. fired_old) then
               if (old_needs_rebuild(diff_old, buffer)) then
                  fired_old = .true.
                  sep_old = sep
               end if
            end if
            if (fired_new .and. fired_old) exit
         end do
!        One increment of slack: the criterion is only consulted between steps.
         if (sep_new < rcut - 2.d0*dx) bad_new = bad_new + 1
         if (sep_old < rcut - 2.d0*dx) bad_old = bad_old + 1
         write (*, '(F11.2,F19.4,F19.4,A)') buffer, sep_new, sep_old, &
            trim(loss_tag(sep_new, sep_old, rcut - 2.d0*dx))
      end do

      write (*, '(A,F5.2,A)') '     a pair is lost below a separation of ', rcut, ' A'
      if (bad_new > 0) then
         write (*, '(A,I0,A)') '     FAIL: this module was late for ', bad_new, ' of the buffers'
         nfail = nfail + 1
      else if (bad_old == 0) then
         write (*, '(A)') '     FAIL: the replaced criterion was never late, so nothing here is being tested'
         nfail = nfail + 1
      else
         write (*, '(A,I0,A)') '     PASS (the replaced criterion was late for ', bad_old, ' of them)'
      end if

   end subroutine check_bound
!**************************************************************************

!**************************************************************************
   function loss_tag(sep_new, sep_old, limit) result(tag)

      implicit none

      real(dp), intent(in) :: sep_new
      real(dp), intent(in) :: sep_old
      real(dp), intent(in) :: limit
      character(len=24) :: tag

      if (sep_new < limit .and. sep_old < limit) then
         tag = '   both'
      else if (sep_new < limit) then
         tag = '   this module'
      else if (sep_old < limit) then
         tag = '   the one replaced'
      else
         tag = '   neither'
      end if

   end function loss_tag
!**************************************************************************

!**************************************************************************
!  The test this module replaces: a squared displacement against half the
!  buffer, unsquared.
!
   pure function old_needs_rebuild(positions_diff, buffer) result(rebuild)

      implicit none

      real(dp), intent(in) :: positions_diff(:, :)
      real(dp), intent(in) :: buffer
      logical :: rebuild
      integer :: i

      rebuild = .false.
      do i = 1, size(positions_diff, 2)
         if (positions_diff(1, i)**2 + positions_diff(2, i)**2 + positions_diff(3, i)**2 > buffer/2.d0) then
            rebuild = .true.
            return
         end if
      end do

   end function old_needs_rebuild
!**************************************************************************

!**************************************************************************
!  All pairs closer than r under the minimum image convention. O(N^2), which
!  is the point: it shares nothing with the code under test.
!
   subroutine pair_set(pos, cell, cell_inv, r, set)

      implicit none

      real(dp), intent(in) :: pos(:, :)
      real(dp), intent(in) :: cell(1:3, 1:3)
      real(dp), intent(in) :: cell_inv(1:3, 1:3)
      real(dp), intent(in) :: r
      logical, intent(out) :: set(:, :)
      integer :: i
      integer :: j

      set = .false.
      do i = 1, size(pos, 2)
         do j = i + 1, size(pos, 2)
            if (mic_dist(pos(1:3, i), pos(1:3, j), cell, cell_inv) < r) then
               set(i, j) = .true.
               set(j, i) = .true.
            end if
         end do
      end do

   end subroutine pair_set
!**************************************************************************

!**************************************************************************
   function contained(pos, cell, cell_inv, r, set) result(ok)

      implicit none

      real(dp), intent(in) :: pos(:, :)
      real(dp), intent(in) :: cell(1:3, 1:3)
      real(dp), intent(in) :: cell_inv(1:3, 1:3)
      real(dp), intent(in) :: r
      logical, intent(in) :: set(:, :)
      logical :: ok
      integer :: i
      integer :: j

      ok = .true.
      do i = 1, size(pos, 2)
         do j = i + 1, size(pos, 2)
            if (mic_dist(pos(1:3, i), pos(1:3, j), cell, cell_inv) < r .and. .not. set(i, j)) then
               ok = .false.
               return
            end if
         end do
      end do

   end function contained
!**************************************************************************

!**************************************************************************
   function mic_dist(pi, pj, cell, cell_inv) result(d)

      implicit none

      real(dp), intent(in) :: pi(1:3)
      real(dp), intent(in) :: pj(1:3)
      real(dp), intent(in) :: cell(1:3, 1:3)
      real(dp), intent(in) :: cell_inv(1:3, 1:3)
      real(dp) :: d
      real(dp) :: s(1:3)
      real(dp) :: v(1:3)
      real(dp) :: best
      real(dp) :: dsq
      integer :: ia
      integer :: ib
      integer :: ic

      s = matmul(cell_inv, pj - pi)
      s = s - dnint(s)
!     The nearest image of a triclinic cell is not always the one nint picks,
!     so the 27 surrounding translations are searched. Slow and correct.
      best = huge(1.d0)
      do ia = -1, 1
         do ib = -1, 1
            do ic = -1, 1
               v = matmul(cell, s + (/real(ia, dp), real(ib, dp), real(ic, dp)/))
               dsq = v(1)**2 + v(2)**2 + v(3)**2
               if (dsq < best) best = dsq
            end do
         end do
      end do
      d = dsqrt(best)

   end function mic_dist
!**************************************************************************

!**************************************************************************
   subroutine wrap(pos, cell, cell_inv)

      implicit none

      real(dp), intent(inout) :: pos(:, :)
      real(dp), intent(in) :: cell(1:3, 1:3)
      real(dp), intent(in) :: cell_inv(1:3, 1:3)
      real(dp) :: s(1:3)
      integer :: i

      do i = 1, size(pos, 2)
         s = matmul(cell_inv, pos(1:3, i))
         s = s - real(floor(s), dp)
         pos(1:3, i) = matmul(cell, s)
      end do

   end subroutine wrap
!**************************************************************************

!**************************************************************************
   subroutine inv3(a, ainv)

      implicit none

      real(dp), intent(in) :: a(1:3, 1:3)
      real(dp), intent(out) :: ainv(1:3, 1:3)
      real(dp) :: det

      ainv(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
      ainv(1, 2) = a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3)
      ainv(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
      ainv(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
      ainv(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
      ainv(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
      ainv(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
      ainv(3, 2) = a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2)
      ainv(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
      det = a(1, 1)*ainv(1, 1) + a(1, 2)*ainv(2, 1) + a(1, 3)*ainv(3, 1)
      ainv = ainv/det

   end subroutine inv3
!**************************************************************************

!**************************************************************************
   subroutine swap(a, b)

      implicit none

      real(dp), intent(inout) :: a
      real(dp), intent(inout) :: b
      real(dp) :: t

      t = a
      a = b
      b = t

   end subroutine swap
!**************************************************************************

end program neighbors_skinverify
