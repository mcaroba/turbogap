! TurboGAP as an i-PI force engine.
!
! In `turbogap ipi` the nuclei belong to i-PI. It owns the beads, the
! thermostat and the integrator; TurboGAP owns the potential. That division is
! what puts path integrals -- and therefore PIGLET, and therefore nuclear
! quantum effects on an OH stretch -- within reach of a GAP without a line of
! ring-polymer code in this tree.
!
! WHERE THIS SITS
!
! src/turbogap.f90 calls compute_md once per pass of the main loop, after the
! forces for the current positions exist. This routine replaces that one call
! and nothing else: same loop, same neighbour lists, same GAP evaluation, and
! in place of an integrator a socket. Everything compute_md does around the
! integration -- the Verlet-skin displacement accounting, the supercell image
! refresh, the broadcast of positions to the other ranks -- has to happen here
! too, and does, in the same order, for the same reasons.
!
! THE ONE-STEP OFFSET
!
! TurboGAP computes forces at the top of an iteration and this routine runs at
! the bottom, so on the FIRST pass there are forces in hand that i-PI never
! asked for: they belong to the configuration in atoms_file. Answering STATUS
! with HAVEDATA there would hand i-PI forces for the wrong geometry, and
! nothing downstream could tell -- the dynamics would simply be wrong by one
! frame's worth of force. So `n_exchanged` gates it: the first pass reports
! READY, takes i-PI's coordinates and returns, and that first GAP evaluation
! is discarded. It costs one force evaluation per run.
!
! UNITS
!
! i-PI is atomic units throughout -- bohr, Hartree, Hartree/bohr -- and
! TurboGAP is eV and Angstrom. Conversions are applied at the socket boundary
! only, so nothing inside TurboGAP sees them.
!
! THE VIRIAL
!
! i-PI wants -dE/d(strain) in Hartree. TurboGAP's `virial` is defined by the
! pressure it forms from it in turbogap_md.f90, P = (N kB T + tr(virial)/3)/V.
! i-PI forms its own as p = tr(vir + kstress)/(3 V) with tr(kstress) = 2 K =
! 3 N kB T, which is the same expression. So the two conventions agree in sign
! and in scale, and the conversion is the eV-to-Hartree factor and nothing
! else. Checked against i-PI 3.3 rather than assumed: a sign error here would
! be invisible in NVT, where the virial reaches nothing but the reported
! pressure, and would silently corrupt any run with a barostat. It is transposed on the way out to match i-PI's index order;
! for every potential in this tree the tensor is symmetric and the transpose is
! a no-op, but a future non-symmetric one would need it.
!
! WHAT IS NOT SUPPORTED
!
! The cell is taken from i-PI on every POSDATA, so a barostat on i-PI's side
! works, but a cell change forces a full neighbour-list rebuild. GETSTRESSES
! (i-PI's per-atom stress extension) is not implemented and is refused rather
! than answered with zeros.
module ipi_driver

   use kinds, only: dp
   use error, only: turbogap_abort
   use ipi_socket, only: ipi_connect, ipi_disconnect, ipi_get_header, ipi_put_header, &
                         ipi_get_int, ipi_put_int, ipi_get_reals, ipi_put_reals, &
                         ipi_put_real, ipi_skip_bytes, IPI_MSGLEN
   use neighbors_skin, only: skin_accumulate, skin_needs_rebuild
   use md, only: wrap_pbc
#ifdef _MPIF90
   use mpi
#endif

   implicit none

   private
   public :: ipi_driver_open
   public :: ipi_driver_exchange
   public :: ipi_driver_close

   real(dp), parameter :: BOHR_ANG = 0.5291772109_dp
   real(dp), parameter :: HARTREE_EV = 27.2113862460_dp

!  Rank 0 alone holds the connection; the others learn what came over it by
!  broadcast, exactly as they do for an MD step.
   integer, save :: fd = -1
!  Wall clock at the end of the last exchange, so the reported cost is the one
!  that matters -- how long i-PI waited for this client, force call included.
   real(dp), save :: t_last = -1.0_dp
   real(dp), save :: t_sum = 0.0_dp
!  Has i-PI sent coordinates yet? Gates the STATUS answer, see the header.
   integer, save :: n_exchanged = 0
!  i-PI asks for INIT once per replica and expects NEEDINIT until it has.
   logical, save :: isinit = .false.

contains

   subroutine ipi_driver_open(address, rank)

      implicit none

      character(len=*), intent(in) :: address
      integer, intent(in) :: rank

      if (rank /= 0) return
      if (len_trim(address) == 0) then
         write (*, *) "ERROR: turbogap ipi needs ipi_address in the input file, e.g."
         write (*, *) "       ipi_address = 'UNIX:turbogap'   or   ipi_address = 'localhost:31415'"
         call turbogap_abort()
      end if
      call ipi_connect(trim(address), fd)
      write (*, *) '                                       |'
      write (*, '(1X,A,A)') 'i-PI driver connected to ', trim(address)
      write (*, *) '                                       |'

   end subroutine ipi_driver_open

   subroutine ipi_driver_close(rank)

      implicit none

      integer, intent(in) :: rank

      if (rank /= 0) return
      call ipi_disconnect(fd)

   end subroutine ipi_driver_close

!  Serve i-PI until it has taken the forces in hand and given back the next
!  set of coordinates, or until it says EXIT.
!
!  On return either exit_loop is set, or positions/a_box/b_box/c_box hold
!  i-PI's next configuration and the caller's next pass evaluates it.
   subroutine ipi_driver_exchange(rank, n_sites, positions, positions_prev, positions_diff, &
                                  velocities, a_box, b_box, c_box, indices, neighbors_buffer, &
                                  forces, energy, virial, exit_loop, rebuild_neighbors_list)

      implicit none

      integer, intent(in) :: rank
      integer, intent(in) :: n_sites
      real(dp), intent(inout) :: positions(:, :)
      real(dp), intent(inout) :: positions_prev(:, :)
      real(dp), intent(inout) :: positions_diff(:, :)
      real(dp), intent(inout) :: velocities(:, :)
      real(dp), intent(inout) :: a_box(1:3)
      real(dp), intent(inout) :: b_box(1:3)
      real(dp), intent(inout) :: c_box(1:3)
      integer, intent(in) :: indices(1:3)
      real(dp), intent(in) :: neighbors_buffer
      real(dp), intent(in) :: forces(:, :)
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: virial(1:3, 1:3)
      logical, intent(inout) :: exit_loop
      logical, intent(inout) :: rebuild_neighbors_list

      character(len=IPI_MSGLEN) :: header
      real(dp) :: cell(1:3, 1:3)
      real(dp) :: cellbuf(1:9)
      real(dp) :: a_new(1:3)
      real(dp) :: b_new(1:3)
      real(dp) :: c_new(1:3)
      real(dp), allocatable :: posbuf(:)
      real(dp), allocatable :: fbuf(:)
      logical :: hasdata
      logical :: cell_changed
      integer :: nat
      integer :: nchar
      integer :: ibead
      integer :: i
      integer :: i2
      integer :: j
      integer :: j2
      integer :: k2
      integer :: n_pos
      integer :: ierr

      a_new = a_box
      b_new = b_box
      c_new = c_box
      cell_changed = .false.

#ifdef _MPIF90
      IF (rank == 0) THEN
#endif
         hasdata = (n_exchanged > 0)
         allocate (posbuf(1:3*n_sites))
         allocate (fbuf(1:3*n_sites))

         do
            call ipi_get_header(fd, header)

            select case (trim(header))

            case ("STATUS")
               if (hasdata) then
                  call ipi_put_header(fd, "HAVEDATA")
               else if (isinit) then
                  call ipi_put_header(fd, "READY")
               else
                  call ipi_put_header(fd, "NEEDINIT")
               end if

            case ("INIT")
               call ipi_get_int(fd, ibead)
               call ipi_get_int(fd, nchar)
               call ipi_skip_bytes(fd, nchar)
               isinit = .true.

            case ("GETFORCE")
               do i = 1, n_sites
                  fbuf(3*(i - 1) + 1:3*i) = forces(1:3, i)*BOHR_ANG/HARTREE_EV
               end do
               call ipi_put_header(fd, "FORCEREADY")
               call ipi_put_real(fd, energy/HARTREE_EV)
               call ipi_put_int(fd, n_sites)
               call ipi_put_reals(fd, fbuf)
               call ipi_put_reals(fd, reshape(transpose(virial)/HARTREE_EV, [9]))
!              i-PI reads a length and then that many characters of free-form
!              "extra" output. Nothing here produces any, but the length must
!              still be on the wire or the next header is read from the middle
!              of this message.
               call ipi_put_int(fd, 0)
               hasdata = .false.

            case ("POSDATA")
               call ipi_get_reals(fd, cellbuf)
               cell = reshape(cellbuf, [3, 3])
!              i-PI stores the lattice vectors as the COLUMNS of h and sends h
!              row-major, so a column-major reshape here puts lattice vector i
!              in cell(i,:). Reading it the other way round is silent for a
!              cubic box and wrong for every other one.
!              i-PI's cell is the PRIMITIVE one. TurboGAP's a_box, b_box and
!              c_box are the SUPERCELL: when the box cannot hold one cutoff
!              sphere read_xyz replicates it by `indices` and scales the
!              lattice vectors to match, and every consumer downstream --
!              build_neighbors_list, wrap_pbc, the image refresh below --
!              divides by `indices` to get back to the primitive cell.
!              Assigning i-PI's cell here without scaling puts a primitive
!              vector where a supercell one belongs: the images are then laid
!              out at half spacing on top of the real atoms, the neighbour
!              lists explode, and the run does not fail, it grinds and then
!              dies. Invisible whenever indices is 1, which is why the socket
!              test carries a case where it is not.
               a_new = cell(1, 1:3)*BOHR_ANG*dfloat(indices(1))
               b_new = cell(2, 1:3)*BOHR_ANG*dfloat(indices(2))
               c_new = cell(3, 1:3)*BOHR_ANG*dfloat(indices(3))
!              The inverse cell, which i-PI sends so a client need not invert
!              it. TurboGAP does not use it.
               call ipi_get_reals(fd, cellbuf)
               call ipi_get_int(fd, nat)
               if (nat /= n_sites) then
                  write (*, *) "ERROR: i-PI sent", nat, "atoms and TurboGAP was set up for", n_sites
                  write (*, *) "       atoms_file and i-PI's structure must be the same system."
                  call turbogap_abort()
               end if
               call ipi_get_reals(fd, posbuf)
               positions_prev(1:3, 1:n_sites) = positions(1:3, 1:n_sites)
               do i = 1, n_sites
                  positions(1:3, i) = posbuf(3*(i - 1) + 1:3*i)*BOHR_ANG
               end do
               cell_changed = (maxval(abs(a_new - a_box)) > 0.d0 .or. &
                               maxval(abs(b_new - b_box)) > 0.d0 .or. &
                               maxval(abs(c_new - c_box)) > 0.d0)
               a_box = a_new
               b_box = b_new
               c_box = c_new
               n_exchanged = n_exchanged + 1
               if (n_exchanged == 1) call report_first_geometry(n_sites, positions, a_box, b_box, c_box, indices)
               call report_rate()
               exit

            case ("EXIT")
               exit_loop = .true.
               exit

            case ("GETSTRESSES")
               write (*, *) "ERROR: i-PI asked for per-atom stresses (GETSTRESSES) and this"
               write (*, *) "       driver does not compute them. Turn off the atomic-stress"
               write (*, *) "       or heat-flux output on the i-PI side."
               call turbogap_abort()

            case default
               write (*, *) "ERROR: unknown message from the i-PI server: '", trim(header), "'"
               write (*, *) "       The stream is out of step; the run cannot continue."
               call turbogap_abort()

            end select
         end do

         deallocate (posbuf)
         deallocate (fbuf)

         if (.not. exit_loop) then
!           i-PI is free to send coordinates that have left the box -- with
!           pbc='false' on its ffsocket it always does -- and every MD step in
!           this tree works with positions wrapped into the primitive cell.
!           get_distance is minimum-image on both its branches so unwrapped
!           input would give the same forces, but the two paths would then hold
!           different numbers for the same atom, and the skin accounting and
!           the supercell images below are written for the wrapped convention.
            call wrap_pbc(positions(1:3, 1:n_sites), a_box/dfloat(indices(1)), &
                          b_box/dfloat(indices(2)), c_box/dfloat(indices(3)))
!           Everything below mirrors the tail of compute_md. The Verlet skin
!           accounts displacement since the last neighbour build under the
!           primitive minimum image, so an atom crossing a boundary contributes
!           its step and not a box length.
            call skin_accumulate(positions(1:3, 1:n_sites), positions_prev(1:3, 1:n_sites), &
                                 reshape([a_box/dfloat(indices(1)), b_box/dfloat(indices(2)), &
                                          c_box/dfloat(indices(3))], [3, 3]), positions_diff)
            rebuild_neighbors_list = .false.
!           A supercell is rebuilt every step for the wrapping reason given in
!           compute_md, and a cell that moved invalidates the list outright.
            if (any(indices > 1)) rebuild_neighbors_list = .true.
            if (cell_changed) rebuild_neighbors_list = .true.
            if (skin_needs_rebuild(positions_diff, neighbors_buffer)) then
               rebuild_neighbors_list = .true.
               positions_diff = 0.d0
            end if

!           Refresh the periodic images of the atoms i-PI just moved.
            j = 0
            do i2 = 1, indices(1)
               do j2 = 1, indices(2)
                  do k2 = 1, indices(3)
                     do i = 1, n_sites
                        j = j + 1
                        if (j > n_sites) then
                           positions(1:3, j) = positions(1:3, i) + dfloat(i2 - 1)/dfloat(indices(1))*a_box &
                                               + dfloat(j2 - 1)/dfloat(indices(2))*b_box &
                                               + dfloat(k2 - 1)/dfloat(indices(3))*c_box
                           velocities(1:3, j) = velocities(1:3, i)
                        end if
                     end do
                  end do
               end do
            end do
         end if

#ifdef _MPIF90
      END IF

      call mpi_bcast(exit_loop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if (.not. exit_loop) then
         n_pos = size(positions, 2)
         call mpi_bcast(positions, 3*n_pos, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call mpi_bcast(c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if
#endif

   end subroutine ipi_driver_exchange

!  What i-PI actually handed over on the first call.
!
!  A units mismatch between i-PI's structure file and this client does not
!  fail: i-PI reads its xyz in whatever units the comment line claims, and if
!  it reads Angstrom as bohr the configuration arrives 1.9x too dense. The GAP
!  still evaluates it. What the user sees is a run that is inexplicably slow --
!  the neighbour lists are enormous -- and a trajectory that is nonsense. The
!  closest approach is the number that says so at once: liquid water is 0.95 A
!  O-H and no condensed phase has a first neighbour below about 0.7 A.
   subroutine report_first_geometry(n_sites, positions, a_box, b_box, c_box, indices)

      implicit none

      integer, intent(in) :: n_sites
      real(dp), intent(in) :: positions(:, :)
      real(dp), intent(in) :: a_box(1:3)
      real(dp), intent(in) :: b_box(1:3)
      real(dp), intent(in) :: c_box(1:3)
      integer, intent(in) :: indices(1:3)
      real(dp) :: d_min
      real(dp) :: d
      real(dp) :: dr(1:3)
      real(dp) :: l(1:3)
      integer :: i
      integer :: j
      integer :: k

      l = [dsqrt(dot_product(a_box, a_box))/dfloat(indices(1)), &
           dsqrt(dot_product(b_box, b_box))/dfloat(indices(2)), &
           dsqrt(dot_product(c_box, c_box))/dfloat(indices(3))]
      d_min = huge(1.0_dp)
!     Orthorhombic minimum image, which is all this diagnostic needs: on a
!     triclinic cell it overestimates the closest approach and so can only
!     fail to warn, never warn wrongly.
      do i = 1, n_sites - 1
         do j = i + 1, n_sites
            dr = positions(1:3, j) - positions(1:3, i)
            do k = 1, 3
               if (l(k) > 0.0_dp) dr(k) = dr(k) - l(k)*dnint(dr(k)/l(k))
            end do
            d = dsqrt(dot_product(dr, dr))
            if (d < d_min) d_min = d
         end do
      end do

      write (*, '(1X,A,3F10.4,A)') "i-PI cell:", l, " A (primitive)"
      write (*, '(1X,A,F10.4,A)') "i-PI closest approach:", d_min, " A"
      if (d_min < 0.5_dp) then
         write (*, *) "WARNING: that is closer than any chemical bond. The most likely"
         write (*, *) "         cause is a units mismatch on i-PI's side -- an xyz whose"
         write (*, *) "         comment line does not say positions{angstrom} is read as"
         write (*, *) "         bohr. The run will be slow and the trajectory meaningless."
      end if
      flush (6)

   end subroutine report_first_geometry

!  How many force calls this client has served, and what they cost.
!
!  The first few and then every hundredth: a driver that is connected but
!  never asked for anything looks exactly like one that is working, and over
!  a run of 40000 steps there is nothing else that says which. The rate is
!  also the only honest way to size a run before committing to it.
   subroutine report_rate()

      implicit none

      real(dp) :: t_now
      integer :: count
      integer :: count_rate

      call system_clock(count, count_rate)
      t_now = dfloat(count)/dfloat(count_rate)
      if (t_last > 0.0_dp) t_sum = t_sum + (t_now - t_last)
      if (n_exchanged <= 5 .or. modulo(n_exchanged, 100) == 0) then
         if (n_exchanged > 1) then
            write (*, '(1X,A,I8,A,F8.3,A,F8.3,A)') "i-PI force call", n_exchanged, &
               ":", t_now - t_last, " s, mean", t_sum/dfloat(n_exchanged - 1), " s"
         else
            write (*, '(1X,A,I8,A)') "i-PI force call", n_exchanged, ": first"
         end if
         flush (6)
      end if
      t_last = t_now

   end subroutine report_rate

end module ipi_driver
