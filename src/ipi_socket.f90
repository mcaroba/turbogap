! The client half of the i-PI socket protocol, in pure Fortran.
!
! i-PI (Kapil et al., Comput. Phys. Commun. 236 (2019) 214) drives a force
! engine over a socket: i-PI owns the nuclei -- ring-polymer beads, the
! thermostat, the integrator -- and the engine answers one question, "what are
! the forces on these coordinates". That is what makes path integrals, and
! therefore PIGLET, available to TurboGAP without a line of path-integral code
! in TurboGAP.
!
! WHY THIS IS FORTRAN AND NOT THE C SHIM EVERYONE ELSE VENDORS
!
! The reference implementation is a ~150-line sockets.c that FHI-aims, CP2K,
! LAMMPS and QUIP each carry a copy of. Linking it here would mean a C
! compiler, a C rule in the Makefile, and a CC in every one of the
! makefiles/Makefile.* architecture files -- a cross-cutting build change for
! four libc calls. iso_c_binding reaches socket(), connect(), read(), write()
! and close() directly, so this module adds a file to SRC and nothing else.
!
! The cost is that the two socket address structures are packed here by hand.
! They are frozen ABI on Linux and are the only platform assumption in the
! file:
!
!   struct sockaddr_un   16-bit family, then 108 bytes of NUL-padded path
!   struct sockaddr_in   16-bit family, 16-bit port and 32-bit address, both
!                        big-endian on the wire, then 8 unused bytes
!
! THE PROTOCOL, as the reference driver implements it
!
! Every message is a 12-byte, space-padded ASCII header, and the server speaks
! first. The client answers STATUS with READY (waiting for coordinates),
! HAVEDATA (forces computed, not yet collected) or NEEDINIT; receives geometry
! on POSDATA; returns energy, forces and virial on GETFORCE; and stops on EXIT.
! Payloads are raw little-endian doubles and 32-bit integers with no padding,
! in ATOMIC UNITS -- bohr, Hartree, Hartree/bohr -- which is where the unit
! conversions in ipi_driver come from.
!
! A stream socket is free to return a short read, and does whenever a message
! spans a TCP segment: 3*648 doubles is 15552 bytes and the loopback MTU is
! 65536, so a short read is rare on one host and routine on two. Every
! transfer here loops until the byte count is met. Getting that wrong gives a
! run that works for hours and then desynchronises the stream, which presents
! as a garbage header rather than as an I/O error.
!
! This module knows nothing about TurboGAP beyond kinds and error, which is
! what lets tests/ipi_socket compile it against a reference server on its own.
module ipi_socket

   use, intrinsic :: iso_c_binding, only: c_int, c_long, c_size_t, c_ptr, c_loc, c_signed_char
   use kinds, only: dp
   use error, only: turbogap_abort

   implicit none

   private
   public :: ipi_connect
   public :: ipi_disconnect
   public :: ipi_get_header
   public :: ipi_put_header
   public :: ipi_get_int
   public :: ipi_put_int
   public :: ipi_get_reals
   public :: ipi_put_reals
   public :: ipi_put_real
   public :: ipi_skip_bytes
   public :: IPI_MSGLEN

!  Every header on the wire is exactly this long, space padded.
   integer, parameter :: IPI_MSGLEN = 12

!  <sys/socket.h>. AF_UNIX is 1 and AF_INET is 2 on Linux; SOCK_STREAM is 1.
   integer, parameter :: AF_UNIX = 1
   integer, parameter :: AF_INET = 2
   integer, parameter :: SOCK_STREAM = 1

!  sizeof(struct sockaddr_un) and sizeof(struct sockaddr_in).
   integer, parameter :: SA_UN_LEN = 110
   integer, parameter :: SA_IN_LEN = 16

!  i-PI's server puts its UNIX-domain sockets here. The name is not
!  configurable in i-PI either, so it is not configurable here.
   character(len=*), parameter :: UNIX_PREFIX = "/tmp/ipi_"

   interface
      integer(c_int) function c_socket(domain, typ, protocol) bind(C, name="socket")
         import :: c_int
         integer(c_int), value :: domain
         integer(c_int), value :: typ
         integer(c_int), value :: protocol
      end function c_socket

      integer(c_int) function c_connect(fd, addr, addrlen) bind(C, name="connect")
         import :: c_int, c_ptr
         integer(c_int), value :: fd
         type(c_ptr), value :: addr
         integer(c_int), value :: addrlen
      end function c_connect

      integer(c_long) function c_read(fd, buf, n) bind(C, name="read")
         import :: c_int, c_long, c_ptr, c_size_t
         integer(c_int), value :: fd
         type(c_ptr), value :: buf
         integer(c_size_t), value :: n
      end function c_read

      integer(c_long) function c_write(fd, buf, n) bind(C, name="write")
         import :: c_int, c_long, c_ptr, c_size_t
         integer(c_int), value :: fd
         type(c_ptr), value :: buf
         integer(c_size_t), value :: n
      end function c_write

      integer(c_int) function c_close(fd) bind(C, name="close")
         import :: c_int
         integer(c_int), value :: fd
      end function c_close
   end interface

contains

!  Signed-char storage holds -128..127, and an IP octet or a high byte of a
!  port does not. Same bit pattern either way; this is the wrap C does for
!  free on assignment to a char.
   pure function to_byte(v) result(b)

      implicit none

      integer, intent(in) :: v
      integer(c_signed_char) :: b

      if (v > 127) then
         b = int(v - 256, c_signed_char)
      else
         b = int(v, c_signed_char)
      end if

   end function to_byte

!  One direction of one transfer, looping until nbytes have moved. A stream
!  socket may satisfy a request partially; a return of zero means the peer
!  closed, and negative means an error, and both are fatal here because there
!  is no way to resynchronise a half-read message.
   subroutine sock_xfer(fd, buf, nbytes, writing)

      implicit none

      integer, intent(in) :: fd
      integer(c_signed_char), intent(inout), target :: buf(:)
      integer, intent(in) :: nbytes
      logical, intent(in) :: writing
      integer(c_long) :: done
      integer(c_long) :: n

      done = 0_c_long
      do while (done < int(nbytes, c_long))
         if (writing) then
            n = c_write(int(fd, c_int), c_loc(buf(done + 1)), int(nbytes - done, c_size_t))
         else
            n = c_read(int(fd, c_int), c_loc(buf(done + 1)), int(nbytes - done, c_size_t))
         end if
         if (n <= 0_c_long) then
            if (writing) then
               write (*, *) "ERROR: the i-PI socket closed while sending. The i-PI server"
            else
               write (*, *) "ERROR: the i-PI socket closed while receiving. The i-PI server"
            end if
            write (*, *) "       exited or was killed; look at its log, not at this one."
            call turbogap_abort()
         end if
         done = done + n
      end do

   end subroutine sock_xfer

!  Connect to an i-PI server. `address` is either
!
!     UNIX:name          a UNIX-domain socket at /tmp/ipi_name
!     host:port          a TCP socket
!
!  matching what i-PI's own <address>/<port> pair and the reference driver's
!  -a/-p options produce.
!
!  A hostname is NOT resolved. gethostbyname would be a sixth libc binding
!  carrying a struct hostent, and every use of this is either a UNIX socket or
!  loopback; "localhost" is special-cased and anything else must be a dotted
!  quad. The error says so.
   subroutine ipi_connect(address, fd)

      implicit none

      character(len=*), intent(in) :: address
      integer, intent(out) :: fd
      integer(c_signed_char), target :: sa(1:SA_UN_LEN)
      character(len=256) :: host
      character(len=256) :: path
      integer :: domain
      integer :: salen
      integer :: port
      integer :: icolon
      integer :: iostatus
      integer :: ip(1:4)
      integer :: i

      sa = 0_c_signed_char
      port = 0

      if (len_trim(address) >= 5 .and. (address(1:5) == "UNIX:" .or. address(1:5) == "unix:")) then
         domain = AF_UNIX
         salen = SA_UN_LEN
         path = UNIX_PREFIX//trim(adjustl(address(6:)))
         if (len_trim(path) > SA_UN_LEN - 3) then
            write (*, *) "ERROR: the i-PI UNIX socket path is longer than 107 characters:"
            write (*, *) "       ", trim(path)
            call turbogap_abort()
         end if
         sa(1) = to_byte(AF_UNIX)
         sa(2) = 0_c_signed_char
         do i = 1, len_trim(path)
            sa(2 + i) = to_byte(iachar(path(i:i)))
         end do
      else
         domain = AF_INET
         salen = SA_IN_LEN
         icolon = index(address, ":", back=.true.)
         if (icolon < 2) then
            write (*, *) "ERROR: ipi_address = ", trim(address)
            write (*, *) "       Expected UNIX:name for a UNIX socket, or host:port for TCP."
            call turbogap_abort()
         end if
         host = adjustl(address(1:icolon - 1))
         read (address(icolon + 1:), *, iostat=iostatus) port
         if (iostatus /= 0 .or. port <= 0 .or. port > 65535) then
            write (*, *) "ERROR: ipi_address = ", trim(address)
            write (*, *) "       The text after the last colon is not a port number."
            call turbogap_abort()
         end if
         if (trim(host) == "localhost") then
            ip = [127, 0, 0, 1]
         else
            call parse_dotted_quad(trim(host), ip)
         end if
         sa(1) = to_byte(AF_INET)
         sa(2) = 0_c_signed_char
!        Port and address travel big-endian regardless of the host's byte order.
         sa(3) = to_byte(port/256)
         sa(4) = to_byte(modulo(port, 256))
         do i = 1, 4
            sa(4 + i) = to_byte(ip(i))
         end do
      end if

      fd = int(c_socket(int(domain, c_int), int(SOCK_STREAM, c_int), 0_c_int))
      if (fd < 0) then
         write (*, *) "ERROR: could not create a socket for the i-PI connection."
         call turbogap_abort()
      end if

      if (c_connect(int(fd, c_int), c_loc(sa(1)), int(salen, c_int)) < 0_c_int) then
         write (*, *) "ERROR: could not connect to the i-PI server at ", trim(address)
         if (domain == AF_UNIX) then
            write (*, *) "       Looked for ", trim(UNIX_PREFIX)//trim(adjustl(address(6:)))
         end if
         write (*, *) "       Start i-PI first and let it create the socket, then start"
         write (*, *) "       TurboGAP. There is no retry: a driver that outlives its"
         write (*, *) "       server would hang rather than fail."
         call turbogap_abort()
      end if

   end subroutine ipi_connect

   subroutine parse_dotted_quad(host, ip)

      implicit none

      character(len=*), intent(in) :: host
      integer, intent(out) :: ip(1:4)
      integer :: i
      integer :: j
      integer :: k
      integer :: iostatus
      character(len=64) :: field

      j = 1
      do i = 1, 4
         if (i < 4) then
            k = index(host(j:), ".")
            if (k == 0) then
               write (*, *) "ERROR: ipi_address host ", trim(host), " is neither localhost nor"
               write (*, *) "       a dotted-quad IP address. Hostnames are not resolved."
               call turbogap_abort()
            end if
            field = host(j:j + k - 2)
            j = j + k
         else
            field = host(j:)
         end if
         read (field, *, iostat=iostatus) ip(i)
         if (iostatus /= 0 .or. ip(i) < 0 .or. ip(i) > 255) then
            write (*, *) "ERROR: ipi_address host ", trim(host), " is not a valid IP address."
            call turbogap_abort()
         end if
      end do

   end subroutine parse_dotted_quad

   subroutine ipi_disconnect(fd)

      implicit none

      integer, intent(inout) :: fd
      integer(c_int) :: ignored

      if (fd >= 0) ignored = c_close(int(fd, c_int))
      fd = -1

   end subroutine ipi_disconnect

   subroutine ipi_get_header(fd, header)

      implicit none

      integer, intent(in) :: fd
      character(len=IPI_MSGLEN), intent(out) :: header
      integer(c_signed_char) :: buf(1:IPI_MSGLEN)
      integer :: i

      call sock_xfer(fd, buf, IPI_MSGLEN, .false.)
      do i = 1, IPI_MSGLEN
         header(i:i) = achar(iand(int(buf(i)), 255))
      end do

   end subroutine ipi_get_header

!  Headers are fixed width and space padded. A header sent short desynchronises
!  the stream in a way that surfaces several messages later.
   subroutine ipi_put_header(fd, header)

      implicit none

      integer, intent(in) :: fd
      character(len=*), intent(in) :: header
      integer(c_signed_char) :: buf(1:IPI_MSGLEN)
      integer :: i

      buf = to_byte(iachar(" "))
      do i = 1, min(len_trim(header), IPI_MSGLEN)
         buf(i) = to_byte(iachar(header(i:i)))
      end do
      call sock_xfer(fd, buf, IPI_MSGLEN, .true.)

   end subroutine ipi_put_header

!  i-PI's integers on the wire are 32 bit, whatever the client's default
!  integer kind happens to be.
   subroutine ipi_get_int(fd, v)

      implicit none

      integer, intent(in) :: fd
      integer, intent(out) :: v
      integer(c_signed_char) :: buf(1:4)
      integer(kind=4) :: tmp

      call sock_xfer(fd, buf, 4, .false.)
      tmp = transfer(buf, tmp)
      v = int(tmp)

   end subroutine ipi_get_int

   subroutine ipi_put_int(fd, v)

      implicit none

      integer, intent(in) :: fd
      integer, intent(in) :: v
      integer(c_signed_char) :: buf(1:4)
      integer(kind=4) :: tmp

      tmp = int(v, kind=4)
      buf = transfer(tmp, buf)
      call sock_xfer(fd, buf, 4, .true.)

   end subroutine ipi_put_int

   subroutine ipi_get_reals(fd, x)

      implicit none

      integer, intent(in) :: fd
      real(dp), intent(out) :: x(:)
      integer(c_signed_char), allocatable :: buf(:)

      allocate (buf(1:8*size(x)))
      call sock_xfer(fd, buf, 8*size(x), .false.)
      x = transfer(buf, x)
      deallocate (buf)

   end subroutine ipi_get_reals

   subroutine ipi_put_reals(fd, x)

      implicit none

      integer, intent(in) :: fd
      real(dp), intent(in) :: x(:)
      integer(c_signed_char), allocatable :: buf(:)

      allocate (buf(1:8*size(x)))
      buf = transfer(x, buf)
      call sock_xfer(fd, buf, 8*size(x), .true.)
      deallocate (buf)

   end subroutine ipi_put_reals

   subroutine ipi_put_real(fd, v)

      implicit none

      integer, intent(in) :: fd
      real(dp), intent(in) :: v
      real(dp) :: x(1:1)

      x(1) = v
      call ipi_put_reals(fd, x)

   end subroutine ipi_put_real

!  Read and throw away. The INIT message carries a per-replica parameter
!  string that electronic-structure clients use and this one does not, but it
!  is on the wire and has to leave it.
   subroutine ipi_skip_bytes(fd, nbytes)

      implicit none

      integer, intent(in) :: fd
      integer, intent(in) :: nbytes
      integer(c_signed_char), allocatable :: buf(:)

      if (nbytes <= 0) return
      allocate (buf(1:nbytes))
      call sock_xfer(fd, buf, nbytes, .false.)
      deallocate (buf)

   end subroutine ipi_skip_bytes

end module ipi_socket
