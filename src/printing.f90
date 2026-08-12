! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, printing.f90, is copyright (c) 2019-2025, Miguel A. Caro and
! HND X   Tigany Zarrouk
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module printing
   use kinds, only: dp

   implicit none

   integer, parameter, private :: max_length = 42
   integer, parameter, private :: max_length_half = 31
   integer, parameter, private :: max_length_line_type = 87

contains

   subroutine print_matrix_dp(values)
      real(dp), intent(in) :: values(:, :)
      character(len=max_length) :: string
      character(len=max_length) :: temp_string
      character(len=max_length) :: fmt
      integer :: n1
      integer :: n2
      integer :: i
      integer :: j

      n1 = size(values, 1)
      n2 = size(values, 2)

      if (n1 < 10) then
         write (fmt, '(A1,I1,A)') "(", n1, "(1X, F12.8))"
      elseif (n1 < 100) then
         write (fmt, '(A1,I2,A)') "(", n1, "(1X, F12.8))"
      elseif (n1 < 1000) then
         write (fmt, '(A1,I3,A)') "(", n1, "(1X, F12.8))"
      elseif (n1 < 10000) then
         write (fmt, '(A1,I4,A)') "(", n1, "(1X, F12.8))"
      elseif (n1 < 100000) then
         write (fmt, '(A1,I5,A)') "(", n1, "(1X, F12.8))"
      elseif (n1 < 1000000) then
         write (fmt, '(A1,I6,A)') "(", n1, "(1X, F12.8))"
      else
         write (fmt, '(A1,I7,A)') "(", n1, "(1X, F12.8))"
      end if

      do i = 1, n2
         write (*, fmt) values(1:n1, i)
      end do

   end subroutine print_matrix_dp

!  The input-file echo.
!
!  These lines are meant to be read back: the point of echoing a deck is that a
!  reader can see what the run actually decided, and copy any of it into an
!  input file of their own. So they carry no `|` end cap, and nothing is
!  truncated -- both of which the fixed-width buffers here used to do. A name
!  longer than 20 characters lost its tail (`write_pair_distribut`,
!  `mc_planes_restrict_t`), and a value longer than the 42-character line lost
!  its tail too, silently. Everything is built in a deferred-length string now,
!  so a long path or a long species list comes out whole.
!
!  Values are written with G0/I0 rather than a fixed F or I width, so they are
!  valid Fortran literals whatever their magnitude, which is what makes the
!  line paste back into a deck.
   subroutine print_parameters(desc, values, unit)
      class(*), intent(in) :: values(:)
      character(len=*) :: desc
      character(len=*), optional :: unit
      character(len=:), allocatable :: line
      integer :: i

!     One line for the whole list, "name = v1 v2 v3", because that is the shape
!     an input file wants it back in. The old form printed the name on its own
!     line and then one value per line, which read back as nothing.
      line = pad_name(desc)//' ='
      do i = 1, size(values)
         line = line//' '//value_string(values(i))
      end do
!     The unit goes after a `!`, which is what the input parser treats as a
!     comment, so the line stays paste-able as it stands.
      if (present(unit)) line = line//'   ! '//trim(adjustl(unit))

      write (*, '(1X,A)') line

   end subroutine print_parameters

   subroutine print_parameter(desc, value, unit)
      class(*), intent(in) :: value
      character(len=*) :: desc
      character(len=*), optional :: unit
      character(len=:), allocatable :: line

      line = pad_name(desc)//' = '//value_string(value)
      if (present(unit)) line = line//'   ! '//trim(adjustl(unit))

      write (*, '(1X,A)') line

   end subroutine print_parameter

!  Keyword names padded to a common column so the `=` line up. The width is the
!  longest keyword there is (mc_planes_restrict_to_polyhedron, 32); a longer one
!  is not truncated, it just pushes its own `=` out.
   function pad_name(desc) result(padded)
      character(len=*), intent(in) :: desc
      character(len=:), allocatable :: padded
      integer, parameter :: name_width = 32
      character(len=name_width) :: buffer

      if (len_trim(desc) > name_width) then
         padded = trim(adjustl(desc))
      else
         buffer = adjustl(desc)
         padded = buffer
      end if

   end function pad_name

!  One value, as a literal an input file would accept back.
   function value_string(value) result(text)
      class(*), intent(in) :: value
      character(len=:), allocatable :: text
      character(len=64) :: buffer

      select type (value)
      type is (character(len=*))
         text = trim(adjustl(value))
      type is (real(dp))
         write (buffer, '(G0.10)') value
         text = trim(adjustl(buffer))
      type is (integer)
         write (buffer, '(I0)') value
         text = trim(adjustl(buffer))
      type is (logical)
         if (value) then
            text = '.true.'
         else
            text = '.false.'
         end if
      class default
         text = '<unprintable>'
      end select

   end function value_string

   subroutine print_line(input_string, line_type)
      character(len=*) :: input_string
      character(len=max_length_line_type) :: string = ''
      character(len=*), optional :: line_type
      character(len=1) :: end_cap = '|'
      integer :: length
      integer :: length_diff

      length = len_trim(input_string)
      length_diff = max_length - 2 - length

      ! Now pad remaining with space
      if (length_diff > 0) then
         if (present(line_type)) then
            if (trim(line_type) == 'warning') then
               string = trim(input_string)//repeat(' ', length_diff - 1)//end_cap//' <-- WARNING     '
            else if (trim(line_type) == 'error') then
               string = trim(input_string)//repeat(' ', length_diff - 1)//end_cap//' <-- !! ERROR !! '
            else if (trim(line_type) == 'note') then
               string = trim(input_string)//repeat(' ', length_diff - 1)//end_cap//' --- Note        '
            else if (trim(line_type) == 'debug') then
               string = trim(input_string)//repeat(' ', length_diff - 1)//end_cap//' @_@ DEBUG @_@   '
            else
               string = trim(input_string)//repeat(' ', length_diff - 1)//end_cap//'                 '
            end if
         else
            string = trim(input_string)//repeat(' ', length_diff - 1)//end_cap//'                 '
         end if
      else
         string = trim(input_string)
      end if

      write (*, *) string

   end subroutine print_line

   subroutine print_unbroken_lines(string, line_type)
      character(len=*), intent(in) :: string
      character(len=*) :: line_type
      character(len=max_length_line_type) :: line
      integer :: i, word_start, word_end, line_length
      logical :: found_word = .false.

      line = ''
      word_start = 1
      line_length = 0

      ! have an index for the start of the word and then we loop over the line
      ! and see what fits

      do i = 1, len(trim(string))
         if (string(i:i) == ' ' .or. i == len(trim(string))) then
            ! We have found a word boundary
            if (i == len(trim(string))) then
               word_end = i
            else
               word_end = i - 1
               found_word = .true.
            end if

            if (line_length + (word_end - word_start + 1) > max_length - 6) then
               call print_line(trim(line), line_type)
               line = ''            ! Reset the line
               line_length = 0
               found_word = .false.
            end if
            if (found_word) then
               line = trim(line)//trim(string(word_start:word_end))
            else
               line = trim(string(word_start:word_end))
            end if
            line_length = len(trim(line))
            word_start = i
         end if
      end do

      if (line_length > 0) then
         call print_line(trim(line), line_type)
      end if
   end subroutine print_unbroken_lines

   subroutine print_separator(separation_character)
      character, intent(in) :: separation_character

      call print_line(repeat(separation_character, max_length - 3), 'normal')

   end subroutine print_separator

   subroutine print_small_message(message)
      character(len=*) :: message
      integer :: length

      length = len_trim(message)
      !call print_line(repeat('-', length))
      !call print_separator('-')
      call print_unbroken_lines(message, 'normal')
      call print_line(repeat('-', length))
      ! call print_separator('-')

   end subroutine print_small_message

   subroutine printf_small_message(message, value, message2)
      character(len=*) :: message
      character(len=*) :: message2
      class(*) :: value
      character(len=max_length*4) :: string

      !call print_separator('-')
      call print_separator(' ')

      select type (value)
      type is (character(len=*))
         write (string, '(a,a,a)') adjustl(message), value, message2
      type is (real(dp))
         write (string, '(a,f12.6,a)') adjustl(message), value, message2
      type is (integer)
         write (string, '(a,i8,a)') adjustl(message), value, message2
      end select
      call print_unbroken_lines(trim(string), 'normal')

   end subroutine printf_small_message

   subroutine printf_message(message, value, message2)
      character(len=*) :: message
      character(len=*) :: message2
      class(*) :: value
      character(len=max_length*4) :: string

      call print_separator('-')
      call print_separator(' ')

      select type (value)
      type is (character(len=*))
         write (string, '(a,a,a)') adjustl(message), value, message2
      type is (real(dp))
         write (string, '(a,f12.6,a)') adjustl(message), value, message2
      type is (integer)
         write (string, '(a,i8,a)') adjustl(message), value, message2
      end select
      call print_unbroken_lines(trim(string), 'normal')

      call print_separator(' ')
      call print_separator('-')

   end subroutine printf_message

   subroutine print_message(message)
      character(len=*) :: message

      call print_separator('-')
      call print_separator(' ')

      call print_unbroken_lines(message, 'normal')

      call print_separator(' ')
      call print_separator('-')

   end subroutine print_message

   subroutine print_debug(message, location, rank)
      character(len=*) :: message, location
      character*32 :: debug = 'debug'
      integer, intent(in), optional :: rank
      character*8 :: rank_str

      if (present(rank)) then
         write (rank_str, '(I4,1X)') rank
      else
         rank_str = ""
      end if

      call print_separator('X')
      call print_separator(' ')
      call print_unbroken_lines(location, debug)
      call print_separator(' ')

      call print_unbroken_lines(rank_str//message, debug)

      call print_separator(' ')
      call print_separator('X')

   end subroutine print_debug

   subroutine print_note(message)
      character(len=*) :: message
      character*32 :: debug = 'note'

      call print_separator('-')
      call print_separator(' ')

      call print_unbroken_lines(message, debug)

      call print_separator(' ')
      call print_separator('-')

   end subroutine print_note

   subroutine print_error(message)
      character(len=*) :: message
      character*32 :: debug = 'error'

      call print_separator('!')
      call print_separator(' ')

      call print_unbroken_lines(message, debug)

      call print_separator(' ')
      call print_separator('!')

   end subroutine print_error

   subroutine print_warning(message)
      character(len=*) :: message
      character*32 :: debug = 'warning'

      call print_separator('~')
      call print_separator(' ')

      call print_unbroken_lines(message, debug)

      call print_separator(' ')
      call print_separator('~')

   end subroutine print_warning

   subroutine print_end_of_execution()
      write (*, *) '                                       |'
      write (*, *) '.......................................|'
      write (*, *) '                                       |'
      write (*, *) 'All GAPs have been Turboed.            |'
      write (*, *) '                                       |'
      write (*, *) 'End of execution                       |'
      write (*, *) '_______________________________________/'
   end subroutine print_end_of_execution

end module printing

! ! Uncomment to test the module
! program test_print
!    use printing
!    implicit none
!    character(len=:), allocatable :: inputString
!    real*8 :: f = 12.66775599
!    integer :: i = 12345
!    logical :: true = .true., false = .false.
!    character(len=30) :: blah
!    inputString = "This is an example supercalifragilisticexpia&
!         & of a thing long string with many words that need&
!         & to be printed on separate lines without breaking&
!         & any of the words."

!    blah = "check?!!?"
!    call print_parameter("vibes? ", blah)
!    call print_parameter("f ", f)
!    call print_parameter("i ", i)
!    call print_parameter("true ", true)
!    call print_parameter("false ", false)

!    call print_line('yadda yadda', 'normal')
!    call print_debug(inputString, 'test_print')
!    call print_warning(inputString)
!    call print_error(inputString)
!    call print_note(inputString)
!    call print_message(inputString)
!    call print_message("Turbogap needs more information&
!         &from the input file for this calculation!! ")
! end program test_print
