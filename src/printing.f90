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

   subroutine print_parameters(desc, values, unit)
      class(*), intent(in) :: values(:)
      character(len=*) :: desc
      character(len=*), optional :: unit
      integer :: i, n

      n = size(values)
      call print_line(desc, 'normal')
      do i = 1, n
         if (present(unit)) then
            call print_parameter("", values(i), unit)
         else
            call print_parameter("", values(i))
         end if

      end do
   end subroutine print_parameters

   subroutine print_parameter(desc, value, unit)
      class(*), intent(in) :: value
      character(len=*) :: desc
      character(len=max_length) :: string
      character(len=max_length_half) :: temp
      character(len=max_length_half) :: temp2
      character(len=3) :: equals = ' = '
      character(len=*), optional :: unit

      write (temp, '(a20)') adjustl(trim(desc))
      if (present(unit)) then

         select type (value)
         type is (character(len=*))
            write (string, '(a17,a3,a,a)') adjustl(temp), equals, trim(value), adjustl(unit)
         type is (real(dp))
            if (value < 0.0_dp) then
               !write (temp2, '(g13.6)') value
               write (temp2, '(f13.6)') value
               write (string, '(a20,a3,a13,1X,a)') adjustl(temp), equals, adjustl(temp2), adjustl(unit)
            else
               !write (temp2, '(g13.7)') value
               write (temp2, '(f13.7)') value
               write (string, '(a20,a3,a13,1X,a)') adjustl(temp), equals, " "//adjustl(temp2), adjustl(unit)
            end if

            !write (string, '(a17,a3,g18.6,a)') adjustl(temp), equals, value, adjustl(unit)
         type is (integer)
            write (string, '(a17,a3,i8,a)') adjustl(temp), equals, value, adjustl(unit)
         type is (logical)
            write (string, '(a17,a3,1x,l2,a)') adjustl(temp), equals, value, adjustl(unit)
         end select

      else

         select type (value)
         type is (character(len=*))
            if (len_trim(value) > max_length - 24) then
               write (temp, '(a20,a3)') adjustl(temp), equals
               call print_line(temp//adjustl(trim(value)), 'normal')
               return
            else

               write (string, '(a20,a3,a)') adjustl(temp), equals, " "//trim(value)
            end if

         type is (real(dp))
            if (value < 0.0_dp) then
               !write (temp2, '(g18.6)') value
               write (temp2, '(f18.6)') value
               write (string, '(a20,a3,a18)') adjustl(temp), equals, adjustl(temp2)
            else
               !jwrite (temp2, '(g19.7)') value
               write (temp2, '(f19.7)') value
               write (string, '(a20,a3,a19)') adjustl(temp), equals, " "//adjustl(temp2)
            end if

         type is (integer)
            if (value < 0) then
               write (temp2, '(i8)') value
               write (string, '(a20,a3,a18)') adjustl(temp), equals, adjustl(temp2)
            else
               write (temp2, '(i9)') value
               write (string, '(a20,a3,a19)') adjustl(temp), equals, " "//adjustl(temp2)
            end if

         type is (logical)
            if (value) then
               write (temp2, '(a6)') '.true.'
            else
               write (temp2, '(a7)') '.false.'
            end if
            write (string, '(a20,a3,1x,a18)') adjustl(temp), equals, adjustl(temp2)
         end select
      end if

      call print_line(string, 'normal')

   end subroutine print_parameter

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
