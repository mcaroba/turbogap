! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_utils.f90, is copyright (c) 2019-2026, Miguel A. Caro and
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

!  Input-file reading helpers, taken from turbogap_2.0's src/read/read_utils.f90.
!
!  Two things the old reading code could not do:
!
!   * say what it read. Every keyword handler was
!     `read (10, *, iostat=iostatus) cjunk, cjunk, params%x` and nothing echoed
!     the value, so a typo'd or defaulted parameter was invisible.
!   * fail when a list is the wrong length. A list read of n elements from a line
!     carrying fewer silently leaves the rest at whatever they held, and a line
!     carrying more silently drops the excess. Both run to completion and give
!     wrong numbers.
!
!  check_read_parameters_count closes the second by counting the fields on the
!  line before consuming it, and read_parameters wraps the pair together.
!
!  Not brought across from 2.0: check_deprecated, print_deprecation_message and
!  upper_to_lower_case, which already exist in read_files.f90 here -- defining
!  them in both would be an ambiguity, not a convenience.
module read_utils

   use kinds, only: dp
   use printing, only: print_error
   use error, only: turbogap_abort

   implicit none

!  read_parameters_char8 exists in 2.0 but is deliberately absent from the
!  generic interface there, and from this one: its `character*8, allocatable`
!  dummy and the `character(len=*), allocatable` one below are not
!  distinguishable for generic resolution. The len=* version handles character*8
!  actual arguments perfectly well, so there is nothing to add.
   interface read_parameters
      module procedure read_parameters_double
      module procedure read_parameters_int
      module procedure read_parameters_char
   end interface read_parameters

contains

   subroutine read_parameters_double(unit, iostatus, n_elements, values)
      real(dp), allocatable, intent(inout) :: values(:)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: cjunk

      if (.not. allocated(values)) then
         call print_error("The real variables one is trying to read to is not&
              & allocated. This is likely due to you not ordering elements in&
              & the input file correctly. The number of elements should be&
              & specified before the values to be read in")
         call turbogap_abort()
      end if

      call check_read_parameters_count(unit, iostatus, n_elements)

      read (unit, *, iostat=iostatus) cjunk, cjunk, values(1:n_elements)
   end subroutine read_parameters_double

   subroutine read_parameters_int(unit, iostatus, n_elements, values)
      integer, allocatable, intent(inout) :: values(:)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: cjunk

      if (.not. allocated(values)) then
         call print_error("The integer variables one is trying to read to is not&
              & allocated. This is likely due to you not ordering elements in&
              & the input file correctly. The number of elements should be&
              & specified before the values to be read in")
         call turbogap_abort()
      end if

      call check_read_parameters_count(unit, iostatus, n_elements)

      read (unit, *, iostat=iostatus) cjunk, cjunk, values(1:n_elements)
   end subroutine read_parameters_int

   subroutine read_parameters_char(unit, iostatus, n_elements, values)
      character(len=*), allocatable, intent(inout) :: values(:)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: cjunk

      if (.not. allocated(values)) then
         call print_error("The character variables one is trying to read to is not&
              & allocated. This is likely due to you not ordering elements in&
              & the input file correctly. The number of elements should be&
              & specified before the values to be read in")
         call turbogap_abort()
      end if

      call check_read_parameters_count(unit, iostatus, n_elements)

      read (unit, *, iostat=iostatus) cjunk, cjunk, values(1:n_elements)
   end subroutine read_parameters_char

!  Count the fields on the next line without consuming it, and refuse to read
!  more elements than are actually there.
!
!  The field count is done by collapsing runs of blanks rather than counting
!  every blank: `species = C  O` has a double space and 2.0's version counts it
!  as an extra field, which makes the check too permissive by one for every
!  repeated separator. Tabs count as separators too.
   subroutine check_read_parameters_count(unit, iostatus, n_elements)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: string
      integer :: i, n_words, n_values, trimmed
      logical :: in_word
      character*1 :: c

      read (unit, '(A)', iostat=iostatus) string
      call check_iostatus(iostatus, "read_parameters")
      backspace (unit)

      n_words = 0
      in_word = .false.
      trimmed = len_trim(string)
      do i = 1, trimmed
         c = string(i:i)
         if (c == ' ' .or. c == char(9)) then
            in_word = .false.
         else if (.not. in_word) then
            in_word = .true.
            n_words = n_words + 1
         end if
      end do

!     The first two fields are the keyword and the '=' sign.
      n_values = n_words - 2

      if (n_words > 0 .and. n_elements > n_values) then
         call print_error("The number of elements to read in the following line is wrong! "//trim(string))
         write (*, '(A,1X,I0,A,1X,I0)') "Expected number of elements", n_elements, &
            ", number of elements found", max(n_values, 0)
         call turbogap_abort()
      end if
   end subroutine check_read_parameters_count

   subroutine check_iostatus(iostatus, keyword)
      integer, intent(in) :: iostatus
      character(len=*) :: keyword

      if (iostatus /= 0) then
         call print_error("TurboGAP has had an issue reading the input file with&
             & the keyword "//trim(keyword)//". Please make sure&
             & the number of variables and the types of variables&
             & are consistent with the documentation and that the&
             & input file ends with a new line character.")
         call turbogap_abort()
      end if
   end subroutine check_iostatus

   function file_exists(filename) result(res)
      character(len=*), intent(in) :: filename
      logical                     :: res
      inquire (file=trim(filename), exist=res)
   end function file_exists

   subroutine check_file_exists(filename)
      character(len=*), intent(in) :: filename
      character(len=256) :: string
      if (.not. file_exists(filename)) then
         write (string, '(A,1X,A,1X,A)') "The file", trim(filename), &
            "does not exist! Please specify the right path in the input file!"
         call print_error(string)
         call turbogap_abort()
      end if
   end subroutine check_file_exists

end module read_utils
