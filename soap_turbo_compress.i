# 1 "src/soap_turbo/src/soap_turbo_compress.f90"
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2021, Miguel A. Caro
! HND X
! HND X   soap_turbo is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   soap_turbo is distributed in the hope that it will be useful for non-commercial
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

# 26
module soap_turbo_compress

  contains

!**************************************************************************
  subroutine get_compress_indices( compress_mode, alpha_max, l_max, dim, indices, what_to_do )

  implicit none

! Input variables
  integer, intent(in) :: alpha_max(:), l_max
  character(*), intent(in) :: compress_mode, what_to_do

! Input-Output variables
  integer, intent(inout) :: dim, indices(:)

! Internal vairables
  integer, allocatable :: pivot(:)
  integer :: n_species, n_max, i, counter, n, m, l, k
  logical :: set_indices

  if( what_to_do == "get_dim" )then
    set_indices = .false.
  else if( what_to_do == "set_indices" )then
    set_indices = .true.
  else
    write(*,*) "ERROR: incorrect what_to_do in get_compress_indices() subroutine!"
    stop
  end if

  n_species = size(alpha_max)
  n_max = sum(alpha_max)

  if( compress_mode == "trivial" .or. compress_mode == "Trivial" )then
    allocate( pivot(1:n_species) )
    pivot = 0
    pivot(1) = 1

    do i = 1, n_species-1
      pivot(i+1) = pivot(i) + alpha_max(i)
    end do

    counter = 0
    k = 1

    if( set_indices )then
      indices = 0
    end if

    do n = 1, n_max    
      do m = n, n_max    
        do l = 0, l_max    
          if( any(n == pivot) .or. any(m == pivot) )then
            counter = counter + 1
            if( set_indices )then
              indices(counter) = k
            end if
          end if
          k = k + 1
        end do
      end do
    end do
    dim = counter
    deallocate( pivot )
  else
    write(*,*) "ERROR: I don't understand compress_mode =", compress_mode
    stop
  end if

  end subroutine
!**************************************************************************

end module
