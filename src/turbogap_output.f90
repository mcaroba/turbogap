! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_output.f90, is copyright (c) 2026, Miguel A. Caro and
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

!  What the driver itself prints: the banner and the options summary.
module turbogap_output

   use turbogap_comm, only: comm_t, comm_with_mpi
   use types, only: input_parameters
   use turbogap_setup, only: model_t

   implicit none

   private
   public :: print_banner
   public :: print_options

contains

   subroutine print_banner(comm)
      type(comm_t), intent(in) :: comm

      if (comm%rank == 0) then
         write (*, *) '_________________________________________________________________ '
         write (*, *) '                             _                                   \'
         write (*, *) ' ___________            __   \\ /\        _____     ___   _____  |'
         write (*, *) '/____  ____/           / / /\|*\|*\/\    / ___ \   /   | |  _  \ |'
         write (*, *) '    / / __  __  __    / /  \********/   / /  /_/  / /| | | / | | |'
         write (*, *) '   / / / / / / / /_  / /__  \**__**/   / / ____  / / | | | |_/ / |'
         write (*, *) '  / / / / / / / __/ / ___ \ /*/  \*\  / / /_  / / /__| | |  __/  |'
         write (*, *) ' / / / /_/ / / /   / /__/ / \ \__/ / / /___/ / / ____  | | |     |'
         write (*, *) '/_/_/_____/_/_/___/______/___\____/__\______/_/_/____|_|_|_|____ |'
         write (*, *) '_____________________________________________________________  / |'
         write (*, *) '*************************************************************|/  |'
         write (*, *) '                  Welcome to the TurboGAP code                   |'
         write (*, *) '                         Maintained by                           |'
         write (*, *) '                                                                 |'
         write (*, *) '               Miguel A. Caro and Tigany Zarrouk                 |'
         write (*, *) '                       mcaroba@gmail.com                         |'
         write (*, *) '                      miguel.caro@aalto.fi                       |'
         write (*, *) '                                                                 |'
         write (*, *) '          Department of Chemistry and Materials Science          |'
         write (*, *) '                     Aalto University, Finland                   |'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) '====================>>>>>  turbogap.fi  <<<<<====================|'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) 'Contributors (code and methodology) in chronological order:      |'
         write (*, *) '                                                                 |'
         write (*, *) 'Miguel A. Caro, Patricia Hernández-León, Suresh Kondati          |'
         write (*, *) 'Natarajan, Albert P. Bartók, Eelis V. Mielonen, Heikki Muhli,    |'
         write (*, *) 'Mikhail Kuklin, Gábor Csányi, Jan Kloppenburg, Richard Jana,     |'
         write (*, *) 'Tigany Zarrouk, Cristian V. Achim                                |'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) '                     Last updated: Jun. 2026                     |'
         write (*, *) '                                        _________________________/'
         write (*, *) '.......................................|'
         if (comm_with_mpi) then
            write (*, *) '                                       |'
            write (*, *) 'Running TurboGAP with MPI support:     |'
            write (*, *) '                                       |'
            write (*, '(A,I6,A)') ' Running TurboGAP on ', comm%size, ' MPI tasks   |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         else
            write (*, *) '                                       |'
            write (*, *) 'Running the serial version of TurboGAP |'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         end if
      end if
   end subroutine print_banner

   subroutine print_options(comm, params, model)
      type(comm_t), intent(in) :: comm
      type(input_parameters), intent(in) :: params
      type(model_t), intent(in) :: model
      character*8 :: i_char
      integer :: i

      if (comm%rank == 0) then
         ! Print out chosen options:
         write (*, *) '                                       |'
         write (*, '(1X,A)') 'You specified the following options:   |'
         write (*, *) '                                       |'
         write (*, *) '---------------------------------      |'
         if (len(trim(params%atoms_file)) > 20) then
            write (*, '(1X,A,A20,A)') 'Atoms file = ', adjustr(trim(params%atoms_file)), '...   |'
         else
            write (*, '(1X,A,A20,A)') 'Atoms file = ', adjustr(trim(params%atoms_file)), '      |'
         end if
         write (*, *) '---------------------------------      |'
         write (i_char, '(I8)') model%n_species
         write (*, '(1X,A,A8,A)') 'No. of species   = ', adjustl(i_char), '            |'
         do i = 1, model%n_species
            write (i_char, '(I8)') i
           write (*, '(1X,A,A2,A,A8,A)') '  *) Species #', adjustl(i_char), ' =       ', adjustr(params%species_types(i)), '      |'
         end do
         write (*, *) '---------------------------------      |'
         write (*, '(1X,A,F15.4,A)') 'rcut_max = ', model%rcut_max, ' Angst.      |'
         write (*, *) '---------------------------------      |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
      end if
   end subroutine print_options

end module turbogap_output
