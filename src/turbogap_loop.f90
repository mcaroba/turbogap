! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2026, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_loop.f90, is copyright (c) 2026, Miguel A. Caro and
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

!  Where the main loop is: which frame, MD step or MC step it is on, whether it
!  should go round again, and the progress bar's state.
module turbogap_loop

   use turbogap_comm, only: comm_t
   use types, only: input_parameters

   implicit none

   private
   public :: loop_init
   public :: loop_continues
   public :: loop_begin_step

   character*1, parameter, public :: creturn = achar(13)

   type, public :: loop_t
      integer :: md_istep
      integer :: mc_istep
      integer :: n_xyz
      logical :: repeat_xyz = .true.
      logical :: exit_loop = .true.
      integer :: n_sites_prev = 0
      integer :: counter = 0
      integer :: update_bar
      integer :: bar_frac = 0
   end type loop_t

contains

   subroutine loop_init(loop, params, comm)
      type(loop_t), intent(inout) :: loop
      type(input_parameters), intent(in) :: params
      type(comm_t), intent(in) :: comm

      loop%md_istep = -1
      loop%mc_istep = -1
      loop%n_xyz = 0

      if (params%do_md) then
         if (comm%rank == 0) then
            write (*, *) '                                       |'
            write (*, *) 'Doing molecular dynamics...            |'
            if (params%print_progress .and. loop%md_istep > 0) then
               write (*, *) '                                       |'
               write (*, *) 'Progress:                              |'
               write (*, *) '                                       |'
               write (*, '(1X,A)', advance='no') '[                                    ] |'
            end if
         end if
         loop%update_bar = params%md_nsteps/36
         if (loop%update_bar < 1) then
            loop%update_bar = 1
         end if
         loop%counter = 1
      end if
   end subroutine loop_init

   logical function loop_continues(loop, params)
      type(loop_t), intent(in) :: loop
      type(input_parameters), intent(in) :: params

      loop_continues = loop%repeat_xyz .or. (params%do_md .and. loop%md_istep < params%md_nsteps) &
                       .or. (params%do_mc .and. loop%mc_istep < params%mc_nsteps)
   end function loop_continues

!  Advance the step counters and redraw the MD progress bar.
   subroutine loop_begin_step(loop, params, comm)
      type(loop_t), intent(inout) :: loop
      type(input_parameters), intent(in) :: params
      type(comm_t), intent(in) :: comm
      integer :: i
      integer :: j

      if (params%do_mc) then
         loop%mc_istep = loop%mc_istep + 1
         ! Undo if the step is md related
         if (loop%md_istep > -1) loop%mc_istep = loop%mc_istep - 1

      end if

      if (params%do_md) then
         loop%md_istep = loop%md_istep + 1
      else
         loop%n_xyz = loop%n_xyz + 1
      end if

      !   Update progress bar
      !
      !   md_nsteps = 0 is a legitimate input -- one configuration, forces, no
      !   dynamics -- and the bar divides by it. Integer division by zero is a
      !   SIGFPE, so the run died on the first step with a backtrace and no
      !   message rather than producing the single frame it was asked for.
      loop%bar_frac = 0
      if (params%md_nsteps > 0) loop%bar_frac = 36*loop%md_istep/params%md_nsteps
      if (loop%bar_frac > 36) loop%bar_frac = 36
      if (params%print_progress .and. loop%counter == loop%update_bar .and. (.not. params%do_mc)) then
         if (comm%rank == 0) then
            do j = 1, 36 + 3
               write (*, "(A)", advance="no") creturn
            end do
            write (*, "(1X,A)", advance="no") "["
            do i = 1, loop%bar_frac
               write (*, "(A)", advance="no") "."
            end do
            do i = loop%bar_frac + 1, 36
               write (*, "(A)", advance="no") " "
            end do
            write (*, "(A)", advance="no") "] |"
            if (loop%md_istep == params%md_nsteps) then
               write (*, *)
            end if
         end if
         loop%counter = 1
      else
         loop%counter = loop%counter + 1
      end if
   end subroutine loop_begin_step

end module turbogap_loop
