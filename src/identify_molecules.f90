! ----------------------------------------------------
! Idea here is that we have some reference database or 
! a built in list of supported structures. Within the
! loop, this code identifies the molecules and updates
! a database that contains the counts. The counts are
! written into a file at every step or every k steps.
! 
! I WILL NOT SPEND TOO MUCH TIME ON THIS, WE CAN JUST 
! DO THE COUNTS IN POST-PROCESSING BUT WE DO NEED THE 
! ABILITY TO REMOVE ATOMS FROM THE SIMULATION
!
! In the mean time, the main loop will identify 
! molecules/particles that have broken off from other
! atoms. The length scale can be specified as a keyword
! The vacuum has to be large enough in the simulation.


module molecule_identifier
    use kinds
    use types
#ifdef _MPIF90
   use mpi
#endif


    contains

    subroutine remove_particles(params, positions, positions_diff, n_sites)
        implicit none
        type(input_parameters), intent(in) :: params
        real(dp), intent(inout) :: positions(:,:), positions_diff(:,:)
        integer, intent(inout) :: n_sites

        integer :: i, j

        if (params%remove_detached) then
            do i = 1, n_sites
                if (positions(3,i) > params%removal_height) then
                    do j = i, n_sites - 1
                        positions(:,j) = positions(:,j+1)
                        positions_diff(:,j) = positions_diff(:,j+1)
                    enddo
                    n_sites = n_sites - 1
                endif
            enddo
        endif
    end subroutine remove_particles


end module molecule_identifier