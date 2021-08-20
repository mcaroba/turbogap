! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2021, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, xyz.f90, is copyright (c) 2021, Miguel A. Caro
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

module xyz_module


  use soap_turbo_functions

  contains


!**************************************************************************
! This subroutine writes the trajectory_out.xyz file with ASE's extended
! XYZ format
!
! Each property is designated with a number:
!
!  1 -> Number of atoms
!  2 -> Lattice vectors
!  3 -> Temperature
!  4 -> Pressure
!  5 -> Time step
!  6 -> Time
!  7 -> Total energy
!  8 -> Virial tensor
!  9 -> Stress tensor
! 10 -> Unit cell volume
! 11 -> Step number
!
! Each array property is designated with a number:
!  1 -> Species
!  2 -> Positions
!  3 -> Velocities
!  4 -> Forces
!  5 -> Local energy
!  6 -> Masses
!  7 -> Hirshfeld volumes
!  8 -> Fix atoms
!
! If the corresponding write_property(i) or write_array_property(i)
! is .true., we write out the corresponding property
!
  subroutine write_extxyz( Nat, md_istep, dt, temperature, pressure, a_cell, b_cell, c_cell, virial, &
                           species, positions, velocities, forces, local_energies, masses, &
                           hirshfeld_v, write_property, write_array_property, fix_atom )

    implicit none

!   In variables:
    real*8, intent(in) :: dt, temperature, pressure, a_cell(1:3), b_cell(1:3), c_cell(1:3), virial(1:3,1:3)
    real*8, intent(in) :: forces(:,:), velocities(:,:), positions(:,:), local_energies(:), masses(:)
    real*8, intent(in) :: hirshfeld_v(:)
    integer, intent(in) :: Nat, md_istep
    character(len=*), intent(in) :: species(:)
    logical, intent(in) :: write_property(:), write_array_property(:), fix_atom(:,:)

!   Internal variables:
    real*8 :: vol
    integer :: n_properties, n_array_properties, i, j
    character*1024 :: properties_string
    character*16 :: lattice_string(1:16), temp_string

    n_properties = 0
    do i = 1, size( write_property )
      if( write_property(i) )then
        n_properties = n_properties + 1
      end if
    end do

    n_array_properties = 0
    do i = 1, size( write_array_property )
      if( write_array_property(i) )then
        n_array_properties = n_array_properties + 1
      end if
    end do

    if( md_istep == 0 .or. md_istep == -1 )then
      open(unit=10, file="trajectory_out.xyz", status="unknown")
    else
      open(unit=10, file="trajectory_out.xyz", status="old", position="append")
    end if

!   We always write the number of atoms on the first line
    write(10, "(I8)") Nat

!   Now we write the properties on the "comment" line
!
!   First check which array properties we will write out
    if( n_array_properties > 0 )then
      properties_string = "Properties="
      i = 0
      if( write_array_property(1) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "species:S:1"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      if( write_array_property(2) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "positions:R:3"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      if( write_array_property(3) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "velocities:R:3"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      if( write_array_property(4) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "forces:R:3"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      if( write_array_property(5) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "local_energy:R:1"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      if( write_array_property(6) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "masses:R:1"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      if( write_array_property(7) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "hirshfeld_v:R:1"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      if( write_array_property(8) )then
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "fix_atoms:S:3"
        i = i + 1
        if( i < n_array_properties )then
          write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
        end if
      end if
      write(10, "(1X,A)", advance="no") trim(adjustl(properties_string))
    end if
!   Now write the NON array properties
!
!   Lattice vectors
    if( write_property(2) )then
      do i = 1, 3
        write(lattice_string(i),'(F16.10)') a_cell(i)
        write(lattice_string(i+3),'(F16.10)') b_cell(i)
        write(lattice_string(i+6),'(F16.10)') c_cell(i)
      end do
      write(10, "(1X,11A)", advance="no") "Lattice=""", adjustl(lattice_string(1)), lattice_string(2:9), """"
    end if
!
!   Temperature
    if( write_property(3) )then
      write(temp_string, "(F16.8)") temperature
      write(10, "(1X,2A)", advance="no") "temperature=", trim(adjustl(temp_string))
    end if
!
!   Pressure
    if( write_property(4) )then
      write(temp_string, "(F16.8)") pressure
      write(10, "(1X,2A)", advance="no") "pressure=", trim(adjustl(temp_string))
    end if
!
!   Time step
    if( write_property(5) )then
      write(temp_string, "(F16.4)") dt
      write(10, "(1X,2A)", advance="no") "time_step=", trim(adjustl(temp_string))
    end if
!
!   Time
    if( write_property(6) )then
      write(temp_string, "(F16.4)") dfloat(md_istep)*dt
      write(10, "(1X,2A)", advance="no") "time=", trim(adjustl(temp_string))
    end if
!
!  Total energy
    if( write_property(7) )then
      write(temp_string, "(F16.6)") sum(local_energies)
      write(10, "(1X,2A)", advance="no") "energy=", trim(adjustl(temp_string))
    end if
!
!   Virial tensor
    if( write_property(8) )then
      do i = 1, 3
        write(lattice_string(i),'(F16.8)') virial(1,i)
        write(lattice_string(i+3),'(F16.8)') virial(2,i)
        write(lattice_string(i+6),'(F16.8)') virial(3,i)
      end do
      write(10, "(1X,11A)", advance="no") "virial=""", adjustl(lattice_string(1)), lattice_string(2:9), """"
    end if
!
!   Stress tensor
    if( write_property(9) )then
      vol = dot_product(cross_product(a_cell, b_cell), c_cell)
      do i = 1, 3
        write(lattice_string(i),'(F16.8)') -virial(1,i)/vol
        write(lattice_string(i+3),'(F16.8)') -virial(2,i)/vol
        write(lattice_string(i+6),'(F16.8)') -virial(3,i)/vol
      end do
      write(10, "(1X,11A)", advance="no") "stress=""", adjustl(lattice_string(1)), lattice_string(2:9), """"
    end if
!
!   Volume
    if( write_property(10) )then
      write(temp_string, "(F16.6)") dot_product(cross_product(a_cell, b_cell), c_cell)
      write(10, "(1X,2A)", advance="no") "volume=", trim(adjustl(temp_string))
    end if
!
!   Step
    if( write_property(11) )then
      if( md_istep >= 0 )then
        write(temp_string, "(I10)") md_istep
        write(10, "(1X,2A)", advance="no") "step=", trim(adjustl(temp_string))
      else
        write(temp_string, "(I10)") -md_istep
        write(10, "(1X,2A)", advance="no") "i_config=", trim(adjustl(temp_string))
      end if
    end if
!
!   Advance
    write(10,*)



!   Write the array properties now
    do i = 1, Nat
!     Species
      if( write_array_property(1) )then
        write(10, "(1X,A8)", advance="no") species(i)
      end if
!     Positions
      if( write_array_property(2) )then
        write(10, "(1X,F16.8,1X,F16.8,1X,F16.8)", advance="no") positions(1:3, i)
      end if
!     Velocities
      if( write_array_property(3) )then
        write(10, "(1X,F16.8,1X,F16.8,1X,F16.8)", advance="no") velocities(1:3, i)
      end if
!     Forces
      if( write_array_property(4) )then
        write(10, "(1X,F16.8,1X,F16.8,1X,F16.8)", advance="no") forces(1:3, i)
      end if
!     Local energy
      if( write_array_property(5) )then
        write(10, "(1X,F16.8)", advance="no") local_energies(i)
      end if
!     Masses
      if( write_array_property(6) )then
        write(10, "(1X,F16.8)", advance="no") masses(i)/103.6426965268d0
      end if
!     Hirshfeld volumes
      if( write_array_property(7) )then
        write(10, "(1X,F16.8)", advance="no") hirshfeld_v(i)
      end if
!     Fix atoms
      if( write_array_property(8) )then
        write(10, "(1X,L1,1X,L1,1X,L1)", advance="no") fix_atom(1:3, i)
      end if
!     Advance
      write(10,*)
    end do



    close(10)

  end subroutine
!**************************************************************************

end module
