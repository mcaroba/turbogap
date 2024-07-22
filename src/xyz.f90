! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, xyz.f90, is copyright (c) 2021-2023, Miguel A. Caro
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
!  7 -> Hirshfeld volumes || Note, this is being replaced by the local properties
!  8 -> Fix atoms
!
! If the corresponding write_property(i) or write_array_property(i)
! is .true., we write out the corresponding property
!

    subroutine get_xyz_energy_string(energies_soap, energies_2b,&
         & energies_3b, energies_core_pot, energies_vdw, energies_exp&
         &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
         & valid_pdf, valid_sf, valid_xrd, valid_nd,  do_pair_distribution,&
         & do_structure_factor, do_xrd, do_nd, string)
      implicit none
      real*8, intent(in), allocatable :: energies_soap(:), energies_2b(:),&
           & energies_3b(:), energies_core_pot(:), energies_vdw(:),&
           & energies_exp(:), energies_lp(:), energies_pdf(:), energies_sf(:),&
           & energies_xrd(:), energies_nd(:)
      logical, intent(in) :: valid_pdf, valid_sf, valid_xrd, valid_nd, do_pair_distribution,&
           & do_structure_factor, do_xrd, do_nd
      character*1024, intent(out) :: string
      character*32 :: temp_string


      write(temp_string, "(F16.8)") sum(energies_soap)
      write(string, "(1X,A)") "energy_soap=" // trim(adjustl(temp_string))

      write(temp_string, "(F16.8)") sum(energies_2b)
      write(string, "(A)") adjustl(trim(string)) // " energy_2b=" // trim(adjustl(temp_string))

      write(temp_string, "(F16.8)") sum(energies_3b)
      write(string, "(A)") adjustl(trim(string)) // " energy_3b=" // trim(adjustl(temp_string))

      write(temp_string, "(F16.8)") sum(energies_core_pot)
      write(string, "(A)") adjustl(trim(string)) // " energy_core_pot=" // trim(adjustl(temp_string))

      write(temp_string, "(F16.8)") sum(energies_vdw)
      write(string, "(A)") adjustl(trim(string)) // " energy_vdw=" // trim(adjustl(temp_string))

      write(temp_string, "(F16.8)") sum(energies_exp)
      write(string, "(A)") adjustl(trim(string)) // " energy_exp=" // trim(adjustl(temp_string))

      write(temp_string, "(F16.8)") sum(energies_lp)
      write(string, "(A)") adjustl(trim(string)) // " energy_xps=" // trim(adjustl(temp_string))

      if ( valid_pdf .and. do_pair_distribution )then
         write(temp_string, "(F16.8)") sum(energies_pdf)
         write(string, "(A)") adjustl(trim(string)) // " energy_pdf=" // trim(adjustl(temp_string))
      end if
      if ( valid_sf .and. do_structure_factor )then
         write(temp_string, "(F16.8)") sum(energies_sf)
         write(string, "(A)") adjustl(trim(string)) // " energy_sf=" // trim(adjustl(temp_string))
      end if
      if ( valid_xrd .and. do_xrd )then
         write(temp_string, "(F16.8)") sum(energies_xrd)
         write(string, "(A)") adjustl(trim(string)) // " energy_xrd=" // trim(adjustl(temp_string))
      end if
      if ( valid_nd .and. do_nd )then
         write(temp_string, "(F16.8)") sum(energies_nd)
         write(string, "(A)") adjustl(trim(string)) // " energy_nd=" // trim(adjustl(temp_string))
      end if

    end subroutine get_xyz_energy_string




    
  subroutine write_extxyz( Nat, md_istep, md_time, dt, temperature, pressure, a_cell, b_cell, c_cell, virial, &

                           species, positions, velocities, forces, local_energies, masses, &
                           write_property,&
                           & write_array_property,&
                           & write_local_properties,&
                           & local_property_labels, local_properties, fix_atom,&
                           & filename , string,  overwrite )

    implicit none

!   In variables:
    real*8, intent(in) :: md_time, dt, temperature, pressure, a_cell(1:3), b_cell(1:3), c_cell(1:3), virial(1:3,1:3)
    real*8, intent(in) :: forces(:,:), velocities(:,:), positions(:,:), local_energies(:), masses(:)
    real*8, intent(in) :: local_properties(:,:)
    integer, intent(in) :: Nat, md_istep
    character(len=*), intent(in) :: species(:), filename, string
    logical, intent(in) :: write_property(:), write_array_property(:), fix_atom(:,:), overwrite
    logical, allocatable, intent(in) :: write_local_properties(:)
    character*1024, allocatable, intent(in) :: local_property_labels(:)
!   Internal variables:
    real*8 :: vol
    integer :: n_properties, n_array_properties, i, j, k
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
    ! Adding in for the local properties
    if(allocated(write_local_properties))then
       do i = 1, size( write_local_properties )
          if( write_local_properties(i) )then
             n_array_properties = n_array_properties + 1
          end if
       end do
    end if


    if( md_istep == 0 .or. md_istep == -1 .or. overwrite)then
      open(unit=10, file=filename, status="unknown")
    else
      open(unit=10, file=filename, status="old", position="append")
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
        write(properties_string, "(A)") trim(adjustl(properties_string)) // "pos:R:3"
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
     ! Now we write in the local properties of those which are passed in
     if(allocated(write_local_properties))then

        do k=1, size(write_local_properties,1)
           if( write_local_properties(k) )then
              write(properties_string, "(A)") trim(adjustl(properties_string)) // trim(local_property_labels(k)) // ":R:1"
              i = i + 1
              if( i < n_array_properties )then
                 write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
              end if
           end if
        end do
     end if

! Not removing yet for compatibility
      ! if( write_array_property(7) )then
      !   write(properties_string, "(A)") trim(adjustl(properties_string)) // "hirshfeld_v:R:1"
      !   i = i + 1
      !   if( i < n_array_properties )then
      !     write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
      !   end if
      ! end if
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
          !******** time is actually the md_time not md_istep*dt since dt changes, so md_time is put in the trajectory_out.xyz file
    
      write(temp_string, "(F16.6)") md_time			!dfloat(md_istep)*dt
      write(10, "(1X,2A)", advance="no") "time=", trim(adjustl(temp_string))
    end if
!
!  Total energy
    if( write_property(7) )then
      write(temp_string, "(F16.6)") sum(local_energies)
      write(10, "(1X,2A)", advance="no") "energy=", trim(adjustl(temp_string))

      write(10, "(1X,A)", advance="no") trim(adjustl(string))

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
     !  if( write_array_property(7) )then
     !    write(10, "(1X,F16.8)", advance="no") hirshfeld_v(i)
     ! end if
! Local properties
     if( allocated(write_local_properties ))then
        do j = 1, size(local_properties, 2)
           write(10, "(1X,F16.8)", advance="no") local_properties(i, j)
        end do
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







!**************************************************************************
  subroutine read_xyz_line( properties, line, species, positions, velocities, fix_atom, has_velocities, &
                            masses, has_masses )

    implicit none

!   Input variables
    character*1024, intent(in) :: properties, line

!   Output variables
    real*8, intent(inout) :: velocities(1:3), positions(1:3), masses
    character*8 :: species
    logical, intent(inout) :: fix_atom(1:3)
    logical, intent(out) :: has_velocities, has_masses

!   Internal variables
    integer :: i, j, k, iostatus
    character*1 :: c, junk
    character*32 :: property

    has_velocities = .false.
    has_masses = .false.

    j = 0
    property = ""
    do i = 1, len(properties)
      c = properties(i:i)

      if( property == "properties" )then
        property = ""
      else if( c == ":" )then
        if( property == "species" )then
          read(line, *) (junk, k = 1, j), species
        else if( property == "pos" .or. property == "positions" )then
          read(line, *) (junk, k = 1, j), positions(1:3)
        else if( property == "vel" .or. property == "velocities" )then
          read(line, *) (junk, k = 1, j), velocities(1:3)
          has_velocities = .true.
        else if( property == "fix_atoms" .or. property == "fix_atom" )then
          read(line, *) (junk, k = 1, j), fix_atom(1:3)
        else if( property == "mass" .or. property == "masses" )then
          read(line, *) (junk, k = 1, j), masses
          has_masses = .true.
        else
!         Advance the pointer by the correct number of fields
          read(property, *, iostat=iostatus) k
          if( iostatus == 0 )then
            j = j + k
          end if
        end if
        property = ""
      else
        property = adjustl(trim(property)) // c
      end if

    end do

  end subroutine
!**************************************************************************

end module
