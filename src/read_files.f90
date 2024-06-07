! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_files.f90, is copyright (c) 2019-2023, Miguel A. Caro and
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

module read_files


  use neighbors
  use types
  use splines
  use vdw
  use soap_turbo_compress_module
  use xyz_module
  use md


  contains


!**************************************************************************
! This subroutine reads in the XYZ file
!
! WE NEED TO WRITE A PROPER EXTXYZ READER THAT CAN IDENTIFY WHICH COLUMN CONTAINS
! EACH PROPERTY. THIS SUBROUTINE CAN ONLY READ IN FILES WITH THE FOLLOWING CONVENTION:
!
! SPECIES X Y Z (VX VY VZ (FIXX FIXY FIXZ))
!
  subroutine read_xyz(filename, ase_format, all_atoms, do_timing, n_species, species_types, &
                      repeat_xyz, rcut_max, which_atom, positions, &
                      do_md, velocities, masses_types, masses, xyz_species, xyz_species_supercell, &
                      species, species_supercell, indices, a_box, b_box, c_box, n_sites, &
                      supercell_check_only, fix_atom, t_beg, write_masses, recalculate_supercell )

    implicit none

!   Input variables
    real*8, intent(in) :: rcut_max, masses_types(:), t_beg
    integer, intent(in) :: which_atom, n_species
    character*8, intent(in) :: species_types(:)
    character*1024, intent(in) :: filename
    logical, intent(in) :: ase_format, all_atoms, do_timing, do_md, supercell_check_only, &
         recalculate_supercell

!   In and out variables
    real*8, allocatable, intent(inout) :: positions(:,:), velocities(:,:), masses(:)
    real*8, intent(inout) :: a_box(1:3), b_box(1:3), c_box(1:3)
    integer, allocatable, intent(inout) :: species(:), species_supercell(:)
    integer, intent(inout) :: n_sites
    integer, intent(inout) :: indices(1:3)
    character*8, allocatable, intent(inout) :: xyz_species(:), xyz_species_supercell(:)
    logical, intent(inout) :: repeat_xyz, write_masses
    logical, allocatable, intent(inout) :: fix_atom(:,:)

!   Internal variables
    real*8, allocatable :: positions_supercell(:,:), velocities_supercell(:,:)
    real*8 :: time1, time2, dist(1:3), read_time, E_kinetic, instant_temp
    real*8 :: kB = 8.6173303d-5, rjunk(1:3), rjunk1d
    integer :: i, iostatus, j, n_sites_supercell, counter, ijunk, k2, i2, j2
    integer :: indices_prev(1:3)
    character*8 :: i_char
    character*128 :: cjunk, cjunk_array(1:100)
    character*1024 :: cjunk1024, properties
    character*12800 :: cjunk_array_flat
    logical :: masses_from_xyz, has_velocities, ljunk(1:3)

    indices_prev = indices

    if( do_timing )then
      call cpu_time(time1)
    end if

if( .not. supercell_check_only )then
    inquire(file=filename, number=i)
    if( i /= 11 )then
      open(unit=11, file=filename, status="old")
    end if
    read(11, *) n_sites
    if( ase_format )then
      read(11, fmt='(A)') cjunk_array_flat
      cjunk_array = ""
      read(cjunk_array_flat, *, iostat=iostatus) cjunk_array(:)
!     Read in lattice vectors
      i = 0
      do
        i = i + 1
        cjunk = cjunk_array(i)
        call upper_to_lower_case(cjunk)
        if( cjunk(1:7) == "lattice" )then
          read(cjunk(10:), *) a_box(1)
          read(cjunk_array(i+1), *) a_box(2)
          read(cjunk_array(i+2), *) a_box(3)
          read(cjunk_array(i+3), *) b_box(1)
          read(cjunk_array(i+4), *) b_box(2)
          read(cjunk_array(i+5), *) b_box(3)
          read(cjunk_array(i+6), *) c_box(1)
          read(cjunk_array(i+7), *) c_box(2)
          cjunk = adjustr(cjunk_array(i+8))
          read(cjunk(1:127), *) c_box(3)
          exit
        end if
      end do
!     Read in properties string
      i = 0
      do
        i = i + 1
        cjunk = cjunk_array_flat(i:i+9)
        call upper_to_lower_case(cjunk)
        if( cjunk == "properties" )then
          j = i+9
          do
            j = j + 1
            if( cjunk_array_flat(j:j) == " " )then
              properties = cjunk_array_flat(i:j-1)
              call upper_to_lower_case(properties)
              exit
            end if
          end do
          exit
        end if
      end do
    else
      write(*,*) "You must use ASE's XYZ format so that I can read the lattice vectors <-- ERROR"
      stop
!    else
!      a_box = 0.d0
!      b_box = 0.d0
!      c_box = 0.d0
!      read(10, *) cjunk, cjunk, a_box(1), junk, junk, junk, b_box(2), junk, junk, junk, c_box(3)
   end if
   if (.not. recalculate_supercell)then
    if( allocated(positions) )deallocate(positions)
    if( allocated(xyz_species) )deallocate(xyz_species)
    if( allocated(species) )deallocate(species)
    allocate( positions(1:3, 1:n_sites) )
    allocate( xyz_species(1:n_sites) )
    allocate( species(1:n_sites) )
    xyz_species = ""
    species = 0
!   We need to comment this out here for nested sampling
!    if( do_md )then
    if( .true. )then
      if( allocated(velocities) )deallocate(velocities)
      if( allocated(masses) )deallocate(masses)
      if( allocated(fix_atom) )deallocate(fix_atom)
      allocate( velocities(1:3, 1:n_sites) )
      velocities = 0.d0
      allocate( masses(1:n_sites) )
      masses = 0.d0
      masses_from_xyz = .false.
      allocate( fix_atom(1:3, 1:n_sites) )
      fix_atom = .false.
    end if
!
!   I should raise a warning if the user is specifying more species than exist in the simulation,
!   since that would lead to arrays that are bigger than necessary                                  <------- Fix this
!    allocate( species(1:max_species_multiplicity, 1:n_sites) )
!    allocate( species_multiplicity(1:n_sites) )
!    species = 0
!    species_multiplicity = 0
    do i = 1, n_sites
      read(11, '(A)') cjunk1024
      if( do_md )then
        call read_xyz_line( properties, cjunk1024, i_char, positions(1:3, i), velocities(1:3, i), fix_atom(1:3, i), &
                            has_velocities, masses(i), masses_from_xyz )
        if( masses_from_xyz )then
          masses(i) = masses(i) * 103.6426965268d0
          write_masses = .true.
        end if
      else
        call read_xyz_line( properties, cjunk1024, i_char, positions(1:3, i), rjunk(1:3), ljunk(1:3), has_velocities, &
                            rjunk1d, masses_from_xyz )
      end if
      do j = 1, n_species
        if( trim(i_char) == trim(species_types(j)) )then
!          species_multiplicity(i) = species_multiplicity(i) + 1
!          species(species_multiplicity(i), i) = j
          xyz_species(i) = species_types(j)
          species(i) = j
!         This is commented out because we also need masses with nested sampling when used in combination with MD
!          if( do_md .and. .not. masses_from_xyz )then
          if( .not. masses_from_xyz )then
            masses(i) = masses_types(j)
          end if
!          exit
        end if
      end do
!      if( species(1, i) == 0 )then
      if( xyz_species(i) == "" )then
        write(*,*)'                                       |'
        write(*,*)'ERROR: atom', i, 'has no known species |  <-- ERROR'
        write(*,*)'                                       |'
        write(*,*)'.......................................|'
        stop
      end if
    end do
!   Randomize velocities if velocities are not provided
    if( do_md .and. .not. has_velocities )then
      write(*,*)'                                       |'
      write(*,*)'WARNING: you have not provided initial |  <-- WARNING'
      write(*,*)'velocities. I am randomizing them so   |'
      write(*,*)'that they match your initial target    |'
      write(*,*)'temperature:                           |'
      write(*,*)'                                       |'
      write(*,'(A, F16.4, A)')' t_beg = ', t_beg, ' K             |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
      call random_number(velocities)
      call remove_cm_vel(velocities(1:3,1:n_sites), masses(1:n_sites))
      E_kinetic = 0.d0
      do i = 1, n_sites
        E_kinetic = E_kinetic + 0.5d0 * masses(i) * dot_product(velocities(1:3, i), velocities(1:3, i))
      end do
      instant_temp = 2.d0/3.d0/dfloat(n_sites-1)/kB*E_kinetic
      velocities = velocities * dsqrt(t_beg/instant_temp)
    end if
!   Check if there are more structures in the xyz file
    read(11, *, iostat=iostatus) cjunk
    if( iostatus == 0 )then
      backspace(11)
      repeat_xyz = .true.
    else
      close(11)
      repeat_xyz = .false.
    end if
    indices_prev = 1
 end if
end if

!   Now we construct a supercell of the required size to accommodate the given rcut_max
!   This needs to be done when the simulation box cannot accommodate one cutoff sphere
    a_box = a_box/dfloat(indices_prev(1))
    b_box = b_box/dfloat(indices_prev(2))
    c_box = c_box/dfloat(indices_prev(3))
    call number_of_unit_cells_for_given_cutoff(a_box, b_box, c_box, rcut_max, [.true., .true., .true.], indices)

    if( .not. supercell_check_only .or. (supercell_check_only .and. any(indices /= indices_prev)) &
         .or. recalculate_supercell )then
    if( indices(1) > 1 .or. indices(2) > 1 .or. indices(3) > 1 )then
      n_sites_supercell = n_sites * indices(1) * indices(2) * indices(3)
      allocate( positions_supercell(1:3, 1:n_sites_supercell) )
!     We need to comment this out here for nested sampling
!      if( do_md )then
      if( .true. )then
        allocate( velocities_supercell(1:3, 1:n_sites_supercell) )
      end if
!      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
!      species_supercell = 0
      if( allocated(xyz_species_supercell) )deallocate(xyz_species_supercell)
      if( allocated(species_supercell) )deallocate(species_supercell)
      allocate( xyz_species_supercell(1:n_sites_supercell) )
      allocate( species_supercell(1:n_sites_supercell) )
      xyz_species_supercell = ""
      species_supercell = 0
      counter = 0
      do i2 = 1, indices(1)
        do j2 = 1, indices(2)
          do k2 = 1, indices(3)
            do i = 1, n_sites
              counter = counter + 1
              positions_supercell(1:3, counter) = positions(1:3, i) + dfloat(i2-1)*a_box(1:3) &
                                                                    + dfloat(j2-1)*b_box(1:3) &
                                                                    + dfloat(k2-1)*c_box(1:3)
!             We need to comment this out here for nested sampling
!              if( do_md )then
              if( .true. )then
                velocities_supercell(1:3, counter) = velocities(1:3, i)
              end if
!              species_supercell(:, counter) = species(:, i)
              xyz_species_supercell(counter) = xyz_species(i)
              species_supercell(counter) = species(i)
            end do
          end do
        end do
      end do
      deallocate( positions )
      allocate( positions(1:3, 1:n_sites_supercell) )
      positions(1:3, 1:n_sites_supercell) = positions_supercell(1:3, 1:n_sites_supercell)
      deallocate( positions_supercell )
!     We need to comment this out here for nested sampling
!      if( do_md )then
      if( .true. )then
        deallocate( velocities )
        allocate( velocities(1:3, 1:n_sites_supercell) )
        velocities(1:3, 1:n_sites_supercell) = velocities_supercell(1:3, 1:n_sites_supercell)
        deallocate( velocities_supercell )
      end if
      a_box = dfloat(indices(1)) * a_box
      b_box = dfloat(indices(2)) * b_box
      c_box = dfloat(indices(3)) * c_box
    else
      n_sites_supercell = n_sites
      if( allocated(xyz_species_supercell) )deallocate(xyz_species_supercell)
      if( allocated(species_supercell) )deallocate(species_supercell)
      allocate( xyz_species_supercell(1:n_sites_supercell) )
      allocate( species_supercell(1:n_sites_supercell) )
      xyz_species_supercell = xyz_species
      species_supercell = species
!      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
!      species_supercell = species
      if( supercell_check_only )then
        allocate( positions_supercell(1:3, 1:n_sites_supercell) )
        positions_supercell = positions(1:3, 1:n_sites_supercell)
        deallocate( positions )
        allocate( positions(1:3, 1:n_sites_supercell) )
        positions(1:3, 1:n_sites_supercell) = positions_supercell(1:3, 1:n_sites_supercell)
        deallocate( positions_supercell )
!       We need to comment this out here for nested sampling
!        if( do_md )then
        if( .true. )then
          allocate( velocities_supercell(1:3, 1:n_sites_supercell) )
          velocities_supercell = velocities(1:3, 1:n_sites_supercell)
          deallocate( velocities )
          allocate( velocities(1:3, 1:n_sites_supercell) )
          velocities(1:3, 1:n_sites_supercell) = velocities_supercell(1:3, 1:n_sites_supercell)
          deallocate( velocities_supercell )
        end if
      end if
    end if

!   This is perhaps not the most efficient way to select only one atom, fix in the future <----- FIX THIS
    if( .not. all_atoms )then
      n_sites = 1
      dist(1:3) = positions(1:3, 1)
      positions(1:3, 1) = positions(1:3, which_atom)
      positions(1:3, which_atom) = dist(1:3)

!      deallocate( species )
!      allocate( species(1:max_species_multiplicity, 1:n_sites) )
!      species(1:max_species_multiplicity, 1) = species_supercell(1:max_species_multiplicity, which_atom)
!      species_supercell(1:max_species_multiplicity, which_atom) = species_supercell(1:max_species_multiplicity, 1)
!      species_supercell(1:max_species_multiplicity, 1) = species(1:max_species_multiplicity, 1)

!      ijunk = species_multiplicity(which_atom)
!      deallocate( species_multiplicity )
!      allocate( species_multiplicity(1:n_sites) )
!      species_multiplicity(1) = ijunk
    end if
else
a_box = a_box*dfloat(indices_prev(1))
b_box = b_box*dfloat(indices_prev(2))
c_box = c_box*dfloat(indices_prev(3))
end if

    if( do_timing )then
      call cpu_time(time2)
      read_time = time2 - time1
      write(*,*)'                                       |'
      write(*,*)'Atoms timings (read):                  |'
      write(*,*)'                                       |'
      write(*,'(A, F11.3, A)') '  *) Read xyz file: ', read_time, ' seconds |'
      write(*,*)'                                       |'
      write(*,*)'.......................................|'
    end if


  end subroutine
!**************************************************************************







!**************************************************************************
  subroutine read_exp_data(file_data, n_points, data)

    implicit none

!   Input variables
    character*1024, intent(in) :: file_data
!   Output variables
    real*8, allocatable, intent(out) :: data(:,:)
    integer, intent(out) :: n_points

!   Internal variables
    integer :: i, j, iostatus, dim, unit_number


    ! if the file_data == none then we allocate and exit
    if ( trim(file_data) == "none" )then
       n_points = 1
       allocate( data(1:2,1:n_points) )
    else

       !   Read data file to figure out data file size
       open(newunit=unit_number, file=file_data, status="old")
       iostatus = 0
       n_points = -1
       do while(iostatus == 0)
          read(unit_number, *, iostat=iostatus)
          n_points = n_points + 1
       end do
       close(unit_number)

       allocate( data(1:2,1:n_points) )
       !     Read local_property data
       open(newunit=unit_number, file=file_data, status="old")
       do i = 1, n_points
          read(unit_number, *)  data(1,i), data(2,i)
       end do
       close(unit_number)
    end if

  end subroutine read_exp_data


  subroutine write_exp_data(x, y, overwrite, filename, label)

    implicit none

!   Input variables
    character(len = *), intent(in) :: filename, label
!   Output variables
    real*8, allocatable, intent(in) :: x(:), y(:)
    logical, intent(in) :: overwrite
!   Internal variables
    integer :: i

    if( overwrite )then
       open(unit=200, file=filename, status="unknown")
       write(200,'(A,1X,A)') '# ', trim(label)
    else
       open(unit=200, file=filename, status="old", position="append")
       write(200,*) ' '
    end if

    do i = 1, size(x)
       write(200, '(1X,F20.8,1X,F20.8)') x(i), y(i)
    end do
    close(200)

  end subroutine write_exp_data

  subroutine write_exp_datan(x, y, overwrite, filename, label)

    implicit none

!   Input variables
    character(len = *), intent(in) :: filename, label
!   Output variables
    real*8, intent(in) :: x(:), y(:)
    logical, intent(in) :: overwrite
!   Internal variables
    integer :: i

    if( overwrite )then
       open(unit=200, file=filename, status="unknown")
       write(200,'(A,1X,A)') '# ', trim(label)
    else
       open(unit=200, file=filename, status="old", position="append")
       write(200,*) ' '
    end if

    do i = 1, size(x)
       write(200, '(1X,F20.8,1X,F20.8)') x(i), y(i)
    end do
    close(200)

  end subroutine write_exp_datan



!**************************************************************************

!**************************************************************************
  subroutine read_alphas_and_descriptors(file_desc, file_alphas, n_sparse, descriptor_type, alphas, Qs, cutoff)

    implicit none

!   Input variables
    character*1024, intent(in) :: file_desc, file_alphas
    character(len=*), intent(in) :: descriptor_type

!   Output variables
    real*8, allocatable, intent(out) :: alphas(:), Qs(:,:), cutoff(:)
    integer, intent(out) :: n_sparse

!   Internal variables
    integer :: i, j, iostatus, dim, unit_number


!   Read alphas to figure out sparse set size
    open(newunit=unit_number, file=file_alphas, status="old")
    iostatus = 0
    n_sparse = -1
    do while(iostatus == 0)
      read(unit_number, *, iostat=iostatus)
      n_sparse = n_sparse + 1
    end do
    close(unit_number)

!   Read descriptor vectors in spare set
    open(newunit=unit_number, file=file_desc, status="old")
    iostatus = 0
    i = -1
    do while(iostatus == 0)
      read(unit_number, *, iostat=iostatus)
      i = i + 1
    end do
    dim = i / n_sparse
    close(unit_number)

!   We do things differently for each descriptor
    if( descriptor_type == "soap_turbo" )then
!     Allocate stuff
      allocate( alphas(1:n_sparse) )
      allocate( Qs(1:dim, 1:n_sparse) )
!     Read alphas SOAP
      open(newunit=unit_number, file=file_alphas, status="old")
      do i = 1, n_sparse
        read(unit_number, *) alphas(i)
      end do
      close(unit_number)
!     Read sparse set descriptors
      open(newunit=unit_number, file=file_desc, status="old")
      do i = 1, n_sparse
        do j = 1, dim
          read(unit_number, *) Qs(j, i)
        end do
      end do
      close(unit_number)
    else if( descriptor_type == "distance_2b" )then
      if( dim /= 1 )then
        write(*,*) "ERROR: Bad 2b descriptor/alphas file(s), dimensions/n_sparse don't match number of data entries"
        stop
      end if
!     Allocate stuff
      allocate( alphas(1:n_sparse) )
      allocate( cutoff(1:n_sparse) )
      allocate( Qs(1:n_sparse, 1:1) )
!     Read alphas 2b and cutoff
      open(newunit=unit_number, file=file_alphas, status="old")
        do i = 1, n_sparse
          read(unit_number, *) alphas(i), cutoff(i)
        end do
      close(unit_number)
!     Read soap vectors in spare set
      open(newunit=unit_number, file=file_desc, status="old")
      do i = 1, n_sparse
        read(unit_number, *) Qs(i, 1)
      end do
      close(unit_number)
    else if( descriptor_type == "angle_3b" )then
      if( dim /= 3 )then
        write(*,*) "ERROR: Bad 3b descriptor/alphas file(s), dimensions/n_sparse don't match number of data entries"
        stop
      end if
!     Allocate stuff. NOTE: the array indices are the opposite of SOAP convention, this makes execution faster
      allocate( alphas(1:n_sparse) )
      allocate( cutoff(1:n_sparse) )
      allocate( Qs(1:n_sparse, 1:3) )
!     Read alphas 3b and cutoff
      open(newunit=unit_number, file=file_alphas, status="old")
      do i = 1, n_sparse
        read(unit_number, *) alphas(i), cutoff(i)
      end do
      close(unit_number)
!     Read soap vectors in spare set
      open(newunit=unit_number, file=file_desc, status="old")
      do i = 1, n_sparse
        do j = 1, 3
          read(unit_number, *) Qs(i,j)
        end do
      end do
      close(unit_number)
    end if

  end subroutine
!**************************************************************************




!**************************************************************************
! This reads the input file
  subroutine read_input_file(n_species, mode, params, rank)

    implicit none

!   Input variables
    integer, intent(in) :: n_species, rank
    character(len=*) :: mode

!   Output variables
    type(input_parameters), intent(out) :: params

!   Internal variables
    real*8 :: c6_ref, r0_ref, alpha0_ref, bsf, k
    integer :: iostatus, i, j, i2,  nw, iostatus2
    character*1024 :: long_line
    character*128, allocatable :: long_line_items(:)
    character*64 :: keyword, cjunk
    character*32 :: implemented_thermostats(1:3)
    character*32 :: implemented_barostats(1:2)
    character*32 :: implemented_mc_types(1:8)
    character*32 :: implemented_exp_observables(1:5)
    character*2 :: element
    character*1 :: keyword_first
    logical :: are_vdw_refs_read(1:3), valid_choice, masses_in_input_file = .false.

    implemented_thermostats(1) = "none"
    implemented_thermostats(2) = "berendsen"
    implemented_thermostats(3) = "bussi"

    implemented_barostats(1) = "none"
    implemented_barostats(2) = "berendsen"

    implemented_mc_types(1) = "none"
    implemented_mc_types(2) = "move"
    implemented_mc_types(3) = "insertion"
    implemented_mc_types(4) = "removal"
    implemented_mc_types(5) = "relax"
    implemented_mc_types(6) = "md"
    implemented_mc_types(7) = "swap"
    implemented_mc_types(8) = "volume"

    implemented_exp_observables(1) = "xps"
    implemented_exp_observables(2) = "xrd"
    implemented_exp_observables(3) = "saxs"
    implemented_exp_observables(4) = "pair_distribution"
    implemented_exp_observables(5) = "structure_factor"


    k = 0.d0

!   Some defaults before reading the input file (the values in the input file will override them)
    if( mode == "md" )then
      params%do_md = .true.
      params%do_prediction = .true.
      params%do_forces = .true.
      params%do_derivatives = .true.
    else if( mode == "mc" )then
      params%do_mc = .true.
      params%do_prediction = .true.
      params%do_forces = .true.
      params%do_derivatives = .true.
    else if( mode == "soap" )then
      params%write_soap = .true.
    else if( mode == "predict" )then
      params%do_prediction = .true.
      params%do_forces = .true.
    end if

!   Let's allocate some arrays:
    allocate( params%species_types(1:n_species) )
    allocate( params%masses_types(1:n_species) )
    allocate( params%radii(1:n_species) )
    allocate( params%e0(1:n_species) )
    allocate( params%vdw_c6_ref(1:n_species) )
    allocate( params%vdw_r0_ref(1:n_species) )
    allocate( params%vdw_alpha0_ref(1:n_species) )
!   Some defaults before reading from file
    params%masses_types = 0.d0
    params%radii = 0.5d0
    params%e0 = 0.d0
    params%vdw_c6_ref = 0.d0
    params%vdw_r0_ref = 0.d0
    params%vdw_alpha0_ref = 0.d0
    are_vdw_refs_read = .false.


!   Read the input file now
    iostatus = 0
    i2 = 0
    do while(iostatus==0)
      read(10, *, iostat=iostatus) keyword
      call upper_to_lower_case(keyword)
      keyword = trim(keyword)
      i2 = len(trim(keyword))
      if(iostatus/=0)then
        exit
      end if
      keyword_first = keyword(1:1)
      if(keyword_first=='#' .or. keyword_first=='!' .or. keyword=='pot_file' .or. keyword=='n_species')then
        continue
      else if(keyword=='do_md')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_md
      else if(keyword=='do_mc')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_mc
      else if(keyword=='verbosity' .or. keyword=='verb')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%verb
      else if(keyword=='do_prediction')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_prediction
      else if(keyword=='do_forces')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_forces
      else if(keyword=='do_derivatives')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_derivatives
      else if(keyword=='do_derivatives_fd')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_derivatives_fd
      else if(keyword=='write_soap')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_soap
      else if(keyword=='write_derivatives')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_derivatives
      else if(keyword=='timing')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_timing
      else if(keyword=='neighbors_buffer')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%neighbors_buffer
      else if(keyword=='t_beg')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%t_beg
      else if(keyword=='t_end')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%t_end
      else if(keyword=='p_beg')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%p_beg
      else if(keyword=='p_end')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%p_end
      else if(keyword=='tau_t')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%tau_t
      else if(keyword=='tau_p')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%tau_p
      else if(keyword=='n_t_hold')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%n_t_hold
        allocate( params%t_hold( 1:params%n_t_hold*3 ) )
      else if(keyword=='t_hold')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, (params%t_hold(nw),nw=1,params%n_t_hold*3)
      else if(keyword=='gamma_p')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%gamma_p
      else if(keyword=='md_step')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%md_step
      else if(keyword=='md_nsteps')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%md_nsteps
      else if(keyword=='mc_nsteps')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_nsteps
      else if(keyword=='n_mc_types')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%n_mc_types
        allocate( params%mc_types(1:params%n_mc_types) )
        allocate( params%mc_acceptance(1:params%n_mc_types) )
        params%mc_acceptance = 1.d0 / params%n_mc_types
      else if(keyword=='n_mc_swaps')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%n_mc_swaps
        allocate( params%mc_swaps(1:2*params%n_mc_swaps) )
        allocate( params%mc_swaps_id(1:2*params%n_mc_swaps) )
      else if(keyword=='mc_swaps')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, (params%mc_swaps(nw),nw=1,2*params%n_mc_swaps)
        !       Need the check the implemented types
        valid_choice = .false.
        do j = 1, 2*params%n_mc_swaps
           valid_choice = .false.
           do i = 1, n_species
              if( trim(params%species_types(i)) == trim(params%mc_swaps(j)) )then
                 params%mc_swaps_id(i) = i
                 valid_choice = .true.
              end if
           end do
           if( .not. valid_choice )then
              if( rank == 0 )then
                 write(*,*) "ERROR -> Invalid mc_swaps species keyword:", params%mc_swaps(j)
                 write(*,*) "This is a list of valid options:"
                 write(*,*) params%species_types
              end if
              stop
           end if
        end do

      else if(keyword=='mc_types')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, (params%mc_types(nw),nw=1,params%n_mc_types)
        !       Need the check the implemented types
        valid_choice = .false.
        do j = 1, params%n_mc_types
           call upper_to_lower_case(params%mc_types(j))
           valid_choice = .false.
           do i = 1, size(implemented_mc_types)
              if( trim(params%mc_types(j)) == trim(implemented_mc_types(i)) )then
                 valid_choice = .true.
              end if
           end do
           if( .not. valid_choice )then
              if( rank == 0 )then
                 write(*,*) "ERROR -> Invalid mc_type keyword:", params%mc_types(j)
                 write(*,*) "This is a list of valid options:"
                 write(*,*) implemented_mc_types
              end if
              stop
           end if
        end do
      else if(keyword=='mc_move_max')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_move_max
      else if(keyword=='mc_min_dist')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_min_dist
      else if(keyword=='mc_lnvol_max')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_lnvol_max
      else if(keyword=='mc_mu')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_mu
      else if(keyword=='mc_species')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_species
      else if(keyword=='mc_write_xyz')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_write_xyz
      else if(keyword=='mc_hamiltonian')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_hamiltonian
      else if(keyword=='mc_relax')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_relax
      else if(keyword=='mc_nrelax')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_nrelax
      else if(keyword=='mc_relax_opt')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_relax_opt
      else if(keyword=='mc_hybrid_opt')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_hybrid_opt
     else if(keyword=='mc_acceptance')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, (params%mc_acceptance(nw),nw=1,params%n_mc_types)
        ! The acceptance probability is based on this sum and normalised
        do i=1, params%n_mc_types
           k = k + params%mc_acceptance(i)
        end do

        do i=1, params%n_mc_types
           params%mc_acceptance(i) = params%mc_acceptance(i) / k
        end do
      else if(keyword=='mc_optimize_exp')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_optimize_exp
      else if(keyword=='mc_reverse')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_reverse
      else if(keyword=='mc_reverse_lambda')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%mc_reverse_lambda
      else if(keyword=='accessible_volume')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%accessible_volume
      else if(keyword=='exp_forces')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%exp_forces
! do experimental
        params%do_exp = .true.

      else if(keyword=='exp_energies')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%exp_energies
! do experimental
        params%do_exp = .true.

     else if(keyword=='exp_energy_scales' .or. keyword&
          &=='exp_energy_scales_initial' .or. keyword&
          &=='exp_energy_scales_beg')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, (params&
             &%exp_energy_scales(nw),nw=1,params&
             &%n_exp)

        ! Set the final gamma to the initial in case
        do nw = 1, params%n_exp
           params%exp_energy_scales_initial(nw) = params%exp_energy_scales(nw)
           params%exp_energy_scales_final(nw) = params%exp_energy_scales(nw)
        end do

      else if(keyword=='exp_energy_scales_final' .or. keyword=='exp_energy_scales_end')then
         backspace(10)
         if (params%n_moments > 0)then
            read(10, *, iostat=iostatus) cjunk, cjunk, (params&
             &%exp_energy_scales_final(nw),nw=1,params&
             &%n_moments)
         else
            read(10, *, iostat=iostatus) cjunk, cjunk, (params&
                 &%exp_energy_scales_final(nw),nw=1,params&
                 &%n_exp)
         end if


      else if(keyword=='exp_input_type')then
         backspace(10)
         read(10, *, iostat=iostatus) cjunk, cjunk, &
              (params%exp_data(nw)%input, nw=1, params%n_exp)

      else if(keyword=='xps_sigma')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xps_sigma
      else if(keyword=='xps_force_type')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xps_force_type
      else if(keyword=='print_lp_forces')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%print_lp_forces
      else if(keyword=='print_vdw_forces')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%print_vdw_forces
      else if(keyword=='exp_similarity_type')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%exp_similarity_type
      else if(keyword=='xrd_alpha')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_alpha
      else if(keyword=='xrd_damping')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_damping
      else if(keyword=='xrd_wavelength')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_wavelength
      else if(keyword=='xrd_method')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_method

      else if(keyword=='nd_wavelength')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%nd_wavelength

     else if(keyword=='xrd_output')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_output

     else if(keyword=='nd_output')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%nd_output

     ! else if(keyword=='xrd_input')then
     !    backspace(10)
     !    read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_output


     else if(keyword=='pair_distribution_output')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_output


      else if(keyword=='xrd_iwasa')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_iwasa
      else if(keyword=='write_xyz')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_xyz
      else if(keyword=='write_thermo')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_thermo
      else if(keyword=='n_nested')then
        if( mode /= "predict" )then
          write(*,*) 'ERROR: the "n_nested" option for nested sampling can only be used with "turbogap predict"'
          stop
        end if
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%n_nested
        if( params%n_nested > 0 )then
          params%do_nested_sampling = .true.
        end if
      else if(keyword=='t_extra')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%t_extra
      else if(keyword=='p_nested')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%p_nested
      else if(keyword=='nested_max_strain')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%nested_max_strain
      else if(keyword=='nested_max_volume_change')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%nested_max_volume_change
      else if(keyword=='scale_box_nested')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%scale_box_nested
      else if(keyword=='n_local_properties')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%n_local_properties
        allocate( params%write_local_properties(1:params%n_local_properties) )
        allocate( params%compute_local_properties(1:params%n_local_properties) )
        params%write_local_properties = .true.
      else if(keyword=='compute_local_properties')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, (params&
             &%compute_local_properties(nw),nw=1 ,params&
            &%n_local_properties)

      else if(keyword=='do_pair_distribution')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_pair_distribution

      else if(keyword=='do_structure_factor')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_structure_factor
        if (params%do_structure_factor)then
           params%do_pair_distribution = .true.
         end if

      else if(keyword=='structure_factor_window')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_window

      else if(keyword=='do_xrd')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_xrd

        if (params%do_xrd)then
           params%do_pair_distribution = .true.
!           params%do_structure_factor = .true.
        end if

      else if(keyword=='do_nd')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_nd

        if (params%do_nd)then
           params%do_pair_distribution = .true.
!           params%do_structure_factor = .true.
        end if

      else if(keyword=='do_exp')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%do_exp

      else if(keyword=='n_exp')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%n_exp
        allocate( params%exp_data(1:params%n_exp) )
        allocate( params%exp_energy_scales(1:params%n_exp) )
        allocate( params%exp_energy_scales_initial(1:params%n_exp) )
        allocate( params%exp_energy_scales_final(1:params%n_exp) )

        ! Turning on exp prediction
        params%do_exp = .true.

      else if( keyword == "exp_labels" )then
         backspace(10)
         read(10, *, iostat=iostatus) cjunk, cjunk, &
              (params%exp_data(nw)%label, nw=1, params%n_exp)
         do nw=1, params%n_exp
            call upper_to_lower_case(params%exp_data(nw)%label)
            if(     trim(params%exp_data(nw)%label) == "xps")then
               params%xps_idx = nw
               if (rank == 0) write(*,*)' - Valid exp. XPS found                |'

            else if(trim(params%exp_data(nw)%label) == "xrd")then
               params%xrd_idx = nw
               params%valid_xrd = .true.
               if (rank == 0) write(*,*)' - Valid exp. XRD found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.
            else if(trim(params%exp_data(nw)%label) == "nd")then
               params%nd_idx = nw
               params%valid_nd = .true.
               if (rank == 0) write(*,*)' - Valid exp. ND found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.

            else if(trim(params%exp_data(nw)%label) == "saxs")then
               params%saxs_idx = nw
               params%valid_xrd = .true.
               if (rank == 0) write(*,*)' - Valid exp. XRD found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.
            else if(trim(params%exp_data(nw)%label) == "pair_distribution")then
               params%pdf_idx = nw
               params%valid_pdf = .true.
               if (rank == 0) write(*,*)' - Valid exp. pair distribution found  |'
            else if(trim(params%exp_data(nw)%label) == "structure_factor")then
               params%sf_idx = nw
               params%valid_sf = .true.
               if (rank == 0) write(*,*)' - Valid exp. structure factor found   |'
            end if
         end do
      else if( keyword == "exp_data_files" )then
         backspace(10)
         read(10, *, iostat=iostatus) cjunk, cjunk, &
              (params%exp_data(nw)%file_data, nw=1, params%n_exp)

         do nw = 1, params%n_exp
            if ( trim( params%exp_data(nw)%file_data ) == "none" )then
               ! Make sure that no type of exp data is written
               params%exp_data(nw)%compute_exp = .false.
               params%exp_data(nw)%compute_similarity = .false.
               ! If the compute exp is false, then a user range must be specified
               params%exp_data(nw)%wrote_exp = .true.
            else

               call read_exp_data(&
                    params%exp_data(nw)%file_data,&
                    params%exp_data(nw)%n_data,&
                    params%exp_data(nw)%data)

               params%exp_data(nw)%compute_exp = .true.
               params%exp_data(nw)%compute_similarity = .true.
               params%exp_data(nw)%range_min = params%exp_data(nw)%data(1,1)
               params%exp_data(nw)%range_max = params&
                    &%exp_data(nw)%data(1,params%exp_data(nw)%n_data)
            end if
         end do

      else if( keyword == "xrd_rcut" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_rcut

      else if( keyword == "nd_rcut" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%nd_rcut

      else if( keyword == "pair_distribution_rcut" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_rcut

      else if( keyword == "pair_distribution_partial" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_partial

      else if( keyword == "structure_factor_from_pdf" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_from_pdf
      else if( keyword == "structure_factor_matrix" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_matrix


      else if( keyword == "pair_distribution_kde_sigma" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_kde_sigma

      else if( keyword == "write_pair_distribution" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_pair_distribution
      else if( keyword == "write_structure_factor" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_structure_factor

      else if( keyword == "write_xrd" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_xrd

      else if( keyword == "write_nd" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_nd

      else if( keyword == "write_exp" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_exp

      else if( keyword == "pair_distribution_n_samples" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_n_samples


      else if( keyword == "structure_factor_n_samples" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_n_samples


      else if( keyword == "xrd_n_samples" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_n_samples
        params%structure_factor_n_samples = params%xrd_n_samples

      else if( keyword == "nd_n_samples" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_n_samples
        params%structure_factor_n_samples = params%nd_n_samples

      else if( keyword == "r_range_min" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%r_range_min

      else if( keyword == "r_range_max" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%r_range_max


      else if( keyword == "q_range_min" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%q_range_min

      else if( keyword == "q_range_max" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%q_range_max

      else if( keyword == "q_units" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%q_units

      else if( keyword == "exp_n_samples" )then
         backspace(10)
         read(10, *, iostat=iostatus) cjunk, cjunk, &
              (params%exp_data(nw)%n_samples, nw=1, params%n_exp)

      else if (keyword(i2-4:i2) == "range" .or.  keyword(i2-8:i2) ==&
           & "file_data" .or. keyword(i2-8:i2) == "n_samples"  )then
         backspace(10)
         ! Check if experimental range or data files are specified
         do nw = 1, params%n_exp
            ! See if the keyword matches any exp observables
            if ( keyword == trim(params%exp_data(nw)%label)//"_range")then
               ! Expect two values which are in order of lower higher for the range to do the prediction
               params%exp_data(nw)%user_range = .true.
               read(10, *, iostat=iostatus) cjunk, cjunk, params%exp_data(nw)%range_min, params%exp_data(nw)%range_max
            elseif ( keyword == trim(params%exp_data(nw)%label)//"_file_data")then

               read(10, *, iostat=iostatus) cjunk, cjunk, params%exp_data(nw)%file_data
               if ( trim( params%exp_data(nw)%file_data ) /= "none" )then

                  call read_exp_data(&
                       params%exp_data(nw)%file_data,&
                       params%exp_data(nw)%n_data,&
                       params%exp_data(nw)%data)

                  params%exp_data(nw)%wrote_exp = .false.
                  params%exp_data(nw)%compute_exp = .true.
                  params%exp_data(nw)%compute_similarity = .true.
                  params%exp_data(nw)%range_min = params%exp_data(nw)%data(1,1)
                  params%exp_data(nw)%range_max = params&
                       &%exp_data(nw)%data(1,params%exp_data(nw)%n_data)
               elseif ( trim( params%exp_data(nw)%file_data ) == "none" )then
                  ! Make sure that no type of exp data is written
                  params%exp_data(nw)%compute_exp = .false.
                  params%exp_data(nw)%compute_similarity = .false.
                  ! If the compute exp is false, then a user range must be specified
                  params%exp_data(nw)%wrote_exp = .true.

               end if
            elseif ( keyword == trim(params%exp_data(nw)%label)//"_n_samples")then
               read(10, *, iostat=iostatus) cjunk, cjunk, params%exp_data(nw)%n_samples
            end if
         end do

      else if(keyword=='write_velocities')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_velocities
      else if(keyword=='write_forces')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_forces
      else if(keyword=='write_fixes')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_fixes
      else if(keyword=='write_stress')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_stress
      else if(keyword=='write_virial')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_virial
      else if(keyword=='write_pressure')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_pressure
      else if(keyword=='write_hirshfeld_v')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_hirshfeld_v
      else if(keyword=='write_local_properties')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, (params&
             &%write_local_properties(nw),nw=1 ,params&
             &%n_local_properties)
      else if(keyword=='write_local_energies')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_local_energies
      else if(keyword=='write_masses')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_masses
      else if(keyword=='target_pos_step')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%target_pos_step
        params%variable_time_step = .true.
      else if(keyword=='tau_dt')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%tau_dt
      else if(keyword=='print_progress')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%print_progress
      else if(keyword=='thermostat')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%thermostat
        call upper_to_lower_case(params%thermostat)
        valid_choice = .false.
        do i = 1, size(implemented_thermostats)
          if( trim(params%thermostat) == trim(implemented_thermostats(i)) )then
            valid_choice = .true.
          end if
        end do
        if( .not. valid_choice )then
          if( rank == 0 )then
            write(*,*) "ERROR -> Invalid thermostat keyword:", params%thermostat
            write(*,*) "This is a list of valid options:"
            write(*,*) implemented_thermostats
          end if
          stop
        end if
      else if(keyword=='barostat')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%barostat
        call upper_to_lower_case(params%barostat)
        valid_choice = .false.
        do i = 1, size(implemented_barostats)
          if( trim(params%barostat) == trim(implemented_barostats(i)) )then
            valid_choice = .true.
          end if
        end do
        if( .not. valid_choice )then
          if( rank == 0 )then
            write(*,*) "ERROR -> Invalid barostat keyword:", params%barostat
            write(*,*) "This is a list of valid options:"
            write(*,*) implemented_barostats
          end if
          stop
        end if
      else if(keyword=='barostat_sym')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%barostat_sym
      else if(keyword=='which_atom')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%which_atom
      else if(keyword=='masses')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%masses_types(1:n_species)
!       We convert the masses in amu to eV*fs^2/A^2
        params%masses_types = params%masses_types * 103.6426965268d0
        masses_in_input_file = .true.
      else if(keyword=='radii')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%radii(1:n_species)
!       We convert the masses in amu to eV*fs^2/A^2
      else if(keyword=='e0')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%e0(1:n_species)
      else if(keyword=='atoms_file' .or. keyword=='input_file')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%atoms_file
      else if(keyword=='max_gbytes_per_process')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%max_Gbytes_per_process
      else if(keyword=='e_tol')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%e_tol
      else if(keyword=='f_tol')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%f_tol
      else if(keyword=='p_tol')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%p_tol
      else if(keyword=='scale_box')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%scale_box
      else if(keyword=='write_lv')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%write_lv
      else if(keyword=='box_scaling_factor')then
        backspace(10)
        read(10, '(A)', iostat=iostatus) long_line
        allocate( long_line_items(1:9) )
        do i = 1, 9
          read(long_line, *, iostat=iostatus2) cjunk, cjunk, long_line_items(1:i)
          if( iostatus2 == -1 )exit
        end do
        i = i - 1
        if( i == 1 )then
          read(long_line_items(1), *) bsf
          params%box_scaling_factor(1,1) = bsf
          params%box_scaling_factor(2,2) = bsf
          params%box_scaling_factor(3,3) = bsf
        else if( i == 3 )then
          read(long_line_items(1), *) bsf
          params%box_scaling_factor(1,1) = bsf
          read(long_line_items(2), *) bsf
          params%box_scaling_factor(2,2) = bsf
          read(long_line_items(3), *) bsf
          params%box_scaling_factor(3,3) = bsf
        else if( i == 9 )then
          read(long_line_items(1), *) bsf
          params%box_scaling_factor(1,1) = bsf
          read(long_line_items(2), *) bsf
          params%box_scaling_factor(1,2) = bsf
          read(long_line_items(3), *) bsf
          params%box_scaling_factor(1,3) = bsf
          read(long_line_items(4), *) bsf
          params%box_scaling_factor(2,1) = bsf
          read(long_line_items(5), *) bsf
          params%box_scaling_factor(2,2) = bsf
          read(long_line_items(6), *) bsf
          params%box_scaling_factor(2,3) = bsf
          read(long_line_items(7), *) bsf
          params%box_scaling_factor(3,1) = bsf
          read(long_line_items(8), *) bsf
          params%box_scaling_factor(3,2) = bsf
          read(long_line_items(9), *) bsf
          params%box_scaling_factor(3,3) = bsf
        else
          write(*,*) "ERROR: the box_scaling_factor must be given by 1, 3 or 9 numbers"
          stop
        end if
        deallocate( long_line_items )
      else if( keyword == "vdw_type" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_type
        call upper_to_lower_case(params%vdw_type)
        if( params%vdw_type == "ts" )then
          continue
        else if( params%vdw_type == "none" )then
          continue
        else
          write(*,*) "ERROR: I do not recognize the vdw_type keyword ", params%vdw_type
          stop
        end if
      else if( keyword == "vdw_sr" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_sr
      else if( keyword == "vdw_d" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_d
      else if( keyword == "vdw_rcut" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_rcut
      else if( keyword == "vdw_buffer" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_buffer
      else if( keyword == "vdw_rcut_inner" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_rcut_inner
      else if( keyword == "vdw_buffer_inner" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_buffer_inner
      else if( keyword == "vdw_c6_ref" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_c6_ref(1:n_species)
        are_vdw_refs_read(1) = .true.
      else if( keyword == "vdw_r0_ref" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_r0_ref(1:n_species)
        are_vdw_refs_read(2) = .true.
      else if( keyword == "vdw_alpha0_ref" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_alpha0_ref(1:n_species)
        are_vdw_refs_read(3) = .true.
      else if( keyword == "vdw_scs_rcut" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_scs_rcut
      else if( keyword == "vdw_mbd_nfreq" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_nfreq
      else if( keyword == "vdw_mbd_grad" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%vdw_mbd_grad
      else if( keyword == "core_pot_cutoff" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%core_pot_cutoff
      else if( keyword == "core_pot_buffer" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%core_pot_buffer
      else if( keyword == "optimize" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%optimize
        if( params%optimize == "vv" .or. params%optimize == "gd" .or. params%optimize == "gd-box" .or. &
            params%optimize == "gd-box-ortho" )then
          continue
        else
          write(*,*) "ERROR: optimize algorithm not implemented:", params%optimize
          stop
        end if
      else if( keyword == "gamma0" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%gamma0
      else if( keyword == "max_opt_step" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%max_opt_step
      else if( keyword == "max_opt_step_eps" )then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%max_opt_step_eps
      else if(keyword=='species')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, cjunk, params%species_types(1:n_species)
        if( iostatus > 0 )then
          write(*,*)'                                       |'
          write(*,*)'ERROR: your "species" keyword is wrong |  <-- ERROR'
          stop
        end if
      else
        write(*,*)"ERROR: I do not recognize the input file keyword", keyword
        stop
     end if
    end do

!   Do some checks
    if( params%write_xyz == 0 )then
      params%write_xyz = params%md_nsteps
    end if

    if( params%do_md )then
      params%do_prediction = .true.
      params%do_forces = .true.
      params%which_atom = 0
    end if

    if( params%do_forces )then
      params%do_derivatives = .true.
    end if

    if( params%which_atom /= 0 )then
      params%all_atoms = .false.
    else
      params%all_atoms = .true.
    end if

    do i = 1, n_species
      call get_vdw_ref_params( params%species_types(i), c6_ref, r0_ref, alpha0_ref, rank )
      if( .not. are_vdw_refs_read(1) )then
        params%vdw_c6_ref(i) = c6_ref
      end if
      if( .not. are_vdw_refs_read(2) )then
        params%vdw_r0_ref(i) = r0_ref
      end if
      if( .not. are_vdw_refs_read(3) )then
        params%vdw_alpha0_ref(i) = alpha0_ref
      end if
    end do

!   If we don't use van der Waals, then unset the default cutoff
    if( params%vdw_type == "none" )then
      params%vdw_rcut = 0.d0
    else
!     If van der Waals is enabled, make sure the inner and outer cutoff regions do not overlap
!     and other sanity checks
      if( params%vdw_rcut - params%vdw_buffer < params%vdw_rcut_inner + params%vdw_buffer_inner )then
        write(*,*) "ERROR: vdW inner and outer cutoff regions can't overlap. Check your vdw_* definitions"
        stop
      end if
    end if

!   Nested sampling checks
    if( params%do_nested_sampling )then
      if( params%thermostat /= "none" )then
        write(*,*)'                                       |'
        write(*,*)'WARNING: Nested sampling only works    |  <-- WARNING'
        write(*,*)'(currently) in combination with total  |'
        write(*,*)'energy MD. The selected thermostat has |'
        write(*,*)'been disabled.                         |'
      end if
!     Prepare directory where we create the latest version of the walkers
      call system("rm -rf walkers/")
      call system("mkdir -p walkers/")
    end if


!   Experimental prediction checks
    if( params%do_exp )then
       if (rank == 0) write(*,*)'                                       |'
       if (rank == 0) write(*,*)' Experimental prediction mode          |'
       do i = 1, params%n_exp
          ! check if a user range has been submitted
          write(*,*)'                                       |'

          if (params%exp_data(i)%user_range)then
             if (rank == 0) write(*,'(A,1X,A,1X,A)')'User exp. range specified for:', trim(params%exp_data(i)%label),'     |'
             if (rank == 0) write(*,*)'                                       |'
             if (rank == 0) write(*,*)' WARNING!! This feature is obselete    |'
             if (rank == 0) write(*,*)'                                       |'
          else
             if (rank == 0) write(*,'(A,1X,A,1X,A)')'Exp data range will be used for:', trim(params%exp_data(i)%label),' |'
             if (rank == 0) write(*,'(A,1X,A,1X,A)')' from the file:', trim(params%exp_data(i)%file_data),' |'

          end if


          if( params%exp_data(i)%range_min == 0.d0 .and. params%exp_data(i)%range_max == 1.d0 )then
             if (rank == 0) write(*,*)'                                       |'
             if (rank == 0) write(*,*)'WARNING: Data range being used for exp.|'
             if (rank == 0) write(*,*)' observable is the default (0.0, 1.0)! |'
             if (rank == 0) write(*,*)'                                       |'
             if (rank == 0) write(*,*)' To modify specify:                    |'
             if (rank == 0) write(*,'(A,1X,A,1X,A)')'  `range_',trim(params%exp_data(i)%label),  ' = {lower_bound} {upper_bound}` |'
             if (rank == 0) write(*,*)' in the input file.                    |'
             if (rank == 0) write(*,*)'                                       |'
          end if

          if ( trim(params%exp_data(i)%label) == 'pair_distribution')then
             if (rank == 0) write(*,'(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                  & ' found, setting r_range_min/max ', '     |'
             ! Note: for consistency with the implementation, we can
             ! change the value of r_min/r_max such that the x_i
             ! generated
             ! by the bin_edges of the pair_distribution function
             ! match those
             ! of the actual experimental data

             params%do_pair_distribution = .true.

             params%r_range_min = params%exp_data(i)%range_min - &
                  & ( params%exp_data(i)%range_max - params%exp_data(i)%range_min ) / &
                  & ( dfloat( 2 * (params%exp_data(i)%n_samples - 1) ) )

             params%r_range_max = params%exp_data(i)%range_min + &
                  & ( dfloat( 2 * params%exp_data(i)%n_samples - 1 ) * &
                  & ( params%exp_data(i)%range_max - params%exp_data(i)%range_min ) / &
                  & ( dfloat( 2 * (params%exp_data(i)%n_samples - 1) ) ) )

             params%pair_distribution_n_samples = params%exp_data(i)%n_samples
          elseif ( trim(params%exp_data(i)%label) == 'xrd')then
             if (rank == 0) write(*,'(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                  & ' found, setting q_range_min/max with q_units = ' // trim(params%q_units) , ' |'

             params%do_pair_distribution = .true.
             params%pair_distribution_partial = .true.
             params%do_structure_factor = .true.
             params%structure_factor_from_pdf = .true.
             params%do_xrd = .true.

             params%q_range_min = params%exp_data(i)%range_min
             params%q_range_max = params%exp_data(i)%range_max
             ! params%q_units = 'twotheta'
             params%xrd_n_samples = params%exp_data(i)%n_samples
             params%structure_factor_n_samples = params%exp_data(i)%n_samples

          elseif ( trim(params%exp_data(i)%label) == 'nd')then
             if (rank == 0) write(*,'(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                  & ' found, setting q_range_min/max with q_units = ' // trim(params%q_units) , ' |'

             params%do_pair_distribution = .true.
             params%pair_distribution_partial = .true.
             params%do_structure_factor = .true.
             params%structure_factor_from_pdf = .true.
             params%do_nd = .true.

             params%q_range_min = params%exp_data(i)%range_min
             params%q_range_max = params%exp_data(i)%range_max
             ! params%q_units = 'twotheta'
             params%nd_n_samples = params%exp_data(i)%n_samples
             params%structure_factor_n_samples = params%exp_data(i)%n_samples

          elseif ( trim(params%exp_data(i)%label) == 'saxs')then
             if (rank == 0) write(*,'(A,1X,A,1X,A)') trim(params&
                  &%exp_data(i)%label), ' found, setting q_range_min&
                  &/max with q_units = "q"', ' |'

             params%do_pair_distribution = .true.
             params%pair_distribution_partial = .true.
             params%do_structure_factor = .true.
             params%structure_factor_from_pdf = .true.
             params%do_xrd = .true.

             params%q_range_min = params%exp_data(i)%range_min
             params%q_range_max = params%exp_data(i)%range_max
             params%q_units = 'q'
             params%xrd_n_samples = params%exp_data(i)%n_samples
             params%structure_factor_n_samples = params%exp_data(i)%n_samples
          elseif ( trim(params%exp_data(i)%label) == 'structure_factor')then
             if (rank == 0) write(*,'(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                  & ' found, setting q_range_min/max with q_units =&
                  & "q"', ' |'

             params%do_pair_distribution = .true.
             params%pair_distribution_partial = .true.
             params%do_structure_factor = .true.
             params%structure_factor_from_pdf = .true.

             params%q_range_min = params%exp_data(i)%range_min
             params%q_range_max = params%exp_data(i)%range_max
             params%q_units = 'q'
             params%structure_factor_n_samples = params%exp_data(i)%n_samples
             params%xrd_n_samples = params%exp_data(i)%n_samples
          end if


          if (rank == 0) write(*,'(A,1X,F12.6,1X,A,F12.6,1X,A)')' min =', params&
               &%exp_data(i)%range_min, ' max =', params%exp_data(i)&
               &%range_max, ' |'

          if (rank == 0) write(*,'(A,1X,I8,1X,A)')' n_samples   =', params%exp_data(i)%n_samples,'                |'
          if (rank == 0) write(*,'(A,1X,L4,1X,A)')' compute_exp =', params%exp_data(i)%compute_exp,'                    |'


          if (.not. allocated(params%exp_energy_scales) .and. ( params%exp_forces .or. params%mc_optimize_exp ) )then
             if (rank == 0) write(*,*)'WARNING: No energy scales set for exp .|'
             if (rank == 0) write(*,*)' optimisation by forces / MC!          |'
             if (rank == 0) write(*,*)'                                       |'
             if (rank == 0) write(*,*)' To modify specify:                    |'
             if (rank == 0) write(*,'(A)')'  `exp_energy_scales = {E1} {E2}`  |'
             if (rank == 0) write(*,*)' In the input file.                    |'
             if (rank == 0) write(*,*)' (example above is for n_exp = 2)      |'
             if (rank == 0) write(*,*)'                                       |'
          end if

       end do
    end if


!   Monte-carlo checks
    if( params%do_mc )then
       do i = 1, params%n_mc_types
          if (params%mc_types(i) == "md")then
             if( params%thermostat == "none" )then
                if (rank == 0) write(*,*)'                                       |'
                if (rank == 0) write(*,*)'WARNING: You need to specify a         |  <-- WARNING'
                if (rank == 0) write(*,*)'thermostat when using md type mc steps!|'
             end if
          end if

          if (params%mc_types(i) == "relax")then
             if( params%optimize == "none" )then
                if (rank == 0) write(*,*)'                                       |'
                if (rank == 0) write(*,*)'WARNING: You need to specify an        |  <-- WARNING'
                if (rank == 0) write(*,*)'optimizer when using relax type mc     |'
                if (rank == 0) write(*,*)'steps!!                                |'
             end if
          end if

          if (params%mc_types(i) == "volume")then
             if( params%p_beg == 1.0d0 )then
                if (rank == 0) write(*,*)'                                       |'
                if (rank == 0) write(*,*)'WARNING: p_beg is the default          |  <-- WARNING'
                if (rank == 0) write(*,*)'value of 1.0 bar. For MC volume moves  |'
                if (rank == 0) write(*,*)'please make sure this is specified!!   |'
             end if
             if( params%mc_lnvol_max == 0.01d0 )then
                if (rank == 0) write(*,*)'                                       |'
                if (rank == 0) write(*,*)'WARNING: mc_lnvol_max is the default   |  <-- WARNING'
                if (rank == 0) write(*,*)'value of 0.01. For MC volume moves     |'
                if (rank == 0) write(*,*)'please make sure this is specified!!   |'
             end if

          end if
       end do

       do i = 1, n_species
          if( params%accessible_volume .and. (params%radii(i) == 0.5d0 ))then
             if (rank == 0) write(*,*)'                                       |'
             if (rank == 0) write(*,*)'WARNING: radii for accessible volume   |  <-- WARNING'
             if (rank == 0) write(*,*)'is the default value of 0.5A.          |'
             if (rank == 0) write(*,*)'please make sure this correct!!        |'
          end if
       end do


    end if



!   Set the writeouts
    if( .not. params%do_md )then
!     Do not write temperature
      params%write_property(3) = .false.
!     Do not write pressure
      params%write_property(4) = .false.
!     Do not write time step
      params%write_property(5) = .false.
!     Do not write time
      params%write_property(6) = .false.
!     Do not write MD step
      params%write_property(11) = .false.
!     Do not write velocities
      params%write_array_property(3) = .false.
!     Do not write masses
      params%write_array_property(6) = .false.
!     Do not write fixes
      params%write_array_property(8) = .false.
    end if
    if( .not. params%do_forces )then
!     Do not write pressure
      params%write_property(4) = .false.
!     Do not write virial
      params%write_property(8) = .false.
!     Do not write stress
      params%write_property(9) = .false.
!     Do not write forces
      params%write_array_property(4) = .false.
    end if
    if( params%vdw_type == "none" )then
!     Do not write Hirshfeld volume
      params%write_array_property(7) = .false.
    end if
!   Now individual flags
    if( .not. params%write_velocities )then
      params%write_array_property(3) = .false.
    end if
    if( .not. params%write_forces )then
      params%write_array_property(4) = .false.
    end if
    if( .not. params%write_local_energies )then
      params%write_array_property(5) = .false.
    end if
    if( .not. params%write_masses )then
      params%write_array_property(6) = .false.
    end if
    if( .not. params%write_fixes )then
      params%write_array_property(8) = .false.
    end if
    if( .not. params%write_pressure )then
      params%write_property(7) = .false.
    end if
    if( .not. params%write_virial )then
      params%write_property(8) = .false.
    end if
    if( .not. params%write_stress )then
      params%write_property(9) = .false.
    end if

!   Get masses from database
    if( params%do_md .and. .not. masses_in_input_file )then
      if( rank == 0 )then
        write(*,*)'                                       |'
        write(*,*)'WARNING: you have not provided masses  |  <-- WARNING'
        write(*,*)'in your input file. I am attempting to |'
        write(*,*)'read them from a database. If you have |'
        write(*,*)'provided masses in your XYZ file these |'
        write(*,*)'values will be overwritten and you can |'
        write(*,*)'safely disregard any further warnings  |'
        write(*,*)'printed below if a given element is not|'
        write(*,*)'in the database (usually because you   |'
        write(*,*)'provided a non-standard name; note that|'
        write(*,*)'element names are case sensitive).     |'
        write(*,*)'                                       |'
        write(*,*)'               Element      Mass (amu) |'
      end if
      do i = 1, n_species
        call get_atomic_mass( params%species_types(i), params%masses_types(i), valid_choice )
        if( rank == 0 )then
          write(*,*)'                                       |'
          if( valid_choice )then
            write(*,'(A, A8, A, F15.6, A)')' ', adjustr(params%species_types(i)), ' (in database) ', params%masses_types(i), ' |'
          else
            write(*,'(A, A8, A, F11.6, A)')' ', adjustr(params%species_types(i)), ' (not in database) ', params%masses_types(i), &
                                           ' |  <-- WARNING'
          end if
        end if
      end do
!     We convert the masses in amu to eV*fs^2/A^2
      params%masses_types = params%masses_types * 103.6426965268d0
      if( rank == 0 )then
        write(*,*)'                                       |'
        write(*,*)'.......................................|'
      end if
    end if

  end subroutine
!**************************************************************************




!**************************************************************************
  subroutine read_gap_hypers(file_gap, &
                             n_soap_turbo, soap_turbo_hypers, &
                             n_distance_2b, distance_2b_hypers, &
                             n_angle_3b, angle_3b_hypers, &
                             n_core_pot, core_pot_hypers, &
                             rcut_max, do_prediction, params )

    implicit none

!   Input variables
    type(input_parameters), intent(in) :: params
    logical, intent(in) :: do_prediction
    character(len=*), intent(in) :: file_gap

!   Output variables
    real*8, intent(out) :: rcut_max
    integer, intent(out) :: n_soap_turbo, n_distance_2b, n_angle_3b, n_core_pot
    integer :: nw
    type(soap_turbo), allocatable, intent(out) :: soap_turbo_hypers(:)
    type(distance_2b), allocatable, intent(out) :: distance_2b_hypers(:)
    type(angle_3b), allocatable, intent(out) :: angle_3b_hypers(:)
    type(core_pot), allocatable, intent(out) :: core_pot_hypers(:)
!   Internal variables
    real*8, allocatable :: u(:), x(:), V(:)
    real*8 :: sig, p, qn, un
    integer :: iostatus, i, counter, n_species, n_sparse, ijunk, n, n_nonzero, j
    character*64 :: keyword, cjunk, compress_string
    character*1 :: keyword_first

    open(unit=10, file=file_gap, status="old", iostat=iostatus)
!   Look for the number of instances of each GAP
    n_soap_turbo = 0
    n_distance_2b = 0
    n_angle_3b = 0
    n_core_pot = 0
    rcut_max = 0.d0
    iostatus = 0
    do while(iostatus==0)
      read(10, *, iostat=iostatus) keyword
      keyword = trim(keyword)
      if(iostatus/=0)then
        exit
      end if
      keyword_first = keyword(1:1)
      if(keyword_first=='#' .or. keyword_first=='!')then
        continue
      else if(keyword=='gap_beg')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, keyword
        if( keyword == "soap_turbo" )then
          n_soap_turbo = n_soap_turbo + 1
        else if( keyword == "distance_2b" )then
          n_distance_2b = n_distance_2b + 1
        else if( keyword == "angle_3b" )then
          n_angle_3b = n_angle_3b + 1
        else if( keyword == "core_pot" )then
          n_core_pot = n_core_pot + 1
        end if
      end if
    end do
!   Allocate the variables
    if( n_soap_turbo > 0 )then
      allocate( soap_turbo_hypers(1:n_soap_turbo) )
    end if
    if( n_distance_2b > 0 )then
      allocate( distance_2b_hypers(1:n_distance_2b) )
    end if
    if( n_angle_3b > 0 )then
      allocate( angle_3b_hypers(1:n_angle_3b) )
    end if
    if( n_core_pot > 0 )then
      allocate( core_pot_hypers(1:n_core_pot) )
    end if

!   Now record the hypers
    rewind(10)
    n_soap_turbo = 0
    n_distance_2b = 0
    n_angle_3b = 0
    n_core_pot = 0
    iostatus = 0
    do while(iostatus==0)
      read(10, *, iostat=iostatus) keyword
      keyword = trim(keyword)
      if(iostatus/=0)then
        exit
      end if
      keyword_first = keyword(1:1)
      if(keyword_first=='#' .or. keyword_first=='!')then
        continue
      else if(keyword=='gap_beg')then
        backspace(10)
        read(10, *, iostat=iostatus) cjunk, keyword
!       soap_turbo definitions here
        if( keyword == "soap_turbo" )then
          n_soap_turbo = n_soap_turbo + 1
          counter = 0
          do while(iostatus == 0)
            read(10, *, iostat=iostatus) keyword
            counter = counter + 1
            if( keyword == "n_species" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, n_species
              soap_turbo_hypers(n_soap_turbo)%n_species = n_species
              allocate( soap_turbo_hypers(n_soap_turbo)%nf(1:n_species) )
              soap_turbo_hypers(n_soap_turbo)%nf = 4.d0
              allocate( soap_turbo_hypers(n_soap_turbo)%rcut_hard(1:n_species) )
              allocate( soap_turbo_hypers(n_soap_turbo)%rcut_soft(1:n_species) )
              soap_turbo_hypers(n_soap_turbo)%rcut_soft = 0.5d0
              allocate( soap_turbo_hypers(n_soap_turbo)%atom_sigma_r(1:n_species) )
              allocate( soap_turbo_hypers(n_soap_turbo)%atom_sigma_t(1:n_species) )
              allocate( soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling(1:n_species) )
              soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling = 0.d0
              allocate( soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling(1:n_species) )
              soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling = 0.d0
              allocate( soap_turbo_hypers(n_soap_turbo)%amplitude_scaling(1:n_species) )
              soap_turbo_hypers(n_soap_turbo)%amplitude_scaling = 1.d0
              allocate( soap_turbo_hypers(n_soap_turbo)%central_weight(1:n_species) )
              soap_turbo_hypers(n_soap_turbo)%central_weight = 1.d0
              allocate( soap_turbo_hypers(n_soap_turbo)%global_scaling(1:n_species) )
              soap_turbo_hypers(n_soap_turbo)%global_scaling = 1.d0
              allocate( soap_turbo_hypers(n_soap_turbo)%alpha_max(1:n_species) )
              allocate( soap_turbo_hypers(n_soap_turbo)%species_types(1:n_species) )
              do i = 1, counter
                backspace(10)
              end do
              exit
            end if
          end do
          iostatus = 0
          do while( keyword /= "gap_end" .and. iostatus == 0 )
            read(10, *, iostat=iostatus) keyword
            if( keyword == "nf" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%nf(1:n_species)
            else if( keyword == "rcut" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%rcut_hard(1:n_species)
              soap_turbo_hypers(n_soap_turbo)%rcut_max = 0.d0
              do i = 1, n_species
                if( soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) > soap_turbo_hypers(n_soap_turbo)%rcut_max )then
                  soap_turbo_hypers(n_soap_turbo)%rcut_max = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
                end if
              end do
            else if( keyword == "buffer" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%rcut_soft(1:n_species)
            else if( keyword == "atom_sigma_r" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_r(1:n_species)
            else if( keyword == "atom_sigma_t" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_t(1:n_species)
            else if( keyword == "atom_sigma_r_scaling" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling(1:n_species)
            else if( keyword == "atom_sigma_t_scaling" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling(1:n_species)
            else if( keyword == "amplitude_scaling" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%amplitude_scaling(1:n_species)
            else if( keyword == "central_weight" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%central_weight(1:n_species)
            else if( keyword == "global_scaling" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%global_scaling(1:n_species)
            else if( keyword == "n_max" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%alpha_max(1:n_species)
!             But we can actually use n_max to referred to the "total" radial basis (the sum of orthogonal bases
!             for different species)
              soap_turbo_hypers(n_soap_turbo)%n_max = 0
              do i = 1, n_species
                soap_turbo_hypers(n_soap_turbo)%n_max = soap_turbo_hypers(n_soap_turbo)%n_max + &
                                                        soap_turbo_hypers(n_soap_turbo)%alpha_max(i)
              end do
            else if( keyword == "species" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%species_types(1:n_species)
            else if( keyword == "l_max" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%l_max
            else if( keyword == "radial_enhancement" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%radial_enhancement
              if( soap_turbo_hypers(n_soap_turbo)%radial_enhancement < 0 .or. &
                soap_turbo_hypers(n_soap_turbo)%radial_enhancement > 2 )then
                write(*,*)'                                       |'
                write(*,*)'WARNING: radial_enhancement must be    |  <-- WARNING'
                write(*,*)'and 0 <= n <= 2. I am defaulting to 0! |'
                write(*,*)'                                       |'
                write(*,*)'.......................................|'
                soap_turbo_hypers(n_soap_turbo)%radial_enhancement = 0
              end if
            else if( keyword == "compress_soap" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%compress_soap
            else if( keyword == "file_compress_soap" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_compress
            else if( keyword == "compress_mode" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%compress_mode
            else if( keyword == "zeta" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%zeta
            else if( keyword == "delta" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%delta
            else if( keyword == "central_species" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%central_species
            else if( keyword == "basis" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%basis
              if( soap_turbo_hypers(n_soap_turbo)%basis /= "poly3" .and. &
                soap_turbo_hypers(n_soap_turbo)%basis /= "poly3gauss" )then
                write(*,*)'                                       |'
                write(*,*)'WARNING: I didn''t understand your      |  <-- WARNING'
                write(*,*)'keywork for basis; defaulting to       |'
                write(*,*)'"poly3"                                |'
                soap_turbo_hypers(n_soap_turbo)%basis = "poly3"
              end if
            else if( keyword == "scaling_mode" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%scaling_mode
              if( soap_turbo_hypers(n_soap_turbo)%scaling_mode /= "polynomial" )then
                write(*,*)'                                       |'
                write(*,*)'WARNING: I didn''t understand your      |  <-- WARNING'
                write(*,*)'keywork for scaling_mode; defaulting   |'
                write(*,*)'to "polynomial"                        |'
                soap_turbo_hypers(n_soap_turbo)%scaling_mode = "polynomial"
              end if
            else if( keyword == "alphas_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_alphas
            else if( keyword == "desc_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_desc
            else if( keyword == "has_vdw" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_vdw
            else if( keyword == "vdw_qs" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_vdw_desc
            else if( keyword == "vdw_alphas" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_vdw_alphas
            else if( keyword == "vdw_zeta" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_zeta
            else if( keyword == "vdw_delta" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_delta
            else if( keyword == "vdw_v0" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_v0
            else if( keyword == "has_local_properties" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_local_properties
            else if( keyword == "n_local_properties" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%n_local_properties
              ! Now allocate the local_property_soap_turbo object in the soap_turbo_hypers
              allocate( soap_turbo_hypers(n_soap_turbo)%local_property_models(&
                   1:soap_turbo_hypers(n_soap_turbo)%n_local_properties) )


            else if( keyword == "local_property_labels" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, &
                   (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                   &%label,nw=1&
                   &,soap_turbo_hypers(n_soap_turbo)%n_local_properties)
              do nw=1, soap_turbo_hypers(n_soap_turbo)%n_local_properties
                 if(trim(soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                      &%label) == "hirshfeld_v")then
                    soap_turbo_hypers(n_soap_turbo)%has_vdw=.true.
                    soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)%do_derivatives=.true.
                    soap_turbo_hypers(n_soap_turbo)%vdw_index=nw
                 end if

                 if(trim(soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                      &%label) == "core_electron_be")then
                    soap_turbo_hypers(n_soap_turbo)%has_core_electron_be=.true.
                    soap_turbo_hypers(n_soap_turbo)%core_electron_be_index=nw
                 end if

              end do

            else if( keyword == "local_property_qs" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, &
                   (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                   &%file_desc,nw=1&
                   &,soap_turbo_hypers(n_soap_turbo)%n_local_properties)
            else if( keyword == "local_property_alphas" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, &
                   (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                   &%file_alphas,nw=1&
                   &,soap_turbo_hypers(n_soap_turbo)%n_local_properties)
            else if( keyword == "local_property_zetas" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk,&
                   & (soap_turbo_hypers(n_soap_turbo)&
                   &%local_property_models(nw)%zeta,nw=1&
                   &,soap_turbo_hypers(n_soap_turbo)%n_local_properties)
            else if( keyword == "local_property_deltas" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, &
                   & (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                   &%delta ,nw=1,soap_turbo_hypers(n_soap_turbo)&
                   &%n_local_properties)
            else if( keyword == "local_property_v0s" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, &
                   & (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                   &%V0 ,nw=1,soap_turbo_hypers(n_soap_turbo)&
                   &%n_local_properties)
            end if
          end do
!         We actually read in the "buffer" zone width, so transform to rcut_soft:
          do i = 1, n_species
            if( soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) == 0.d0 )then
              soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
            else
              soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) - &
                                                             soap_turbo_hypers(n_soap_turbo)%rcut_soft(i)
            end if
          end do
!         Read the sparse set information
          if( do_prediction )then
            call read_alphas_and_descriptors(soap_turbo_hypers(n_soap_turbo)%file_desc, &
                                             soap_turbo_hypers(n_soap_turbo)%file_alphas, &
                                             soap_turbo_hypers(n_soap_turbo)%n_sparse, &
                                             "soap_turbo", soap_turbo_hypers(n_soap_turbo)%alphas, &
                                             soap_turbo_hypers(n_soap_turbo)%Qs, &
                                             soap_turbo_hypers(n_soap_turbo)%cutoff)
            ! Commenting this out as it will be subsumed into local property prediction
            ! if( soap_turbo_hypers(n_soap_turbo)%has_vdw )then
            !    call read_alphas_and_descriptors(soap_turbo_hypers(n_soap_turbo)%file_vdw_desc, &
            !         soap_turbo_hypers(n_soap_turbo)%file_vdw_alphas, &
            !         soap_turbo_hypers(n_soap_turbo)%vdw_n_sparse, &
            !         "soap_turbo", soap_turbo_hypers(n_soap_turbo)%vdw_alphas, &
            !         soap_turbo_hypers(n_soap_turbo)%vdw_Qs, &
            !         soap_turbo_hypers(n_soap_turbo)%vdw_cutoff)

            ! end if

            if( soap_turbo_hypers(n_soap_turbo)%has_local_properties )then
               do j=1, soap_turbo_hypers(n_soap_turbo)%n_local_properties

                  call read_alphas_and_descriptors(&
                       soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%file_desc, &
                       soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%file_alphas, &
                       soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%n_sparse, &
                       "soap_turbo", &
                       soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%alphas, &
                       soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%Qs, &
                       soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%cutoff)

                  ! Really, this could actually just not be associated
                  ! with the soap turbo type as the same data might be
                  ! reread into separate soap turbo descriptors when
                  ! only one is needed, and further this is
                  ! broadcasted. But this way, all the files are
                  ! specified in the .gap file rather than in the
                  ! input file.


                  ! soap_turbo_hypers(n_soap_turbo)&
                  !      &%local_property_models(j)%dim =&
                  !      & size(soap_turbo_hypers(n_soap_turbo)&
                  !      &%local_property_models(j)%Qs,1)



               end do
            end if

          end if
          do i = 1, n_species
            if( soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) > rcut_max )then
              rcut_max = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
            end if
          end do
!         Handle SOAP compression here
!         Here we read in the compression information from a file (compress_file) or rely on a keyword provided
!         by the user (compress_mode) which leads to a predefined recipe to compress the soap_turbo descriptor
!         The file always takes precedence over the keyword.
          if( soap_turbo_hypers(n_soap_turbo)%compress_soap )then
!           A compress file takes priority over compress mode
            if( soap_turbo_hypers(n_soap_turbo)%file_compress /= "none" )then
              open(unit=20, file=soap_turbo_hypers(n_soap_turbo)%file_compress, status="old")
              read(20, *) (ijunk, i=1,n_species), ijunk, soap_turbo_hypers(n_soap_turbo)%dim
!             This enables definition of arbitrary compression transformations via a file
              read(20, '(A)') compress_string
              if( compress_string == "P_transformation" )then
                n_nonzero = -1
                do while( compress_string /= "end_transformation" )
                  read(20, '(A)') compress_string
                  n_nonzero = n_nonzero + 1
                end do
                soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero = n_nonzero
                allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:n_nonzero) )
                allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:n_nonzero) )
                allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:n_nonzero) )
                do i = 1, n_nonzero+1
                  backspace(20)
                end do
                do i = 1, n_nonzero
                  read(20,*) soap_turbo_hypers(n_soap_turbo)%compress_P_i(i), &
                             soap_turbo_hypers(n_soap_turbo)%compress_P_j(i), &
                             soap_turbo_hypers(n_soap_turbo)%compress_P_el(i)
                end do
              else
!               Old way to handle compression for backcompatibility
                backspace(20)
                soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero = soap_turbo_hypers(n_soap_turbo)%dim
                allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:soap_turbo_hypers(n_soap_turbo)%dim) )
                allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:soap_turbo_hypers(n_soap_turbo)%dim) )
                allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:soap_turbo_hypers(n_soap_turbo)%dim) )
                do i = 1, soap_turbo_hypers(n_soap_turbo)%dim
                  read(20, *) soap_turbo_hypers(n_soap_turbo)%compress_P_j(i)
                  soap_turbo_hypers(n_soap_turbo)%compress_P_i(i) = i
                  soap_turbo_hypers(n_soap_turbo)%compress_P_el(i) = 1.d0
                end do
              end if
              close(20)
            else if( soap_turbo_hypers(n_soap_turbo)%compress_mode /= "none" )then
              call get_compress_indices( soap_turbo_hypers(n_soap_turbo)%compress_mode, &
                                         soap_turbo_hypers(n_soap_turbo)%alpha_max, &
                                         soap_turbo_hypers(n_soap_turbo)%l_max, &
                                         soap_turbo_hypers(n_soap_turbo)%dim, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_i, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_j, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_el, &
                                         "get_dim" )
              allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero) )
              allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero) )
              allocate( soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero) )
              call get_compress_indices( soap_turbo_hypers(n_soap_turbo)%compress_mode, &
                                         soap_turbo_hypers(n_soap_turbo)%alpha_max, &
                                         soap_turbo_hypers(n_soap_turbo)%l_max, &
                                         soap_turbo_hypers(n_soap_turbo)%dim, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_i, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_j, &
                                         soap_turbo_hypers(n_soap_turbo)%compress_P_el, &
                                         "set_indices" )
            else
              write(*,*) "ERROR: you're trying to use compression but neither a file_compress_soap nor", &
                         "compress_mode are defined!"
              stop
            end if
         else
            soap_turbo_hypers(n_soap_turbo)%dim = soap_turbo_hypers(n_soap_turbo)%n_max * &
                                                  (soap_turbo_hypers(n_soap_turbo)%n_max+1)/2 * &
                                                  (soap_turbo_hypers(n_soap_turbo)%l_max+1)
         end if
!       distance_2b definitions here
        else if( keyword == "distance_2b" )then
          n_distance_2b = n_distance_2b + 1
          iostatus = 0
          do while( keyword /= "gap_end" .and. iostatus == 0 )
            read(10, *, iostat=iostatus) keyword
            if( keyword == "delta" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%delta
            else if( keyword == "sigma" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%sigma
            else if( keyword == "rcut" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%rcut
            else if( keyword == "Z1"  .or. keyword == "z1" .or. keyword == "species1" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%species1
            else if( keyword == "Z2"  .or. keyword == "z2" .or. keyword == "species2" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%species2
            else if( keyword == "desc_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%file_desc
            else if( keyword == "alphas_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%file_alphas
            end if
          end do
!         Read the sparse set information
          call read_alphas_and_descriptors(distance_2b_hypers(n_distance_2b)%file_desc, &
                                           distance_2b_hypers(n_distance_2b)%file_alphas, &
                                           distance_2b_hypers(n_distance_2b)%n_sparse, &
                                           "distance_2b", distance_2b_hypers(n_distance_2b)%alphas, &
                                           distance_2b_hypers(n_distance_2b)%Qs, &
                                           distance_2b_hypers(n_distance_2b)%cutoff)
          if( distance_2b_hypers(n_distance_2b)%rcut > rcut_max )then
            rcut_max = distance_2b_hypers(n_distance_2b)%rcut
          end if
!       angle_3b definitions here
        else if( keyword == "angle_3b" )then
          n_angle_3b = n_angle_3b + 1
          iostatus = 0
          do while( keyword /= "gap_end" .and. iostatus == 0 )
            read(10, *, iostat=iostatus) keyword
            if( keyword == "delta" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%delta
            else if( keyword == "sigma" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%sigma(1:3)
            else if( keyword == "rcut" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%rcut
            else if( keyword == "Z1"  .or. keyword == "z1" .or. keyword == "species1" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species1
            else if( keyword == "Z2"  .or. keyword == "z2" .or. keyword == "species2" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species2
            else if( keyword == "Z_center"  .or. keyword == "z_center" .or. keyword == "species_center" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species_center
            else if( keyword == "kernel_type" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%kernel_type
            else if( keyword == "desc_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%file_desc
            else if( keyword == "alphas_sparse" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%file_alphas
            end if
          end do
!         Read the sparse set information
          call read_alphas_and_descriptors(angle_3b_hypers(n_angle_3b)%file_desc, &
                                           angle_3b_hypers(n_angle_3b)%file_alphas, &
                                           angle_3b_hypers(n_angle_3b)%n_sparse, &
                                           "angle_3b", angle_3b_hypers(n_angle_3b)%alphas, &
                                           angle_3b_hypers(n_angle_3b)%Qs, &
                                           angle_3b_hypers(n_angle_3b)%cutoff)
          if( angle_3b_hypers(n_angle_3b)%rcut > rcut_max )then
            rcut_max = angle_3b_hypers(n_angle_3b)%rcut
          end if
!       core_pot definitions here
        else if( keyword == "core_pot" )then
          n_core_pot = n_core_pot + 1
          iostatus = 0
          do while( keyword /= "gap_end" .and. iostatus == 0 )
            read(10, *, iostat=iostatus) keyword
!            if( keyword == "n" .or. keyword == "N" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%n
!            else if( keyword == "yp1" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%yp1
!            else if( keyword == "ypn" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%ypn
!            else if( keyword == "Z1"  .or. keyword == "z1" .or. keyword == "species1" )then
            if( keyword == "Z1"  .or. keyword == "z1" .or. keyword == "species1" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%species1
            else if( keyword == "Z2"  .or. keyword == "z2" .or. keyword == "species2" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%species2
            else if( keyword == "core_pot_file" )then
              backspace(10)
              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%core_pot_file
            end if
          end do
!         Allocate some arrays, read in potential, etc.
          open(20, file=core_pot_hypers(n_core_pot)%core_pot_file, status="unknown")
          read(20, *) core_pot_hypers(n_core_pot)%n, core_pot_hypers(n_core_pot)%yp1, core_pot_hypers(n_core_pot)%ypn
          n = core_pot_hypers(n_core_pot)%n

          allocate( V(1:n) )
          allocate( x(1:n) )
          counter = 0
          do i = 1, n
            read(20,*) x(i), V(i)
            if( x(i) <= params%core_pot_cutoff )then
              counter = counter + 1
              x(counter) = x(i)
              if( x(i) <= params%core_pot_cutoff - params%core_pot_buffer )then
                V(counter) = V(i)
              else
                V(counter) = V(i) * 0.5d0 * ( dcos(dacos(-1.d0)/params%core_pot_buffer*(x(i) - params%core_pot_cutoff &
                             + params%core_pot_buffer)) + 1.d0 )
              end if
            end if
          end do
          close(20)
          n = counter
          core_pot_hypers(n_core_pot)%n = n
          allocate( core_pot_hypers(n_core_pot)%V(1:n) )
          core_pot_hypers(n_core_pot)%V(1:n) = V(1:n)
          deallocate( V )
          allocate( core_pot_hypers(n_core_pot)%x(1:n) )
          core_pot_hypers(n_core_pot)%x(1:n) = x(1:n)
          deallocate( x )
          allocate( core_pot_hypers(n_core_pot)%dVdx2(1:n) )
!         This code below for spline second derivative is more or less copy-pasted from QUIP. It's the
!         easiest way to make sure both interpolations give the same numbers
          allocate( u(1:n) )
          if( core_pot_hypers(n_core_pot)%yp1 > 0.99d30 )then
            core_pot_hypers(n_core_pot)%dVdx2(1) = 0.d0
            u(1) = 0.d0
          else
            core_pot_hypers(n_core_pot)%dVdx2(1) = -0.5d0

            u(1) =  (3.d0 / (core_pot_hypers(n_core_pot)%x(2) - core_pot_hypers(n_core_pot)%x(1) ) ) * &
                    ( (core_pot_hypers(n_core_pot)%V(2) - core_pot_hypers(n_core_pot)%V(1)) / &
                      (core_pot_hypers(n_core_pot)%x(2) - core_pot_hypers(n_core_pot)%x(1)) - core_pot_hypers(n_core_pot)%yp1)
          end if
          do i = 2, n-1
            sig = ( core_pot_hypers(n_core_pot)%x(i) - core_pot_hypers(n_core_pot)%x(i-1) ) / &
                  ( core_pot_hypers(n_core_pot)%x(i+1) - core_pot_hypers(n_core_pot)%x(i-1) )
            p = sig * core_pot_hypers(n_core_pot)%dVdx2(i-1) + 2.d0
            core_pot_hypers(n_core_pot)%dVdx2(i) = (sig-1.d0)/p
            u(i) = ( core_pot_hypers(n_core_pot)%V(i+1) - core_pot_hypers(n_core_pot)%V(i) ) / &
                   ( core_pot_hypers(n_core_pot)%x(i+1) - core_pot_hypers(n_core_pot)%x(i)) - &
                   ( core_pot_hypers(n_core_pot)%V(i) - core_pot_hypers(n_core_pot)%V(i-1) ) / &
                   ( core_pot_hypers(n_core_pot)%x(i) - core_pot_hypers(n_core_pot)%x(i-1))
            u(i) = ( 6.d0 * u(i) / ( core_pot_hypers(n_core_pot)%x(i+1) - core_pot_hypers(n_core_pot)%x(i-1) ) &
                     - sig * u(i-1) ) / p
          end do
          if( core_pot_hypers(n_core_pot)%ypn > 0.99d30 )then
            qn = 0.d0
            un = 0.d0
          else
            qn = 0.5d0
            un = ( 3.d0 / (core_pot_hypers(n_core_pot)%x(n) - core_pot_hypers(n_core_pot)%x(n-1)) ) * &
                 ( core_pot_hypers(n_core_pot)%ypn - (core_pot_hypers(n_core_pot)%V(n) - core_pot_hypers(n_core_pot)%V(n-1)) &
                 / (core_pot_hypers(n_core_pot)%x(n) - core_pot_hypers(n_core_pot)%x(n-1)) )
          end if
          core_pot_hypers(n_core_pot)%dVdx2(n) = ( un - qn*u(n-1) ) / ( qn*core_pot_hypers(n_core_pot)%dVdx2(n-1) + 1.d0 )
          do i = n-1, 1, -1
            core_pot_hypers(n_core_pot)%dVdx2(i) = core_pot_hypers(n_core_pot)%dVdx2(i) * &
                                                  core_pot_hypers(n_core_pot)%dVdx2(i+1) + u(i)
          end do
          if( core_pot_hypers(n_core_pot)%yp1 > 0.99d30 )then
            u(1:1) = spline_der(core_pot_hypers(n_core_pot)%x, core_pot_hypers(n_core_pot)%V, &
                                core_pot_hypers(n_core_pot)%dVdx2, core_pot_hypers(n_core_pot)%yp1, &
                                core_pot_hypers(n_core_pot)%ypn, core_pot_hypers(n_core_pot)%x(1:1), 1.d100)
            core_pot_hypers(n_core_pot)%yp1 = u(1)
          end if
          if( core_pot_hypers(n_core_pot)%ypn > 0.99d30 )then
            u(1:1) = spline_der(core_pot_hypers(n_core_pot)%x, core_pot_hypers(n_core_pot)%V, &
                                core_pot_hypers(n_core_pot)%dVdx2, core_pot_hypers(n_core_pot)%yp1, &
                                core_pot_hypers(n_core_pot)%ypn, core_pot_hypers(n_core_pot)%x(n:n), 1.d100)
            core_pot_hypers(n_core_pot)%ypn = u(1)
          end if
          deallocate( u )
          if( maxval(core_pot_hypers(n_core_pot)%x) > rcut_max )then
            rcut_max = maxval(core_pot_hypers(n_core_pot)%x)
          end if
        end if
      end if
    end do
    close(10)

  end subroutine
!**************************************************************************




!**************************************************************************
  subroutine upper_to_lower_case(string)

    implicit none

    character(len=*), intent(inout) :: string
    character*1 :: upper_case_dict(1:26), lower_case_dict(1:26)
    integer :: i, j, n

    upper_case_dict = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", &
                       "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", &
                       "U", "V", "W", "X", "Y", "Z" ]
    lower_case_dict = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", &
                       "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", &
                       "u", "v", "w", "x", "y", "z" ]


    do i = 1, len(string)
      do j = 1, size(upper_case_dict)
        if( string(i:i) == upper_case_dict(j) )then
          string(i:i) = lower_case_dict(j)
        end if
      end do
    end do

  end subroutine
!**************************************************************************




end module
