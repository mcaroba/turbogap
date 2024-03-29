! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, groups.f90, is copyright (c) 2023, Uttiyoarnab Saha
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

!**************************************************************************
! This program will make the different groups of atoms as required in the 
! model simulation using input data.
! Note:
! This program allows for repeated use of the same input keyword 'make_group'. 
! From the input information, groups are first registered. Then the registered 
! data are used when the main program calls to make the groups of atoms specific
! to the physical processes. After that, if dynamic update of group(s) is not
! invoked, the arrays with registered data are freed from memory.
! -------  
! 				
!********* by Uttiyoarnab Saha

!**************************************************************************

module groups

type groups_class
	
	integer :: group_num_atoms
	integer :: count_groups = 0, max_num_groups = 17 
	!! made group IDs and group types
	character*16, allocatable :: group_IDs(:), group_styles(:)
	!! sizes of each group operation
	integer, allocatable :: group_block_values(:), group_add_values(:), &
		group_subtract_values(:), group_ID_values(:), group_sphere_values(:), &
		group_atomtypes_values(:)
	!! parameters of each group operation
	real*8, allocatable :: group_block_limits(:), group_sphere_limits(:)
	character*16, allocatable :: group_add_IDs(:), group_subtract_IDs(:)
	character*8, allocatable :: group_atom_types(:)
	integer, allocatable :: group_atom_IDs(:)
	!! atom-group map
	integer, allocatable :: group_atoms_tallies(:,:)
	!! for dynamic groups 
	integer, allocatable :: dynamic_group_tallies(:), dynamic_update_steps(:)
	logical :: presence_of_dynamic_groups = .false.

	logical :: recorded_group_for_all_atoms = .false.

	!! to be obtained on call from main
	integer :: n_sites

	contains
	procedure :: defaultInitialGroups, recordGroups, checkValidGroupStyle
	procedure :: checkValidGroupID, makeGroups, getAtomsInGroups, freeGroupsRegister
	procedure :: updateDynamicGroups

end type groups_class

contains


!! Initialization and default case
subroutine defaultInitialGroups (this)

	implicit none
	class (groups_class) :: this

	if (.not. allocated(this%group_IDs)) allocate(this%group_IDs(0))
	if (.not. allocated(this%group_styles)) allocate(this%group_styles(0))
	if (.not. allocated(this%dynamic_group_tallies)) allocate(this%dynamic_group_tallies(0))
	if (.not. allocated(this%dynamic_update_steps)) allocate(this%dynamic_update_steps(0))

	!! register group for 'all' atoms with its type and group id both as 'all'

	if ( .not. this%recorded_group_for_all_atoms) then
		this%group_IDs = [this%group_IDs, 'all']
		this%group_styles = [this%group_styles, 'all']
		this%count_groups = size(this%group_IDs)
		this%recorded_group_for_all_atoms = .true.
	end if

end subroutine defaultInitialGroups


!! To check if invoked type of group is in the implemented list
subroutine checkValidGroupStyle (this, rank, ierr, used_group_style)

	implicit none
	class (groups_class) :: this
	character*16, intent(in) :: used_group_style
	integer, intent(in) :: rank
	integer :: ierr
	character*16, allocatable :: implemented_group_styles(:)
	logical :: groupstyle_validity = .true.

	allocate(implemented_group_styles(0))
	implemented_group_styles = [implemented_group_styles, &
	'all', 'block', 'add', 'subtract', 'id', 'sphere', 'atomtype', 'dynamic']

	if ( findloc(implemented_group_styles(2:size(implemented_group_styles)), &
			used_group_style, dim = 1) == 0 ) groupstyle_validity = .false.
	if (.not. groupstyle_validity) then
		if (rank == 0) then
			write(*,*) 'ERROR: Group type ', trim(used_group_style), ' is not possible.'
			write(*,*) 'See list of possible group styles.'
		end if
		call mpi_finalize(ierr)
		stop
	end if

end subroutine checkValidGroupStyle


!! To checck if given group name is valid for use in any calculation
subroutine checkValidGroupID (this, rank, ierr, used_groupID)

	implicit none
	class (groups_class) :: this
	character*16, intent(in) :: used_groupID
	integer, intent(in) :: rank
	integer :: ierr
	logical :: groupID_validity = .true.

	if ( findloc(this%group_IDs, used_groupID, dim = 1) == 0 ) groupID_validity = .false.
	if (.not. groupID_validity) then
		if (rank == 0) write(*,*) 'ERROR: checkValidGroupID: Group ID ', trim(used_groupID), ' is not known.'
		call mpi_finalize(ierr)
		stop
	end if 

end subroutine checkValidGroupID


!! Records all the groups made when called after reading input file
subroutine recordGroups (this, rank, ierr, group_ID, group_style, group_style_value, &
	group_style_IDs, group_atom_IDs, group_style_limits, group_atom_types)

	implicit none
	
	class (groups_class) :: this
	character*16, intent(in) :: group_ID, group_style, group_style_IDs(:)
	character*8, intent(in) :: group_atom_types(:)
	real*8, intent(in) :: group_style_limits(:)
	integer, intent(in) :: rank, group_style_value, group_atom_IDs(:)
	integer :: ierr
	logical :: is_group_dynamic = .false.
	integer :: i, group_label

	!! the following arrays will contain the corresponding sizes and dimensions 
	!! respective to the group types of different
	!! group ids that have been made at their own definite positions 
	
	if (group_style == 'block') then
		if (.not. allocated(this%group_block_values)) allocate(this%group_block_values(0))
		if (.not. allocated(this%group_block_limits)) allocate(this%group_block_limits(0))

	else if (group_style == 'add') then
		if (.not. allocated(this%group_add_values)) allocate(this%group_add_values(0))
		if (.not. allocated(this%group_add_IDs)) allocate(this%group_add_IDs(0))

	else if (group_style == 'subtract') then
		if (.not. allocated(this%group_subtract_values)) allocate(this%group_subtract_values(0))
		if (.not. allocated(this%group_subtract_IDs)) allocate(this%group_subtract_IDs(0))

	else if (group_style == 'id') then
		if (.not. allocated(this%group_ID_values)) allocate(this%group_ID_values(0))
		if (.not. allocated(this%group_atom_IDs)) allocate(this%group_atom_IDs(0))

	else if (group_style == 'sphere') then
		if (.not. allocated(this%group_sphere_values)) allocate(this%group_sphere_values(0))
		if (.not. allocated(this%group_sphere_limits)) allocate(this%group_sphere_limits(0))

	else if (group_style == 'atomtype') then
		if (.not. allocated(this%group_atomtypes_values)) allocate(this%group_atomtypes_values(0))
		if (.not. allocated(this%group_atom_types)) allocate(this%group_atom_types(0))

	else if (group_style == 'dynamic') then
		if (.not. this%presence_of_dynamic_groups) this%presence_of_dynamic_groups = .true.
		is_group_dynamic = .true.
	end if

	!! check if the group type is not dynamic, consider adding its record
	if (.not. is_group_dynamic) then
		!! check if maximum number of groups are made
		!! check if group ID already exists
		!! if not , then register group with ID, type and also count it in.
		if ( this%count_groups <= this%max_num_groups ) then
			group_label = findloc(this%group_IDs, group_ID, dim = 1)
			if ((group_label == 0) .and. group_ID /= '') then
				this%group_IDs = [this%group_IDs, group_ID]
				this%group_styles = [this%group_styles, group_style]
				this%count_groups = size(this%group_IDs)
			else
				if (rank == 0)  write(*,*) 'ERROR: Group ID cannot be repeated or blank.'
				call mpi_finalize(ierr)
				stop
			end if
		else
			if (rank == 0) write(*,*) 'ERROR: Too many groups are already made.'
			call mpi_finalize(ierr)
			stop
		end if

		group_label = findloc(this%group_IDs, group_ID, dim = 1)
	
		if (group_style == 'block') then
			do i = 1, (group_label - size(this%group_block_values) - 1)
				this%group_block_values = [this%group_block_values, 0]
			end do
			this%group_block_values = [this%group_block_values, group_style_value]
			this%group_block_limits = [this%group_block_limits, group_style_limits]
	
		else if (group_style == 'add') then
			do i = 1, (group_label - size(this%group_add_values) - 1)
				this%group_add_values = [this%group_add_values, 0]
			end do
			this%group_add_values = [this%group_add_values, group_style_value]
			do i = 1, group_style_value
				call checkValidGroupID (this, rank, ierr, group_style_IDs(i))
			end do
			this%group_add_IDs = [this%group_add_IDs, group_style_IDs]
	
		else if (group_style == 'subtract') then
			do i = 1, (group_label - size(this%group_subtract_values) - 1)
				this%group_subtract_values = [this%group_subtract_values, 0]
			end do
			this%group_subtract_values = [this%group_subtract_values, group_style_value]
			do i = 1, group_style_value
				call checkValidGroupID (this, rank, ierr, group_style_IDs(i))
			end do
			this%group_subtract_IDs = [this%group_subtract_IDs, group_style_IDs]
	
		else if (group_style == 'id') then
			do i = 1, (group_label - size(this%group_ID_values) - 1)
				this%group_ID_values = [this%group_ID_values, 0]
			end do
			this%group_ID_values = [this%group_ID_values, group_style_value]
			this%group_atom_IDs = [this%group_atom_IDs, group_atom_IDs]
	
		else if (group_style == 'sphere') then
			do i = 1, (group_label - size(this%group_sphere_values) - 1)
				this%group_sphere_values = [this%group_sphere_values, 0]
			end do
			this%group_sphere_values = [this%group_sphere_values, group_style_value]
			this%group_sphere_limits = [this%group_sphere_limits, group_style_limits]

		else if (group_style == 'atomtype') then
			do i = 1, (group_label - size(this%group_atomtypes_values) - 1)
				this%group_atomtypes_values = [this%group_atomtypes_values, 0]
			end do 
			this%group_atomtypes_values = [this%group_atomtypes_values, group_style_value]
			this%group_atom_types = [this%group_atom_types, group_atom_types]

		end if
	else
		!! consider recording the given group as dynamic
		group_label = findloc(this%group_IDs, group_ID, dim = 1)
		this%dynamic_group_tallies = [this%dynamic_group_tallies, group_label]
		this%dynamic_update_steps = [this%dynamic_update_steps, group_style_value]
	end if 

end subroutine recordGroups


!! This subroutine has to be called from turbogap.f90 only once before starting MD.
!! It uses the registered data of its class and form all the groups at one call.

!! Makes all the groups at once when called from main turbogap.f90 
subroutine makeGroups(this, rank, ierr, n_sites, positions, masses, species_types, species_mass)

	implicit none
	class (groups_class) :: this
	integer, intent(in) :: rank, n_sites
	real*8, intent(in) :: positions(:,:), masses(:), species_mass(:)
	character*8, intent(in) :: species_types(:)
	integer :: i, lower_index, higher_index
	integer :: ierr

	this%n_sites = n_sites

	if (.not. allocated(this%group_atoms_tallies)) then
		allocate(this%group_atoms_tallies(this%n_sites, this%count_groups))
		this%group_atoms_tallies = 0
	end if

	!! default group for all atoms
	call makeGroupAll(this)

	!! for other groups
	do i = 2, this%count_groups
		if ( this%group_styles(i) == 'block' ) then
			lower_index = sum( this%group_block_values(1:i-1) ) + 1
			higher_index = sum( this%group_block_values(1:i) )
			call makeGroupBlock ( this, rank, ierr, positions, i, this%group_block_values(i), &
								this%group_block_limits(lower_index:higher_index) )
		else if ( this%group_styles(i) == 'add' ) then
			lower_index = sum( this%group_add_values(1:i-1) ) + 1
			higher_index = sum( this%group_add_values(1:i) )
			call makeGroupAdd ( this, i, this%group_add_values(i), this%group_add_IDs(lower_index:higher_index) )
		else if ( this%group_styles(i) == 'subtract' ) then
			lower_index = sum( this%group_subtract_values(1:i-1) ) + 1
			higher_index = sum( this%group_subtract_values(1:i) )
			call makeGroupSubtract ( this, i, this%group_subtract_values(i), this%group_subtract_IDs(lower_index:higher_index) )
		else if ( this%group_styles(i) == 'id' ) then
			lower_index = sum( this%group_ID_values(1:i-1) ) + 1
			higher_index = sum( this%group_ID_values(1:i) )
			call makeGroupID ( this, rank, ierr, i, this%group_ID_values(i), this%group_atom_IDs(lower_index:higher_index) )
		else if ( this%group_styles(i) == 'sphere' ) then
			lower_index = sum( this%group_sphere_values(1:i-1) ) + 1
			higher_index = lower_index + 3
			call makeGroupSphere ( this, positions, i, this%group_sphere_limits(lower_index:higher_index) )
		else if ( this%group_styles(i) == 'atomtype' ) then
			lower_index = sum( this%group_atomtypes_values(1:i-1) ) + 1
			higher_index = sum( this%group_atomtypes_values(1:i) )
			call makeGroupAtomTypes ( this, i, this%group_atom_types(lower_index:higher_index), &
					masses, species_types, species_mass)
		end if
	end do

	if (this%count_groups > 1) call printNumberOfAtomsInGroups(this, rank)

end subroutine makeGroups


!! Update dynamic groups
subroutine updateDynamicGroups (this, rank, ierr, at_index, positions)

	implicit none
	class (groups_class) :: this
	integer, intent(in) :: rank, at_index
	real*8, intent(in) :: positions(:,:)
	integer :: ierr
	integer :: group_label, lower_index, higher_index

	group_label = this%dynamic_group_tallies(at_index)

	!! first uncheck dynamic group label
	this%group_atoms_tallies(:,group_label) = 0

	!! update corresponding dynamic group
	if ( this%group_styles(group_label) == 'block' ) then
		lower_index = sum( this%group_block_values(1:group_label-1) ) + 1
		higher_index = sum( this%group_block_values(1:group_label) )
		call makeGroupBlock ( this, rank, ierr, positions, group_label, this%group_block_values(group_label), &
				this%group_block_limits(lower_index:higher_index) )

	else if ( this%group_styles(group_label) == 'add' ) then
		lower_index = sum( this%group_add_values(1:group_label-1) ) + 1
		higher_index = sum( this%group_add_values(1:group_label) )
		call makeGroupAdd ( this, group_label, this%group_add_values(group_label), &
				this%group_add_IDs(lower_index:higher_index) )

	else if ( this%group_styles(group_label) == 'subtract' ) then
		lower_index = sum( this%group_subtract_values(1:group_label-1) ) + 1
		higher_index = sum( this%group_subtract_values(1:group_label) )
		call makeGroupSubtract ( this, group_label, this%group_subtract_values(group_label), &
				this%group_subtract_IDs(lower_index:higher_index) )

	else if ( this%group_styles(group_label) == 'sphere' ) then
		lower_index = sum( this%group_sphere_values(1:group_label-1) ) + 1
		higher_index = lower_index + 3
		call makeGroupSphere ( this, positions, group_label, this%group_sphere_limits(lower_index:higher_index) )
	
	!! atom types and ids do not change dynamically, so no need to update groups made based on them
	end if

	if ( rank == 0 ) then
		write(*,*)sum(this%group_atoms_tallies(:,group_label)), ' atoms in dynamic group ',this%group_IDs(group_label)
	end if

end subroutine updateDynamicGroups


!! Routines specific to types of groups
!! Makes a group of all atoms
subroutine makeGroupAll (this)

	implicit none
	class (groups_class) :: this
	
	this%group_atoms_tallies(:,1) = 1

end subroutine makeGroupAll


!! Makes block type group of atoms
subroutine makeGroupBlock (this, rank, ierr, positions, group_label, num_values, group_dimensions)

	implicit none
	class (groups_class) :: this
	real*8, intent(in) :: positions(:,:), group_dimensions(:)
	integer, intent(in) :: rank, group_label, num_values
	real*8 :: group_xlow, group_xhigh, group_ylow, group_yhigh, group_zlow, group_zhigh
	integer :: i, j, ierr

	do i = 1, num_values, 6
		group_xlow = group_dimensions(i)
		group_xhigh = group_dimensions(i+1)
		group_ylow = group_dimensions(i+2)
		group_yhigh = group_dimensions(i+3)
		group_zlow = group_dimensions(i+4)
		group_zhigh = group_dimensions(i+5)

		if ( group_xlow > group_xhigh .or. &
			group_ylow > group_yhigh .or. group_zlow > group_zhigh) then
			if ( rank == 0 ) write(*,*) 'ERROR: Lower limit is larger than higher limit in group dimensions.'
			call mpi_finalize(ierr)
			stop
		end if
		
		do j = 1, this%n_sites
			if ((group_xlow <= positions(1,j) .and. positions(1,j) <= group_xhigh) &
			.and. (group_ylow <= positions(2,j) .and. positions(2,j) <= group_yhigh) &
			.and. (group_zlow <= positions(3,j) .and. positions(3,j) <= group_zhigh)) then

				this%group_atoms_tallies(j, group_label) = 1
			end if
		end do
	end do

end subroutine makeGroupBlock


!! Makes group by adding groups of atoms
subroutine makeGroupAdd (this, group_label, num_values, ids_to_add)
	
	implicit none
	class (groups_class) :: this
	character*16, intent(in) :: ids_to_add(:)
	integer, intent(in) :: group_label, num_values
	integer :: i, j, dynamic_labels
	
	do i = 1, num_values
		dynamic_labels = findloc( this%group_IDs, ids_to_add(i), dim = 1 )
		do j = 1, this%n_sites
			if ( this%group_atoms_tallies(j, dynamic_labels) == 1 ) then
				this%group_atoms_tallies(j, group_label) = 1
			end if
		end do
	end do

end subroutine makeGroupAdd


!! Makes group by subtracting groups of atoms from a group
subroutine makeGroupSubtract (this, group_label, num_values, ids_to_subtract)
	
	implicit none
	class (groups_class) :: this
	character*16, intent(in) :: ids_to_subtract(:)
	integer, intent(in) :: group_label, num_values
	integer :: i, j, fixed_label, dynamic_labels, flag

	fixed_label = findloc( this%group_IDs, ids_to_subtract(1), dim = 1 )
	do j = 1, this%n_sites
		if ( this%group_atoms_tallies(j, fixed_label) == 1 ) then
			flag = 0
			do i = 2, num_values
				dynamic_labels = findloc( this%group_IDs, ids_to_subtract(i), dim = 1 )
				if ( this%group_atoms_tallies(j, dynamic_labels) == 1 ) then
					flag = 1
					exit
				end if
			end do
			if (flag == 0) this%group_atoms_tallies(j, group_label) = 1
		end if
	end do

end subroutine makeGroupSubtract


!! Makes group by the ids of atoms
subroutine makeGroupID (this, rank, ierr, group_label, num_values, ids_of_atoms)
	
	implicit none
	class (groups_class) :: this
	integer, intent(in) :: rank, ids_of_atoms(:), group_label, num_values
	integer :: i, ierr
	
	if (minval(ids_of_atoms) < 1 .or. maxval(ids_of_atoms) > this%n_sites) then
		if ( rank == 0 ) write(*,*) 'ERROR: Atom id not found for grouping.'
		call mpi_finalize(ierr)
		stop
	end if

	do i = 1, num_values
		this%group_atoms_tallies(ids_of_atoms(i), group_label) = 1
	end do

end subroutine makeGroupID


!! Makes sphere type group of atoms
subroutine makeGroupSphere (this, positions, group_label, group_dimensions)

	implicit none
	class (groups_class) :: this
	real*8, intent(in) :: positions(:,:), group_dimensions(:)
	integer, intent(in) :: group_label
	real*8 :: centre_x, centre_y, centre_z, sphere_radius, dist
	integer :: j

	sphere_radius = group_dimensions(1)
	centre_x = group_dimensions(2)
	centre_y = group_dimensions(3)
	centre_z = group_dimensions(4)

	do j = 1, this%n_sites
		dist = get_distance ( centre_x, centre_y, centre_z, &
					positions(1,j), positions(2,j), positions(3,j) )
		if ( dist <= sphere_radius ) this%group_atoms_tallies(j, group_label) = 1
	end do

end subroutine makeGroupSphere


!! Makes group atoms by specified types
subroutine makeGroupAtomTypes ( this, group_label, species_in_group, &
			masses, species_types, species_mass)
	
	implicit none
	class (groups_class) :: this
	integer, intent(in) :: group_label
	character*8, intent(in) :: species_in_group(:), species_types(:)
	real*8, intent(in) :: masses(:), species_mass(:)
	real*8, allocatable :: group_species_mass(:)
	integer :: i, j

	allocate(group_species_mass(size(species_in_group)))
	do i = 1, size(species_in_group)
		group_species_mass(i) = species_mass(findloc(species_types, species_in_group(i), dim = 1))
	end do
	do i = 1, this%n_sites
		if ( findloc(group_species_mass, masses(i), dim = 1) /= 0 ) &
				this%group_atoms_tallies(i, group_label) = 1
	end do
	
end subroutine makeGroupAtomTypes


!! Return the atom ids in specific groups
subroutine getAtomsInGroups (this, used_groupID, is_group_dynamic, update_interval, &
		md_step, atoms_in_group)

	implicit none
	class (groups_class) :: this
	character*16, intent(in) :: used_groupID
	integer, allocatable, intent(out) :: atoms_in_group(:)
	integer, intent(out) :: update_interval
	integer, intent(in) :: md_step
	logical, intent(inout) :: is_group_dynamic
	integer :: i, group_label, at_index

	group_label = findloc( this%group_IDs, used_groupID, dim = 1 )
	allocate(atoms_in_group(0))

	do i = 1, this%n_sites
		if ( this%group_atoms_tallies(i,group_label) == 1 ) atoms_in_group = [atoms_in_group, i]
	end do

	if ( md_step == 0 ) then
		at_index = findloc(this%dynamic_group_tallies, group_label, dim = 1)
		if ( at_index /= 0 ) then
			is_group_dynamic = .true.
			update_interval = this%dynamic_update_steps(at_index)
		end if	
	end if

end subroutine getAtomsInGroups


!! After the call of getAtomsInGroups before starting MD, this subroutine 
!! is called from main turbogap.f90 to deallocate all/one (depending on static/dynamic) 
!! of the arrays registered with groups because atoms specific to the processes
!! have already been obtained for doing the MD.
!! Mostly simulations are with static groups only.
subroutine freeGroupsRegister (this)
	
	implicit none
	class (groups_class) :: this
	
	if (.not. this%presence_of_dynamic_groups) then
		!! none of these arrays will be used later
		if ( allocated(this%group_IDs) ) deallocate (this%group_IDs)
		if ( allocated(this%group_styles) ) deallocate (this%group_styles)
		if ( allocated(this%group_block_values) ) deallocate (this%group_block_values)
		if ( allocated(this%group_add_values) ) deallocate (this%group_add_values)
		if ( allocated(this%group_subtract_values) ) deallocate (this%group_subtract_values)
		if ( allocated(this%group_ID_values) ) deallocate (this%group_ID_values)
		if ( allocated(this%group_block_limits) ) deallocate (this%group_block_limits)
		if ( allocated(this%group_add_IDs) ) deallocate (this%group_add_IDs)
		if ( allocated(this%group_subtract_IDs) ) deallocate(this%group_subtract_IDs)
		if ( allocated(this%group_atom_IDs) ) deallocate(this%group_atom_IDs)
		if ( allocated(this%group_atoms_tallies) ) deallocate(this%group_atoms_tallies)
	end if

end subroutine freeGroupsRegister


!! count and print out number of atoms in groups after (static) groups are made
subroutine printNumberOfAtomsInGroups (this, rank)

	implicit none
	class (groups_class) :: this
	integer, intent(in) :: rank
	integer :: i

	write(*,*)
	do i = 2, size(this%group_IDs)
		if ( rank == 0 ) write(*,*) sum(this%group_atoms_tallies(:,i)),' atoms in group ',this%group_IDs(i) 
	end do

end subroutine printNumberOfAtomsInGroups


!! Find distance between two atoms.
real function get_distance (xi, yi, zi, xj, yj, zj) result (r_ij)
	implicit none
	real*8, intent(in) :: xi, yi, zi, xj, yj, zj
	r_ij = sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
end function get_distance

end module groups
