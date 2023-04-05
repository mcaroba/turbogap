> force_timing.dat
for d in $(seq 2 1 6); do
echo $d
> c60_forces_$d.dat
cat <<EOF > input
atoms_file = 'c60.xyz'
pot_file = 'gap_files/carbon.gap'

! Species info
n_species = 1 !2
species = C !H
!e0 = -.16138053
masses = 12.01 1.00784

! van der Waals info
vdw_type = mbd
vdw_sr = 0.83 !0.97
vdw_d = 6. !20.            ! Use d = 20 for TS(SCS) and d = 6 for MBD
vdw_rcut = 20.
vdw_r0_ref = 1.900 !1.64
vdw_alpha0_ref = 1.778 !0.667
vdw_c6_ref = 27.8 !3.88
! The order of the cutoffs is vdw_scs_rcut < vdw_mbd_rcut < vdw_2b_rcut
! If two cutoffs are equal, there is no buffer region between them!
vdw_buffer = 2.          ! Buffer for transitions between cut-off regions. Type: REAL
vdw_loc_rcut = 4.5
vdw_scs_rcut = 4.5        ! vdw_scs_rcut > vdw_buffer. Type: REAL
vdw_mbd_rcut = 7.        ! Cut-off for atoms to include for local MBD energy (vdw_mbd_rcut >= vdw_scs_rcut + vdw_buffer). Type: REAL
vdw_2b_rcut = 20.          ! Cut-off for local TS-SCS (vdw_2b_rcut >= vdw_mbd_rcut + vdw_buffer), Type: REAL
vdw_mbd_nfreq = 12        ! Number of frequency values for MBD integration. Type: INT
vdw_mbd_norder = $d        ! Contributions up to n-body interactions (i.e. cut-off degree for Taylor expansion of ln(I-AT)). Type: INT
vdw_mbd_grad = .true.     ! Calculate MBD forces. Type: LOGICAL
vdw_hirsh_grad = .true.   ! Include Hirshfeld gradients in the forces. Type: LOGICAL
vdw_polynomial = .false.  ! Use polynomial approximation for inverse matrices. Type: LOGICAL
vdw_omega_ref = 4.d0
EOF
> output
mpirun -np 1 ../bin/turbogap predict >> output
grep "MBD force" output | awk '{print $5}' >> c60_forces_$d.dat
grep "MBD timing" output | awk '{print $3}' >> force_timing.dat
done
