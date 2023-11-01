> AA_pol.dat
for i in $(seq 1.0 0.1 6.0); do
echo $i
cat <<EOF > input
atoms_file = 'AA_$i.xyz'
pot_file = 'gap_files/carbon.gap'

! Species info
n_species = 1 !2
species = C !H
e0 = -.16138053
masses = 12.01 !1.00784

! van der Waals info
vdw_type = mbd
vdw_sr = 0.83
vdw_d = 6.            ! Use d = 20 for TS(SCS) and d = 6 for MBD
vdw_rcut = 20.
vdw_r0_ref = 1.900 !1.64
vdw_alpha0_ref = 1.778 !0.667
vdw_c6_ref = 27.8 !3.88
! The order of the cutoffs is vdw_scs_rcut < vdw_mbd_rcut < vdw_2b_rcut
! If two cutoffs are equal, there is no buffer region between them!
vdw_buffer = 0.5          ! Buffer for transitions between cut-off regions. Type: REAL
vdw_scs_rcut = 6.        ! vdw_scs_rcut > vdw_buffer. Type: REAL
vdw_loc_rcut = 0.
vdw_mbd_rcut = 5.        ! Cut-off for atoms to include for local MBD energy (vdw_mbd_rcut >= vdw_scs_rcut + vdw_buffer). Type: REAL
vdw_mbd_rcut2 = 5.
vdw_2b_rcut = 0.          ! Cut-off for local TS-SCS (vdw_2b_rcut >= vdw_mbd_rcut + vdw_buffer), Type: REAL
vdw_mbd_nfreq = 11        ! Number of frequency values for MBD integration. Type: INT
vdw_mbd_norder = 6        ! Contributions up to n-body interactions (i.e. cut-off degree for Taylor expansion of ln(I-AT)). Type: INT
vdw_mbd_grad = .false.     ! Calculate MBD forces. Type: LOGICAL
vdw_hirsh_grad = .false.   ! Include Hirshfeld gradients in the forces. Type: LOGICAL
vdw_polynomial = .false.  ! Use polynomial approximation for inverse matrices. Type: LOGICAL
vdw_omega_ref = 1.3d0
do_nnls = .false.
EOF
../bin/turbogap predict | grep "a_iso" | awk '{print $2 $3}' >> AA_pol.dat
done
