> energies_p_2829.dat

for i in {1..9}; do

for c in {1..3}; do

echo $i $c

cat <<EOF > input
atoms_file = '2829_${i}_p${c}.xyz'
pot_file = 'gap_files/carbon.gap'

! Species info
n_species = 1 !2
species = C !H
!e0 = -.16138053
masses = 12.01 1.00784

! van der Waals info
vdw_type = ts
vdw_sr = 0.83 !0.97
vdw_d = 6. !20.            ! Use d = 20 for TS(SCS) and d = 6 for MBD
vdw_rcut = 9.
vdw_buffer = 0.
vdw_r0_ref = 1.900 !1.64
vdw_alpha0_ref = 1.778 !0.667
vdw_c6_ref = 27.8 !3.88
vdw_scs_rcut = 4.5
vdw_mbd_nfreq = 3
vdw_mbd_grad = .true.
!vdw_mbd_grad = .false.
EOF

../bin/turbogap predict | grep "sphere           $i" | awk '{print $7}' >> energies_p_2829.dat

done

done
