> cutoff_scs.dat
> pol_scs2.dat

for i in $(seq 1.0 0.2 8.0); do

cat <<EOF > input
atoms_file = 'atoms.xyz'
!atoms_file = 'atoms_p.xyz'
!atoms_file = 'atoms_m.xyz'
pot_file = 'gap_files/carbon.gap'

! Species info
n_species = 1
species = C
e0 = -.16138053
masses = 12.01

! van der Waals info
vdw_type = ts
vdw_sr = 0.83
vdw_d = 6.
vdw_rcut = $i
vdw_buffer = 0.
vdw_r0_ref = 1.900
EOF

echo $i >> cutoff_scs.dat
../bin/turbogap predict | grep "alpha_SCS" | awk '{print $2}' >> pol_scs2.dat

done
