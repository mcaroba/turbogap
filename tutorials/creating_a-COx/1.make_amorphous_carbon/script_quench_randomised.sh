
sim_name="md_quench_randomise"

randomised_atoms="atoms_randomise.xyz"

if [ ! -f $randomised_atoms ]; then
    echo "!! Randomised atoms don't exist! run `bash 1.run_randomise.sh` "
    exit 1
fi


input_atoms="atoms.xyz"

cp $randomised_atoms $input_atoms

output_atoms="atoms_quench.xyz"

ln -sf ../gap_files ./

# 1. Run TurboGAP to quench the structure, this time using a barostat
cat <<EOF > input
! Species-specific info
atoms_file = '${input_atoms}'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99

! MD options
md_nsteps = 5500               ! 5.5ps quench
md_step = 1.
thermostat = berendsen
t_beg = 9000                   ! Quenching from 9000K
t_end = 1000                   !             to 1000K

tau_t = 100.
write_thermo = 1
write_xyz = 200
EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap md >  turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"


sleep 1

python3 get_final_config.py $output_atoms

mv thermo.log thermo_2.log
mv trajectory_out.xyz traj_2.xyz
