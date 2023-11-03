
sim_name="md_anneal_quench"

quenched_atoms="atoms_quench.xyz"

if [ ! -f $quenched_atoms ]; then
    echo "!! Quenched atoms don't exist! run `bash 2.run_quench.sh` "
    exit 1
fi

input_atoms="atoms.xyz"

cp $quenched_atoms $input_atoms

output_atoms="atoms_anneal.xyz"

ln -sf ../gap_files ./

# 1. Run TurboGAP to anneal the structure, this time at constant pressure, using a barostat
cat <<EOF > input
! Species-specific info
atoms_file = '${input_atoms}'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99

! MD options
md_nsteps = 10000                 ! Graphitization (actual time in the paper is 200ps)
md_step = 1.
barostat = berendsen              ! Using barostat (for illustration purposes) 
t_beg = 1000
t_end = 1000
p_beg = 1.0                       ! Initial pressure [bar]
p_end = 1.0                       ! Final pressure   [bar]
tau_t = 100.

write_thermo = 1
write_xyz = 200
EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap md >  turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"


sleep 1

python3 get_final_config.py $output_atoms

mv thermo.log thermo_3.log
mv trajectory_out.xyz traj_3.xyz
