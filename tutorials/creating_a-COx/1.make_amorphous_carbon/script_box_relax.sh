
sim_name="md_box_relax"

annealed_atoms="atoms_anneal.xyz"

if [ ! -f $annealed_atoms ]; then
    echo "!! Annealed atoms don't exist! run `bash 3.run_anneal.sh` "
    exit 1
fi


input_atoms="atoms.xyz"

cp $annealed_atoms $input_atoms

output_atoms="atoms_relax.xyz"

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
optimize = gd                 ! optimize option allows us to specify the type of relaxation
                              !  Here we keep the orthogonality of the box 

! e_tol = 1.d-6               ! Default energy/force tolerances used 
! f_tol = 0.010
md_nsteps = 2000

write_xyz = 200

core_pot_cutoff = 1.
neighbors_buffer = 0.5

EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap md >  turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"


sleep 1

python3 get_final_config.py $output_atoms

mv thermo.log thermo_4.log
mv trajectory_out.xyz traj_4.xyz

