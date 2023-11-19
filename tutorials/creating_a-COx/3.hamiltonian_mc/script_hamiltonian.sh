
sim_name="hamiltonian_mc"

prev_dir="../2.standard_gcmc"
gcmc_atoms="${prev_dir}/mc_current.xyz"

if [ ! -f $gcmc_atoms ]; then
    echo "!! Gcmc atoms don't exist! Go to $prev_dir and run `bash 3.run_anneal.sh` "
    exit 1
fi

input_atoms="atoms.xyz"

cp $gcmc_atoms $input_atoms

output_atoms="atoms_gcmc.xyz"

ln -sf ../gap_files ./

# 1. Run TurboGAP with GCMC to add oxygen to amorphous carbon structure.
cat <<EOF > input
! Species-specific info
atoms_file = '${input_atoms}'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99

md_nsteps = 20                           ! specifying the number of steps for velocity verlet
md_step = 0.1                            ! 0.1fs timestep

mc_nsteps = 50                           ! Number of mc trial steps to be performed
n_mc_types = 1
mc_types = 'md'
mc_acceptance = 1                        ! Relative rate of choosing the respective trial moves
mc_hamiltonian = .true.                  ! For Hamiltonian Monte-Carlo
                                         !   (NVE ensemble used to increase trial acceptance) for md
					 !   move type
write_xyz = 10
EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap mc >  turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"

cp mc_current.xyz $output_atoms
