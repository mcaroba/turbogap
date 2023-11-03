
sim_name="gcmc_standard_oxygen"

prev_dir="../1.make_amorphous_carbon"
relaxed_atoms="${prev_dir}/atoms_relax.xyz"

if [ ! -f $relaxed_atoms ]; then
    echo "!! Relaxed atoms don't exist! Go to $prev_dir and run `bash 4.run_relax.sh` "
    exit 1
fi

input_atoms="atoms.xyz"

cp $relaxed_atoms $input_atoms

output_atoms="atoms_gcmc.xyz"

ln -sf ../gap_files ./

# The chemical potential used for the simulation:
#   mu = -5.16 eV is related to half the binding energy of O_2 at 1atm, 300K.
# We will use a much larger one to see saturation of O
mu=0.0

# 1. Run TurboGAP with GCMC to add oxygen to amorphous carbon structure. 
cat <<EOF > input
! Species-specific info
atoms_file = '${input_atoms}'
pot_file = 'gap_files/CO.gap'
n_species = 2
species = C O
masses = 12.01 15.99


mc_nsteps = 10000                         ! Number of mc trial steps to be performed

n_mc_types = 3                           ! Number of mc trial types
mc_types = 'move' 'insertion' 'removal'  ! MC types can be: 'insertion' 'removal' 'md' 'swap' 'move'
                                         !                  'volume'
mc_acceptance = 1 1 1                    ! Ratios for choosing the respective trial moves (all equally likely here)

mc_move_max = 0.2                        ! Maximum distance for the move of a particular particle

mc_mu = ${mu}                              ! gcmc: Chemical potential [eV]
mc_species = 'O'                         ! gcmc: species to insert / remove 
mc_min_dist = 0.1                        ! gcmc: minimum distance between particles for insertion

t_beg = 300.0                            ! Temperature

! Note: The following files are written to every write_xyz steps
!       mc.log: self-explanatory,
!          format: mc_step  mc_move  accepted  E_trial  E_current  N_sites(trial)
!                                                                     N_gcmc_species(trial)
!       mc_all.xyz:     an appended file which contains all accepted moves
!       mc_trial.xyz:   a single configuration which contains the trial move
!       mc_current.xyz: a single configuration which contains the current accepted



write_xyz = 200
EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap mc >  turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"

cp mc_current.xyz $output_atoms
cp mc.log mc_1.log
