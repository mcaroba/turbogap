
sim_name="md_diamond_randomise"

input_atoms="atoms.xyz"
output_atoms="atoms_randomise.xyz"

ln -sf ../gap_files ./

# 1. Create diamond structure (1000 atoms in atoms.xyz file)
echo "> Running: python create_diamond.py"
python3 create_diamond.py

cp diamond.xyz $input_atoms

# 2. Run TurboGAP to randomize the structure.
cat <<EOF > input
! Species-specific info
atoms_file = '${input_atoms}'  ! Input file
pot_file = 'gap_files/CO.gap'  ! path to gap_files
n_species = 2	               ! > Actually the number of species in atoms.xyz is 1 (C), 
                               !   but we will add oxygen in future simulations
species = C O
masses = 12.01 15.99	       
			       
! MD options		       
md_nsteps = 5000	       ! Number of MD steps (5ps randomise - actual time in paper is 20ps)
md_step = 1		       ! MD timestep [fs]
thermostat = bussi	       ! Either bussi / berendsen

t_beg = 9000		       ! Initial temperature [K]
t_end = 9000		       ! Final temperature   [K]
tau_t = 100.		       ! Time constant [fs]

! Output
write_thermo = 1	       ! Write thermodynamic information every step 
                               !       (Step, Time, Temp, Kin_E, Pot_E, Pres)

write_xyz = 200                ! Write extended xyz trajectory (trajectory_out.xyz) every 200 steps

! Core potential cutoff
core_pot_cutoff = 1.

! Neighbors buffer for rebuilding the neighbor list
neighbors_buffer = 0.5

EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap md > turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"


sleep 1

python3 get_final_config.py $output_atoms

mv thermo.log thermo_1.log
mv trajectory_out.xyz traj_1.xyz
