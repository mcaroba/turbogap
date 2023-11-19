#!/bin/bash -l
#SBATCH --job-name=turbogap
#SBATCH --account=project_2008666
##SBATCH -p test
#SBATCH -p medium
##SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH -t 0-00:10:00
##SBATCH -t 0-00:10:00
##SBATCH --mem-per-cpu=4000
##SBATCH -o out_files/out_%a
##SBATCH --array=0-10


module reset
module load gcc/11.2.0
module load openblas/0.3.18-omp
module load openmpi/4.1.2
module load python-data/3.9-22.04

cpus=$SLURM_NTASKS

###--- Turbogap 
export PATH="/projappl/project_2008666/turbogap/bin:${PATH}"


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
thermostat = berendsen	       ! Either berendsen / bussi

t_beg = 9000		       ! Initial temperature [K]
t_end = 9000		       ! Final temperature   [K]
tau_t = 100.		       ! Time constant for berendsen [fs]

! Output
write_thermo = 1	       ! Write thermodynamic information every step 
                               !       (Step, Time, Temp, Kin_E, Pot_E, Pres)

write_xyz = 200                ! Write extended xyz trajectory (trajectory_out.xyz) every 200 steps
EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap md > turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"


sleep 1

python3 get_final_config.py $output_atoms

mv trajectory_out.xyz traj_1.xyz
