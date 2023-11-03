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


sim_name="volume_mc_oxygen"

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


mc_nsteps = 200

n_mc_types = 3        
mc_types = 'move' 'volume' 'md'                  ! We now specify an MD move type for hamiltonian MC
                                
mc_acceptance = 1 1 1
mc_move_max = 0.2    

mc_lnvol_max = 0.01                              ! Maximum lnvol to modify the volume

! Specify MD configuration
t_beg = 300
t_end = 300
p_beg = 1.0 
p_end = 1.0 

md_nsteps = 30
md_step = 0.1                                    ! Reducing the timestep  
barostat = berendsen
tau_t = 100.

write_xyz = 10
EOF

echo "> Running: Turbogap for $sim_name"

srun turbogap mc >  turbogap.out_$sim_name

echo ">>> Completed: Turbogap for $sim_name <<<"

cp mc_current.xyz $output_atoms

