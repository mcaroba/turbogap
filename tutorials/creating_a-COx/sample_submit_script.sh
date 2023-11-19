#!/bin/bash -l
#SBATCH --job-name=turbogap
#SBATCH --account=project_
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

