#!/bin/bash -l
#SBATCH --job-name=GPU_test
#SBATCH --account=project_2002883
#SBATCH -p gputest
#SBATCH --gres=gpu:v100:4
#SBATCH --ntasks-per-node=4
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -t 0-0:05:00
#SBATCH -o out_files/out_%a
#SBATCH --array=1


module --force purge

module load gcc
# module load openblas
module load openmpi
module load cuda
module load intel-oneapi-mkl/2022.1.0

export HOP_ROOT=/projappl/project_2006384/hop/hop


srun ../bin/turbogap predict
