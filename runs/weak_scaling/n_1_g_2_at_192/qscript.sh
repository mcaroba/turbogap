#!/bin/bash -l
#SBATCH --job-name=n_1_g_2_at_192
#SBATCH --account=project_2002883
#SBATCH -p gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --ntasks-per-node=2
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -t 0-2:00:00
#SBATCH -o out_files/out_n_1_g_2_at_192
#SBATCH --array=1


module --force purge

module load gcc
# module load openblas
module load openmpi
module load cuda
module load intel-oneapi-mkl/2022.1.0

export HOP_ROOT=/projappl/project_2006384/hop/hop


srun /projappl/project_2002883/turbogap_mbd/bin/turbogap predict
