#!/usr/bin/bash

cat <<EOF > sample_submit_script.sh
#!/bin/bash -l
#SBATCH --job-name=turbogap
#SBATCH --account=project_2008666
#SBATCH -p small
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH -t 0-00:10:00

module reset

module load gcc/11.3.0
module load openmpi/4.1.4
module load intel-oneapi-mkl/2022.1.0
module load python-data/3.9-22.04

cpus=$SLURM_NTASKS

###--- Turbogap
export PATH="/projappl/project_2008666/turbogap/bin:${PATH}"

EOF

cd 1.make_amorphous_carbon
# Change the number of atoms to something smaller
sed -i 's/(5,5,5)/(3,3,3)/g' create_diamond.py
cd ../


module load python-data/3.9-22.04
pip install ase --user
