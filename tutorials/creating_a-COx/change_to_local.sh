#!/usr/bin/bash

# If you want to use this tutorial locally we can change the scripts to reflect this

change_run(){
    run_command=$1
    scripts=$(ls --color=never *.sh)
    # replace srun with mpirun
    for script in $scripts; do
	sed -i "s/srun/${run_command}/g" $script
	sed -i "s/sbatch/bash/g" $script
    done
}

# specify the number of processors
n=4
run_command="mpirun -np $n "

# Remove module commands
sed -i 's/module/# module/g' sample_submit_script.sh

# NOTE: You must change your path in sample_submit_script.sh to where you've installed turbogap/bin

cd 1.make_amorphous_carbon
change_run
# Change the number of atoms to something smaller
sed -i 's/(5,5,5)/(2,2,2)/g' create_diamond.py
cd ../

cd 2.standard_gcmc
change_run
cd ../

cd 3.volume_mc
change_run
cd ../
