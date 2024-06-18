#!/usr/bin/bash

calc_type=$1
valid_calculation_types=(prediction standard_gcmc xps_optimization_gcmc)

which turbogap

if [ ! $? ]; then
    echo "> Binary for turbogap has not been installed! Check that you have run 'make clean all' in turbogap/ and specified the correct makefile for your system in Makefile and exported the turbogap/bin PATH "
    exit 1
fi

binary=$(which turbogap)

run_turbogap="mpirun -np 4 $binary "

run_dir="runs/${calc_type}"
cwd=$(pwd)
if [ ! -f input_files/input_${calc_type} ]; then
    echo -e "> $calc_type is not a valid calculation type! \n>> -- Use one of:"
    for i in "${!valid_calculation_types[@]}"; do
	echo "     '${valid_calculation_types[$i]}'"
    done
else

    if [ -d $run_dir ]; then
	rm -rf $run_dir
    fi

    mkdir -p $run_dir
    cp input_files/input_${calc_type} ${run_dir}/input
    cp xps_spectra_interp.dat ${run_dir}/
    cp plot_files/plot_data.py ${run_dir}/
    cp atoms.xyz ${run_dir}/
    cd $run_dir

    ln -sf ${cwd}/gap_files ./

    if [ "${calc_type}" == "prediction" ]; then
	# Using mc mode!
	echo "> Running '$run_turbogap predict > output.${calc_type}' in ${cwd}/${run_dir}"

	$run_turbogap predict > output.${calc_type}
	echo ">> Finished '$run_turbogap predict > output.${calc_type}' in ${cwd}/${run_dir}"

    else
	echo "> Running '$run_turbogap mc > output.${calc_type}' in ${cwd}/${run_dir}"

	$run_turbogap mc > output.${calc_type}
	echo ">> Finished '$run_turbogap mc > output.${calc_type}' in ${cwd}/${run_dir}"

	python3 plot_data.py $calc_type

    fi
    cd -
fi
