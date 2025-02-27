module --force purge

module load gcc
# module load openblas
module load openmpi
module load cuda
module load intel-oneapi-mkl/2022.1.0

export HOP_ROOT=/projappl/project_2006384/hop/hop

# Change the Makefile to make sure that we are loading the correct makefile

#cp ./Makefile_debug ./Makefile
# cp ./Makefile_omp_debug ./Makefile

make deepclean
make

