module --force purge
module load PrgEnv-gnu-amd/8.4.0
module load LUMI/23.09; 
module load partition/G;
module load rocm/5.4.6;

# Change the Makefile to make sure that we are loading the correct makefile
cp Makefile_gnu Makefile

make clean all 

cp bin/* bin_debug/