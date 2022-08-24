# Compile with gpu
Load the modules:


``` 
module load gcc; module load openblas; module load openmpi ; 
module load cuda; 
```  
Run with:

``` 
srun  --time=00:15:00 --partition=gputest --account=project_2000634 --nodes=1 --ntasks-per-node=4  --cpus-per-task=32 --gres=gpu:a100:4 ../turbogap_dev/bin/turbogap  predict
``` 


Temporary fix:
```
rm cuda_wrappers.o;
nvcc -c -O2 -arch=sm_80 src/cuda_wrappers.cu; 
make -B;
```

