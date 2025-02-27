# Notes for scaling tests 

## Weak Scaling 
Increase the number of cpus/gpus _with_ the system size 
- 1 node, 1 gpu, ntasks-per-process=1  , n_atoms = 100
- 1 node, 2 gpu, ntasks-per-process=2  , n_atoms = 200
- 1 node, 4 gpu, ntasks-per-process=4  , n_atoms = 400
- 2 node, 4 gpu, ntasks-per-process=4  , n_atoms = 800
- 4 node, 4 gpu, ntasks-per-process=4  , n_atoms = 1600
- 8 node, 4 gpu, ntasks-per-process=4  , n_atoms = 3200
- 16 node, 4 gpu, ntasks-per-process=4 , n_atoms = 6400



## Strong Scaling 
Increase the number of cpus/gpus for the same system size 
- 1 node, 1 gpu, ntasks-per-process=1
- 1 node, 2 gpu, ntasks-per-process=2
- 1 node, 4 gpu, ntasks-per-process=4
- 2 node, 4 gpu, ntasks-per-process=4
- 4 node, 4 gpu, ntasks-per-process=4
- 8 node, 4 gpu, ntasks-per-process=4
- 16 node, 4 gpu, ntasks-per-process=4
