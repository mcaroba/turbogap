# Compile with gpu

Temporary fix:
```
nvcc -c -O2 --arch=sm_80 src/cuda_wrappers.cu;
make -B;
```

