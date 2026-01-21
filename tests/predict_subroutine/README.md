# Description

Test program with the carbon potential and diamond structure, which splits the mpi_comm_world into two groups, and launches two instances of turbogap_predict at once, one with each mpi group. The input filenames are different, each instance is computing a different atomic structure.



# Compile and run

compile the test:

```bash
mpif90 -I../../include -o main.x main.f90 ../../lib/libturbogap.a -llapack -lblas
```

potential files (from tutorial "Simple molecular dynamics"):

```bash
wget https://zenodo.org/record/4000211/files/carbon.gap.tar.gz
tar -xvf carbon.gap.tar.gz
```

run:

```bash
mpirun main.x
```
