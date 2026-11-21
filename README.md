# The TurboGAP code

**TurboGAP** (c) 2018-2026 by **Miguel A. Caro** and others (see "contributors" below
for detailed authorship info).

[www.turbogap.fi](http://www.turbogap.fi)

## Contributors, copyright and license

The following people, listed chronologically, have contributed code or ideas to the **TurboGAP**
code. Those whose names are in bold have contributed code _to the master branch_ (relevant for
purpose of copyright; each file in the **TurboGAP** repo that contains code has a copyright
statement at the beginning).

- **Miguel A. Caro** (Aalto University)
- Patricia Hernández-León (Aalto University)
- Suresh Kondati Natarajan (formerly @ Aalto University)
- **Albert P. Bartók-Pártay** (Warwick University)
- Eelis V. Mielonen (formerly @ Aalto University)
- **Heikki Muhli** (Aalto University)
- **Mikhail Kuklin** (formerly @ Aalto University)
- Gábor Csányi (University of Cambridge)
- **Jan Kloppenburg** (Aalto University)
- **Richard Jana** (Aalto University)
- **Tigany Zarrouk** (Aalto University)
- **Uttiyoarnab Saha** (Aalto University)

**TurboGAP** is licensed under the Academic Software License (ASL), an "available source"
non-commercial license. This means that you are free to use and distribute the code for
non-commercial academic research (or teaching) under the terms of the license. See
`LICENSE.md` for details. If you want to obtain a commercial license for **TurboGAP**, please
contact Miguel Caro (mcaroba@gmail.com).

Some third-party code is included with **TurboGAP** for convenience, under the
`src/third_party` directory. These codes are licensed independently from **TurboGAP** and their
respective licenses have been verified to be compatible for redistribution with **TurboGAP**.
They may be redistributed separately from **TurboGAP** under their respective licenses.
Refer to each piece of software in that subdirectory for further information.

The **soap_turbo** submodule, under `src/soap_turbo`, is a separate distribution from
**TurboGAP**, but it is required
for running **TurboGAP**, since it contains the `soap_turbo` routines. These routines are
copyright (c) of Miguel A. Caro and they are also distributed under the ASL. Therefore, you
can freely use this code for non-commercial academic research or teaching. If you want to
obtain a commercial license for **soap_turbo** please contact Miguel Caro (mcaroba@gmail.com).

## Overview of the code

The **TurboGAP** code is a specialized sampling/prediction engine which is made
for use with Gaussian Approximation Potentials that use many-body `soap_turbo`
descriptors.

Unique features are included in the code (e.g., Molecular Augmented Dynamics and experimental spectra prediction and
full support for van der Waals corrections).
Feedback can be sent to Tigany Zarrouk (tiganyzarrouk@gmail.com) or Miguel Caro (mcaroba@gmail.com)
or left on the Issues section of the Github page.

## Installation Quick Start

More detailed information on installation is found further on.

_**tl;dr, CPU**_ (standard linux, prerequisites installed):

```sh
git clone --recursive https://github.com/mcaroba/turbogap.git
cd turbogap
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
./compile_cpu.sh
export PATH="$(realpath bin):$PATH"
```

_**tl;dr, GPU**_ (CUDA or ROCm, on top of the above):

```sh
git clone --recursive https://github.com/mcaroba/turbogap.git
cd turbogap
export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda
export CUDA_ARCH=sm_80    # your card: sm_80 A100, sm_90 H100/GH200, ...
./compile_gpu.sh
export PATH="$(realpath bin-gpu):$PATH"
```

_**tl;dr, GPU through Kokkos**_ (optional; needed to run TurboGAP from LAMMPS):

```sh
# If you haven't installed Kokkos you can use the script in turbogap/tools
tools/install_kokkos.sh              # builds Kokkos, ~15 min, once
export KOKKOS_ROOT=$HOME/.local/kokkos-4.7.04 # Change this to where your kokkos installation is
KOKKOS=1 ./compile_gpu.sh
export PATH="$(realpath bin-kokkos-gpu):$PATH"
```

_**tl;dr, developer tooling**_ (only needed to commit, or to run the python
test drivers):

```sh
tools/setup_dev_env.sh
export PATH="$HOME/.venvs/turbogap-tools/bin:$PATH"
tools/setup_dev_env.sh --check
```

### Prerequisites

To use TurboGAP for CPU architectures, you must have the libraries `openmpi`, `lapack` + `blas` (or `openblas`)
and the Fortran compilers necessary for your system (e.g. `gcc`).

For GPU installations, you can use the CUDA/ROCm which is done via the `hop`
library that is included, or via Kokkos which is used for LAMMPS integration.

_For more detailed info on download, installation, etc., you can visit the
[TurboGAP wiki](https://turbogap.fi/wiki/index.php/Installation)._

### Getting the code

To get the **TurboGAP** code and the necessary **soap_turbo** routines, do a recursive
`git clone` with `--depth=1` to speed up things:

```sh
git clone --recursive --depth=1 http://github.com/mcaroba/turbogap.git
```

### What gets installed

A `--recursive` clone brings three submodules:

| submodule                                                 | what it is                                                          |
| --------------------------------------------------------- | ------------------------------------------------------------------- |
| `src/soap_turbo`                                          | CPU SOAP turbo routines uses                                        |
| `src/soap_turbo_gpu`                                      | GPU SOAP routines                                                   |
| `src/hop`                                                 | Header Only Porting for interoperability on CUDA/ROCm architectures |
| `tools/setup_dev_env.sh` creates a virtual environment at |
| `~/.venvs/turbogap-tools` holding, all pinned:            |

| tool                        | why                                                                                                                                 |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `fprettify`, `clang-format` | the formatting the pre-commit hook applies. Pinned because an unpinned formatter rewrites files nobody touched                      |
| `pre-commit`                | the hook itself: formatting, the generated keyword reference, whitespace                                                            |
| `numpy`                     | what the python test drivers compare against. Without it `tests/ir_fft`, `tests/mad_ir` and `tests/ipi_socket` skip rather than run |
| `i-PI`                      | drives TurboGAP over a socket in `turbogap ipi` mode, for path-integral and other advanced dynamics.                                |
| `cmake`                     | only installed by `tools/install_kokkos.sh`, and only when the machine has none of its own. Kokkos has no other build system        |

### Building Manually

To build the **TurboGAP** binary and library, you need to select the options
that best match your architecture by exporting the environment variable
`TURBOGAP_ARCH=<your_architecture>` which is used in the `Makefile` e.g.:

```sh
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
```

This will include the architecture specific examples found in `makefiles`, e.g.
the above will include the specification for that architecture,

```sh
include makefiles/Makefile.Ubuntu_gfortran_mpi
```

A list of example makefiles is provided under the `makefiles/` directory for different systems.

Once you are happy with your `Makefile`, to build the code just type

```sh
make j<N_processes>
```

where `<N_processes>` is the number of processes you wish to build the code with.

Then add `turbogap/bin` to your path.

```sh
turbogap_dir=$(realpath bin)
export PATH="$turbogap_dir:$PATH"
echo "export PATH=\"${turbogap_dir}:\$PATH\"" >> ~/.bashrc
```

If you need to rebuild the code,
you can `make clean; make` or `make deepclean; make`.

## Running TurboGAP

TurboGAP reads the input file `input` from the directory it runs in, and it can run in several modes:

```sh
turbogap predict    # single points: energy, forces, stress, observables
turbogap md         # molecular dynamics, or a relaxation (optimize = ...)
turbogap mc         # Monte Carlo
turbogap ipi        # a force engine for i-PI: path integrals and more
```

`turbogap --help <mode>` lists the keywords that act in that mode and what
each needs; `turbogap --help gap` lists those of the potential file. All of
them are in `docs/keywords.html` (searchable) and `docs/keywords.md`.

### Which setup

| run | binary | ranks x threads | in `input` |
| --- | --- | --- | --- |
| the GAP on CPUs | `bin/turbogap` | one rank per core, 1 thread | `dd_type = rcb`, `dd_integrate = domain` |
| the GAP on GPUs | `bin-gpu/turbogap` | one rank per GPU, its cores as threads | the same, with more than one GPU |
| MBD on CPUs | `bin/turbogap` | 4-16 threads per rank | decomposition is ignored |
| MBD on GPUs | `bin-gpu/turbogap` | one rank per GPU, its cores as threads | decomposition is ignored |
| experimental biasing | either, by observable | as the GAP | `dd_type = rcb`, `dd_integrate = domain` |

On large cells a 384-core CPU node does the work of about 1.9 GH200s. Give
each CPU core a few hundred atoms or more, and each GPU 100k or more.
`docs/PERFORMANCE.md` has the measurements behind this section.

### Domain decomposition

By default (`dd_type = none`) every rank holds every atom and rank 0
integrates: each step the structure is broadcast and the forces reduced
through rank 0, which stops scaling at about two CPU nodes. Two lines in
`input` give each rank a region of the cell, with a halo `rcut_max` wide, and
let each region integrate its own atoms:

```
dd_type = rcb             # none (default), legacy or rcb
dd_integrate = domain     # rank0 (default) or domain
```

`rcb` cuts the cell by recursive bisection weighted by pair count, so it takes
any number of ranks and balances pores, clusters and surfaces; `legacy` is the
older grid balanced by atom count. Results differ only in the last digits of
sums. 308k-atom CO, seconds per MD step on 384-core CPU nodes:

| nodes | 1 | 2 | 4 | 6 |
| --- | ---: | ---: | ---: | ---: |
| `none` | 1.50 | 1.10 | 1.18 | 1.08 |
| `rcb` + `domain` | 0.46 | 0.25 | 0.13 | 0.10 |

What the domains do not cover falls back by itself, and the top of the output
says so:

- the whole run as `none`: MBD, a run without `do_prediction`, a build
  without MPI;
- integrated on rank 0, still evaluated on the domains: optimisers,
  thermostats other than none, Berendsen and Bussi, a barostat or box
  scaling, a variable time step, electronic stopping, MC, nested sampling,
  MAD-IR, and the pdf, structure factor, XRD (Debye included) and ND
  observables.

### Environment

| variable | value | why |
| --- | --- | --- |
| `OMP_NUM_THREADS` | threads per rank; always set it | unset, a rank divides the cores it can see by the ranks on its node, which is 1 once each rank is bound to its cores |
| `OMP_PLACES`, `OMP_PROC_BIND` | `cores`, `close` | each rank's threads stay on its own cores |
| `OPENBLAS_NUM_THREADS` | `1` | OpenBLAS starts no threads of its own inside each rank |
| `CUDA_VISIBLE_DEVICES` (ROCm: `HIP_VISIBLE_DEVICES`) | one card per rank | only under `mpirun`; `srun --gpu-bind=closest` does it under Slurm |

With more than one thread the top of the output prints `OpenMP threads/rank:`;
check it. In the scripts below, add your site's `--account` and `--partition`
and load the modules TurboGAP was built with. Submit from the directory that
holds `input`: `sbatch job.sh`.

### CPU nodes

One rank per physical core and one thread, with `rcb` and `dd_integrate =
domain`. Threads help only when memory forces fewer ranks: 16 threads per rank
run at about 65% of the speed of 16 ranks. A core is within 15% of its best
from about 150 atoms.

```sh
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=384       # one per physical core
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

module load gcc openmpi openblas
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
srun --cpu-bind=cores turbogap md
```

On a workstation: `OMP_NUM_THREADS=1 mpirun -np <cores> --bind-to core turbogap md`.

### GPUs

One rank per GPU, with that GPU's share of the node's cores as threads, which
prepare the next batch while the GPU runs. Give each GPU at least 50k atoms,
and 100k or more when spreading over several; there, `rcb` and `dd_integrate =
domain` take a 308k-atom step on four GH200s from 0.35 s to 0.27 s. Batch
sizes follow from the card's memory, so leave `gpu_batch_gbytes` unset unless
another program shares the card.

```sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4         # one per GPU
#SBATCH --cpus-per-task=72          # cores per GPU: 288 / 4 on a GH200 node
#SBATCH --gpus-per-node=4
#SBATCH --time=01:00:00

module load gcc openmpi openblas cuda
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1
srun --cpus-per-task=$SLURM_CPUS_PER_TASK --cpu-bind=cores --gpu-bind=closest turbogap md
```

Without `--gpu-bind=closest` every rank can land on GPU 0: each `GPUdev:` line
(stderr) must say its rank sees 1 device. With one GPU, no MPI is needed:
`OMP_NUM_THREADS=<cores> turbogap md`. Under `mpirun`, give each rank its own
card through a wrapper:

```sh
cat > one_gpu.sh << 'EOF'
#!/bin/bash
export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK    # ROCm: HIP_VISIBLE_DEVICES
exec "$@"
EOF
chmod +x one_gpu.sh
OMP_NUM_THREADS=<cores per GPU> mpirun -np <GPUs> ./one_gpu.sh turbogap md
```

### MBD van der Waals

`vdw_type = mbd` makes the whole run replicated, whatever `dd_type` says. Each
rank takes a share of the sites, one per atom, and each site is one iteration
of a threaded loop, so here threads help. A site costs according to the atoms
within `vdw_rcut` (25 A) of it, so dense and dilute systems differ about a
hundredfold. vdW time of one single point:

| | one CPU node (384 cores) | one GH200 node (4 x 72 cores) |
| --- | --- | --- |
| dense: 576-atom black phosphorus, 3,457 atoms around each site | 0.43 s at 96 x 4, 0.48 s at 24 x 16, 0.57 s at 384 x 1 | 0.39 s on the GPUs, 0.54 s on the host cores |
| dilute: 27 P4 clusters, 2,916 atoms, 108 around each site | 0.027 s at 24 x 16, 0.081 s at 384 x 1 | 0.14 s on the host cores |

- On CPUs, start at 4 threads per rank for dense solids and 16 for molecules,
  clusters and porous systems. Each thread holds about 100 MB on a dense
  system.
- On GPUs (`vdw_mbd_device = .true.`, the default) the host threads list each
  site's pairs and solve its polarisabilities, and the device does the rest.
  A site with fewer than `vdw_mbd_device_min_pairs` (100,000) pairs stays on
  the host threads, which are faster there. Dilute systems run best on CPU
  nodes.
- In MD, `vdw_type = ts+mbd` with `mbd_correction_freq = N` recomputes MBD
  every N steps (100 by default) and reuses it in between.
- The potential has to provide Hirshfeld volumes.

`input`, with the cutoffs of `tests/regression/cases/vdw_mbd_cell_mpi2`:

```
vdw_type = mbd
vdw_r0_ref = 2.12             # one value per species
vdw_alpha0_ref = 3.7046
vdw_c6_ref = 110.54
vdw_rcut = 25.
vdw_mbd_rcut = 15.
vdw_mbd_rcut2 = 10.
vdw_2b_rcut = 15.
vdw_2b_rcut2 = 10.
```

```sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=4           # dense; 24 x 16 for dilute systems
#SBATCH --time=01:00:00

module load gcc openmpi openblas
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1
srun --cpus-per-task=$SLURM_CPUS_PER_TASK --cpu-bind=cores turbogap md
```

On GPUs, the GPU script above runs MBD as it is.

### Experimental biasing (MAD)

Each observable named in `exp_labels` adds its mismatch with its data file to
the energy, weighted by `exp_energy_scales`, and with `exp_forces` the
mismatch's gradient to the forces. Keep `rcb` and `dd_integrate = domain`:
every observable is evaluated on the domains and XPS also integrates there,
while the pdf, structure factor, XRD, ND and MAD-IR integrate on rank 0.

The observable decides the hardware. 307k-atom CO, microseconds per
atom-step; the CPU column was measured before the CPU potential became three
times faster, so its numbers are now lower:

| | 1 GH200 | 4 GH200 | 1 CPU node |
| --- | ---: | ---: | ---: |
| plain MD | 2.7 | 0.90 | 4.9 (1.44 now) |
| XRD from the pdf | 6.1 | 2.6 | 26 |
| Debye XRD (`xrd_debye`) | 11 | 3.9 | 9.2 |
| XPS | 6.6 | 2.6 | 5.3 |

XPS and the Debye sum split their work over the ranks like the potential;
XRD from the pdf does not yet, and belongs on GPUs. A CPU node runs XPS faster
than one GH200, and the GPU build does the Debye sum on its host cores.
Memory is modest: on CPUs XPS adds nothing per rank and Debye XRD 0.2-0.4 GB;
on four GH200s either takes a rank from 1.8 to 2.5 GB.

`input`, biasing towards an XPS spectrum (from
`tests/regression/cases/xps_mad_diamond`):

```
n_exp = 1
exp_labels = 'xps'
exp_data_files = 'xps_target.dat'
exp_n_samples = 601
exp_energy_scales = 2.48
exp_energies = .true.
exp_forces = .true.
xps_sigma = 0.2

dd_type = rcb
dd_integrate = domain
```

For XRD, `xrd_debye_mad` (Debye sum) and `xrd_mad` (from the pdf) in the same
directory have the keyword blocks. Run with the CPU or GPU script above.

## TurboGAP Tutorials

For simple test cases of various simulation modes, one can consult `tests/regression/cases`.

For practical use cases (i.e. useful simulations for understanding a particular
research question) one can follow the tutorials which were ran for the TurboGAP
School, and adapt them for their machine: [TurboGAP
School](https://github.com/mcaroba/TurboGAP_School).

Simple tutorials can be found in the [TurboGAP Tutorials](https://github.com/TiganyZ/turbogap_tutorials) repository `git clone https://github.com/TiganyZ/turbogap_tutorials.git`. Further tutorials and all documentation can be found on the [Turbogap Website](https://turbogap.fi/wiki/index.php/Tutorials).

## Testing TurboGAP

The data for tests can be by cloning the [TurboGAP Tests](https://github.com/TiganyZ/turbogap_tests). This can be done using the script `tests/fetch_test_data.sh`.

Instead of doing a shallow clone (with depth=1) you can do the full clone so you have access to the baseline commit for regression tests.

```sh
git clone --recursive http://github.com/mcaroba/turbogap.git
cd turbogap
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
make -j4
turbogap_dir=$(realpath bin)
export PATH="$turbogap_dir:$PATH"
cd tests/regression
./make_baseline.sh
../fetch_test_data.sh
TURBOGAP_KEEP=1 ./run.sh
```

and the tests will be found in `$TMPDIR/turbogap_regression.xxxxx`.

### Testing a GPU build

The same case list runs against the device binary. Build both, then:

```sh
tests/regression/run.sh --gpu      # every case, on bin-gpu/turbogap
tests/regression/run.sh --both     # the host pass, then the device pass
```

`--gpu` compares the device binary against the **host build of the same
source**, not against the frozen baseline: the baseline is a snapshot of an
older commit, and comparing against it would fold every intended change since
into the same number as the host/device difference.

The comparison is numerical rather than bit-exact, through
`tests/regression/compare_tol.py`. The device sums the SOAP batches and the
cuBLAS reductions in a different order, so the last digits move on a run that is
entirely correct; a bit-exact diff would call every device run a failure.
Non-numeric tokens -- species columns, `Properties=` strings -- must still match
exactly, so a structural change is still caught. The default tolerance is
`rtol = atol = 1e-6`, overridable per case with `GPU_RTOL` and `GPU_ATOL` in its
`case.conf`, or globally with `TURBOGAP_GPU_RTOL` and `TURBOGAP_GPU_ATOL`.

### Testing a Kokkos build

```sh
tools/verify_kokkos.sh             # every case; add case names to narrow it
```

Other tests can be done by running the scripts in the `tests/<test_name>/run.sh` directories respectively.
Each honours `TURBOGAP_BIN`, so the same suite can be pointed at either build:

```sh
TURBOGAP_BIN=$(realpath bin-gpu/turbogap) tests/ipi_pimd/run.sh
```

## Developing TurboGAP

Follow the details of how the code is structured in `docs/DEVELOPMENT.md`. Install the developer tooling as stated above.

To develop for TurboGAP, one can install the development tools python environment (which has fprettify and pre-commit and so on which is installed through uv) such that formatting is preserved. Make a new branch or fork and then once ready submit a pull request. Please add tests for your new feature in the ` tests/` folder such that the CI interface can test upon pushing.

If you're including new source files, make sure to run the script which
regenerates the `makefiles/Makefile.deps`.

```sh
python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps
```

This may require you to first install a python environment if you don't have a
reasonable version of python. This can be done using the `source
tools/setup_dev_env.sh` which will install one for you using `uv` which will be
placed in `$HOME/.venvs`

To add new keywords in to TurboGAP, please follow "Documenting a keyword" in
`docs/DEVELOPMENT.md`. Following this format allows for automated
documentation and help information for usage.

### Debugging

Debug flags can be enabled by exporting

```sh
DEBUG=1
```

This will create a `bin-dbg` folder which contains the (now slower) binary.

## Attribution

When using **TurboGAP**, you should give attribution to the
**TurboGAP** author(s). The appropriate way to do that is to provide a link to the
[**TurboGAP** website](https://www.turbogap.fi) and, if you publish results obtained
using **TurboGAP** or the **soap_turbo** library,
even if it is through one of its external interfaces, you should cite:

> **Miguel A. Caro**. _Optimizing many-body atomic descriptors for enhanced computational
> performance of machine learning based interatomic potentials_. [Phys. Rev. B 100, 024112
> (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.024112).

In addition, you should cite any other relevant literature and code websites (e.g., the
original SOAP/GAP papers) as appropriate.
