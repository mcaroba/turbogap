# Performance: speed, sizing and where the time goes

How fast TurboGAP runs, how many atoms to give each CPU core or GPU, what
experimental biasing (MAD) costs on top of the potential, and what a profile
of each build shows. For how to profile, see `docs/PROFILING.md`; for how to
launch a run, the README's "Running in parallel".

## What was measured, and on what

**System.** Amorphous carbon monoxide, the `CO` set of
[turbogap_tests](https://github.com/TiganyZ/turbogap_tests): the 7,176-,
24,219- and 57,408-atom cells as they are there, a 307,671-atom bulk cell, and
that cell tiled 2x1x1 (615,342 atoms) and 2x2x1 (1,230,684). The potential is
`CO/gap_files/CO.gap`: two `soap_turbo` descriptors (C- and O-centred,
n_max 8+8, l_max 8, zeta 4, trivially compressed), `distance_2b`, `angle_3b`
and `core_pot` terms, and two local properties (`core_electron_be`,
`hirshfeld_v`) whose gradients are formed only when XPS or van der Waals
forces read them.

**Run.** NVE MD from velocities drawn at 300 K, `dd_type = rcb` with
`dd_integrate = domain`, no trajectory output. Each case ran for 2 and for 42
steps (12 for the MAD runs); the cost of a step is the difference of the two
wall times, less each run's setup and structure read, divided by the steps
between them. One run per case; repeats on Roihu vary by about 5%, and
partial-node GPU jobs that share a node with other users by up to 25%, so
the GPU runs used whole nodes.

**Machines.** CSC Roihu.

| | GPU node | CPU node |
| --- | --- | --- |
| hardware | 4 x NVIDIA GH200 (Grace: 72 Arm cores and 95 GB HBM3 per GPU) | 2 x AMD EPYC 9965, 384 cores, 760 GB |
| software | GCC 13.4, CUDA 12.6.3, Open MPI 5.0.10, OpenBLAS 0.3.30 | GCC 15.2, Open MPI 5.0.10, OpenBLAS 0.3.30 |
| build | `compile_gpu.sh`, `CSC-Roihu_gfortran_openblas_cuda` | `make OPENMP=0`, `CSC-Roihu_gfortran_openblas_mpi` |
| layout | one rank per GPU, 72 threads, `--gpu-bind=closest` | one rank per core, one thread |

**Code.** Branch `perf/soap-der-one` as of 2026-10-31; the Debye and XPS
rows of the MAD table, `perf/mad-rank-local` (`b156426`).

## Speed of the potential

Microseconds of wall time per MD step per atom. In brackets, that times the
GPUs or cores used: the cost in device- or core-microseconds per atom-step,
which stays constant when a run scales perfectly.

| atoms | 1 GH200 | 4 GH200 | 1 CPU node (384) | 2 CPU nodes (768) |
| ---: | --- | --- | --- | --- |
| 7,176 | 5.3 | 4.2 (17) | 6.7 (2,580) | |
| 24,219 | 3.5 | 2.3 (9.2) | 5.4 (2,070) | |
| 57,408 | 3.0 | 1.4 (5.8) | 5.2 (1,980) | 2.7 (2,090) |
| 307,671 | 2.7 | 0.89 (3.5) | 4.9 (1,870) | 2.6 (1,990) |
| 615,342 | 2.7 | 0.83 (3.3) | 4.9 (1,880) | 2.5 (1,900) |
| 1,230,684 | 2.7 | 0.78 (3.1) | 5.0 (1,910) | 2.4 (1,880) |

One GH200 does the work of about two of these 384-core CPU nodes.

**Atoms per GPU.** One GPU reaches its best, 2.7 microseconds per atom-step,
from about 50k atoms; below that the fixed costs of a step dominate, and at 7k
atoms each atom costs twice as much. Over four GPUs the efficiency against
one is 77% at 77k atoms per GPU and 85% at 300k per GPU. Give each GPU 100k
atoms or more.

**Atoms per CPU core.** The cost per core is flat from 150 atoms per core up,
on one node and on two (94-100% efficient from one node to two), and is still
73% of that at 19 atoms per core. 100 or more per core is the comfortable
range. Threads do not help on the CPU: 16 per rank reach about 60% of the
speed of 16 ranks, so run one rank per core.

## What experimental biasing costs

The same 307,671 atoms, microseconds per step per atom, with the factor over
plain MD run in the same job:

| | 1 GH200 | 4 GH200 | 1 CPU node | 2 CPU nodes |
| --- | --- | --- | --- | --- |
| plain MD | 2.7 | 0.90 | 4.9 | 2.6 |
| XRD from the pdf | 6.1 (2.3x) | 2.6 (2.9x) | 26 (5.3x) | 14 (5.3x) |
| Debye XRD | 11 (4.4x) | 3.9 (4.1x) | 9.2 (1.9x) | 5.4 (2.1x) |
| XPS | 6.6 (2.6x) | 2.6 (2.7x) | 5.3 (1.09x) | 2.9 (1.12x) |

The pdf, the Debye sum and the XPS spectrum are properties of the whole
system, computed under a decomposition from a gathered copy of the positions
or binding energies. Each rank broadens only its own atoms' XPS lines and
sums only its own atoms' Debye pairs, so those two scale like the potential.
XRD from the pdf does not yet, and is now the costliest on CPUs. The GPU
build runs the Debye sum as host code, and a CPU node now runs XPS faster
than one GH200; both are where the GPU build has most to gain.

Until `b156426` every rank also broadened every atom's XPS line and held the
Debye gradient of every atom:

| per step, 307k atoms | 1 GH200 | 4 GH200 | 1 CPU node | 2 CPU nodes |
| --- | --- | --- | --- | --- |
| XPS, s | 2.20 -> 2.02 | 0.88 -> 0.79 | 7.73 -> 1.64 | 7.26 -> 0.89 |
| Debye XRD, s | 3.43 -> 3.43 | 1.25 -> 1.21 | 5.98 -> 2.83 | 4.92 -> 1.66 |
| XPS, peak GB per rank | 7.3 -> 5.9 | 3.8 -> 2.4 | 1.8 -> 0.6 | 1.8 -> 0.6 |
| Debye XRD, peak GB per rank | 6.0 -> 6.0 | 3.6 -> 2.5 | 1.9 -> 0.8 | 2.2 -> 1.0 |

**Memory.** Peak resident memory per rank (sacct MaxRSS) with XPS is now that
of plain MD on CPUs, 0.6 GB at 307k atoms on 384 ranks, and Debye XRD adds
0.2-0.4 GB. On four GH200s plain MD peaks at 1.8 GB per rank and either term
at 2.4-2.5 GB.

## Where the time goes

**GPU** (four GH200, 307k atoms, nsys on one rank). The device is busy about
half of each 0.27 s step. Kernel time per step, 133 ms, is led by the
per-pair contraction of the SOAP forces (23 ms), the radial basis (16 ms), the
two-body term (12 ms), the angular expansion derivatives (11 ms) and the
three-body term (10 ms). The other half is host code with the device idle,
spread thinly: zeroing and copying result arrays, accumulating each batch's
results, assembling contributions over the decomposition, the neighbour
geometry, preparing the next batch, and Open MPI/UCX bookkeeping around
device allocations.

**CPU** (384 ranks, 307k atoms, perf on four ranks). Three SOAP routines take
most of a rank's compute: the expansion (`get_soap`, 22%), forming the
descriptor's derivatives with respect to every neighbour
(`soap_derivatives_from_expansion`, 21%) and contracting them with the
model's gradient (`get_soap_energy_and_forces`, 12%). The device build no
longer forms those derivatives -- it contracts first, which took its device
time per step from 0.50 to 0.19 s on four GPUs -- and the same change would
carry over to the host build.
Waiting in MPI is the next largest item on ranks 1-383, and rank 0 spends a
few percent formatting output.

## The checks behind these numbers

What the measured code was tested with, and where:

| check | what it compares | machine |
| --- | --- | --- |
| `tests/regression --cpu`, all 44 cases | outputs against the previous binary, bit for bit; only the XPS and Debye cases differ, as intended | alt: x86, 12 cores |
| `tests/regression --gpu`, all 41 cases with data | device against host within tolerance; the same 36 agree as before | Roihu GH200 |
| `tests/soap_descriptor_gpu` | the device descriptor kernels bit for bit against the originals, and the force contraction against the explicit derivative chain (to 4e-15) | alt: RTX 5070 |
| `tests/xps_fd` | XPS forces and virial against finite differences of the XPS energy, serial and on two ranks | alt |
| `tests/xrd_debye` | the Debye pattern against an independent numpy sum, forces and virial against finite differences | alt |
| XPS and Debye on 3 ranks under `dd_type` none, rcb and legacy | forces and energies against the binary before `b156426`: identical | alt |
| `tests/pdf_virial`, `tests/xrd_lp_pdf` (gradient) | pdf, structure factor and XRD forces and virial against finite differences, with the device build | alt: RTX 5070 |
