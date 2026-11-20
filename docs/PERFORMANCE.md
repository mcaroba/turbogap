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
rows of the MAD table, `perf/mad-rank-local` (`b156426`); the rows marked as
measured since, `perf/host-openmp` (`91acaa1`); the MBD section,
`perf/host-openmp` (`9a871ab`; CPU nodes `0548fb4`, the same host code);
2 and 22 steps.

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

One GH200 did the work of about two of these 384-core CPU nodes. The CPU
columns predate the host build's contract-first SOAP forces (`perf/host-openmp`,
2026-11-19), which made it about 3.3 times faster. Measured since, 307,671
atoms:

| | 1 GH200 | 1 CPU node | 2 CPU nodes | 4 CPU nodes | 6 CPU nodes |
| --- | --- | --- | --- | --- | --- |
| microseconds per atom-step | 2.4 | 1.49 (571) | 0.81 (622) | 0.43 (660) | 0.33 (749) |

A CPU node now does the work of about 1.6 GH200s.

**Atoms per GPU.** One GPU reaches its best, 2.7 microseconds per atom-step,
from about 50k atoms; below that the fixed costs of a step dominate, and at 7k
atoms each atom costs twice as much. Over four GPUs the efficiency against
one is 77% at 77k atoms per GPU and 85% at 300k per GPU. Give each GPU 100k
atoms or more.

**Atoms per CPU core.** The cost per core is flat from 150 atoms per core up,
on one node and on two (94-100% efficient from one node to two), and is still
73% of that at 19 atoms per core. 100 or more per core is the comfortable
range. Threads do not help on the CPU: 16 per rank reach about 65% of the
speed of 16 ranks, so run one rank per core -- except for MBD (below).

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

**CPU** (384 ranks, 307k atoms, perf on ranks 0 and 191). The host build now
contracts first too (`soap_contract`, on by default): the model's gradient is
contracted with the expansion before differentiating, so the descriptor's
derivatives with respect to every neighbour are never formed. It took the
step on one node from 1.47 to 0.42 s. What a rank computes is now led by the
contraction (`soap_contract`, 14%), the prediction's gradient
(`get_soap_energy_and_forces`, 6%), the expansion (`get_soap`, 5%) and its
angular part (5%). Waiting in MPI is the largest item on rank 191 (38%);
rank 0 waits little but spends ~8% formatting output.

## MBD van der Waals

MBD runs replicated under every `dd_type`: each rank takes a share of the
sites, and each site -- its SCS polarisabilities, then the MBD energy and
forces of its cluster -- is one iteration of a threaded loop. Measured with
the `vdw_mbd_cell_mpi2` settings (`vdw_rcut` 25 A, MBD cutoffs 15/10 A, 13
frequencies, order 6), vdW time of one single point:

| system | where | master | now |
| --- | --- | ---: | ---: |
| P4 gas, 108 atoms | alt, 1 rank x 1 / 4 / 12 threads | 6.1 s | 0.29 / 0.079 / 0.039 s |
| P4 gas, 2,916 atoms | one CPU node, 384 x 1 | 0.68 s | 0.088 s |
| | one CPU node, 24 x 16 | 7.7 s | 0.044 s |
| | one GH200 node, 4 x 72 | 67 s | 0.18 s |
| black phosphorus, 16 atoms | alt, 1 rank x 4 threads (master: 1) | 217 s | 1.10 s |
| black phosphorus, 576 atoms | one CPU node, 96 x 4 | -- | 0.83 s |
| | one GH200 node, 4 x 72, MBD on the GPUs | -- | 0.53 s |
| black phosphorus, 144 atoms | one GH200, 1 x 72, MBD on the GPU | -- | 0.36 s |

Master cannot run the dense cells in these layouts: each of its sites
allocated several 3N x 3N matrices (0.86 GB each for a cluster of 3,500
atoms), and with them 72 threads ran a GH200 node out of host memory. Through
`f28397a` the CPU build's output was master's byte for byte. From `62ea743`
the host's default path sums in its own order (below); the printed outputs
of the MBD regression and golden cases, at 1, 2, 4 and 12 threads and on 4
and 6 ranks, are still master's, and so is the 2,916-atom cell's at 384 x 1
on a CPU node (with threads, the GAP's SOAP energy moves in the last printed
digit there, as master's own does). MBD on the device and on the host
threads print the same digits on every cell above.

The P4 cell holds about a hundred atoms in each MBD cluster; black
phosphorus, at 0.053 atoms/A^3, about 3,500 (580k pairs), and there the sites
are what costs. One CPU node does the 576 sites of the dense cell in 0.83 s
at 96 x 4 (0.91 s at 384 x 1, 0.87 s at 24 x 16; 5.2, 6.7 and 5.9 s before
`9963898`). On a loaded node a thread spends about 0.43 s on a site: 0.25 s
in the energy series, 0.11 s building the cluster, 0.02 s counting its
pairs, 0.03 s setting up and 0.03 s in the force passes; one idle alt core
takes 0.25 s.

What made the difference, in order of effect:

- `cpu_time` was called for every pair to feed timing prints that are
  switched off, and that clock is a system call summing over every thread of
  the process, so threads serialised on it.
- The sites run in threads, SCS included.
- Only the first frequency's dipole matrix is kept: every frequency's matrix
  is A_w = D_w A_1 D_w with per-atom factors, so the series multiplies A_1
  with all 13 frequencies' columns at once, with 13 times less memory
  traffic (the dense cell on a CPU node at 384 x 1, 6.7 -> 2.7 s).
- The force passes' 2b terms are formed per pair from two sums over its nine
  elements, and only for the pairs that move with the atom whose force is
  computed: pairs with an atom within `vdw_loc_rcut`, whose SCS
  polarisability moves, or with that atom itself. The other ~99% of a dense
  cluster's pairs add exactly nothing. One core: 0.38 -> 0.012 s a site.
- Per-pair arrays are kept per thread instead of allocated, cleared and
  freed for every site. glibc returned their pages after each site and the
  next one faulted them in again, and on a loaded node the clearing alone
  took a fifth of a site.
- The matrix is kept as one 3 x 3 block per pair, so a block reads its
  column atom's rows once for three output rows.
- What depends on one atom (its polarisability, its frequency factors, its
  damping radius) is formed once per atom, not once per pair, and the
  never-read 3N x 3N arrays are gone (`b9ecedb`).

**Domains.** Decomposing MBD itself would not make it faster: it is already
split by site, and its communication -- the SCS polarisabilities of every
atom (a reduce and a broadcast) and the Hirshfeld-volume derivatives (an
all-to-all) -- is a handful of collectives a step. What a decomposition would
buy is the rest of the step on domains (`docs/LEGACY_DD.md`, "Not yet done").

**GPU.** In the GPU build the host threads list each site's pairs and solve
its SCS polarisabilities; the device forms every pair's values, as one 3 x 3
block per pair, runs the energy series over the blocks and forms the 2b force
terms, each host thread on a stream and a device allocation of its own
(`vdw_mbd_device`, on by default). A site with fewer than
`vdw_mbd_device_min_pairs` pairs (default 100,000) stays on the host
threads, which are faster for small clusters: the P4 gas has a few thousand
pairs a site and ran slower on a GH200 node's GPUs than on its host threads
(0.19 against 0.13 s). vdW time of one single point on GH200s, one rank per
GPU with 72 threads:

| system | GPUs | MBD on the host threads | on the device |
| --- | ---: | ---: | ---: |
| black phosphorus, 144 atoms | 1 | 0.86 s | 0.36 s |
| black phosphorus, 576 atoms | 4 | 0.88 s | 0.53 s |

For dense systems a GH200 node with MBD on the device is 1.6 times as fast as
a CPU node (0.53 s against 0.83 s). Of the device's work, the series' four
products take 0.36 ms each per dense site on a GH200. Between jobs the same
binary's GH200 times vary by up to a third, so compare builds within one
job.

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
| `tests/mbd_gpu` | MBD on the device against MBD on the host threads of the same GPU build, on the P4 cell (every site sent to the device) and a dense black-phosphorus cell, to the printed precision | alt: RTX 5070 |
| MBD regression and golden cases, P4 cell at 1-12 threads and 4 and 6 ranks | the CPU build against master: printed outputs identical | alt |
| two dense black-phosphorus cells | each MBD commit against the one before it: printed outputs identical, CPU and device | alt |
| `tests/regression --gpu`, the vdW cases | vdw_mbd and vdw_mbd_cell_mpi2 pass; the seven others fail with MBD on the device and on the host alike | alt: RTX 5070 |
