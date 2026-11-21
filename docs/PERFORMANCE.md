# Performance: speed, sizing and where the time goes

How fast TurboGAP runs, how many atoms to give each CPU core or GPU, what
experimental biasing (MAD) costs on top of the potential, and what a profile
of each build shows. For how to launch a run, see the README's "Running
TurboGAP".

## What was measured, and on what

**System.** Amorphous oxygenated carbon, the `CO` set of
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

**Code.** The speed of the potential: master at `8581780`, 2026-11-20. The
MAD tables: branch `perf/soap-der-one` as of 2026-10-31, and the Debye and
XPS rows `perf/mad-rank-local` (`b156426`). The MBD section: `perf/mbd-cpu`
(`42f2141`), single points.

## Speed of the potential

Microseconds of wall time per MD step per atom. In brackets, that times the
GPUs or cores used: the cost in device- or core-microseconds per atom-step,
which stays constant when a run scales perfectly.

| atoms | 1 GH200 | 4 GH200 | 1 CPU node (384) | 2 CPU nodes (768) |
| ---: | --- | --- | --- | --- |
| 7,176 | 5.5 | 4.1 (17) | 3.4 (1,310) | 3.9 (2,960) |
| 24,219 | 3.5 | 2.5 (10) | 1.9 (730) | 2.1 (1,640) |
| 57,408 | 3.1 | 1.6 (6.3) | 1.6 (630) | 1.2 (890) |
| 307,671 | 2.8 | 0.95 (3.8) | 1.44 (550) | 0.79 (610) |
| 615,342 | 2.7 | 0.87 (3.5) | 1.43 (550) | 0.74 (570) |
| 1,230,684 | 2.7 | 0.83 (3.3) | 1.41 (540) | 0.73 (560) |

On the large cells a CPU node does the work of about 1.9 GH200s, and a GH200
node, four GPUs, that of about 1.7 CPU nodes. The two largest cells on more
CPU nodes:

| atoms | 1 CPU node | 2 CPU nodes | 4 CPU nodes | 6 CPU nodes |
| ---: | --- | --- | --- | --- |
| 307,671 | 1.44 (550) | 0.79 (610) | 0.43 (660) | 0.35 (810) |
| 1,230,684 | 1.41 (540) | 0.73 (560) | 0.37 (570) | 0.26 (590) |

From one node to six the largest cell keeps 92% of its per-core speed and the
307,671-atom cell 68%: on six nodes that is 134 atoms per core, where the
fixed costs of a step weigh in again (below).

**Atoms per GPU.** One GPU is within 15% of its best, 2.7 microseconds per
atom-step, from about 50k atoms and at it from about 300k; below that the
fixed costs of a step dominate, and at 7k atoms each atom costs twice as much.
Over four GPUs the efficiency against one is 74% at 77k atoms per GPU, 77% at
154k and 81% at 308k. Give each GPU 100k atoms or more.

**Atoms per CPU core.** The cost per core is within 15% of its best from about
150 atoms per core and within 5% from 800, and going from one node to two
loses 3-9% from 300k atoms up (30% at 57k). Small systems pay more than they
did: the SOAP forces became three times cheaper, the fixed costs of a step did
not, and at 19 atoms per core a core now runs at 41% of its best. A few
hundred atoms per core or more is the comfortable range. Threads do not help
on the CPU: 16 per rank reach about 65% of the speed of 16 ranks, so run one
rank per core -- except for MBD (below).

## What experimental biasing costs

The same 307,671 atoms, microseconds per step per atom, with the factor over
plain MD run in the same job. Measured at `b156426`, before the CPU build's
contract-first SOAP forces, so its CPU plain-MD row is over three times the
current one (above), and the CPU factors are against that:

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
frequencies, order 6), vdW time of one single point.

A site's cost is set by the atoms within 25 A of it, its cluster, not by the
size of the cell, so the table keeps two kinds of system apart:

- **dilute**: `p4_cell.xyz` is one cluster of 27 P4 molecules (108 atoms) in
  a 64.5 A box, and the 2,916-atom system is 27 copies of it in a 193 A box,
  too far apart to see each other. Every site's cluster is its own 108 atoms,
  with on average 57 of them within 10 A.
- **dense**: black phosphorus, 0.053 atoms/A^3. Every site's cluster holds
  3,457 atoms, periodic images included, each with about 210 others within
  10 A: 580k pairs a site, roughly a hundred times a dilute site's.

| system | atoms in a site's cluster | where | master | now | speed-up |
| --- | ---: | --- | ---: | ---: | ---: |
| dilute: one P4 cluster, 108 atoms | 108 | alt, 1 rank x 1 / 4 / 12 threads | 6.1 s | 0.17 / 0.045 / 0.024 s | 37 / 136 / 254 |
| dilute: 27 P4 clusters, 2,916 atoms | 108 | one CPU node, 384 x 1 | 0.68 s | 0.081 s | 8.4 |
| | | one CPU node, 24 x 16 | 7.6 s | 0.027 s | 280 |
| | | one GH200 node, 4 x 72 | 67 s | 0.17 s | 400 |
| dense: black phosphorus, 16 atoms | 3,457 | alt, 1 rank x 4 threads (master: 1) | 217 s | 0.41 s | 540 |
| dense: black phosphorus, 576 atoms | 3,457 | one CPU node, 96 x 4 (master: 96 x 1, the most ranks its memory allows) | 62 s | 0.47 s | 130 |
| | | one CPU node, 384 x 1 | out of memory | 0.55 s | |
| | | one GH200 node, 4 x 72, MBD on the GPUs | -- | 0.39 s | |
| | | one GH200 node, 4 x 72, MBD on the host cores | -- | 0.54 s | |
| dense: black phosphorus, 144 atoms | 3,457 | one GH200, 1 x 72, MBD on the GPU | -- | 0.32 s | |

Master's MBD is not threaded, and on the dense cell each of its ranks peaks
at 5.5 GB (several 3N x 3N matrices, 0.86 GB each for a cluster of 3,500
atoms): at 384 x 1 it ran a CPU node out of memory, so its row is 96 ranks of
one thread with four cores each, the most that fit, run in the same job as
the current code's two dense rows (on other nodes the 96 x 4 run took 0.43
s). Master's dense cells were not run on a GH200 node; its dilute run took 67 s
there.

Through `f28397a` the CPU build's output was master's byte for byte.
`62ea743`, `950b66a` and `0548fb4` sum the host's default path in their own
order (below): the printed outputs of the MBD regression cases and of five of
the seven vdW golden cases, at 1, 2, 4 and 12 threads and on 4 and 6 ranks,
are still master's; the two MD goldens, `vdw_tsmbd_md` and
`vdw_tsmbd_md_mpi2`, move by up to 2e-11 in the per-atom MBD scaling after
four steps (their first frame is master's). Every commit since prints what
`0548fb4` printed, byte for byte, on alt and on a Roihu CPU node. Against
master, both CPU-node cells' forces differ by at most 1e-8 eV/A, the last
printed digit, and their vdW energies print the same. MBD on the device and
on the host threads print the same digits on every cell above.

In the dense cells the sites are what costs. One CPU node does the 576 sites
of the dense cell in 0.43 s at 96 x 4 (0.57 s at 384 x 1, 0.48 s at 24 x 16;
5.2, 6.7 and 5.9 s before `9963898`, 0.83, 0.85 and 0.87 s at `b514ca9`). On a
loaded node a thread spends about 0.22 s on a site: 0.09 s in the energy
series, 0.09 s building the cluster's pair blocks, 0.02 s counting the pairs,
0.01 s setting up and 0.02 s in the force passes; one idle alt core takes
0.074 s (0.25 s at `b514ca9`).

At 96 x 4 that cell gives each thread 1.5 sites, so most sites are their
thread's first and fault its per-thread buffers in (about 100 MB in dense
black phosphorus). A 2,304-atom cell of the same density, 6 sites a thread,
measured at `d0ad564`: setting up 0.031 -> 0.008 s a site and the build
0.091 -> 0.065 s. The buffers live for one call, so every MD step pays this
again; keeping them across calls is not done yet.

What made the difference, in order of effect:

- `cpu_time` was called for every pair to feed timing prints that are
  switched off, and that clock is a system call summing over every thread of
  the process, so threads serialised on it.
- The sites run in threads, SCS included.
- Only the first frequency's dipole matrix is kept: every frequency's matrix
  is A_w = D_w A_1 D_w with per-atom factors, so the series multiplies A_1
  with all 13 frequencies' columns at once, with 13 times less memory
  traffic (the dense cell on a CPU node at 384 x 1, 6.7 -> 2.7 s).
- The energy series applies each pair's 3 x 3 block to an atom's three
  columns in one pass over the 13 frequencies' columns, packed 40 wide (39
  padded to a multiple of 8), instead of nine 39-long updates. The first
  power takes only the central atom's partners, where the first column is
  nonzero, and the last power only their rows, the only rows the integrand
  and the forces read: about a third fewer block products. Loaded node:
  0.24 -> 0.09 s a site (`496dd81`).
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
- The cluster's pairs are found by one scan of each pair, not two scans of
  every atom against every other, and the build walks per-row partner
  lists; a pair with neither atom near the central one has its block formed
  once for both directions (`edea662`). One core: 0.069 -> 0.025 s a site
  counting and building.
- G_sym, which the force passes dot with the series' result, is nonzero only
  on the central atom's partners (a fifth of the cluster) and is now formed
  there only; AT_sym and pol_sym live in per-thread buffers (`d0ad564`), and
  seven per-pair arrays that only ever held 1 or 0 are constants
  (`42f2141`). Loaded node: force passes 0.026 -> 0.016 s and setting up
  0.031 -> 0.011 s a site.

**Domains.** Decomposing MBD itself would not make it faster: it is already
split by site, and its communication -- the SCS polarisabilities of every
atom (a reduce and a broadcast) and the Hirshfeld-volume derivatives (an
all-to-all) -- is a handful of collectives a step. What a decomposition would
buy is the rest of the step on domains.

**GPU.** In the GPU build the host threads list each site's pairs and solve
its SCS polarisabilities; the device forms every pair's values, as one 3 x 3
block per pair, runs the energy series over the blocks and forms the 2b force
terms, each host thread on a stream and a device allocation of its own
(`vdw_mbd_device`, on by default). A site with fewer than
`vdw_mbd_device_min_pairs` pairs (default 100,000) stays on the host threads,
which are faster for small clusters: the dilute P4 clusters have a few
thousand pairs a site and ran slower on a GH200 node's GPUs than on its host
threads (0.19 against 0.13 s). vdW time of one single point on GH200s, one
rank per GPU with 72 threads:

| system | atoms in a site's cluster | GH200s | MBD on the host threads | on the device | at `b514ca9`: host / device |
| --- | ---: | ---: | ---: | ---: | ---: |
| dense: black phosphorus, 144 atoms | 3,457 | 1 | 0.48 s | 0.32 s | 0.86 / 0.36 s |
| dense: black phosphorus, 576 atoms | 3,457 | 4 | 0.54 s | 0.39 s | 0.88 / 0.53 s |
| dilute: 27 P4 clusters, 2,916 atoms | 108 | 4 | 0.14 s | 0.17 s | 0.20 / 0.18 s |

The host and device columns (`42f2141`) ran in one job and the `b514ca9`
column in another (the same code at `9a871ab`); run in one job, the
dense 576-atom cell on four GH200s went from 0.87 to 0.54 s on the host
threads and from 0.50 to 0.40 s on the device. The dilute clusters have far
fewer pairs a site than `vdw_mbd_device_min_pairs`, so with the device on
their sites still run on the host threads: the device column there differs
from the host column only by run-to-run noise.

For dense systems a GH200 node with MBD on the device now runs about as fast
as a CPU node (0.39 s against 0.43 s), and its host cores alone are not far
behind (0.54 s). The host's share of a device site -- counting and listing
its pairs, setting up -- went from 0.055 to 0.016 s, and the device, fed by
72 threads per GPU, is what limits it now: the series and force phases,
mostly device work, take 0.042 s a site with the old host code and 0.058 s
with the new, as more threads queue on it. Splitting the sites between the
device and the host threads would use both; that is not done yet. Of the
device's work, the series' four products take 0.36 ms each per dense site
on a GH200. Between jobs the same binary's GH200 times vary by up to a
third, so compare builds within one job.

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
| `tests/mbd_gpu` | MBD on the device against MBD on the host threads of the same GPU build, on the dilute P4 cell (every site sent to the device) and a dense black-phosphorus cell, to the printed precision | alt: RTX 5070 |
| MBD regression and golden cases, P4 cell at 1-12 threads and 4 and 6 ranks | the CPU build against master: printed outputs identical but for the two MD goldens (above); each commit from `496dd81` against `b514ca9`: identical | alt |
| two dense black-phosphorus cells | each MBD commit against the one before it: printed outputs identical, CPU and device | alt |
| 576-atom black phosphorus at 384 x 1, 96 x 4, 24 x 16 | `42f2141` against `0548fb4`, in one job: identical | Roihu CPU node |
| `tests/regression --gpu`, the vdW cases | vdw_mbd and vdw_mbd_cell_mpi2 pass; the seven others fail with MBD on the device and on the host alike | alt: RTX 5070 |
