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

**Code.** The speed of the potential: master at `8581780`, 2026-11-20. XRD and
ND from the pair distribution: `e28f1d1`, 2026-11-25, on master since
`038151b`. The Debye and XPS rows: `perf/mad-rank-local` (`b156426`). The MBD
section: `perf/mbd-cpu` (`42f2141`), single points.

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

### XRD and ND from the pair distribution

Time per atom-step, with the factor over plain MD run in the same job and
layout, then seconds per step and the wall time per picosecond and per
nanosecond at the run's time step. Measured at `74d636b`. Steps as above, from
runs of 2 and 42 steps on the CPU node and 2 and 12 on the GPUs (CO), and of
20 and 220 (GST), with each run's Setup time taken off first: it varies by
~0.4 s between runs, which is 2 ms a step over 200 GST steps.

**CO**, the 307,671-atom cell, `dd_type = rcb`, NVE at 1 fs. XRD (`q*F(q)`) from
partial pdfs to 6 A (201 samples, sigma 0.1 A), fitted with `exp_forces`:

| | us / atom-step | s / step | per ps | per ns |
| --- | --- | --- | --- | --- |
| 1 GH200, plain MD | 2.41 | 0.741 | 12 min | 206 h |
| 1 GH200, XRD | 3.83 (1.6x) | 1.18 | 20 min | 328 h |
| 4 GH200, plain MD | 0.83 | 0.255 | 4.2 min | 71 h |
| 4 GH200, XRD | 2.02 (2.4x) | 0.623 | 10 min | 173 h |
| 1 CPU node, 384 x 1, plain MD | 1.56 | 0.481 | 8.0 min | 134 h |
| 1 CPU node, 384 x 1, XRD | 4.52 (2.9x) | 1.39 | 23 min | 386 h |
| 1 CPU node, 48 x 8, plain MD | 1.59 | 0.490 | 8.2 min | 136 h |
| 1 CPU node, 48 x 8, XRD | 3.19 (2.0x) | 0.981 | 16 min | 272 h |

**GST**, 7,200 atoms of Ge2Sb2Te5 from a 1507 K liquid with a three-species
`soap_turbo` potential (2b, core_pot), NPT at 1507 K and 2 fs. ND (`q*F(q)`, 61
points) from partial pdfs to 12 A (12.6 A cutoff, 301 samples, sigma 0.1 A),
fitted with `exp_forces`:

| | us / atom-step | s / step | per ps | per ns |
| --- | --- | --- | --- | --- |
| 1 GH200, plain MD | 1.24 | 0.0089 | 4.5 s | 1.2 h |
| 1 GH200, ND | 2.60 (2.1x) | 0.0187 | 9.4 s | 2.6 h |
| 1 CPU node, 96 x 4, plain MD | 1.74 | 0.0125 | 6.2 s | 1.7 h |
| 1 CPU node, 96 x 4, ND | 2.38 (1.4x) | 0.0171 | 8.6 s | 2.4 h |
| 1 CPU node, 48 x 8, plain MD | 1.99 | 0.0143 | 7.2 s | 2.0 h |
| 1 CPU node, 48 x 8, ND | 3.14 (1.6x) | 0.0226 | 11 s | 3.1 h |

Over the same GST cell a 140,000-step production stage at 96 x 4 and 0.5 fs
ran at 0.0177 s/step, with a frame and the patterns written every 1 ps.

The pattern and its forces themselves (the Exp. pred timer) are now a small
part of a step: 0.09 s on one GH200, 0.03 on four and 0.08-0.11 on a CPU node
for CO, and 4-6 ms for GST. The potential no longer walks the longer list the
pdf cutoff needs (below). Against `e28f1d1` in the same jobs that took 46%
off a GST nd MAD step on one GH200 (0.0344 to 0.0187 s), 29% off its plain MD,
12% off CO XRD on one GH200 and 8% on a CPU node at 48 x 8; GST on a CPU node
moved by less than its run-to-run spread.

What changed for this: on the host the pdf kernel and the forces are threaded,
and each pair's d g / d r is rebuilt in the force pass instead of being stored
as an n_samples x n_pairs array, which every rank allocated and cleared every
step. On the device each pair's kernel is evaluated only where the host
evaluates it, within 9 sigma, for every species pair in one pass and with no
n_samples x n_pairs array, about 150 bytes per pair. The GAP terms walk only
the pairs within their own cutoffs: the list is cut once a step to the largest
GAP cutoff for SOAP, and that again to the 2b, core-potential and 3b cutoffs
for each of those. On the device the 2b and core-pot terms run a thread per
pair, every descriptor in one launch.

**Layout.** Threads now pay on the CPU when a pattern is fitted from a pdf
with a long list: 48 x 8 beats one rank per core for CO XRD by 29%, while
plain MD is fastest at one rank per core. For 7,200-atom GST, 96 x 4 is
fastest. The GH200 is under-filled at 7,200 atoms.

### Debye XRD and XPS

The same 307,671 atoms, microseconds per step per atom, with the factor over
plain MD run in the same job. Measured at `b156426`, before the CPU build's
contract-first SOAP forces, so its CPU plain-MD row is over three times the
current one (above), and the CPU factors are against that:

| | 1 GH200 | 4 GH200 | 1 CPU node | 2 CPU nodes |
| --- | --- | --- | --- | --- |
| plain MD | 2.7 | 0.90 | 4.9 | 2.6 |
| Debye XRD | 11 (4.4x) | 3.9 (4.1x) | 9.2 (1.9x) | 5.4 (2.1x) |
| XPS | 6.6 (2.6x) | 2.6 (2.7x) | 5.3 (1.09x) | 2.9 (1.12x) |

The Debye sum and the XPS spectrum are properties of the whole system,
computed under a decomposition from a gathered copy of the positions or
binding energies. Each rank broadens only its own atoms' XPS lines and sums
only its own atoms' Debye pairs, so both scale like the potential. The GPU
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

**No domain decomposition.** MBD runs replicated under every `dd_type`, and so
does the rest of any step that includes it. Each rank takes a share of the
sites, and each site -- its SCS polarisabilities, then the MBD energy and forces
of its cluster -- is one iteration of a threaded loop. Decomposing MBD itself
would not make it faster: it is already split by site and exchanges only a
handful of collectives a step. A decomposition would put the GAP part of the
step on domains; that is not done yet.

Measured with the `vdw_mbd_cell_mpi2` settings (`vdw_rcut` 25 A, MBD cutoffs
15/10 A, 13 frequencies, order 6): vdW time of one single point at `42f2141`.
A site's cost is set by the atoms within 25 A of it, its cluster, not by the
size of the cell:

- **dilute**: 27 P4 molecules (108 atoms) in a 64.5 A box, and 27 copies of
  that, 2,916 atoms in a 193 A box, too far apart to see each other. Each
  site's cluster is its own 108 atoms.
- **dense**: black phosphorus, 0.053 atoms/A^3. Each site's cluster holds
  3,457 atoms, periodic images included: 580k pairs, about a hundred times a
  dilute site's.

| system | atoms in a site's cluster | where, ranks x threads | vdW time |
| --- | ---: | --- | ---: |
| dilute, 108 atoms | 108 | alt, 1 x 1 / 1 x 4 / 1 x 12 | 0.17 / 0.045 / 0.024 s |
| dilute, 2,916 atoms | 108 | one CPU node, 384 x 1 / 24 x 16 | 0.081 / 0.027 s |
| | | one GH200 node, 4 x 72 | 0.14 s |
| dense, 16 atoms | 3,457 | alt, 1 x 4 | 0.41 s |
| dense, 576 atoms | 3,457 | one CPU node, 384 x 1 / 96 x 4 / 24 x 16 | 0.57 / 0.43 / 0.48 s |
| | | one GH200 node, 4 x 72: on the GPUs / on the host cores | 0.39 / 0.54 s |
| dense, 144 atoms | 3,457 | one GH200, 1 x 72: on the GPU / on the host cores | 0.32 / 0.48 s |

On CPUs use threads: 4 per rank for dense systems, 16 for dilute ones. Each
thread holds about 100 MB on a dense system. On GPUs (`vdw_mbd_device`, on by
default) the host threads list each site's pairs and solve its SCS
polarisabilities, and the device does the rest; a site with fewer than
`vdw_mbd_device_min_pairs` pairs (100,000) stays on the host threads, which
are faster for small clusters, so dilute systems run on the host cores either
way.

What made it fast:

- **Threads**: the sites, SCS included, are a threaded loop, with nothing in
  it that serialises the threads.
- **One dipole matrix**: every frequency's matrix is the first one scaled per
  atom, so the series multiplies the first with all 13 frequencies' columns
  at once.
- **3 x 3 blocks**: the matrix is one block per pair, applied to all 13
  frequencies' columns in one pass, skipping the rows and columns the result
  never reads.
- **Sparse forces**: the 2b force terms are formed only for the pairs that
  move with the atom whose force is computed, about 1% of a dense cluster's.
- **Memory**: per-pair arrays live in per-thread buffers instead of being
  allocated and cleared for every site; per-atom quantities are formed once
  per atom; the unread 3N x 3N arrays are gone.
- **Pair search**: a cluster's pairs are found in one scan, and a pair away
  from the central atom is formed once for both directions.
- **GPU**: the pair blocks, the series and the 2b force terms run on the
  device, each host thread on a stream of its own.

## The checks behind these numbers

What the measured code was tested with, and where:

| check | what it compares | machine |
| --- | --- | --- |
| pdf, sf, xrd and nd on the host, before and after `65737af` and `88658ec` | every output file byte for byte: nd, xrd and sf MAD on the matrix and loop forces, pdf MAD, a histogram pdf and 7,200-atom GST nd MAD at 1 x 1, 1 x 4 and 2 x 2; on Roihu, step 0 of GST MAD at 96 x 4, 48 x 8 and 24 x 16 | alt, Roihu |
| `tests/exp_omp` | pdf, sf, xrd and nd outputs at 1 and 4 threads and 2 ranks x 2, and the xrd and pdf MAD forces | alt |
| the windowed device route (`299685b`, `e28f1d1`) against the host build | forces at step 0 to 1e-8, patterns and thermo to the printed digits: nd and xrd MAD, three batches, two ranks, 7,200-atom GST; the Kokkos build too | alt: RTX 5070 |
| `tests/regression --gpu`, all cases | the same 23 fail on this card with and without the windowed route; the pdf-route cases pass | alt: RTX 5070 |
| the GAP terms on cut lists (`fe81322`, `74d636b`): `tests/regression --cpu`, all 44 cases | every output against a master (`8a7df2b`) build, each run as the binary under test: 107 files identical | alt |
| `tests/gap_neighbors` | the CO cell with the GAP alone, with a pdf to 9 A at 1 and 3 ranks, and with its 3b descriptors alone, against master's output to 1e-6; a GAP list 0.05 A short or a 3b list 5% short fails it | alt |
| the device 2b and core-pot a thread per pair (`74d636b`) against the host build | CO MD and CO and GST nd MAD at 1 and 2 ranks, as close as master's device build is to its host (forces 4e-8); `tests/regression --gpu` fails 22 cases against master's 23; the Kokkos build too | alt: RTX 5070 |
| `tests/xrd_species_order`, absent-species leg | a species listed and absent from the cell: totals, partials and every energy term as when it is not listed, no NaN | alt, host and device |
| `tests/regression --cpu`, all 44 cases | outputs against the previous binary, bit for bit; only the XPS and Debye cases differ, as intended | alt: x86, 12 cores |
| `tests/regression --gpu`, all 41 cases with data | device against host within tolerance; the same 36 agree as before | Roihu GH200 |
| `tests/soap_descriptor_gpu` | the device descriptor kernels bit for bit against the originals, and the force contraction against the explicit derivative chain (to 4e-15) | alt: RTX 5070 |
| `tests/xps_fd` | XPS forces and virial against finite differences of the XPS energy, serial and on two ranks | alt |
| `tests/xrd_debye` | the Debye pattern against an independent numpy sum, forces and virial against finite differences | alt |
| XPS and Debye on 3 ranks under `dd_type` none, rcb and legacy | forces and energies against the binary before `b156426`: identical | alt |
| `tests/pdf_virial`, `tests/xrd_lp_pdf` (gradient) | pdf, structure factor and XRD forces and virial against finite differences, with the device build | alt: RTX 5070 |
| `tests/mbd_gpu` | MBD on the device against MBD on the host threads of the same GPU build, on the dilute P4 cell (every site sent to the device) and a dense black-phosphorus cell, to the printed precision | alt: RTX 5070 |
| MBD regression and golden cases, P4 cell at 1-12 threads and 4 and 6 ranks | printed outputs against the MBD before these changes: identical, but for the per-atom MBD scaling of the two MD goldens, `vdw_tsmbd_md` and `vdw_tsmbd_md_mpi2`, which moves by up to 2e-11 after four steps; each commit from `496dd81` against `b514ca9`: identical | alt |
| two dense black-phosphorus cells | each MBD commit against the one before it: printed outputs identical, CPU and device | alt |
| 576-atom black phosphorus at 384 x 1, 96 x 4, 24 x 16 | `42f2141` against `0548fb4`, in one job: identical | Roihu CPU node |
| `tests/regression --gpu`, the vdW cases | vdw_mbd and vdw_mbd_cell_mpi2 pass; the seven others fail with MBD on the device and on the host alike | alt: RTX 5070 |
