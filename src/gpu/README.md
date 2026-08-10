# `src/gpu` — the device backend

Everything TurboGAP runs on a GPU lives here. It used to live in three files in
`src/`: `cuda_wrappers.cu` (2684 lines), `gpu_exp.cu` (1709) and `3b_final.cc`
(886). The split is by **what is computed**, in three groups.

| prefix | what it is | who calls it |
| --- | --- | --- |
| `gpu_*` | infrastructure — memory, BLAS, scan primitives | both groups below |
| `gap_*` | the interatomic potential: 2b, 3b, soap_turbo | `gap_backend_gpu.f90`, `gap_interface.f90` |
| `mad_*` | the experimental-data path: pdf, structure factor, electrostatics | `exp_interface.f90`, `turbogap_exp.f90`, `electrostatics.f90` |

## The files

**Infrastructure**

- `gpu_common.h` — the HIP/CUDA includes, `gpuErrchk`, `TG_CHECK_LAUNCH`.
  Included by every `.cu` here and nothing else.
- `gpu_memory.cu` / `.h` — the device memory ledger and its OOM report, the
  `hipMallocAsync`/free/memcpy wrappers, the pinned host allocator, device
  enumeration and selection, `gpu_meminfo`, the debug pointer printers.
- `gpu_blas.cu` / `.h` — hipBLAS handle lifetime and the dense products
  (`mmul_t_n`, `mvmul_n`, `mmul_n_t`, `dgemm_n_n`).
- `gpu_scan.cu` / `.h` — block reduction, recursive reduction and inclusive
  scan over per-pair flag arrays. Used by both `mad_pdf` and
  `mad_electrostatics` to compact their neighbour lists.

**GAP** (declarations in `gap_gpu.h`)

- `gap_predict.cu` — kernel power `k^zeta`, the `e0` offset, and the two
  mat-vec products over the sparse set.
- `gap_2b.cu` — the two-body descriptor and the core repulsive potential; both
  evaluate a cubic spline, so `gpu_spline`/`gpu_spline_der` live here.
- `gap_3b.cc` — the three-body descriptor. **Stays `.cc`** so the Makefile
  compiles it with `$(CC)`, which passes `-std=c++20`; it needs `<numbers>` and
  `<bit>`.
- `gap_soap_radial.cu` — the poly3 and poly3gauss radial expansions and their
  derivatives, plus the global scaling. Two ~200-line kernels; on a GH200 this
  is where `ncu` measures 5.5% achieved occupancy.
- `gap_soap_angular.cu` — associated Legendre arrays, the `exp(i m phi)`
  coefficients and their derivatives, and the `cnk` contraction.
- `gap_soap_descriptor.cu` — the power spectrum from `cnk`, its normalisation,
  the radial/azimuthal/polar derivatives and their conversion to Cartesian, and
  the transposes those need.
- `gap_soap_forces.cu` — contraction of the Cartesian derivatives into forces
  and the virial, and the same for local properties.

**MAD** (declarations in `mad_gpu.h`)

- `mad_pdf.cu` — the pair distribution function: pair counting inside the pdf
  cutoff, compaction into `k`-indexed buffers, the kernel density estimate of
  `g(r)` and its derivative, and the reduction over pairs.
- `mad_xrd.cu` — the structure factor: `Gk` and its position derivative, the
  Hadamard product with the scattering factors, the dgemv that turns that into
  per-pair forces, and their collection into forces and the virial.
- `mad_electrostatics.cu` — shifted-force Coulomb (GSF). Grouped here rather
  than with GAP because it is reached from the same experimental-prediction
  driver and shares the pair compaction, but it is a separate physical model
  from pdf and xrd.

`orthonormalization_kernels.cc` is disabled in the Makefile, not dead; it is
kept for possible reintegration.

## Things the split had to decide, and why

**Kernel launches cannot cross a translation unit.** Without `-rdc=true` a
`<<<>>>` needs the kernel's definition in the same object file. Only one kernel
was launched from two places that are now separate files —
`kernel_multiply_flags`, from the pdf and electrostatics pair counters — so it
gained a host launcher, `gpu_multiply_flags` in `gpu_scan.cu`, computing the
geometry both call sites computed inline. Every `__device__` helper turned out
to be used by exactly one of the new files and simply moved with its callers.
If a future kernel is needed from two files, add a launcher; do not turn on
`-rdc=true` for it, which costs register pressure everywhere.

**`tpb` is per-file and deliberately not unified.** `cuda_wrappers.cu` used 64
and `gpu_exp.cu` used 512, and every kernel was tuned against its own file's
value. Each new file carries the value of the file its code came from. A single
shared `tpb` would silently retune a dozen kernels at once.

**The two `gpuAssert`s had drifted** — the `cuda_wrappers.cu` copy dumped a host
backtrace, the `gpu_exp.cu` copy did not. `gpu_common.h` keeps the backtrace
version, because behind a thin C shim the wrapper's own file and line name the
shim, and only the backtrace lets `addr2line` name the Fortran call site.

**The headers are generated from the definitions**, by `tools/split_gpu.py`,
not typed. Their only job is to fail the build when a signature drifts, which a
hand-written copy would not reliably do.

Dropped in the move, both dead: a file-scope `int counter = 0` that no
uncommented line read, and a commented-out `atomicDoubleAdd` superseded by the
native double `atomicAdd`.

## Adding a kernel

1. Put it in the file named after the physics. If it needs a new group, add a
   `gap_`/`mad_`-prefixed file and a line in `SRC_CUDA` in the Makefile.
2. If it is reachable from Fortran, add the `extern "C"` entry point and
   declare it in `gap_gpu.h` or `mad_gpu.h`, matching the definition exactly.
3. The Makefile lists `src/gpu/*.h` as a prerequisite of every object here, so
   a header edit rebuilds all of them. Nothing else needs touching.
