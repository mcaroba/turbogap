# Profiling the GPU build

Companion to `BUILD_AND_TEST.md`. Everything here is copy-pasteable on `alt`.

```sh
export HOP_ROOT=/u/74/zarrout1/unix/work/hop
export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda
export PATH="$HOME/.local/bin:$PATH"          # uv
```

---

## 1. Once per machine

```sh
tools/setup_profiling_env.sh            # uv venv + analysis packages, and a report
tools/setup_profiling_env.sh --check    # verify only; exits non-zero on a MISS
```

It checks the things that actually stop a profile and names the fix for each.
Two of its findings on `alt` are worth knowing in advance:

* **`nvidia-smi` fails and this does not matter.** `Failed to initialize NVML:
  Driver/library version mismatch` is a userspace/kernel NVML version skew. The
  CUDA runtime is unaffected, and neither nsys nor ncu uses NVML. The check
  script therefore probes the device by compiling a `cudaGetDeviceCount` rather
  than by asking `nvidia-smi`, because taking `nvidia-smi`'s word for it would
  mean refusing to profile a machine that profiles fine.

* **`ncu` cannot run as a normal user here.** `/proc/driver/nvidia/params` has
  `RmProfilingAdminOnly: 1`, so any counter collection fails with
  `ERR_NVGPUCTRPERM`. Fixing it needs root and a module reload:

  ```sh
  echo 'options nvidia NVreg_RestrictProfilingToAdminUsers=0' \
    | sudo tee /etc/modprobe.d/nvidia-profiling.conf
  ```

  Until then `--tool nsys` is the whole workflow. That is less of a loss than it
  sounds: everything in section 4 came out of nsys.

---

## 2. Profile a case

```sh
tools/profile_gpu.sh --list
tools/profile_gpu.sh CO_predict
tools/profile_gpu.sh --no-build XRD_mad estat_gsf
tools/profile_gpu.sh --openmp 4 --input 'gpu_n_batches = 4' --suffix .omp XRD_mad
```

The cases are the regression suite's, read from `tests/gpu/cases.sh` by both
scripts. There is no second copy of an input to drift.

Results land in `profiling/<case>/`:

| | |
|---|---|
| `run.log` | TurboGAP's own output, including its timer table |
| `<case>.nsys-rep` | the timeline; open in the Nsight Systems GUI |
| `<case>.sqlite` | the same data, queryable |
| `stats/*.csv` | nsys's own summaries, including `nvtx_gpu_proj_sum` |
| `report.txt` | what `tools/nsys_report.py` prints |

Re-run the analysis on a capture without re-running TurboGAP:

```sh
python3 tools/nsys_report.py profiling/CO_predict/CO_predict.sqlite --top 30
python3 tools/nsys_report.py profiling/CO_predict/CO_predict.sqlite --plot phases.png
```

### Options that change what you measure

| | |
|---|---|
| `--openmp [N]` | build `OPENMP=1`, run with N host threads |
| `--input 'k = v'` | override one input keyword (repeatable) |
| `--suffix S` | write to `<case>S/`, to compare two variants |
| `--tool ncu` | per-kernel counters (needs the permission fix above) |
| `--ranks N` | run under `mpirun -np N` |

**`--openmp` alone usually measures nothing.** `gpu_n_batches` defaults to 1 and
`gpu_context_init` clamps `n_omp` to it, so one batch means one stream with or
without threads. Pair it with `--input 'gpu_n_batches = 4'`.

### Several ranks on one GPU

This works as of 2026-08-09 and is usually the better lever than `--openmp`:

```sh
tools/run_rank_sweep.sh --atoms atoms_124959.xyz --ranks '1 2 4 6'
tools/run_rank_sweep.sh --mps --ranks '1 4'      # kernels genuinely concurrent
TURBOGAP_MPI_RANKS=2 tests/gpu/run_regression.sh # the suite, multi-rank
```

Measured on 125k atoms: ~1.16x at two ranks, ~1.22x at four with MPS, and flat
after that -- the run is ~98% GPU-bound, so there is little host work to overlap.
**The gain is memory**: device peak per rank falls 6.7x from one rank to six.
A workload with a big cutoff, and so an expensive neighbour build, would gain
much more in time as well.

---

## 3. What the build does

**Before any measurement: check the object tree is not stale.**

```sh
ls -la build/*.o | head          # all from the same build, or you are measuring nothing
```

`make DEBUG=0` recompiles nothing when no source changed, and the architecture
makefiles default `DEBUG ?= 1`. A tree built once with the default and then
built with `DEBUG=0` links the debug objects and reports success. On
2026-08-09 that made the three-body kernel 4.7x slower than the same source
rebuilt clean, and shifted the SOAP energy by 1e-10 relative -- which looks
exactly like a numerical regression from whatever you changed last.

`DEBUG=1` now gets its own tree (`build-dbg/`, `bin-dbg/`), so this cannot
recur. `make deepclean` is the cure for a tree that predates the fix.


```sh
make PROFILE=1 DEBUG=0        # -> build-profile/  bin-profile/turbogap
make OPENMP=1                 # -> build-omp/      bin-omp/turbogap
make OPENMP=1 PROFILE=1 DEBUG=0
```

* `PROFILE=1` adds `-lineinfo` and host `-g`, and defines `_NVTX` so the phase
  markers in `src/nvtx.f90` compile in. It **refuses `DEBUG=1`**: `-G` disables
  device optimisation (2.1x, per `BUILD_AND_TEST.md`) and reorders the kernel
  ranking, so a DEBUG profile does not describe the build you run.
* `OPENMP=1` adds `-fopenmp`. Off by default -- see section 5.
* Each combination gets **its own object tree**. make rebuilds on timestamps,
  not on flags: without the split, `make PROFILE=1` over an existing `build/`
  relinks the old objects into a binary with no line info and no NVTX ranges,
  and reports success. `profile_gpu.sh` additionally checks the binary for an
  NVTX symbol before trusting a single number in the report.

### The phase markers

`time_start`/`time_end` take an optional label. Given one, they also push and
pop an NVTX range, so nsys shows TurboGAP's phases on the timeline beside the
kernels they launched. Labels are attached to the timing buckets rather than
placed independently because `times_t` already names every phase and already
has exactly one start/end pair for each; a second set of markers would be a
second thing to keep in sync.

**A label must be given to both calls of a pair or to neither.** NVTX ranges are
a stack; labelling only the start leaks a range and nests every later phase
inside it.

---

## 4. What the report answers that `nsys stats` does not

All four are questions about time when *nothing* is happening.

**Device occupancy.** Busy time is the measure of the *union* of the kernel and
memcpy intervals, never the sum -- summing double-counts concurrent streams and
can exceed the wall clock. Idle is the budget available for overlap: a run that
is mostly idle is not short of FLOPs.

**Stream concurrency.** Peak concurrent operations has a yes/no answer in the
data. If it is 1, nothing on the device ever overlapped, however many streams
were created.

**Per-phase GPU busy.** `gpu%` per NVTX range is the fraction of that phase
during which the device was doing anything. A phase at 40% is host-bound or
launch-bound and tuning its kernels is wasted effort; a phase at 99% is not
going to be helped by another stream.

**Gaps.** Every interval between GPU operations, attributed to the phase it
falls in.

---

## 5. Measured state of the tree (2026-08-08, RTX A2000, 26 SMs)

Captured with the scripts above; `report.txt` in each `profiling/<case>/`
directory reproduces them.

### The single biggest cost was an allocator barrier

`gpu_malloc_async` ended in `hipDeviceSynchronize()` followed by a
`hipGetLastError()` whose result was assigned to a local and discarded, with the
test commented out. Every device buffer in TurboGAP is allocated through it.

| case | `cudaDeviceSynchronize` calls | % of span |
|---|---|---|
| `CO_predict` | 3982 | 60.5% |
| `estat_gsf` | 1025 | 49.7% |
| `XRD_mad` | 1838 | 31.0% |
| `vdw_ts` | 8929 | 23.8% |

Removing it: `CO_predict` 3.888 s -> **3.094 s of GPU span (-20%)**, with GPU
busy time unchanged at 2.03 s -- so it was pure host-side stall. It was also a
correctness ceiling for the stream work: `hipDeviceSynchronize` waits on every
stream, so no two streams could overlap while it was there.

### Streams are dormant, in three independent ways

Peak device concurrency was **1 on all four cases**, including the two that
create more than one stream. Three separate reasons, all of which had to be
fixed before the fourth question ("would streams help?") could even be asked:

1. **No architecture makefile has ever passed `-fopenmp`**, so `_OPENMP` was
   undefined, `gpu_context_init`'s thread-count block was dead, `n_omp` came out
   1, and the `!$OMP PARALLEL DO` over the batched loops was a comment.
2. **The pdf loop computed a stream index and then assigned `omp_task = 1` over
   the top**, unconditionally.
3. **The xrd loop used `mod(i-1, n_omp)+1`** -- a property of the iteration, not
   of the thread. Under `static` scheduling two concurrent threads can draw the
   same index, which serialises them against each other and lets one thread's
   `hipFreeAsync` land on a buffer the other is still using.

`gpu_omp_task()` in `gpu_context.f90` is now the single answer, and it is
`omp_get_thread_num() + 1`.

### Where the idle actually is

`gpu%` per phase, after the allocator fix:

| case | phase | wall | gpu% |
|---|---|---|---|
| `CO_predict` | `gap:soap` | 2.696 s (87% of span) | **62%** |
| `CO_predict` | `2b` | 0.239 s | 93% |
| `CO_predict` | `3b` | 0.159 s | 85% |
| `estat_gsf` | `neigh` | **4.620 s** | **0%** |
| `estat_gsf` | `gap:soap` | 3.193 s | 42% |
| `XRD_mad` | `exp_batched` | 0.281 s | 94% |

Two things follow, and they point in opposite directions from the obvious plan:

* **The batched loops that already have the stream machinery are the ones that
  least need it.** `exp_batched` is 94% device-busy. Even perfect overlap there
  is worth under 7% of that phase.
* **The idle is in `gap:soap` and in the neighbour build.** `estat_gsf` spends
  4.62 s in `neigh` with the GPU completely idle -- longer than its entire GPU
  span.

### Launch geometry, from `cuda_gpu_kern_gb_sum`

Occupancy cannot be measured directly here (ncu is blocked), but the launch
geometry bounds it. On 26 SMs:

| kernel | grid x block | threads | note |
|---|---|---|---|
| `kernel_get_2b` | 113 x 64 | 7,232 | 11% of `CO_predict`; badly underfills |
| `cuda_get_soap_der_one` | 78 x 128 | 9,984 | 8% of `CO_predict` |
| `kernel_get_radial_poly3gauss` | 155 x 64 | 9,920 | |
| `kernel_2nd_try` (3b) | 7176 x 32 | 229,632 | fills the grid, but 1 warp/block |
| `cuda_get_soap_der_thr_one` | ~9900 x 64 | ~634,000 | saturates |
| `cuda_get_derivatives_new_new` | 155x16x45 x 64 | plenty | saturates |

So "GPU busy" overstates utilisation for `kernel_get_2b` in particular: it
occupies the timeline for 74 ms per launch while offering the device 7,232
threads.

`kernel_get_2b` also recomputes its own neighbour offset by walking the whole
prefix -- `for (i = i_beg-1; i < i_site; i++) k += n_neigh_d[i];` -- so each
thread does up to `n_sites` extra loads. The 3b path already precomputes exactly
this array (`kappas`). That is a larger and simpler win than streams for this
kernel.

---

## 6. Device memory, and large systems

### The ledger

Every wrapper allocation and free goes through a ledger in `cuda_wrappers.cu`,
so at any point TurboGAP can say what **it** is holding, next to what the
**device** reports free. The two answer different questions and you need both:
`hipMemGetInfo` includes every other process on the card, and by the time an
allocation fails the interesting quantity is the high-water mark, not the
current one.

Every run now ends with:

```
GPUmem [end of run]
GPUmem   device         3.158 GB free of    5.653 GB total (44.1% in use)
GPUmem   this rank      0.000 GB live in 0 buffers, peak    1.587 GB
GPUmem   cumulative    41.203 GB requested over 9558 allocations, 9558 frees (0 untracked)
```

`peak` is the number that sizes the next run. `live` at the end should be zero —
anything else is a leak, and the buffer count says how many.

Call `gpu_memory_report("label")` from anywhere in the Fortran to get the same
block mid-run.

### When an allocation fails

It no longer aborts with a bare "out of memory". Two things happen first:

1. **One retry, after draining the stream and trimming the pool.**
   `hipFreeAsync` is stream-ordered: it returns a buffer to the pool when the
   stream *reaches* the free, not when the call is made, so a loop that
   allocates and frees per batch is always holding memory that is logically free
   and physically not. Draining converts it. This is on the failure path only —
   the normal path is one `hipMallocAsync` and nothing else.

2. **A report that names the fix**, including how much was wanted, how much the
   device had, what this rank was holding, and the current value of the keyword
   that controls it.

### Sizing from the device: `gpu_mem_fraction`

```
gpu_mem_fraction = 0.8
```

At start-up TurboGAP reads the device's free memory and sets
`max_Gbytes_per_process` to `fraction × free / ranks_per_device`. That keyword is
what splits the SOAP descriptor loop, and its default of **1.0 GB was chosen
with no card in mind** — on a 5.7 GB A2000 it asks for roughly four times more
batches than needed, and on an 80 GB card eighty times.

Default is 0 (off), so nothing about an existing input changes until you set it.

Divide-by-ranks matters: `cuda_set_device` hands out cards round-robin, so when
there are more ranks than devices they share one card's memory, and each one
budgeting the whole thing is a guaranteed failure that looks like a code bug.

### `gpu_n_batches`

Auto-sized for the **pdf/xrd** path only, from
`estimate_max_exp_forces_device_memory_usage` — which enumerates those buffers
one by one and whose call site had been commented out, so nothing had ever
compared `gpu_n_batches` to the card. Electrostatics stays manual: there is no
validated estimate for it, and a wrong automatic one is worse than none.

**A value in the input is a floor, not a suggestion.** The automatic sizing can
raise it and never lowers it. Someone who set it because the estimate was wrong
for their system is right, and silently overriding them would be worse than not
estimating at all.

### The scaling ladder

`../large_systems/diamond_1M/` holds a diamond ladder from 13,824 to 1,000,000
atoms, cut from the MAD weak-scaling snapshot with `tools/make_subcell.py`.
See its `README.md`. They are profiling cases, not tests — there is no reference
to compare a half-million-atom result against.

```sh
tools/run_scaling_ladder.sh                # every rung, stop at the first failure
tools/run_scaling_ladder.sh --max 262202   # up to a size
tools/run_scaling_ladder.sh --keep-going
tools/profile_gpu.sh diamond_125k          # one rung, under nsys
```

The ladder table reports wall time, **peak host RSS** and **peak device
memory** per size. Both matter and they fail differently: the host runs out
during the neighbour build, the device during the descriptor loop.

---

## 7. Reading the numbers honestly

The GPU span is measured from the first launch to the last completion, so it
excludes start-up and any trailing host work. A phase can therefore show more
than 100% of the span (`estat_gsf`'s `neigh` does) -- that means it began before
the first kernel, not that the arithmetic is wrong.

`--tool ncu` serialises the context and replays each kernel once per metric
pass. `--launch-count` defaults to 3 for that reason: on `CO_predict`, which
launches thousands of kernels, an unbounded ncu run does not finish in a useful
time, and profiling a few launches of each distinct kernel gives the same
per-kernel numbers.
