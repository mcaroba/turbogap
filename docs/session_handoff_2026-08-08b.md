# Session handoff — 2026-08-08 (b): profiling, and what it found

Second session of 2026-08-08. The first cleaned up the driver's merge surface;
this one is about performance. It builds the profiling workflow the tree did not
have, uses it, and acts on the two findings that were both cheap and certain.

| | |
|---|---|
| `CO_predict` GPU span | 3.888 s → **3.094 s** (−20%) |
| GPU idle on `CO_predict` | 47.7% → 34.3% |
| `cudaDeviceSynchronize` calls, `CO_predict` | 3982 → **8** |
| CPU suite | unchanged (not touched) |
| GPU suite | 3 passed, 0 failed, 2 xfail — the documented baseline |

New: `docs/PROFILING.md` is the workflow. This file is why.

---

## 1. What was built

| | |
|---|---|
| `tools/setup_profiling_env.sh` | uv-managed venv + a report on what would stop a profile |
| `tools/profile_gpu.sh` | build `PROFILE=1`, run a suite case under nsys/ncu, analyse |
| `tools/nsys_report.py` | the analysis `nsys stats` does not do |
| `src/nvtx.f90` | phase markers, driven from the timing buckets |
| `tests/gpu/cases.sh` | the case table, now shared with `run_regression.sh` |
| `Makefile` | `PROFILE=1`, `OPENMP=1`, and one object tree per combination |

Three decisions in there are load-bearing:

**The cases are not duplicated.** `run_regression.sh` was rewritten to read
`tests/gpu/cases.sh`, and `profile_gpu.sh` reads the same file. A profile of an
input that is not the one the suite checks is worth very little, and a pasted
copy of an input starts ageing the moment it is pasted. The suite's behaviour is
unchanged; it also takes case names as arguments now.

**NVTX ranges hang off `time_start`/`time_end`, not off new calls.** `times_t`
already names every phase and already has exactly one start/end pair for each.
The label is an optional argument, so unlabelled sites are untouched — but it
must be passed to *both* calls of a pair or neither, because NVTX ranges are a
stack.

**Each flag combination gets its own object tree.** make rebuilds on timestamps,
not on flags. Without the split, `make PROFILE=1` over an existing `build/`
relinks the old objects into a binary with no line info and no NVTX ranges, and
reports success. `profile_gpu.sh` also checks the binary for an NVTX symbol
before believing a number in the report.

On `alt`, `nvidia-smi` fails (NVML userspace/kernel skew) while CUDA works
perfectly, so the env check probes the device by compiling `cudaGetDeviceCount`
instead. `ncu` cannot run at all — `RmProfilingAdminOnly: 1`, needs root. Every
number below is from nsys.

---

## 2. The fix that mattered

`gpu_malloc_async` in `cuda_wrappers.cu` ended with

```c
  hipError_t err;
  hipDeviceSynchronize();
  err = hipGetLastError();      // and nothing ever read err
```

Every device buffer in TurboGAP is allocated through it, so a **whole-device
barrier ran on every allocation**:

| case | `cudaDeviceSynchronize` calls | % of GPU span |
|---|---|---|
| `CO_predict` | 3982 | **60.5%** |
| `estat_gsf` | 1025 | 49.7% |
| `XRD_mad` | 1838 | 31.0% |
| `vdw_ts` | 8929 | 23.8% |

It bought nothing: the error it collected was discarded with the test commented
out, and the `gpuErrchk` on the line above already checks the only status
`hipMallocAsync` reports.

Removing it took `CO_predict` from 3.888 s to 3.094 s of GPU span with **GPU
busy time unchanged at 2.03 s** — so all of it was host stall.

It was also a hard ceiling on everything in §4: `hipDeviceSynchronize` waits on
every stream, so while it was there no two streams could overlap, whatever the
batched loops did with them.

Removed with it: the `hipStreamSynchronize` at the end of
`gpu_get_2b_forces_energies`, commented "temporary, to measure timings", which
blocked the host on every 2b launch. `time%gap_2b` and the `2b` NVTX range moved
from per-descriptor to whole-routine so they still enclose the kernels.

---

## 3. Streams were dormant in three independent ways

Measured peak device concurrency was **1 on all four cases**, including the two
that create more than one stream. Three separate causes, all now fixed, none of
which is sufficient alone:

1. **No architecture makefile has ever passed `-fopenmp`.** `_OPENMP` was
   undefined, so `gpu_context_init`'s thread-count block was dead code, `n_omp`
   came out 1, one stream was created, and the `!$OMP PARALLEL DO` over the
   batched loops was a comment. Everything the tree says about "n_omp streams
   for the batched gpu calculation" described an intention, not a build.

2. **The pdf loop assigned `omp_task = 1` over the computed index**,
   unconditionally, on the line after computing it.

3. **The xrd loop used `mod(i-1, n_omp)+1`** — a property of the iteration, not
   of the thread. Under `static` scheduling two concurrent threads can draw the
   same stream, which serialises them and lets one thread's `hipFreeAsync` land
   on a buffer the other is still using.

`gpu_omp_task()` in `gpu_context.f90` is now the single answer for all sites,
and it is `omp_get_thread_num() + 1`. `OPENMP=1` is a make variable, **off by
default** — see §6.

There is a fourth gate that is not a bug: `gpu_n_batches` defaults to **1**, and
`gpu_context_init` clamps `n_omp` to it. With one batch there is nothing for a
second stream to do, so `--openmp` without `--input 'gpu_n_batches = N'`
measures a serial loop with idle threads.

---

## 4. Would streams help? Mostly no, and the data says where instead

`gpu%` is the fraction of a phase during which the device was doing anything.

| case | phase | wall | gpu% |
|---|---|---|---|
| `CO_predict` | `gap:soap` | 2.696 s (87% of span) | **62%** |
| `CO_predict` | `2b` | 0.239 s | 93% |
| `CO_predict` | `3b` | 0.159 s | 85% |
| `XRD_mad` | `exp_batched` | 0.281 s | **94%** |
| `estat_gsf` | `neigh` | **4.620 s** | **0%** |
| `estat_gsf` | `gap:soap` | 3.193 s | 42% |

**The loops that already have the stream machinery are the ones that least need
it.** `exp_batched` — the pdf/xrd batched path, the one place streams were ever
wired — is 94% device-busy. Perfect overlap there is worth under 7% of a phase
that is 21% of that run.

**SOAP batches will not benefit either, for a different reason.** From
`cuda_gpu_kern_gb_sum`, the SOAP kernels already fill the device on their own:
`cuda_get_soap_der_thr_one` launches ~9900 blocks × 64 threads and
`cuda_get_derivatives_new_new` 155×16×45 blocks. Two of those concurrently
cannot go faster than one after the other on 26 SMs. The 38% idle inside
`gap:soap` is **not** device idle waiting for more device work — it is the host
blocking, see §5.

**2b is the one case where the structure argues for it**, and it is worth
recording why it was not done here. `kernel_get_2b` launches **113 blocks × 64
threads = 7,232 threads** on a device that wants tens of thousands, and takes 74
ms per launch — so it occupies the timeline while badly underfilling the SMs,
and `gpu%` of 93% overstates its utilisation. Three descriptors on three streams
would put 3× the work in flight. What that needs first:

* `alphas_d`, `cutoff_d`, `qs_d` are **module-level** in `gap_backend_gpu.f90`
  and reallocated each iteration — they must become per-iteration locals;
* `kernel_get_2b` does `energies_d[i_site] += ...` and
  `forces_d[3*i_site+i] += ...` **non-atomically**. Safe today because
  descriptors run one after another; a read-modify-write race the moment they
  do not. `virial_d` already uses `atomicAdd`.

But before that, note what the kernel actually spends its time on: every thread
recomputes its own neighbour offset by walking the whole prefix,

```c
  for (i = i_beg - 1; i < i_site; i++)
    k += n_neigh_d[i];
```

so thread `i_site` does up to `n_sites` extra loads. **The 3b path already
precomputes exactly this array** (`kappas`). Doing the same for 2b is simpler
than the stream work, needs no atomics, and is very likely larger.

---

## 5. Where the idle in `gap:soap` actually is

After the allocator fix, the top host-side cost on `CO_predict` is

```
     total s  %span    count    mean us  call
      1.5760  50.9%     1840      856.5  cudaMemcpyAsync
```

856 µs of *host* time per "async" copy. The transfer table says why:

```
      count       volume     GB/s  kind
       1411    55.1 MB      6.7    memcpy HtoD
        120   196.9 MB    102.2    memcpy HtoD(p)
```

**1411 of the transfers are from pageable host memory.** `cudaMemcpyAsync` on
pageable memory is not async: the driver stages through an internal pinned
buffer and the call blocks until prior work in the stream drains. The tree
already has a pinned path — that is what the 102 GB/s row is — and it is 15×
faster on the device side as well.

So the highest-value remaining GPU-side work is **pinning the host buffers the
SOAP path uploads and downloads**, not adding streams. That is what turns
`gap:soap`'s 38% idle into overlap.

---

## 6. Why `OPENMP=1` is off by default

Turning it on makes two loops in `turbogap_exp.f90` run concurrently, and the
only suite case that reaches them (`XRD_mad`) is an **xfail** — so the suite
cannot distinguish a regression there from the drift it already tolerates.
Verify by comparing `trajectory_out.xyz` between an `OPENMP=0` and an
`OPENMP=1` build before trusting it.

Both builds compile and link clean, and `bin-omp/turbogap` links `libgomp`.

---

## 7. A blocker found on the way: the batched path has never run with >1 batch

`tools/profile_gpu.sh --input 'gpu_n_batches = 4' XRD_mad` aborts with

```
GPUassert: invalid configuration argument src/gpu_exp.cu
```

**Not a stream problem and not a memory problem.** The batched experimental
kernels compute their launch extents from per-batch, per-species-pair counts.
Those are never zero when one batch covers the whole cell; they go to zero as
soon as the cell is split, and CUDA rejects a `dim3(0,1,1)` launch. The
zero-byte allocation that goes with it returns `NULL`, which is what the
diagnostic prints show.

Nothing exercised this before: `gpu_n_batches` defaults to 1, and `estat_gsf` —
the only case in the suite that sets it higher — exhausts device memory before
reaching this code.

Fixed here: `gpu_get_pair_distribution_only` and
`gpu_get_pair_distribution_only_falloc` return early for an empty batch (the
answer for no pairs is zero, and the caller has already memset the output).
That gets `gpu_n_batches = 4` from batch 1 to batch 2, where **the next
instance of the same class of bug** stops it. This is systemic across the path,
not a single site.

Two things make it expensive to chase, and one is now fixed:

* `gpuErrchk(hipPeekAtLastError())` — `hipPeekAtLastError` does **not** clear
  the error, so once any launch fails every later check reports it, and the
  reported line is not the failing line. Switching these to `hipGetLastError`
  would make each site report only its own launch. Not done here: it is a
  mechanical change across dozens of sites and this session could not verify it.
* `TG_CHECK_LAUNCH` (new, in `gpu_exp.cu`) prints the kernel name and the
  offending extents before a zero-extent launch, instead of leaving
  "invalid configuration argument" to name neither.

**Until this is finished, `gpu_n_batches > 1` does not work on the experimental
observables — which means the multi-stream path cannot be measured end to end.**
That is the next thing to do if the streams matter.

---

## 8. Two defects in the timing report, noticed while wiring NVTX

* `time%exp_batched` was closed by a hand-written copy of `time_end`
  (`get_time` into (2), then the accumulation into (3)) with the xrd close in
  between. It is now `call time_end(...)`, which is also what pops its NVTX
  range — the hand-written form would have leaked the range and nested every
  later phase inside it.
* **`time%pdf` is printed but never accumulated.** `print_times` reports
  `- pdf:` and adds it into `Exp. pred.`, and no code anywhere calls
  `time_start(time%pdf)`. That line has always read 0.

---

## 9. Ranked, with the evidence

1. **Pin the host buffers on the SOAP upload/download path.** 1.576 s of host
   block on `CO_predict`, 50.9% of the span; the pinned path already exists and
   is 15× faster. (§5)
2. **Precompute the 2b neighbour offsets** the way 3b already does. Removes an
   O(n_sites) walk per thread from a kernel that is 11% of `CO_predict`. (§4)
3. **Finish the empty-batch audit** in the batched experimental path. Blocks
   `gpu_n_batches > 1` entirely, and with it any end-to-end stream measurement.
   Switch the `hipPeekAtLastError` sites to `hipGetLastError` first. (§7)
4. **Parallelise the neighbour build.** `estat_gsf` spends 4.62 s there with the
   GPU idle — longer than its entire GPU span. It is not a pragma: the
   `n_atom_pairs` running counter that indexes `neighbors_list` makes it a
   two-pass (count → prefix sum → fill) refactor. (§4)
5. **Then** streams for the 2b descriptors, which needs the module-level
   `alphas_d`/`cutoff_d`/`qs_d` privatised and `atomicAdd` in `kernel_get_2b`.
   (§4)
