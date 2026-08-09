# Session handoff — 2026-08-08 (c): device memory, and a million atoms

Third session of 2026-08-08. The previous one built the profiling workflow and
found that the dominant cost was a barrier in the allocator. This one is about
the other half of the same problem: **not knowing how much device memory a run
needs until it dies**, and whether the batching keywords that exist for it
actually do anything.

| | |
|---|---|
| new: device memory ledger | every wrapper allocation and free, with a peak |
| new: `gpu_mem_fraction` | sizes `max_Gbytes_per_process` from the card |
| new: scaling ladder | diamond, 13,824 → 1,000,000 atoms |
| GPU suite | 3 passed, 0 failed, 2 xfail — the documented baseline |

---

## 1. What a memory failure used to look like, and what it looks like now

Before, any device allocation that did not fit produced

```
GPUassert: out of memory src/cuda_wrappers.cu 133
GPUassert: host backtrace (10 frames):
  ... ten hex addresses ...
```

No size, no budget, no keyword. On a million-atom cell that is an afternoon.

Now the same failure prints how much was wanted, how much the device had, what
this rank was holding, the high-water mark, and the name and current value of
the keyword that controls it — and it tries once more before giving up.

### The ledger

`cuda_wrappers.cu` keeps a `pointer -> size` map under every wrapper allocation
and free. Cost is one hash operation per allocation, around 100 ns against the
2–3 µs `hipMallocAsync` itself takes.

It exists because `hipMemGetInfo` cannot answer the question. That reports what
the **device** has left, including every other process on the card; the question
the batching keywords turn on is what **this rank** is holding, and by the time
an allocation fails the useful figure is the peak, not the current value. Both
are printed side by side, at the end of every run:

```
GPUmem [end of run]
GPUmem   device         3.158 GB free of    5.653 GB total (44.1% in use)
GPUmem   this rank      0.000 GB live in 0 buffers, peak    1.587 GB
GPUmem   cumulative    41.203 GB requested over 9558 allocations, 9558 frees (0 untracked)
```

`live` at the end should be zero; anything else is a leak and the buffer count
says how big. This is the runtime counterpart to the static
`tools/gpu_check_alloc_pairs.py`, which can only see pairing in the source.

There is already a `total_gpu_memory()` on the Fortran side. It is not this: it
is hand maintained, each call site adds its own byte count, only the
experimental paths call it, and nothing ever compares it to the device. The
ledger sits underneath the wrappers instead, so it cannot disagree with what was
actually allocated.

### The retry

`hipFreeAsync` is **stream-ordered**: it returns a buffer to the pool when the
stream *reaches* the free, not when the call is made. A loop that allocates and
frees per batch is therefore always holding memory that is logically free and
physically not. On failure the allocator now drains the stream, trims the pool,
and tries once more — and says so, because needing it means the run is close to
the limit:

```
GPUmem: note -- a 0.412 GB allocation needed a pool drain to succeed.
```

This is on the failure path only. The normal path is one `hipMallocAsync` and
nothing else, which is the point of the previous session's removal of the
whole-device barrier that used to run there on **every** allocation.

---

## 2. `gpu_mem_fraction`: sizing the batching from the card

```
gpu_mem_fraction = 0.8
```

At start-up TurboGAP reads the device's free memory and sets

    max_Gbytes_per_process = fraction x free / ranks_per_device

That keyword is what splits the SOAP descriptor loop, via
`get_number_of_atom_pairs`. Its default is **1.0 GB, a number chosen with no
reference to any device** — on the 5.7 GB A2000 here it asks for roughly four
times more batches than needed; on an 80 GB card, eighty times.

Default is 0, meaning off, so no existing input changes behaviour until it is
set.

Three details that are not obvious and each of which was got wrong once while
writing this:

* **Where it is called from.** It has to be after `read_input_and_gap_files`,
  because it reads `gpu_mem_fraction`; placed next to `gpu_context_init` — which
  is where it belongs by subject — it ran before the input existed, saw the
  default zero, and reported budgeting off no matter what the input said. It
  also has to be before the first data allocation, or the "free" it budgets from
  already has TurboGAP's own buffers taken out of it. There is exactly one point
  that satisfies both.

* **Dividing by ranks.** `cuda_set_device` hands out cards round-robin
  (`rank % n_devices`), so when there are more ranks than devices they share one
  card's memory. Each budgeting the whole card is a guaranteed failure that
  looks like a code bug. `gpu_device_count()` was added for this.

* **A fraction, not all of it.** cuBLAS/cutlass workspaces, driver context
  growth and fragmentation inside the stream-ordered pool are all invisible
  here and all consume device memory. Sizing to 100% of free reliably produces
  a run that dies near the end.

### `gpu_n_batches`

Auto-sized for the **pdf/xrd** path, from
`estimate_max_exp_forces_device_memory_usage` — which enumerates those device
buffers one by one and **whose call site has been commented out all along**, so
nothing anywhere had ever compared `gpu_n_batches` to the card.

Electrostatics stays manual. There is no validated estimate for that path, and
an unvalidated one is worse than none.

**A value in the input is a floor, not a suggestion.** The automatic sizing can
only raise it. Someone who set it because the estimate was wrong for their
system is right, and silently overriding them is the failure mode to avoid —
this is the concrete answer to "it is likely that a manual setting is
necessary".

If the estimate says one site per batch still does not fit, it says so and
predicts the failure rather than letting the batching look like it coped.

---

## 3. The 1,000,000-atom case

`cpu_vs_gpu_tests/large_systems/diamond_1M/`, outside both repos like
`cpu_vs_gpu_tests/input/`. `atoms_1000000.xyz` symlinks the original MAD
weak-scaling snapshot; the smaller rungs are cut from it by
`tools/make_subcell.py`.

### The input is not the one the system came with

The MAD input is kept as `input.mad_reference` and **will not run here**. It
sets `pair_distribution_rcut = 12.6`, which at diamond's 0.1866 atoms/Å³ is
~1470 neighbours per atom. At a million atoms that is 1.5×10⁹ pairs, and the
**host** arrays alone — `rjs`, `thetas`, `phis`, `xyz` — come to ~75 GB before
the GPU is asked for anything. `alt` has 31 GB. It was a 128-node LUMI run; that
memory is meant to be divided over ranks.

`input.gap` uses the potential only. `carbon.gap`'s own cutoffs are 4.5 Å (SOAP)
and 2.8 Å (2b), giving ~72 neighbours per atom, and that is what makes the
ladder runnable on one card.

### Cutting sub-cells correctly

`tools/make_subcell.py` takes a sub-**box**, not the first N atoms. The first N
lines of a supercell are a slab with two open faces; run periodic, it puts atoms
at bonding distance across a boundary that is not there, and every neighbour
count you were trying to measure is wrong.

Two ways the cut still lands wrong, both of which showed up as a density
*below* the parent's — i.e. atoms lost, and lost one-sidedly:

1. This is a 300 K snapshot, so an atom that started at x = 0 can be at
   x = −0.05. A raw `0 <= x` test drops it, and its image at x ≈ L is outside
   the sub-box too, so it is lost from both ends. Wrap first.
2. Cutting at a whole number of repeats puts the boundary exactly through a
   plane of atoms, where thermal displacement decides each one's fate. Offset
   by an eighth of a repeat and the boundary falls in the gap.

With both, the sub-cell densities match the parent to better than 0.1% — which
is the check that the cut landed, and the reason the script prints it.

---

## 4. Measured: the ladder on one RTX A2000 (5.65 GB, 26 SMs)

`tools/run_scaling_ladder.sh`, `predict`, `gpu_mem_fraction = 0.8`
(budget 3.80 GB/rank), one MPI rank.

<!-- LADDER TABLE -->

Two things this says:

* **Device peak is flat.** It does not grow with system size, because the SOAP
  batching splits the loop to fit the budget — which is the whole design working
  as intended, and the first direct evidence that it does. The card is not the
  limit here.
* **The host is what grows.** Peak RSS is linear in atom count, and it is the
  neighbour list that dominates it. That is the ceiling on this machine, not the
  GPU.

---

## 4b. Where the time goes at scale

TurboGAP's own timer table, per rung, on a **clean DEBUG=0 build** (see §4c --
the first version of this table was not, and was wrong):

| atoms | total s | 3b | soap | neighbours | 2b |
|---|---|---|---|---|---|
|  13,824 |   8.9 |  4.45 (50%) |  4.18 (47%) | 0.31 | 0.24 |
|  32,785 |  19.7 |  9.87 (50%) |  9.28 (47%) | 0.74 | 0.44 |
| 124,959 |  75.7 | 36.99 (49%) | 36.85 (49%) | 2.73 | 1.35 |
| 262,202 | 168.6 | 78.48 (47%) | 86.61 (51%) | 5.88 | 2.36 |

**Three-body and SOAP are co-equal, about half the run each**, and stay that way
across a 19x range of system size. The neighbour build is 3.5% and the two-body
term 1.4%.

That is a different potential's balance from `CO_predict`, where SOAP was 90% of
the wall clock and 3b 4%. Both are true; they are different potentials on
different systems, which is the argument for having a large case at all. What it
means practically is that there is no single hot spot here: halving either term
buys 25% of the run, and halving both buys 50%.

The neighbour build is 3.5% here, not the bottleneck it was on `estat_gsf`.
That case uses a 10 A cutoff; this one uses carbon.gap's 4.5 A. Neighbour cost
is set by the cutoff, not by the atom count.

---

## 4c. A correction, and the defect behind it

**Everything in the first version of §4 and §4b was measured against a binary
that was partly compiled with DEBUG=1, and was wrong by a factor of 2.4 to 2.9.**

`make DEBUG=0` recompiles nothing if no source changed. The architecture
makefiles default `DEBUG ?= 1`, so a tree built once with the default and then
built again with `DEBUG=0` **links the debug objects and reports success**. In
`build/`, `3b_final.o` was three days old and DEBUG=1; so were the four
`soap_turbo_*` objects and ten others.

What that cost, all of it invisible:

| | stale tree | clean DEBUG=0 |
|---|---|---|
| 3b, 125k atoms | 173.0 s | **37.0 s** (4.7x) |
| whole ladder | 23 -> 2029 s | **9.5 -> 750 s** |
| 3b share of the run | 81% | **49%** |
| SOAP energy, 125k | -1245959.20601018 | -1245959.20614126 |

The energy difference is the tell that makes this worth writing down. It is
1e-10 relative -- far too small to fail any test, far too large to be nothing --
and it is a *deterministic* difference between two builds of identical source,
because `-O0` Fortran and `-O2` Fortran contract floating-point differently.
Chasing it as if it were a numerical regression in a GPU kernel is exactly the
wrong afternoon, and that is what nearly happened here.

**Fixed**: `DEBUG=1` now gets its own object tree (`build-dbg/`, `bin-dbg/`),
the same treatment `PROFILE` and `OPENMP` already had. The two builds are
physically separate, so the flag cannot lie about what is in the binary. This
was the one flag that had been left out.

The device-memory numbers in §4 are unaffected -- peak device memory is
identical on both trees, because it is set by the batching, not by optimisation.

---

## 4d. The three-body occupancy change: implemented, measured, reverted

`kernel_2nd_try` launches one 32-thread block per atom. On Ampere the hardware
caps resident blocks at 16 per SM, so that is 16 of 48 warp slots -- a 33%
ceiling from the launch geometry alone.

It was rewritten to put several warps in a block, one atom per warp: the
per-atom shared scratch indexed by warp, every `__syncthreads()` replaced by
`__syncwarp()` (warps in one block are on different atoms with different
neighbour counts, so a block-wide barrier there is a hang, not a slowdown), and
the launch reshaped. `TG_3B_WARPS_PER_BLOCK=1` reproduced the old geometry
exactly, as a self-check.

**Measured, on clean builds, 125k atoms:**

| build | 3b time |
|---|---|
| original | 37.03 s |
| warps/block = 1 | 37.95 s |
| warps/block = 4 | 37.77 s |
| warps/block = 8 | 37.83 s |

No benefit, and consistently ~2% slower. `ptxas` says why the ceiling was not
the ceiling: the kernel uses **82 registers**, which caps occupancy at ~24 warps
(50%) regardless of block size. The block-count limit of 33% and the register
limit of 50% are close enough that lifting the first buys almost nothing, and
the kernel is evidently not occupancy-bound anyway.

**Reverted.** The energies were bit-identical throughout (3b energy agreed to
all printed digits on both systems), so this is not a correctness retreat -- the
change simply did not pay, and it added indirection to a numerically sensitive
kernel that the suite checks bit-exactly.

Worth keeping from it: the register figure, and the fact that 3b is not
occupancy-limited. Whatever makes it faster, it is not more warps in flight.

---

## 4e. Several MPI ranks on one GPU: now works, and what it buys

The question was whether each GPU is tied to one rank, and whether the MPI cores
can drive the device or whether that needs OpenMP and streams.

**Nothing ties a GPU to a rank.** `src/cuda_wrappers.cu` does
`hipSetDevice(my_rank % num_gpus)`, so on a one-GPU node every rank binds to
device 0. The MPI decomposition splits SITES (`i_beg`/`i_end`), and
`build_neighbors_list` takes a `do_list` mask, so each rank already builds only
its own atoms' neighbours.

**But it crashed, and had done since before this work** -- verified against a
pristine `git worktree` of `e70d896`, while the master CPU branch runs two ranks
fine. Four defects, all one bug wearing four hats.

### The pattern

`gap_backend_gpu.f90` mixes two indexing conventions, and it is not arbitrary
which applies where:

* **Sites are GLOBAL.** `i_beg`/`i_end` are absolute site numbers, the host
  arrays (`energies_2b`, `this_forces`, ...) are allocated `(1:n_sites)` over
  the whole system on every rank, and the kernels index them with the global
  site number.
* **Pairs are RANK-LOCAL.** `j_beg = 1` always (turbogap.f90:945); each rank's
  neighbour build produced only its own pairs, numbered from 1.

Every one of the four bugs was a buffer sized or read with the rank-local count
while being indexed globally. **All four are invisible on rank 0**, because
`i_beg == 1` makes the two conventions coincide there -- which is why
single-rank testing never saw any of it.

| where | what was wrong | fix |
|---|---|---|
| `gap_backend_begin` | `n_neigh_d`, `species_d` sized `(i_end-i_beg+1)` | `i_end` |
| `add_2b_contribution`, `add_core_pot_contribution` | `energies_*_d`, `forces_*_d` sized rank-locally | `size(energies_*)` |
| `add_3b_contribution` | `kappas` allocated rank-locally, indexed and uploaded globally | `size(n_neigh)`, zeroed |
| `add_3b_contribution` | READ-BACK length rank-local, from offset zero | `size_energy3b` / `size_forces3b` |

The last one is the instructive one. It did not crash: it copied a slice the
kernel had never written, so every rank but the first silently contributed
**zero** 3b energy. The tell was a 3b energy of exactly 1/2 on two ranks and
exactly 1/4 on four -- but *not* 1/3 on three. The test cell is a 2x2x2
replication, so "rank 0 only" produces exact halves and quarters and a ragged
third; a division would have produced an exact third too. That distinction is
what separated "only rank 0 counted" from "something divided by nranks".

Also fixed, in the suite itself: `mpirun` reads stdin, and the case loop feeds
names to `while read` on the shell's stdin, so under `TURBOGAP_MPI_RANKS` the
first case swallowed the rest of the list and the suite reported
`passed: 1` **as a success**. `</dev/null` on the run.

### Verified

| | |
|---|---|
| 1 / 2 / 3 / 4 ranks, CO_predict | identical energies, to every printed digit |
| trajectories, 2/3/4 vs 1 rank | `compare_xyz.py` PASS |
| full GPU suite at `TURBOGAP_MPI_RANKS=2` | 3 passed, 0 failed, 2 xfail |
| against the CPU branch | 3b 2.13601985 eV, rank-invariant on both |

### Measured: `tools/run_rank_sweep.sh`, 124,959 atoms, one RTX A2000

| ranks | wall s | speedup | wall s (MPS) | speedup (MPS) | device peak/rank |
|---|---|---|---|---|---|
| 1 | 76.2 | 1.00x | 76.9 | 1.00x | 1.687 GB |
| 2 | 65.7 | 1.16x | 69.4 | 1.11x | 0.826 GB |
| 4 | 66.5 | 1.15x | 63.0 | **1.22x** | 0.397 GB |
| 6 | 68.7 | 1.11x | 64.5 | 1.19x | 0.251 GB |

**Extra ranks are worth about 15-20% here, and that is the honest ceiling.**
MPS -- which makes kernels from different ranks genuinely concurrent rather than
time-sliced -- moves it from 1.15x to 1.22x and no further. Neither number is
disappointing once §4b is taken into account: 3b and SOAP together are ~98% of
this run and both are GPU-bound, so there is very little host work left to
overlap and very little device idle to fill. The host-side neighbour build,
which multi-rank parallelises perfectly, is 3.5%.

**The real payoff is memory, not speed.** Device peak per rank falls 6.7x from
one rank to six, and host RSS with it. Multi-rank on one GPU is how a system
that does not fit is made to fit -- which is exactly the question the
million-atom case was about. A workload with a larger host-side fraction (a big
cutoff, so an expensive neighbour build) would also see much more than 20%.

### OpenMP + streams, the backup

`make OPENMP=1` builds and is **verified numerically identical**: the full suite
at `OMP_NUM_THREADS=4` gives 3 passed, 0 failed, 2 xfail, the same as the serial
build. `gpu_omp_task()` gives each thread its own stream.

It remains the narrower option, and deliberately so:

* It covers only the two batched loops in `turbogap_exp.f90` (pdf and xrd).
* It is gated behind `gpu_n_batches > 1`, and that path still hits the
  empty-batch launch bug (§7).
* **The electrostatics loop must not be threaded as it stands.** Its batches are
  independent in device state and in `energies_estat`, but `forces_estat` and
  `virial_estat` are passed whole to every batch and accumulated into -- two
  threads would read-modify-write the same array. The directives are commented
  out for that reason and now say so. Threading it means making forces and
  virial per-batch and summing after the loop, which is what the pdf path
  already does.

So: **use ranks first**, now that they work; OpenMP is for the case where a rank
per core is not wanted, and it needs the §7 work before it can be measured.

---

## 5. What to do next

1. **Give `kernel_2nd_try` more than one warp per block** (§4b). 81% of a large
   diamond run, capped at ~33% occupancy by its launch geometry. Nothing else
   measured here comes close.
2. **The host is the size ceiling, not the device.** Device peak is flat in
   system size — the batching works. Peak RSS is linear, and it is the
   neighbour list. At 4.5 Å that is ~72 neighbours per atom; at the MAD cutoff
   of 12.6 Å it is ~1470, a 20x memory factor and the reason the MAD input
   needs many ranks. That is an MPI decomposition problem, not a batching one.
3. **`build_neighbors_list` reallocates by copying.** `n_neigh_max` starts at
   100 and grows by 10, copying the whole `neighbors_list` each time. Harmless
   at 72 neighbours per atom; quadratic at 1470, i.e. exactly in the case that
   already does not fit.
4. The electrostatics batch estimate is still missing (§2), which is why
   `estat_gsf` remains the only suite case that has ever run out of memory.
5. Everything in `session_handoff_2026-08-08b.md` §9 still stands. Note that it
   ranked SOAP work first on `CO_predict` evidence; §4b above says that ranking
   is system-dependent, and for large carbon systems 3b comes first.
