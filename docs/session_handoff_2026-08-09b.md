# Session handoff, 2026-08-09b: the 1M diamond, with and without MAD

Profiling the million-atom diamond on `alt` (RTX A2000, 31 GB host), for MD with
the potential only and MD with the experimental-observable machinery on, and
acting on what it said. Companion to `docs/PROFILING.md`.

Nothing here is committed. `src/soap_turbo/src/soap_turbo.f90`, `src/turbogap.f90`
and `src/gap_interface.f90` carry the changes in section 3.

---

## 1. The baseline, no MAD

`profiling/diamond_1M_md/`, `input.md` with `md_nsteps = 3`.

```
3 MD steps:     3583.322 s
*  GAP desc/pred:  3423.302 s   95.5%
   - soap_turbo:   2027.837 s   56.6%
   - 3b:           1358.675 s   37.9%
   - lin__turbo:    496.289 s   13.8%
*  Neighbor lists:   100.051 s    2.8%
peak RSS 8.3 GB, 99% CPU, 234,926,305 minor page faults
```

Two cautions on reading this. `lolo__soap` prints 7562 s against a 3583 s total
-- it is double-counted, ignore it. And 99% CPU does **not** mean host-bound:
CUDA busy-waits, so a run blocked in `cudaStreamSynchronize` looks identical to
one doing work. Separating those needs nsys, which is section 2.

## 2. Where the time actually is (nsys, 125k)

A capture at 1M would be huge and slow; the phase structure is the same at 125k
and it costs minutes.

```
GPU busy 82.0%   idle 18.0%   peak concurrent ops 1
gap:2b_3b_corepot  117.5 s  gpu 99.9%
gap:soap            85.2 s  gpu 63.7%
neigh                8.3 s  gpu  0.0%
kernel_2nd_try (3b)  113.2 s = 54.1% of span, 3 launches, 37.7 s each
cudaStreamSynchronize 113.3 s = 54.1%   <- the 99% CPU, explained
```

So the run is device-bound overall, `3b` is genuinely so, and the recoverable
host time is the 36% idle inside `gap:soap`.

## 3. The scaling result, and the experiment that explains it

Per MD step, 125k -> 1M (8x the atoms):

| phase | 125k | 1M | ratio |
|---|---|---|---|
| `3b` | 56.2 s | 452.9 s | **8.06x -- linear** |
| `soap_turbo` | 43.0 s | 675.9 s | 15.7x -- **N^1.32** |

The superlinear excess is ~332 s/step, about **28% of the whole run**.

Hold the system at 125k and vary only `max_Gbytes_per_process`, so only the
batch count changes:

| budget | batches | `soap_turbo` | `3b` |
|---|---|---|---|
| 1.0 GB | ~4x | 77.8 s | 74.7 s |
| 2.0 GB | ~2x | 60.0 s | 75.2 s |
| 4.0 GB | ~1x | 54.0 s | 75.4 s |

`soap_turbo` falls 31% for a 4x cut in batches. `3b` is flat to 1%, which is the
control that makes this a measurement rather than a coincidence.

**Per-batch overhead is O(n_sites) per batch, so at a fixed memory budget it
grows as N^2.** Scaling the 125k per-batch cost (~6 s/step) by 64 predicts ~384
s/step at 1M; ~332 was measured.

Two things that look like the cause and are not, so they do not get re-derived:

* `VmPeak` 33.5 GB against 8.3 GB RSS is **CUDA's address-space reservation**.
  It is not evidence about TurboGAP's allocations.
* `assign_species_multiplicity` rescans the whole supercell once per batch
  (`xyz_species_supercell` is full-system while `xyz_species` is the batch
  slice, so it always takes the `n_sites_supercell > n_sites` branch). It is
  ~25 s/step at 1M. Worth hoisting eventually; not the headline.

## 4. What was changed

Both are behaviour-preserving and both attack the O(n_sites)-per-batch term.

**`get_soap` no longer allocates the CPU implementation's working set.**
`cnk_{rad,azi,pol}_der`, `soap_{rad,azi,pol}_der` and the angular derivative
arrays -- about **28 kB per atom pair** -- were allocated at full size on every
batch and three of them memset, though nothing on this branch reads any of it.
The only consumer is `get_derivatives`, which master calls and this branch does
not; the file already carried a comment saying so. They are shrunk to dummy
extents rather than deleted, because every `st_*` that sizes a **device** buffer
is `sizeof(x(1,1))` -- the element size -- so the declarations and leading
extents have to stay. That is the idiom the file already used for
`radial_exp_coeff_der`.

**The driver stops touching the whole system once per batch.** It zeroed
full-system `this_energies`/`this_forces`, then added the whole array into the
accumulator, to move the contribution of ~1200 sites: 128 MB of traffic per
batch, ~108 GB per step at 1M. `get_gap_soap` now **adds** through
`in_to_out_site` directly into `energies_soap`/`forces_soap`, which the driver
already zeroes once before the descriptor loop. Bit-identical, not merely
equivalent: the per-element sequence of adds is unchanged, and `in_to_out_site`
is injective within a batch so `+` cannot double-count.

### Measured

| | pre | post |
|---|---|---|
| 1M, 3 steps | 3583.3 s | **3321.9 s (-7.3%)** |
| `soap_turbo` | 2027.8 s | 1791.1 s (-11.7%) |
| `get___soap` | 86.9 s | 26.9 s (-69%) |
| `3b` (untouched control) | 1358.7 s | 1336.6 s (-1.6%) |
| minor page faults | 234.9M | 184.6M (-21%) |
| peak RSS | 8306272 kB | 8306232 kB (flat) |
| 125k, 2 steps (clean) | 218.0 s | 214.7 s (-1.5%) |

**The 1M pre-patch baseline overlapped the CPU regression suite for ~16 of its
60 minutes**, so -7.3% is an upper bound. Untouched controls moved -1.6% and
+1.2%, putting noise near +-2%. The page-fault count is a count, not a time, so
-21% is contention-free and is the cleanest confirmation of the mechanism. RSS
being flat is expected: the removed arrays were never touched, so they were
never resident.

### Verification

Built the pre-change sources to a saved binary and passed it as
`TURBOGAP_REF_BIN`, making the suite a direct A/B rather than a GPU-vs-CPU
comparison. `CO_predict`, `CO_md`, `vdw_ts`: **bit-identical**. See section 6
for why the other two prove nothing.

## 5. MAD at a million atoms: it runs, but only just

`input.mad` and the `diamond_*_mad` cases are new (`tests/gpu/cases.sh`). Getting
a million-atom MAD step to complete on this box took two cuts to the cutoff, and
the cutoff is the whole story:

| `pair_distribution_rcut` | outcome | peak RSS |
|---|---|---|
| 12.6 A (the original MAD input) | not attempted -- the ladder README costs it at ~75 GB | -- |
| 6.0 A (what `XRD_mad` uses) | **OOM-killed**, signal 9, on xrd batch 254 of 350 | 28.0 GB |
| 5.0 A | completed | **24.6 GB** |

`alt` has 31 GB total, ~24 GB available. The 5.0 A run finished with essentially
no margin, and it is a truncated PDF -- `r_range_max` 4.5 A covers only the first
few coordination shells of diamond. **Treat 1M MAD as past this machine's
ceiling, not as a supported configuration.**

### One MD step, 1M atoms, `pair_distribution_rcut = 5.0`

```
1 MD step:        2564.112 s        peak RSS 24.6 GB
*  GAP desc/pred:  1821.661 s   71.0%
   - soap_turbo:   1105.634 s
   - 3b:            687.956 s
*  Exp. pred.:      573.263 s   22.4%
   - xrd:           573.252 s
   - sf:              0.011 s
   - pdf:             0.000 s
*  Neighbor lists:    57.349 s    2.2%
```

**The MAD machinery is 22% of the step and it is entirely `xrd`.** `pdf` and
`sf` read as zero not because they do not run but because
`structure_factor_from_pdf = .true.` routes S(q) out of the pair distribution and
XRD out of S(q), so the whole chain lands in the `xrd` bucket. If that phase is
ever optimised, the timers will need splitting before the result can be
attributed.

Note the GAP phase is ~15% more expensive per force evaluation here than in the
no-MAD run, despite being the same potential. It is the same effect as section 3:
the global neighbour list is built at 5.0 A rather than 4.5 A, so every SOAP
batch gathers and filters more input pairs to produce the same descriptors.

### What actually runs out, and why ranks will not save it

The device side sized itself correctly and stayed inside its budget --
`peak 3.743 GB` against `3.796 GB`, with the estimator raising the count as it
went:

```
pdf/xrd batched forces needs ~ 1590.164 GB, budget 3.796 GB/rank -> batches 210 -> 419
```

**What runs out is host memory, and nothing sizes it.** Worse, on this branch
sites are global and the host arrays are full-system on every rank, so adding
MPI ranks does not divide it -- the usual answer to "it does not fit" is not
available. That is the real ceiling for MAD at scale, and it is a different
problem from the device budgeting that `gpu_mem_fraction` solves.

The host high-water mark also **grows between force evaluations** -- ~18 GB on
the first, ~24.6 GB on the third. Whether that is a leak in the exp path or just
allocator fragmentation from 419 allocate/free cycles per evaluation is not
established, and it is worth establishing: it is the difference between "1M MAD
needs 25 GB" and "1M MAD needs 25 GB and rising".

### A bug found on the way

`this_j_beg = ******** this_j_end = ********` in the xrd batch printout. The
format was `'(A,I8,...)'` with six edit descriptors for eight integers, so it
wrapped onto a second line, and `I8` cannot hold a pair index -- there are ~168M
pairs at a million atoms with a 6 A cutoff. Fixed to `I0` with eight descriptors
(`src/turbogap_exp.f90`). The printout itself is unconditional and emits four
lines per batch, which is 890 KB of log for a single step at this size; that is
worth gating behind a verbosity flag, but removing output someone may rely on is
not a call to make in passing.

### Early signal, for contrast

At 13,824 atoms and one step, `xrd` was 28.7 s against 32.5 s for the whole GAP
prediction -- MAD roughly doubled the step. At a million atoms it is 22%. The
MAD phases scale better than the potential does.

## 5b. Recalibrating the bytes-per-pair estimate

`get_number_of_atom_pairs` split the descriptor loop until an estimated batch fit
in `max_Gbytes_per_process`, estimating

    n_atom_pairs_in * n_max * k_max * 150

Two things are wrong with that *shape*, independently of the value. Only one term
is genuinely proportional to `n_max*k_max` -- the three `cnk` derivative arrays,
`48*k_max*n_max`. Everything else goes as `n_soap`, as `n_max`, as `n_species`, or
is flat, so no choice of constant tracks them as the potential changes. It is in
particular blind to **SOAP compression**, which for `carbon.gap` takes `n_soap`
from 324 to 72 and leaves the estimate charging ~5 kB per pair for memory nobody
allocates. And it counted nothing **per site**, though `cnk_d` is `16*k_max*n_max`
per site -- the same order as the per-pair total divided by the neighbour count,
and the term that decides the answer for a short cutoff.

So the constant is gone, replaced by `soap_batch_memory_model` in
`src/neighbors.f90`, which enumerates the buffers `get_soap` and `get_gap_soap`
actually allocate, one line per allocation, checkable against the source:

```
PER PAIR                                    PER SITE
  48*k_max*n_max  cnk_{rad,azi,pol}_der_d     16*k_max*n_max  cnk_d
  64*k_max        angular_exp_coeff + 3 der   16*n_soap       soap_d + host copy
  72*n_soap       soap_*_der_d, soap_cart_der 24              sqrt_dot_p_d, n_neigh_d,
  40*n_max        radial_exp_coeff + 3 temps                  k2_start_d, i_k2_start_d
   8*n_species    mask_d + host copy
  88              rjs/thetas/phis/xyz, indices
```

For `carbon.gap` (n_max 8, l_max 8 so k_max 45, compressed n_soap 72) that is
**25,760 bytes per pair against the 54,000 the constant assumed** -- the loop was
being split about 2.1x more finely than the memory required. `n_soap` and
`n_species` are now passed in from the hypers; `SOAP_BATCH_SAFETY = 1.10` covers
what the enumeration cannot see (the per-descriptor sparse set, allocator
fragmentation, the stream-ordered pool holding logically-free blocks).

### Measured

1M, 3 MD steps, both runs alone on the box with the same binary -- no contention
caveat on this one:

| | patches 1&2 | + recalibration | |
|---|---|---|---|
| **total** | 3321.9 s | **2826.5 s** | **-14.9%** |
| `soap_turbo` | 1791.1 s | **1299.5 s** | **-27.4%** |
| `3b` (control) | 1336.6 s | 1341.5 s | +0.4% |
| `lin__turbo` | 502.1 s | 511.6 s | +1.9% |
| device peak | -- | **3.220 GB** of 3.796 budget | 85% |
| peak RSS | 8306 MB | 8460 MB | +1.8% |

125k, 2 steps: 214.7 -> 207.3 s (-3.4%), `soap_turbo` 81.6 -> 76.0 s (-6.9%),
device peak 3.186 GB. Small there because 125k is already near its floor at this
budget -- the sweep in section 3 shows 2.0 -> 4.0 GB buying only 10%. The gain is
at 1M, where the batch count is 8x higher, which is the whole point.

**The superlinear term is mostly gone.** Section 3 measured a superlinear excess
of ~332 s/step in soap; it is now ~89 s/step, about 73% recovered. `soap_turbo`
has dropped below `3b`, so **3b is now the dominant phase at a million atoms**.

### What it costs

The batch count changes, so the partial sums regroup and the last digits move.
On `CO_predict`, against the same binary without the recalibration: energies
**exactly** equal, forces and virial `1.0e-8` absolute, `5.19e-13` relative on
the virial -- the last ULP of the printed format. Unavoidable and the same class
as `mem_fraction` on the CPU branch.

Verified byte-for-byte while establishing that: **patches 1 and 2 (section 4) are
byte-identical** to the pre-change binary, same 3974 allocations, same output
bytes. The recalibration alone takes that to 2674 allocations.

**A coverage gap this exposed.** The suite reports PASS for the recalibration,
but its comparison has a tolerance wide enough to absorb a batch-count change.
So it cannot detect a batching regression -- only an out-of-memory abort would.
If the estimate is ever changed again, compare allocation counts from the ledger
and the device peak against the budget, not just the suite result.

## 6. Two GPU cases are non-deterministic, which is what their xfails record

Running the **same** binary twice on the same input:

* `XRD_mad`: `local_energy` diverges to **5.4e-6 by frame 5** (exact through
  frame 3). Energies, forces and virial agree to 0.000e+00 in every frame.
* `estat_gsf`: the **virial** varies at rel **6e-5**, and comes out
  **asymmetric** within a single run -- (1,3) = -71.17 against (3,1) = -71.07.
  A symmetric tensor emerging asymmetric is a racy accumulation. This is a
  single-point `predict`, so MD chaos cannot explain it.

Both are marked xfail with reasons attributing the failure to a **GPU-vs-CPU
disagreement**, at the same quantities and magnitudes the GPU build produces
against itself. So part of each xfail is irreproducibility, and **the suite
cannot currently certify those two paths either way**. Anyone A/B-testing a
change should rely on the three deterministic cases and run the same binary
twice before believing a difference.

Not diagnosed. `estat_gsf`'s asymmetric virial points at the accumulation in
`calculate_batched_electrostatics`.

## 7. Streams for the SOAP loop: the measurement argues against it

Asked directly, and the answer is no, for a quantitative reason rather than a
plumbing one. At 1M the host portion of `gap:soap` is roughly 1600 s against
~430 s of GPU work in that phase. Overlap is bounded by the smaller of the two,
so there is nowhere near enough device work to hide the host work behind. And
concurrent batches each need their own device memory, so the budget divides by
the number of streams, the batch count rises, and section 3 says that makes the
quadratic term **worse**. Removing the host work beats hiding it, and it helps
the CPU branch too.

Revisit only once `gap:soap` is near its floor. The plumbing is real but is not
the objection: `get_soap` already takes the stream and cuBLAS handle as
arguments, but `gap_interface.f90` and `gap.f90` take them from module state
(13 and 52 references), and `l_index_d` and `recompute_basis` are shared mutable
state that would race.

## 8. Next, in order of measured value

1. **`3b` is now the largest cost at 1M** (1341 s of 2826, 47%), exactly linear
   and 100% GPU-busy. The warps-per-block rewrite was tried and reverted
   (registers cap occupancy at 50%). Hardest remaining target, and now the one
   that matters most.
1b. **What is left of the soap superlinear term** is ~89 s/step, down from ~332.
   The remaining O(n_sites)-per-batch work is `assign_species_multiplicity`'s
   supercell rescan (~25 s/step) plus `is_atom_seen` / `out_to_in_site` in
   `get_gap_soap`.
2. **Host memory for MAD.** Section 5. Nothing sizes host allocation, ranks do
   not divide it, and the high-water mark grows between force evaluations.
   Settle the leak-vs-fragmentation question first -- it decides whether this is
   a ceiling or a slope.
3. Hoist the supercell scan out of `assign_species_multiplicity` (~25 s/step).
4. Diagnose the two non-deterministic cases (section 6) -- until then two of
   five suite cases certify nothing.
5. Widen the suite's comparison, or add a batching check -- see 5b. Two of five
   cases already certify nothing (section 6); a third blind spot is one too
   many.
