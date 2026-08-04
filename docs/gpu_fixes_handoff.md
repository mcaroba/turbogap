# Handoff: `gpu_fixes` branch (GPU correctness, benchmarking, tests)

Written 2026-08-04. Everything below was verified on `alt` (`d22-0052`); nothing
has been pushed to any remote.

---

## 1. Current state

The GPU build agrees with the CPU reference. On the CO 7176-atom example,
`turbogap_gpu_commit_mahti` vs `turbogap_master_2026`:

| quantity | agreement |
|---|---|
| forces | max abs diff 4.0e-08 (\|F\|max 6.84) |
| virial | max abs diff 2.5e-06 |
| energy, energy_soap | 1.3e-05 |
| energy_2b / energy_3b | 4e-08 / exact |
| per-atom local_energy | 1e-08 |
| 5-step MD | temperature identical to all printed digits, pressure ~1e-7 rel |

Both `turbogap predict` and `turbogap md` run clean — no `GPUassert`, no
`compute-sanitizer` errors in the SOAP/2b/3b path.

**Branch:** `gpu_fixes`, HEAD `23fe556`, in both the parent repo and the
`src/soap_turbo` submodule (submodule HEAD `c605ddf`, and the parent records
it — see §5.2 for why that matters).

## 2. Where things are

| what | path on `alt` |
|---|---|
| GPU code under development | `/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_gpu_commit_mahti` |
| CPU reference | `/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_master_2026` |
| test inputs | `/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/input/CO` |
| scratch comparison dirs | `/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/cmp/` |
| nightly test logs | `<gpu repo>/tests/gpu/cron-logs/` |

Build the GPU code with:

```sh
export HOP_ROOT=/u/74/zarrout1/unix/work/hop
export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda
make DEBUG=0            # see §3, DEBUG defaults to 1
```

## 3. Environment gotchas

These cost real time to discover. Read before debugging anything.

* **`make` defaults to `DEBUG=1`**, which compiles CUDA with `-g -G` (device
  debug, optimisations off). That is a **2.1x** slowdown. Always build
  `DEBUG=0` for anything timing-related.
* **`nvidia-smi` fails** with an NVML driver/library version mismatch. This is
  cosmetic: the CUDA *runtime* works fine. Do not chase it.
* **`/usr/bin/compute-sanitizer` is broken** ("Unable to find injection library
  libsanitizer-collection.so"). Use
  `/usr/local/cuda-13.3/bin/compute-sanitizer` instead — that one works.
* **`ptrace_scope = 2`**, so `gdb` and `cuda-gdb` cannot attach. To locate a GPU
  fault, use the host backtrace that `gpuAssert` now prints (added on this
  branch) and resolve it with `addr2line -f -C -i -e bin/turbogap <addr>`.
* **The build has no `-fopenmp`**, so `_OPENMP` is undefined and `n_omp` is
  hardwired to 1 (`turbogap.f90`, the OPENMP STARTUP block). The GPU-stream
  machinery is therefore inert, and in any case it only wraps the batched
  xrd/pdf paths, not the main SOAP/MD path.
* **`alt` has one RTX A2000 (6 GB) and 12 cores.** No scheduler, no modules.
* **`numpy` is not installed on `alt`.** `tools/compare_xyz.py` is deliberately
  standard-library only for this reason; keep it that way.

## 4. What was fixed

Five distinct bugs. Commits are on `gpu_fixes`.

| # | commit | bug |
|---|---|---|
| 1 | `d06cefc` (submodule) | `st_soap` declared in `get_soap` but **never assigned**, then used as the length of the memset zeroing `soap_d`. Garbage length corrupted the device heap; surfaced later as an illegal access at the next allocation. This was the reported crash. |
| 2 | `f969c90` | Use-after-free: the 3b device buffers are allocated once *before* the descriptor loop but were freed *inside* it, so iteration 2 copied into dangling pointers. |
| 3 | `cff6092` | 2b species indices taken from the descriptor loop counter `i` instead of the species position `j` in `params%species_types`. |
| 4 | `cff6092` | `kernel_get_2b` / `kernel_get_core_pot` *assigned* to `energies_d` while *accumulating* into `forces_d`/`virial_d`; callers then read the buffers back inside the descriptor loop, re-counting each descriptor for every later iteration. |
| 5 | `2be0e10` | **The force bug.** In `kernel_get_radial_poly3gauss`, the amplitude derivative ended `- tmp2 - ampli_tude` where the CPU subtracts the *product* `atom_sigma_scaling/atom_sigma_scaled * amplitude`. A `*` typed as `-`. |

Bug 5 is worth understanding because of how well it hid:

* it is in the `else` branch, reached only when `amplitude_scaling` is neither
  0 nor 1 — and the adjacent `== 1.0` branch is written **correctly**, so the
  code reads as fine at a glance. The CO potential uses `amplitude_scaling = 2.0`.
* it affects only the **derivative**, never the coefficient, so **every energy
  stayed correct while every force was wrong**.
* the error grows as a pair approaches the cutoff, so MD *heated* (300 K to
  28,000 K over 100 steps) rather than merely drifting.

## 5. Testing

### 5.1 Running the suite

```sh
cd <gpu repo>
tests/gpu/run_regression.sh          # predict + 5-step MD vs the CPU build
tests/gpu/test_compare_xyz.py        # self-test of the comparator, no GPU needed
```

`run_regression.sh` honours `TURBOGAP_BIN`, `TURBOGAP_REF_BIN`,
`TURBOGAP_TEST_DATA` and `TURBOGAP_MPI_RANKS`. It fails on a non-zero exit, on
any `GPUassert` in the log, and on any energy/force/virial deviation beyond
tolerance.

**A nightly cron job is installed** (`crontab -l`):

```
30 3 * * * <gpu repo>/tests/gpu/cron_regression.sh >/dev/null 2>&1
```

It tests the binary currently in `bin/` rather than rebuilding — this is a live
working tree and a scheduled job should not silently replace the binary you are
using. `REBUILD=1` opts in. It skips (rather than competes) if another run holds
the lock or a `turbogap` process is already on the GPU. Results:
`tests/gpu/cron-logs/history.log` (one PASS/FAIL line per run) and
`latest.log`.

`.github/workflows/gpu-regression.yml` exists but the GPU job needs a
**self-hosted runner labelled `gpu`** — GitHub-hosted runners have no GPU, no
CUDA and not the test data. The static checks and comparator self-test run
anywhere.

### 5.2 Reproducibility

MD comparisons need identical starting velocities. Two mechanisms exist:

* `random_seed = <n>` input keyword (added to **both** repos) seeds the
  intrinsic PRNG deterministically, folding in the MPI rank.
* `input/CO/atoms_7176_vel.xyz` carries **explicit velocities**, which skips
  randomization entirely. This is the more robust basis and is what the MD
  regression case uses. Regenerate with `tools/add_velocities_to_xyz.py`.

**Submodule pointer:** the parent repo must record the `soap_turbo` commit that
contains the `st_soap` fix. It previously pointed at `c508abb` (pre-fix), so a
fresh clone plus `git submodule update` would have reproduced the original
crash. `23fe556` records `c605ddf`. Check this after any submodule work.

## 6. Benchmarks

Optimised builds (GPU `DEBUG=0` giving Fortran `-O2`; CPU `-fPIC -O3`), 7176
atoms, 1x RTX A2000 vs 12 cores:

| workload | CPU 1 rank | GPU | speedup |
|---|---|---|---|
| `predict` | 17.55 s | 5.56 s | 3.2x |
| `md`, 20 steps | 348.89 s | 86.84 s | 4.0x |
| `predict`, CPU 4 MPI ranks | 6.46 s | — | GPU 1.2x faster |
| `md`, GPU `mpirun -np 2` | — | 60.33 s | 1.44x over 1 GPU rank |

`mpirun` works and two ranks sharing the single physical GPU give correct
results (step-20 thermo matches the CPU and the single-rank GPU to ~1e-11).

## 6b. Experimental observables abort on this branch — OPEN

Found 2026-08-04 while adding regression coverage for the exp-spectra paths.

`turbogap md` with `do_pair_distribution` / `do_structure_factor` / `do_xrd`
**exits after ~1.3 s with status 0, having written no trajectory.** The CPU
build runs the identical input to completion. It is deterministic, and it is
not a hang or a timeout — the process returns success.

The log stops partway through the pair-distribution device setup, after

```
 allocing pdf
 allocing pdf to reduce
 gpu_exp%pair_distribution_partial_d(n_dim_idx)
 pdf_to_reduce_d
 x_d
 dV_d
 gpu_exp%rjs_index_d(n_dim_idx)
```

which is the `gpu_malloc_all` sequence in `electrostatics.f90` around line 440.
No `GPUassert`, no `stop`, no error message. Status 0 with truncated output
means something on that path is terminating the process while reporting
success, so `run_regression.sh` cannot detect it by exit code — the new case
detects it by the absence of `trajectory_out.xyz`.

Reproduce:

```sh
mkdir /tmp/gpuxrd && cd /tmp/gpuxrd
ln -s <data>/xrd_mad/{atoms.xyz,gap_files,xrd_glassy_carbon_zeng_2017.fq} .
cp <cpu repo>/tests/regression/cases/xrd_mad/input .
<gpu repo>/bin/turbogap md          # exits 0, no trajectory_out.xyz
```

**Why this matters beyond the bug.** These paths had no coverage at all here,
so nothing had ever run them. The `XRD_mad` case is marked xfail in
`run_regression.sh` rather than removed, so it stays visible and will report
XPASS the moment this is fixed. Until then the exp-spectra blocks in
`turbogap.f90` cannot be refactored on this branch with any verification —
which is why the CPU branch's Phase 4 extraction has not been brought across.

---

## 7. What is left

Roughly in priority order.

1. **Rename `tmp1`/`tmp2`/`tmp3` in `kernel_get_radial_poly3gauss`.** `tmp1`
   means both the amplitude polynomial `1 + rj^2(2rj-3)` and `dr^2/nf^2` in
   different sections of the same kernel. This reuse is precisely what hid bug 5
   for four rounds of investigation. It is the highest-value cleanup left, but
   it is the most numerically delicate kernel in the tree and deserves its own
   verify cycle against `tests/gpu/run_regression.sh`, section by section.
2. **Failure notification for the cron job.** Output goes to `/dev/null`, so a
   broken nightly is only visible in `history.log`. Needs a decision on the
   channel (email relay, Slack webhook, a watched file).
3. **`turbogap_master_2026` submodule is dirty.** `random_seed` support is
   committed there (`e6eb1aa`), but the `TURBOGAP_DUMP_CNK` instrumentation in
   its `src/soap_turbo/src/soap_turbo.f90` is **not**. Commit or revert
   deliberately — the nightly cron uses that binary as its reference.
4. **Continue the `turbogap.f90` refactor.** `0dafc9a` extracted the GPU 2b /
   core_pot / 3b blocks into internal subroutines (`add_2b_contribution_gpu`
   etc.). The file previously had *no* internal procedures at all; there is
   plenty more inline material worth the same treatment. Internal procedures
   are the cheap route here — host association means no argument lists.
5. **`-fopenmp` and GPU streams**, only if multi-GPU or batched xrd/pdf work
   becomes relevant. See §3.

## 8. Dead ends — already ruled out, do not re-derive

Chasing bug 5 eliminated a lot. Recorded so nobody repeats it.

**Verified to match the CPU, term by term:**

* `cuda_get_soap_der_one` — assembly and multiplicity indexing (`counter2`
  correctly advances only inside the non-skipped branch)
* `cuda_get_soap_der_two_one` / `_two_two` — normalization
  `der/sqrt_p - soap/p^3 * dot(soap,der)`; the shared-memory reduction is sound
  (`tpb = 64`, power of two, and `nthreads == tpb`)
* `cuda_get_soap_der_thr_one` / `_thr_two` — Cartesian transform and the
  central-atom sum; correctly launched over `n_atom_pairs` and `n_sites`
  respectively, and race-free
* `naive_transpose_soap_rad_azi_pol` — parameter *names* are misleading (the
  dimensions are passed swapped) but the arithmetic is correct
* `cuda_soap_forces_virial_two` — contraction and the symmetric virial formula;
  the index transposition is harmless because the expression is symmetric
* the entire soft-cutoff derivative block (`pref_f`, `der_pref_f`,
  `der_rjf_rj`, `der_sjf_rj`, and the `tmp5/tmp6/tmp7` rolling window over
  `exp_coeff_temp2(n, n+1, n+2)` — `tmp7` *is* refreshed each iteration)
* `angular_exp_coeff_rad_der` — agrees with the CPU to 2.8e-14

**Tested and disproved:**

* Making `cuda_global_scaling` apply the same `sqrt(rcut_hard)` factor to the
  coefficient and its derivative. This makes `cnk_rad_der` far worse
  (max deviation 1.1 -> 245). The opposite powers are **correct**: the CPU does
  `exp_coeff * sqrt(rcut_hard)` and `exp_coeff_der / sqrt(rcut_hard)`
  (`soap_turbo_radial.f90`), because the derivative is with respect to the
  reduced coordinate `r/rcut_hard`.
* The `if( .false. .or. ... )` guards in the radial kernel. These look like
  hand-neutered conditions but the **identical idiom is in the CPU source**
  (`soap_turbo_radial.f90` lines 611 and 677). Faithful port, not a bug.
* A "species-dependent" reading of the error. Early dumps appeared to show the
  error concentrated in the second species block; this was an **artifact** of
  each pair only writing its own species block (8 of 16 `n` values). Checking
  `nonzero` against `differing` showed 100% of computed entries were wrong.

**Removed as dead** (`07f5c93`), so do not go looking for them:
`cuda_get_derivatives`, `cuda_get_derivatives_new` (never launched, were
themselves inside a block comment, and wrote `cnk_*_der_d` in a layout
contradicting the live consumer), the commented-out timing harness in
`gpu_get_derivatives`, and `cuda_poly3gauss_one` (body commented out — it read
nothing, wrote nothing, and ran an empty loop nest).

## 9. Debugging technique that worked

The productive loop for numerical disagreement was **dump the same intermediate
from both codes and diff it**, narrowing one stage at a time:

`forces` -> `soap_cart_der` -> `cnk_rad_der` -> `{radial, angular}_exp_coeff_der`
-> the kernel that produces the radial one.

`TURBOGAP_DUMP_CNK` (set the env var; off by default) dumps `cnk_rad_der`,
`radial_exp_coeff_der` and `angular_exp_coeff_rad_der` for the first 20 pairs,
in a matching order, from **both** repos. Always dump `rjs` alongside — that is
what confirms the two codes agree on neighbour-pair ordering and so that the
comparison means anything.

Reading the code alone did not find this bug; every suspicious-looking thing
turned out to be a faithful port. The measurement did.
