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

## 6b. Experimental observables — abort FIXED, virial still wrong

Found 2026-08-04 while adding regression coverage for the exp-spectra paths,
which had none.

### The abort — fixed

`turbogap md` with `do_pair_distribution` / `do_structure_factor` / `do_xrd`
exited after ~1.3 s with **status 0** and no trajectory. The cause was a
leftover debugging `exit( 0 )` in `gpu_exp.cu`, inside
`gpu_get_pair_distribution_only_falloc`, sitting between the KDE evaluation
kernel and the reduction kernel:

```c
  printf("exiting from pdf falloc kernel \n");
  fflush( stdout );
  exit( 0 ) ;
  nblocks=dim3(n_samples,1,1);            // never reached
  kernel_reduce_pair_distribution<<<...>>>(...);
```

So `kernel_reduce_pair_distribution` never ran and the process reported
success. Removed. The GPU build now runs the case to completion.

Two things made this survive: the paths had no test at all, and the exit
status was 0, so nothing that checks a return code would have noticed. This
is why `run_regression.sh` checks for the *absence of a trajectory* and not
just the exit code.

**Not a missing input keyword.** `gpu_batched` (default `.true.`) and
`gpu_n_batches` (default `1`) are the only GPU-specific keywords `read_files`
parses, and all three combinations were tried before the cause was found.
Worth recording separately: with **`gpu_batched = .false.` the run
segfaults** (exit 139) at the same point. That path is untested and is
presumably rotted; the batched path is the one that works.

### What it was hiding — the exp virial is wrong, OPEN (owner: Tigany)

With the abort gone, GPU and CPU now agree on this case for everything
except the virial:

| quantity | agreement |
|---|---|
| `energy` (includes the experimental term) | exact, all printed digits |
| `energy_soap` | 2.4e-07 abs, 3.7e-10 rel |
| `energy_2b` / `energy_3b` / `core_pot` / `vdw` | exact |
| forces | **exact**, maxabsdiff 0.0 |
| `local_energy` | 3.0e-08 |
| **virial** | **maxabsdiff 1.67e+04, rel 90** |

```
ref  diag = [ 184.53,  164.51,  127.92]
gpu  diag = [13885.43, 13985.59, 16792.62]
```

The ratios are not constant (75, 85, 131) and the off-diagonals differ in
sign, so this looks like uninitialised or unreduced memory rather than a
missing scale factor.

#### Narrowed 2026-08-05

**Ruled out.** The virial formula in `get_pair_distribution_forces` is
character-identical on the two branches:

```fortran
virial(k1,k2) = virial(k1,k2) + 0.5d0*(this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
```

`this_virial_pdf`, `_sf`, `_xrd` and `_nd` are not zeroed before accumulation on
*either* branch — only `this_virial_xrd` is, and only here. Zeroing all four on
this branch was tried and changes nothing, so they are already effectively zero.
(It is still worth fixing on both branches as a latent hazard: static `real*8`
arrays relied on to start at zero, which will not hold across snapshots.)

**Where it actually is.** The CPU has no
`get_structure_factor_forces_matrix_original`. The two branches compute the
xrd/sf virial by *entirely different code paths* — this is not a port of the
CPU routine, it is a second implementation. On this branch it runs:

```
get_structure_factor_forces_matrix_original
  -> gpu_exp_force_virial_collection        (gpu_exp.cu:1903)
     -> kernel_exp_force_virial_collection  (gpu_exp.cu:1843)
  -> virial_h  ->  collect_batched_forces   (exp_interface.f90:511)
```

The kernel is:

```c
j2 = j2_list[tid]-1;
this_force = energy_scale * fi[tid ...];
atomicAdd(&forces0[j2].x, this_force.x);     // site index
tmp_xyz = xyz[tid];                          // pair index
loc_viri = 0.5*(f[k1]*xyz[k2] + f[k2]*xyz[k1]);
atomicAdd(&virial[k2+3*k1], loc_viri);
```

Forces at the site index, virial with the pair vector — the same convention the
CPU uses, so the indexing is right. The `k2+3*k1` store is transposed relative
to Fortran's column-major order, but `loc_viri` is symmetric in `k1`/`k2`, so
that is harmless.

**The deduction that should drive the next attempt.** The virial has exactly two
inputs: `this_force` and `xyz`. **`this_force` is proved correct** — the forces
this same kernel accumulates from it match the CPU to `maxabsdiff = 0.0`. So the
fault has to be in **`xyz_k_d`**, the compacted per-batch pair-vector array that
the virial reads and the forces do not.

`xyz_k_d` is built per batch from `k_index_d`, so its ordering differs from the
CPU's `xyz(1:3, j_beg:j_end)` and it cannot be compared element-wise directly.
Compare it against `xyz` gathered through the same `k_index` mapping, or check
whether it holds pair separations at all rather than absolute positions — the
magnitude ratio (13885 against 184) is about what substituting positions for
separations would give in a 20 A cell. Both branches declare and pass `virial` identically
through `calculate_pair_distribution`, so the divergence is below that, in
the GPU kernel or its reduction.

Energies and forces being exact while the virial is garbage is the same
signature as bug 5 in section 4 — a quantity that nothing in the main test
path reads staying wrong indefinitely. `XRD_mad` stays xfail until this is
fixed, and will report XPASS when it is. **Tigany is taking this one** -- it is
not blocking the refactor, which is proceeding around it.

## 6c. The GPU build is not bit-reproducible run to run — OPEN

Found 2026-08-04 while verifying an extraction, and it changes what
verification is available on this branch.

Three runs of the **same unchanged binary** on the `XRD_mad` case produce
three different `trajectory_out.xyz` files. `tools/compare_xyz.py` reports
two differing quantities between two such runs; the divergence is absent in
frame 0 and appears once MD has accumulated a step or two, at ~1e-12
relative. The CPU build on the identical input is reproducible.

Almost certainly non-deterministic reduction order in the device kernels
(atomics, or a block-order-dependent sum). It has not been chased.

**Consequence for refactoring.** The CPU branch's contract -- output must be
bit-identical to a pre-refactor baseline, compared with `diff` -- **cannot be
used here**. There is no baseline to be identical to. A refactor on this
branch can only be checked to a tolerance, which is what `compare_xyz.py`
already does and is why `tests/gpu` was built that way rather than around
`diff`.

That is weaker, and it is worth knowing exactly how much weaker: a change
that shifts results by less than the run-to-run spread is invisible here.
When bringing extractions across from the CPU branch, do them there first,
where the bit-exact contract holds, and port the verified result.

---

## 6d. The gap_backend GPU side — attempted, not landed

The CPU implementation of the seam is on the master branch (`85c1fdb`,
`33c377d`): `module gap_backend` with `gap_backend_begin`, `gap_backend_end`,
`add_2b_contribution`, `add_core_pot_contribution`, `add_3b_contribution`, and
an interface that is physics only. The GPU implementation was written against
it, builds, and is **not correct**. It is parked at
`../phase0_backup/phase2_gpu_wip/` rather than committed.

### What it found, which is worth keeping

**The three procedures had an undocumented ordering dependency.**
`add_2b_contribution_gpu` computed `st_n_sites_double` and `st_virial`, and
`add_core_pot_contribution_gpu` and `add_3b_contribution_gpu` read them through
host association without ever setting them. Nothing said so, and nothing
enforced it: calling 3b without having called 2b first would have copied back
the wrong number of bytes. Making the three independent procedures exposed it
immediately — 3b returned an energy of exactly zero until each computed its own
sizes. Any future attempt must set them per procedure.

This is the concrete cost of "internal procedures are the cheap route, host
association means no argument lists" (section 7.4). The argument list is not
the cost; the hidden coupling is.

**`alphas_d` is reused scratch, not shared state.** The driver allocates,
uploads and frees one `alphas_d` for the SOAP path and separately for the
2b/3b path — different data through one name, the pattern that hid bug 5. A
backend module owning its own removes one instance.

### What is still wrong — narrowed by dump-and-diff

The section 9 technique was applied to this: instrument both builds, dump every
scalar and every buffer checksum the 3b path consumes, and diff.

**Everything on the host is identical.** Verified equal between the working
build and the module build, on `CO_predict`:

| dumped | both builds |
|---|---|
| `i_beg`, `i_end`, `n_sites` | 1, 7176, 7176 |
| `max_np`, `size(rjs)`, `size(n_neigh)`, `size(forces,2)` | 200, 455384, 7176, 7176 |
| `size_maxnp_bytes`, `size_alphas_bytes`, `size_maxnp_qs_bytes` | 1600, 24, 4800 |
| `st_n_sites_double`, `st_virial` | 57408, 72 |
| `sp0_3b`, `sp1_3b`, `sp2_3b`, `c_do_forces` | identical per descriptor |
| checksums of `rjs`, `xyz`, `n_neigh`, `species`, `neighbor_species`, `neighbors_list` | **bit-identical** |
| per-descriptor `n_sparse`, `sum(alphas)`, `sum(cutoff)`, `sum(qs)`, `sum(sigma)`, `delta` | **bit-identical, all six descriptors** |
| `gpu_stream` | associated, same pointer value |

**It is a race, not a wrong value.** Three runs of the module build give
`sum(this_energies)` = 2.2302184918, 2.2299983140, 2.2299886355 — a spread of
~1e-4, against the working build's stable 2.1360198474. And
`sum(this_forces)` is -0.81 where it should be ~0 by Newton's third law
(the working build gives 6.4e-15).

So: identical inputs, identical parameters, identical stream, non-deterministic
wrong output. The fault is in device memory lifetime or synchronisation, not in
any value the host computes. Note this is a *different* non-determinism from
section 6c -- that one is ~1e-12 and present in the working build too; this is
1e-4 and appears only with the module.

The suspects that remain, in order:

1. **Buffer lifetime across the begin/call/end boundary.** `gap_backend_begin`
   issues async uploads and `gap_backend_end` issues async frees; in the driver
   these bracketed the three calls within one scope. If any allocation is not
   ordered against the kernels on the same stream, this is exactly what it
   looks like.
2. **`this_energies`/`this_forces`/`this_virial` as dummies.** `c_loc` is taken
   of all three. If the compiler is materialising anything for the dummy, the
   async `cpy_dtoh` targets one address and the accumulation reads another.
   Passing them as explicit-shape rather than allocatable dummies would rule
   this out.
3. Buffers that were driver variables and are now procedure locals
   (`energies_3b_d`, `forces_3b_d`, `kappas_array_d`, `sigma_d`), if anything
   frees them asynchronously after the procedure returns.

The instrumented module is kept at
`../phase0_backup/phase2_gpu_wip/gap_backend_gpu_instrumented.f90`, so the next
attempt starts from the dumps rather than rebuilding them.

### Original note on what is still wrong

With the ordering dependency fixed, `CO_predict` gives `energy_3b = 2.23018157`
against the CPU's `2.13601985` — 4.4% out. `energy_soap`, `energy_2b`,
`energy_core_pot` and forces all agree. So the fault is confined to the 3b path.

Ruled out: the lifted body is faithful (every one of the 92 code lines of the
original is present, verified line by line); `n_sites`, `st_n_sites_double`,
`st_virial` and `c_do_forces` are computed exactly as the originals did; and
`max_np`, `size_maxnp_bytes`, `size_maxnp_qs_bytes` and `size_alphas_bytes` are
all written before they are read inside the body.

Not yet ruled out: something else the 3b body inherited from the driver's scope
that a read-before-write scan does not catch — most likely a value that is
merely *stale rather than undefined*, so it looks initialised and is simply
wrong. `forces`, passed only so the body can take `size(forces,2)`, is the most
suspicious remaining argument.

### Retried after 6e was fixed — still fails, but now visible to memcheck

The use-after-free fix did not resolve this. `energy_3b` moved from 2.2302 to
0.5099 against the CPU's 2.1360 — a different wrong answer, which is consistent
with the fault being layout-sensitive rather than caused by 6e.

What changed is that it is now **diagnosable**. With the SOAP noise gone,
memcheck on the backend build reports:

```
Invalid __global__ atomic of size 8 bytes
    at kernel_2nd_try<(bool)1, (kern_type)1>+0x4cb0
    Access to 0x7da2c6fe8 is out of bounds
    and is 58,392 bytes before the nearest allocation
    Host Frame: gpu_3b_cc_2nd_try -> gpu_3b -> __gap_backend_MOD_add_3b_contribution
```

68 errors, 9 of them naming the 3b kernel. **Neither the original build nor the
6e-fixed original has a single one** — checked by grepping all three sanitizer
runs for `2nd_try`: 0, 0, 9. Full output at
`../phase0_backup/memcheck_backend_3b.txt`.

So the extraction feeds `kernel_2nd_try` an index that sends its atomics below
the target buffer. The atomics are how the kernel accumulates energy, which is
why the energy comes out low rather than merely different. The indices come
from `kappas_array_d`, so that is where to look — note the host `kappas` array
was already proved identical, so it is the device copy or the offsets derived
from it, not the values themselves.

Ruled out on the retry: a race between the async `kappas` upload and
`deallocate(kappas)` — adding `gpu_stream_sync` before the deallocate changes
nothing (0.5099 vs 0.5199, within the run-to-run spread).

The version used for this retry is kept at
`../phase0_backup/phase2_gpu_wip/gap_backend_gpu_v2.f90`.

#### The indexing was then checked and is not the cause

`kernel_2nd_try` derives everything from the site index and `kappas_array`:

```c
const int i = blockIdx.x + i_beg - 1;      // absolute site index
if (species[i] != sp0) return;
int k = kappas_array[i];
auto i3 = ((neighbors_list[k]-1) % n_sites0);
```

so a wrong `k` gives a garbage `neighbors_list[k]`, a negative `i3`, and an
atomic below the buffer — which is what memcheck reports. **`kappas_array_d`
was therefore read back from the device in both builds and is identical**:
`sum = 1633571860`, `[0] = 0`, `[1] = 61`, `[n_sites-1] = 455312`, matching the
host array exactly.

So the exclusion list is now: host scalars, upload byte counts, host array
checksums, per-descriptor parameters, the stream, **and the device copy of
`kappas_array_d`**. Every input the kernel reads has been proved equal, and
only the module build's atomics leave the buffer.

What that leaves is the *output* pointers — `energies_3b_d`, `forces_3b_d`,
`virial_3b_d`, allocated inside `add_3b_contribution` from `size(n_neigh)` and
`size(forces,2)` (both 7176, both verified). Either one of those allocations is
not what the kernel receives, or the allocation itself is landing somewhere the
pool later reuses. The next step is to read back the three device pointer values
and their allocated sizes on both sides, rather than their contents — that is
the one thing not yet compared.

#### Every argument to gpu_3b is now proved identical — stop comparing inputs

The argument-by-argument comparison is **complete**. Between the working build
and the module build, on `CO_predict`, all of the following are equal:

| compared | result |
|---|---|
| `n_sparse`, `n_sites`, `n_atom_pairs`, `n_sites0` | identical |
| `sp0_3b`, `sp1_3b`, `sp2_3b`, `c_do_forces`, `i_beg`, `i_end` | identical |
| `delta`, `rcut`, `max_np`, all four `size_*` byte counts | identical |
| `c_name_3b` (selects the kernel template) | identical, `"pol "` |
| host checksums: `rjs`, `xyz`, `n_neigh`, `species`, `neighbor_species`, `neighbors_list` | identical |
| per-descriptor host `alphas`, `cutoff`, `qs`, `sigma` | identical, all six |
| **all 11 input device pointers** | identical *addresses* |
| **all 3 output device pointers** and their sizes | identical addresses, 57408 / 172224 / 72 |
| **device contents** of `kappas_array_d` | identical |
| **device contents** of `alphas_d`, `cutoff_d`, `qs_d`, `sigma_d` | identical |
| `gpu_stream` | identical pointer |

The kernel receives byte-identical arguments pointing at byte-identical device
memory, and still writes 59 out-of-bounds atomics that the working build does
not. **Do not re-check any of the above** — it has been measured, not reasoned
about.

That means the cause is not in what the kernel is passed, so input comparison
has been exhausted. The productive next step is a **structural bisect**: start
from the module version and move it back toward the original one change at a
time — fold `gap_backend_begin`'s uploads into `add_2b_contribution` to remove
the begin/end split; make the device pointers driver variables passed as
arguments instead of module state; make the three procedures internal again but
in a module. Whichever step restores 2.1360198474 names the cause.

The instrumented versions used for all of the above are kept in
`../phase0_backup/phase2_gpu_wip/`.

Worth noting for whoever picks this up: the commented-out line directly above
the live one,

```c
//int k = i_beg != 1 ?  kappas_array[i-i_beg+1] : kappas_array[i] ;
int k =  kappas_array[i];
```

means the kernel only indexes `kappas_array` correctly when `i_beg == 1`. Every
test here is single-rank, so it has never been exercised otherwise. That is a
separate latent bug and will bite any multi-rank 3b run.

### Earlier notes

The obvious next move is a bisect rather than more reading: instrument the
working build to dump every scalar the 3b body reads before writing, run both,
and diff. That is the technique from section 9, and it is what found bug 5.

---

## 6e. Use-after-free in the SOAP derivative path — FIXED

**`compute-sanitizer memcheck` reports 181 errors on `CO_predict` with the
unmodified build at HEAD.** Every one is a *use-after-free*, on an allocation
of 1,797,120 bytes, reached from two kernels:

```
Use-after-free on allocation of size 1,797,120 bytes at 0x7dc3a6200
    Address 0x7dc3a6200 is potentially accessed after it is free'd
    Host Frame: cuda_get_soap_der_one(...)
      gpu_get_soap_der -> get_soap -> get_gap_soap -> MAIN__
```

and `cuda_local_property_derivatives`. Nothing in the 3b path appears at all.
Full output kept at `../phase0_backup/memcheck_CO_predict.txt`.

Reproduce:

```sh
/usr/local/cuda-13.3/bin/compute-sanitizer --tool memcheck bin/turbogap predict
```

**Why this is not cosmetic.** The kernels read device memory that has been
freed. What that memory *contains* depends entirely on what has been allocated
into it since — so the answer the code produces is a function of the device
allocation layout, not only of the physics. The build passes its regression
cases today because the layout happens to leave the right bytes there.

This is the same class as bug 2 in section 4, which was a 3b device buffer
freed inside the descriptor loop instead of after it. That one was found because
it crashed. This one does not crash; it silently returns the right answer for
the current layout.

**Cause and fix.** `cuda_malloc_all` is `hipMallocAsync` on `gpu_stream`, so
every allocation in this code is *stream-ordered* and must be released with
`hipFreeAsync` on the same stream. Three sites released theirs with the plain
`gpu_free`, which is `hipFree`:

| freed in | pointer |
|---|---|
| `soap_turbo.f90` `get_soap` | `cnk_d` |
| `gap.f90` `get_soap_energy_and_forces` | `energies_d` |
| `gap_interface.f90` `get_gap_soap` | `soap_cart_der_d` |

Each matches one of the memcheck backtraces exactly, and each sits directly
below lines that already use `gpu_free_async(..., gpu_stream)` — so these were
inconsistencies rather than deliberate choices. Changing the three to
`gpu_free_async` takes memcheck from **181 errors to 1**, and the one that
remains is inside cuBLASLt loading its own kernels, not TurboGAP code. Both
regression cases still pass.

Note the `soap_turbo.f90` fix is in the submodule (`d02eeea`), so the parent's
submodule pointer moves with it.

**It was also what defeated the gap_backend GPU side (6d).** That work changes which
buffers are allocated, in what order and how large, and therefore what sits in
the freed region. Every kernel input was proved bit-identical between the
working and module builds -- scalars, upload byte counts, array checksums,
per-descriptor parameters, `kappas`, the stream -- and the device still returned
a different 3b energy, non-deterministically. Reading freed memory explains
exactly that: identical inputs, layout-dependent output.

So 6d is not a defect in the extraction. **Fix this first, then retry 6d**;
until then any change to device allocation on this branch can move results for
reasons that have nothing to do with the change.

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
