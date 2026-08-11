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

### What it was hiding — the exp virial, FIXED, and it was the CPU that was wrong

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

**Resolved 2026-08-05, and the fault was on the CPU side.**

`get_structure_factor_forces_matrix` in the CPU branch's `exp_utils.f90`
accumulates the virial over the `n_k` selected pairs but indexed `xyz` with
`k` rather than the loop variable `j`. `k` is not live in that loop — it is
left over from the selection loop that ends about sixty lines earlier — so
every one of the `n_k` virial terms used the same fixed pair vector. The
forces in the same loop are correct, using `fi(j,:)` and `j2_list(j)`, which
is why forces agreed exactly while the virial did not.

The routine already builds the compacted per-pair array for exactly this
purpose: line 586 fills `xyz_k(1:3, n_k) = xyz(1:3, k)` inside the selection
loop and line 598 uses it correctly for `Gka`. Only the virial reached past it
to the raw array.

**This branch was right all along.** The device kernel indexes the compacted
array by the selection index, which is what the CPU should have done. Fixed on
the CPU branch (`769c815`); on frame 0, at identical positions, the two virials
now agree to 1.6e-07 absolute and 9.5e-12 relative, against 1.67e+04 before.

The two CPU cases with `exp_forces = .true.` converted to `REFERENCE=golden`
there, since the baseline binary carries the bug.

#### How it was found — the readbacks that led there

`xyz_k_d` checked and correct, which disproved the first hypothesis

`xyz_k_d` was read back from the device during an `XRD_mad` run:

```
XYZK  nk=7348   min -5.99055597   max +5.99055597   mean|v| 2.26   max|v| 5.99
XYZK  nk=3840   min -5.96122864   max +5.96122864   mean|v| 2.25   max|v| 5.96
XYZK  nk=498    min -5.84727751   max +5.84727751   mean|v| 2.26   max|v| 5.85
```

Signed, symmetric about zero, and bounded by the pdf `rcut` of 6.0 — these are
**pair separations, correctly compacted**, not absolute positions (which would
span 0 to 20 in this cell and be mostly positive). The magnitude coincidence
that motivated the guess was just that.

Also checked and identical: master computes the xrd/sf virial with the *same*
formula, in `get_structure_factor_forces_matrix`:

```fortran
virial(k1,k2) = virial(k1,k2) + 0.5d0*(this_force(k1)*xyz(k2,k) + this_force(k2)*xyz(k1,k))
```

**So both of the virial's two inputs are now proved correct on this branch**:
`this_force`, because the per-atom forces this same kernel accumulates from it
match the CPU to `maxabsdiff = 0.0`; and `xyz_k_d`, by the readback above. The
formula matches. The kernel's index convention matches (forces by site, virial
by pair) and its transposed store is harmless by symmetry.

What is left, and it is a different shape of question from anything tried so
far: the two implementations must be accumulating over **different sets of
pairs, or a different number of times**. Per-atom forces would not reveal that
if the extra contributions cancel in the force sum but not in the virial, which
is exactly what an over-count of symmetric pair terms looks like. There are
three call sites of `gpu_exp_force_virial_collection` in `exp_utils.f90` (lines
600, 887, 1008); the next step is to count how many pair terms each branch
accumulates into the virial for one descriptor and compare the counts, not the
values. Both branches declare and pass `virial` identically
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

## 6f. energy_core_pot diverges under MAD — FIXED

Found while confirming the virial fix. With `XRD_mad` run from identical
starting positions, frame 0 agrees on everything including `energy_core_pot`
(both 0.0). From frame 1 the CPU reports `energy_core_pot = 8.53` and this
branch reports `0.00`, and the trajectories then separate.

Frame 0 forces agree to `maxabsdiff = 0.0`, so frame 1 positions must agree
too — which means this is a real disagreement in the core-potential term at
that geometry, not accumulated drift.

`CO_predict` and `CO_md` do not catch it because `energy_core_pot` is 0.0
throughout those runs: the CO structure never brings a pair close enough to
enter the core potential. `XRD_mad` does, because the MAD forces push atoms
together.

**Cause and fix (`a7fe298`).** `add_core_pot_contribution` passed `sp1` and
`sp2` to the kernel without ever setting them. In the driver they were
host-associated, so core_pot silently used whatever species indices
`add_2b_contribution_gpu` had last written — the final 2b descriptor's, not
its own. **That is bug 3 of section 4, fixed for 2b and never for core_pot.**
Moving the procedures into a module turned the stale value into an undefined
one, which is how it finally surfaced.

The kernel filters neighbours by species index, so with the wrong indices no
pair matches and the term comes out as exactly 0. It now derives them from
`core_pot_hypers(i)%species1` and `%species2`, as the CPU does inside
`get_core_pot_energy_and_forces`. Verified against the CPU on `XRD_mad`:

| frame | CPU | GPU |
|---|---|---|
| 2 | 8.53223629 | 8.53223629 (exact) |
| 3 | 82.60515055 | 82.60515056 |
| 4 | 2512.07875798 | 2512.07875785 |

---

## 6g. The exp virial diverged from frame 1 — FIXED, again on the CPU side

What remains after 6b and 6f. Running `XRD_mad` from identical positions with
both fixes in place:

* **frame 0** — everything agrees, including the virial (1.6e-07, rel 9.5e-12)
* **frame 1 onward** — energies, `energy_core_pot`, `local_energy` and
  **forces to `maxabsdiff = 0.0`** all agree; **only the virial differs**

```
frame 1  CPU [14976.64, -14449.08,  2522.62, ...]
         GPU [ 1275.74,  -8742.90,  -550.61, ...]
```

**Cause and fix.** Comparing the pdf, sf and xrd virials *separately* rather
than the total settled it immediately:

| | frame 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| CPU `sum(abs(virial_xrd))` | 72549 | 103851 | 168265 | 192035 | 289353 | 226333 |
| this branch | 72549 | 35239 | 70441 | 45745 | 114767 | 104029 |

The CPU's grows; this branch's fluctuates per frame. The CPU driver zeroes
`virial_pdf/sf/xrd/nd` each prediction block, but under MPI the exp routines
accumulate into the `this_` prefixed copies, **which were never zeroed
anywhere**, so the structure-factor virial accumulated across every MD step of
a run. This branch zeroes `this_virial_xrd` before its batched collection and
was right. Fixed on the CPU branch (`4f99e92`); the two now agree frame by
frame to ~1e-10.

With 6b, 6f and this fixed, `XRD_mad` agrees on energies, forces, core_pot and
the virial in every frame. The only remaining difference is `local_energy`
drift in frames 4 and 5 — 2.4e-06 and 1.5e-05 — which is accumulated
divergence over five MD steps of a chaotic driver, given this build is not
bit-reproducible run to run (§6c). It is not a defect, but it is above the
comparator's tolerance, so `XRD_mad` stays xfail on that alone.

**Recorded as its own item so it is not lost behind the xfail marker:** the
`local_energy` drift is 2.4e-06 at frame 4 and 1.5e-05 at frame 5, growing
with step count. Everything else — energies, forces, `energy_core_pot`, the
virial — agrees in every frame. Before treating this as tolerable, note that
§6c says this build is not bit-reproducible run to run at ~1e-12, and five MAD
steps amplifying that to 1e-05 is plausible but has not been demonstrated. The
way to settle it is to run the same GPU binary twice and see whether frame 5
`local_energy` moves by a comparable amount; if it does, raise the comparator
tolerance for this case and the xfail can go.

---

## 6h. Spurious pdf, sf and nd virial contributions — OPEN, small

Noticed while separating the contributions. With `n_exp = 1` and
`exp_labels = 'xrd'`, only the xrd term is fitted, so only it should carry
forces and a virial. The CPU reports `sum(abs(...))` of exactly 0.0 for pdf,
sf and nd in every frame. This branch reports:

```
pdf 6.126893924420E-03    sf 1.472498760758E+02    nd 1.007955953266E+00
```

**identical to all thirteen digits in all six frames**, while the atoms move.
A geometry-independent value is not physics. Zeroing `this_virial_pdf/sf/nd`
per snapshot does not change them, so they are written after that point, every
frame, to the same value.

**What the correct behaviour is** (confirmed with Tigany, 2026-08-05): when
only `xrd` is given as an experimental observable, the pdf, sf and nd virials
**should not be computed at all, and must not enter the total virial**. Only
the fitted observable carries forces and a virial. So this is not a matter of
getting the numbers to agree with the CPU — the contributions should not be
there in the first place, and the fix is to stop computing and adding them
rather than to correct their values.

`sf` at 147 against an xrd virial of order 1e5 is small but not negligible,
and it is in the reported total today, so any pressure from the exp path on
this branch is wrong by that amount.


## 6i. Out-of-bounds write in `get_gap_soap` — FIXED

`gap_interface.f90` sizes the SOAP neighbour arrays in a counting pass and fills
them in a second pass. The two passes disagreed:

```
line 265 (count):  if (rjs0(k) <  rcut_max .and. species_multiplicity_supercell(j2) > 0)
line 340 (fill) :  if (rjs0(k) <= rcut_max .and. species_multiplicity_supercell(j3) > 0)
```

A pair sitting exactly on `rcut_max` was not counted but was written, so the
fill ran one element past the end of `rjs`, `phis`, `thetas`, `xyz`,
`in_to_out_pairs` and `mask`. Master has the strict comparison at both sites;
the `<=` was a divergence on this branch.

**Why it survived this long, which is the part worth keeping.** On the
7176-atom CO system the overrun lands inside memory the allocator has already
handed out. It corrupts the heap silently instead of crashing, and every test
on this branch used that one system. It only becomes visible on a small cell,
where the overrun crosses the top chunk: the P4 dimer from `vdw_P` aborts with

```
malloc(): corrupted top size
```

before printing a single line, while the CPU build runs the same input fine. A
`-fcheck=bounds` build named it exactly — *"Index '41' of dimension 1 of array
'rjs' above upper bound of 40"*. That is the tool of choice here, per §3:
`ptrace_scope = 2` means gdb cannot attach.

Two lessons:

* **A suite of one system size cannot see this class of bug.** The CO cases pass
  identically before and after the fix — they never place a pair exactly on the
  cutoff. Coverage needs a small cell as well as a large one.
* It is the same shape as `refactor_handoff.md` §5.1 and §7.1: **one condition
  written at more than one site**, with nothing forcing the copies to agree.
  Three separate defects on these branches now have that shape.

## 6j. `local_properties` is unallocated on the vdW path — FIXED

With 6i fixed, the P4 dimer gets further and then dies at `turbogap.f90:1853`:

```
Index '-1101217936' of dimension 2 of array 'local_properties'
   below lower bound of 4391641174017424937
```

Nonsense on both sides of the comparison, which is what an unallocated array
descriptor reads as. The line is

```fortran
v_neigh_vdw(k) = local_properties(j2, vdw_lp_index)
```

in the TS block, and that block is entered on `any(soap_turbo_hypers(:)%has_vdw)`
while `local_properties` is allocated a few hundred lines earlier under
`any(soap_turbo_hypers(:)%has_local_properties)`. The two predicates are not the
same, and nothing checks the second before the first is acted on.

What has been ruled out: the allocation block itself is character-for-character
identical to master's, as is the `vdw_lp_index` assignment in `read_files.f90`.
So this is not a transplanted-code problem, and it predates the vdW adoption —
the same abort happens with the branch's own old `vdw.f90`.

**Consequence: the vdW path on this branch has never run.** Master's `vdw.f90`
was adopted wholesale and is byte-identical to master's copy, so it is correct
by construction, but nothing here exercises it. Fixing this is what makes a vdW
regression case possible, and that case is what would let the branch claim vdW
works at all.

Owner decision needed on the fix: guard the TS block on `has_local_properties`
as well, or make the allocation unconditional when `has_vdw`. Master gets away
without either because its `compute_vdw` is reached on a path where the local
properties have already been computed.

**Fixed (`389cf76`).** The cause was in the GAP-file parsing, not the driver.

A GAP file may declare the Hirshfeld model in the deprecated inline form
(`has_vdw = .true.` plus `vdw_qs`, `vdw_alphas`, `vdw_zeta`, `vdw_delta`,
`vdw_v0`), and `vdw_P/gap_files/phosphorus.gap` does. Master reads `has_vdw`
into `%has_local_properties` and synthesises a `local_property_models` entry
labelled `hirshfeld_v`; this branch read it into the legacy `%has_vdw` field
and never created the model. So `has_vdw` was true, `has_local_properties` was
false, and the TS block ran against an array the allocation guard had skipped.

Master's migration is ported verbatim, with `check_deprecated` and
`print_deprecation_message`, so a deck using the old form is told what to
replace instead of silently taking a different path.

The P4 dimer now runs and agrees with the CPU build: `energy`, `energy_vdw` and
`energy_2b` identical, `energy_soap` differing at 1e-8. It is now the `vdw_ts`
case in `tests/gpu/run_regression.sh`.

Two follow-ons landed with it:

* **Thirteen vdW/MBD input parameters** master has and this branch did not
  (`vdw_mbd_rcut`, `vdw_2b_rcut`, `vdw_omega_ref`, `vdw_sr_mbd`, `do_nnls`, …).
  Master's `vdw.f90` is compiled here now, so those paths exist; without the
  keywords a deck setting them was *silently ignored*, which is worse than
  failing. The two trees now recognise the same keyword set, with the eleven
  `estat_*` keywords the deliberate exception — that is this branch's feature,
  owed to master.
* The suite gained its first **small cell**. Everything else here is the
  7176-atom CO system, which is exactly why §6i went unnoticed.

---

## 6c revisited. The run-to-run irreproducibility, diagnosed — two atomicAdds

§6c recorded that the same binary gives different output run to run and that it
"has not been chased". §6g then hung the `XRD_mad` xfail on a `local_energy`
drift and asked, as the way to settle it, whether the same binary moves frame 5
by a comparable amount on its own. It does. This is that chase, finished.

**The drift is not a GPU-vs-CPU disagreement.** It is the GPU build disagreeing
with itself, and the two places it comes from are now known.

### Method

Three ablations, each run three times, comparing `md5sum` of
`trajectory_out.xyz` over the five `XRD_mad` MD steps. The potential ablations
are `delta = 0.0` on whole `gap_beg` blocks, which zeroes those contributions
exactly and so removes their summation-order noise without touching any code.
The kernel ablation is a build in which `cuda_soap_forces_virial_two` is
launched as a single block looping the pairs in index order, which makes its
scatter sequential.

| build | deck | reproducible |
|---|---|---|
| unmodified | soap only, no exp | **no** |
| serialised soap scatter | soap only, no exp | yes |
| serialised soap scatter | soap + 2b, no exp | yes |
| serialised soap scatter | soap + 2b + 3b, no exp | yes |
| serialised soap scatter | soap + 2b + 3b, with exp | **no** |

A single-point `predict` is bit-identical in every configuration. The
divergence needs an MD step to become visible: the noise is at the 1e-16 level
per force component, which the `F16.8` writer cannot show, and one integration
step turns it into ~1e-12 on the velocities, which it can.

### The two sources

Both are the same shape — a per-pair contribution scattered into a per-atom
force with a double `atomicAdd`, whose completion order is not fixed:

1. `src/gpu/gap_soap_forces.cu:63-65`, `cuda_soap_forces_virial_two`. One block
   per pair; thread 0 adds the block's reduced 3-vector into `forces_d[j2]`.
   Many blocks share a `j2`. The virial at `:92` is the same problem with all
   blocks on nine addresses, worse but harmless here — nothing integrates the
   virial without a barostat.
2. `src/gpu/mad_xrd.cu:178-180` and `:203`,
   `kernel_exp_force_virial_collection`. Identical construction on the exp
   forces, which is why the last row of the table is still red with the soap
   scatter serialised.

`gap_2b.cu` and `gap_3b.cc` also scatter with `atomicAdd`, `gap_3b.cc:243` and
`:315` doing it on `energies[i]` directly. They do not show above: 2b and 3b
carry `delta = 0.5` and `0.01` against soap's `0.1` on a much larger descriptor,
so their last-bit noise never reaches the printed digits here. They are not
deterministic either, and a deck weighted towards them would show it.

### What a fix looks like, and what it costs

Serialising is a diagnosis, not a fix — it is O(n_pairs) sequential in the
hottest kernel on the branch. Two real options:

* **Deterministic gather.** Have each block *store* its 3-vector at
  `pair_force_d[l]` — a plain store, no ordering question — and add a second
  pass that sums each site's incoming pairs in increasing pair index. That
  needs the inverse map (site → its incoming pairs), which is a counting sort
  of `j2_index_d`: histogram with integer atomics (exact), exclusive scan
  (`gpu_scan.cu` is already there), scatter into buckets, then rank each bucket
  by pair index so its order stops depending on which thread got there first.
  Exact — no change to any number — and one reusable helper serves both call
  sites, which have the same `j2_index`/`pair` shape. Perhaps 120 lines of
  CUDA. The extra pass is small next to the `n_soap` contraction each block
  already does, but it wants measuring before it lands in the SOAP path.
* **Fixed-point accumulation.** `atomicAdd` on `unsigned long long` at a fixed
  scale is associative, so any order gives the same answer. About thirty lines.
  It changes every force by the quantisation — at a scale of 2^34 that is
  ~6e-13 relative, comparable to the noise it removes — and it introduces a
  silent overflow mode at ~5e8 eV/A, which is a bad failure to add to a force
  path even if no physical configuration reaches it.

The first is the right one; neither is landed. Until one is, `tests/gpu`
cannot certify anything at a tolerance below the run-to-run spread, which is
§6c's point and still stands.
