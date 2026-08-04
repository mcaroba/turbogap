# Handoff: `refactor/modularize` branch (turbogap.f90 de-monolithing)

Written 2026-08-04. Everything below was verified on `alt` (`d22-0052`);
nothing has been pushed to any remote.

Companion document: `docs/gpu_fixes_handoff.md` in the GPU repo, which covers
the GPU correctness work this branch is eventually meant to merge with.

---

## 1. Current state

`src/turbogap.f90` is **3936 lines, down from 4992** at the branch point
(−21%). Two phases have been lifted into modules, and a bit-exact regression
suite exists to make further lifting safe.

| | lines |
|---|---|
| `src/turbogap.f90` | 3936 (was 4992) |
| `src/turbogap_vdw.f90` | 691 (new) |
| `src/turbogap_setup.f90` | 593 (new) |

**Branch:** `refactor/modularize`, 6 commits on top of `e6eb1aa` (master HEAD).

```
b7781bb Extract the input/GAP-file reading phase out of turbogap.f90
113d0b6 Apply the ts+mbd correction to forces as well as energies
535b1d9 Fix ts+mbd outside MD: recompute when there is no previous step
a4e4ca0 Extract the vdW/MBD phase out of turbogap.f90
852d41f Cover the ts+mbd correction state before refactoring it
27091a7 Add bit-exact regression suite for the turbogap.f90 refactor
```

All 15 regression cases pass. Four of the six commits are pure refactors and
are bit-identical to the pre-refactor binary; two deliberately change `ts+mbd`
results and say so (§5).

---

## 2. The strategy, and why it is this one

The goal is not tidiness. It is to turn **one intractable merge into many
tractable ones**. Three findings drive everything:

1. **The GPU branch's `turbogap.f90` has no CUDA `#ifdef` at all.** 73
   `#ifdef _MPIF90`, zero for CUDA, against 288 uses of
   `gpu_*`/`cuda_*`/`c_loc`. That file cannot compile without CUDA, so a single
   merged `turbogap.f90` is impossible by sprinkling `#ifdef`. It needs a
   **backend interface with build-time selection** — one module name, two
   implementation files, chosen in the Makefile. The seam is already latent:
   `get_gap_soap` exists in both trees with the same name and different bodies.

2. **Both `turbogap.f90` files have identical top-level section banners in
   identical order.** The GPU branch never restructured the program; it rewrote
   block interiors in place. So a common skeleton is achievable, and it is what
   makes per-phase merging possible. `git diff` between the two files is
   currently 9591 lines over 27 hunks spanning the whole file — hopeless as a
   textual merge, tractable as ~10 per-phase merges.

3. **The branches diverged 2024-07-09** (merge-base `d3e06de`); master is +42
   commits, GPU +22. Master gained the entire MBD rewrite (`vdw.f90` 495 →
   6538 lines). The GPU branch touched `vdw.f90`, `md.f90`, `mc.f90` and
   `splines.f90` **not at all**, so those merge for free. The real conflict
   surface is `turbogap.f90` plus `exp_interface`, `exp_utils`, `gap`,
   `gap_interface`, `read_files`, `types`, `neighbors`, `local_properties`.

The consequence for sequencing: refactor master into a skeleton **designed so
the GPU interiors drop into it**, not into whatever reads nicest. Leave
`gpu_fixes` untouched until the skeleton exists.

---

## 3. Where things are

| what | path on `alt` |
|---|---|
| this branch | `/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_master_2026` |
| GPU branch (untouched) | `/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/turbogap_gpu_commit_mahti` |
| test data root | `/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/input` |
| regression suite | `<repo>/tests/regression/` |

Test data, ~35 MB, deliberately outside the repo:

- `input/CO/` — pre-existing; 7176-atom CO, also used by the GPU comparison, so
  timing numbers stay comparable with `gpu_fixes_handoff.md`
- `input/xps_opt/` — copied from `tutorials/xps_optimization`; 512-atom CO with
  core-electron binding energies. Covers XPS, GCMC and local properties.
- `input/vdw_P/` — copied from `TurboGAP_School/vdw_interactions`; P4 dimer and
  108-atom cell. The only data available that exercises TS/MBD, because the CO
  and C potentials are the only ones with a `hirshfeld_v` local property and
  neither ships an MBD-ready input deck.

Build:

```sh
cd <repo>
make clean all          # mpif90 -fPIC -O3, ~44 s
```

---

## 4. Testing

### 4.1 Running

```sh
TURBOGAP_DATA_ROOT=/u/74/zarrout1/unix/work/cpu_vs_gpu_tests/input \
  tests/regression/run.sh            # all 15 cases, ~4 min
tests/regression/run.sh vdw_mbd      # one case
tests/regression/run.sh --list
```

Read `tests/regression/README.md` and `KNOWN_ISSUES.md` before trusting a green
run.

### 4.2 The contract

**Output must be bit-identical to the pre-refactor baseline.** Extracting a
block into a procedure does not reorder floating-point operations, so not one
digit should move. That is why the comparator is `diff` and not a tolerance
check: exact, no tuning, and impossible to quietly loosen as the refactor
proceeds. A change that legitimately alters results is not a refactor and does
not belong on a commit checked this way — it gets its own commit, and its cases
convert to `REFERENCE=golden`.

The baseline binary is **gitignored, not committed** — it is
architecture-specific, and a stale one committed by accident would make the
suite pass vacuously. `tests/regression/make_baseline.sh` rebuilds it from the
SHAs in `baseline/COMMIT` via a detached worktree.

### 4.3 Two kinds of case

Ten cases compare against the baseline binary. Five `ts+mbd` cases use
`REFERENCE=golden` and compare against a checked-in `expected/` directory,
because the baseline cannot serve as a reference for them — it either crashes
on the input or the branch deliberately changed the result. Those are
**characterization** tests: green means "unchanged", not "correct". Regenerate
with `TURBOGAP_BLESS=1`, deliberately, never to make a red suite go green.

Expected outputs are stored **per rank count**: MPI reductions are not
associative, and 1 vs 2 ranks legitimately differ by ~1e-8 in forces.

### 4.4 Coverage

| case | ranks | covers |
|---|---|---|
| `vdw_ts` | 1 | TS pairwise vdW, Hirshfeld local properties |
| `vdw_mbd` | 1 | many-body dispersion |
| `vdw_mbd_cell_mpi2` | 2 | MBD + cross-rank Hirshfeld gradient exchange |
| `vdw_tsmbd_md` / `_mpi2` | 1, 2 | ts+mbd, both the recompute and reuse branch (golden) |
| `vdw_tsmbd_predict` / `_mpi2` | 1, 2 | ts+mbd outside MD (golden) |
| `vdw_tsmbd_mc` | 1 | ts+mbd under mc (golden) |
| `xps_predict` / `_mpi4` | 1, 4 | SOAP loop, local properties, 2b/3b/core_pot, XPS |
| `gcmc_xps` / `_mpi2` | 1, 2 | MC insert/remove: `n_sites` changes, reallocation, rebuilds |
| `co_predict` / `_mpi4` | 1, 4 | 7176-atom single point; **timing reference** |
| `co_md` | 1 | velocity Verlet, Berendsen thermostat, rebuild decision |

Every case fixes `random_seed`; `co_md` starts from explicit velocities. All
were verified reproducible run-to-run *before* being relied on: baseline vs an
identical rebuild passed all of them.

**Known gaps.** No coverage of nested sampling, electronic stopping / eph, box
scaling or NPT. Adding a nested-sampling case is worth doing before touching
that block (§7).

---

## 5. Two vdW bugs found and fixed

> **Correction (2026-08-04, after this document was written).** Only §5.1 is a
> bug. §5.2 was not: the uncorrected forces are correct, and `113d0b6` has been
> reverted in `59997e3`. Read §5.2 as a record of a wrong conclusion, and see
> `tests/regression/KNOWN_ISSUES.md` issue 2 for why the measurement it argues
> from does not establish what it claims. Everything in §5.1 stands.

Both were found while building the safety net, not by reading the code.

### 5.1 `ts+mbd` segfaulted outside MD (`535b1d9`)

`predict` and `mc` died in `this_energies_vdw_corr`. The scheme runs full MBD
every `mbd_correction_freq` steps and reuses a stored correction in between,
and which half of the cycle a call took was decided by
`modulo(md_istep, mbd_correction_freq) == 0` — **written out at six separate
sites that all had to agree**. `md_istep` stays `-1` outside MD, and
`modulo(-1, 100)` is 99, so all six chose the *reuse* half, including the site
that allocates the buffers, which only the *recompute* half reaches. A quieter
second instance: `mbd_ts_scaling` was initialised under `md_istep == 0`, so
outside MD it was read uninitialised.

Fix: the reuse half only means anything during MD, where there was a previous
step to have stored something. Outside MD every call recomputes. One named
predicate now replaces the six copies:

```fortran
is_correction_step = ( md_istep < 0 ) .or. &
                     ( modulo(md_istep, params%mbd_correction_freq) == 0 )
```

For `md_istep >= 0` this is exactly the old expression, so MD was provably
unchanged — verified at the time by the MD cases staying bit-identical.

### 5.2 `ts+mbd` corrected energies but not forces (`113d0b6`) — REVERTED, NOT A BUG

It returned MBD energies with TS-scaled forces. Uncommenting the dormant line
would not have fixed it: `this_forces_vdw_corr` was snapshotted as the *TS*
result and never updated to the MBD one, and the step forming `MBD - TS_scaled`
did not exist for forces at all. The fix mirrors the energy path at all three
points.

Measured on the 108-atom P4 cell against `mbd_correction_freq = 1`:

| | first reuse step | after 10 reuse steps |
|---|---|---|
| without correction | 1.9e-05 (0.33% of \|F\|max) | — |
| with correction | 9.0e-08 (0.00%) | 8.8e-06 (0.15%) |

~200x better at the first reuse step; corrected forces after ten reuse steps
still beat uncorrected forces after one. The cost argument the scheme exists
for, same system, 10 MD steps: **102.0 s at freq=1 vs 8.6 s at freq=1000, 11.9x**.

**Price paid:** this changes `ts+mbd` results, so `vdw_tsmbd_md` and
`vdw_tsmbd_md_mpi2` lost their baseline comparison and became golden. That is a
real reduction in safety for those two, called out in their `case.conf`. It is
contained — those five cases were the only ones that moved.

---

## 6. Environment gotchas

Each of these cost real time.

* **`get_ts_energy_and_forces` assigns, it does not accumulate**, despite
  `energies`/`forces0` being `intent(inout)`. It zeroes them internally
  (`vdw.f90:318`, `343-345`). The whole `ts+mbd` algebra is unreadable until
  you know this — it is what makes the final result land exactly on MBD.
* **Overriding `F90_OPTS` on the `make` command line drops `-J include`**,
  because the Makefile *appends* the module-dir flag to `F90_OPTS`. The result
  is `.mod` files littered across the repo root. Use a makefile fragment, or
  clean up after.
* **One-shot timing is noisy on the short cases.** `vdw_mbd` gave 2.01 / 2.01 /
  3.34 s across repeats of the *same binary*. Use `co_predict` (~18 s) as the
  timing reference and ignore sub-5-second ratios.
* **`ptrace_scope = 2`** on `alt`, so `gdb` cannot attach. `-fcheck=bounds`
  found both vdW bugs immediately; that is the tool of choice here.
* **The negative "Miscellaneous" time** in the end-of-run summary is a
  double-counting bug in the timing buckets. It makes the printed totals
  useless for judging refactor performance — use the wall-clock `run.sh`
  reports. Still open, cosmetic.
* `numpy` is not installed on `alt`. `run.sh` is bash + `diff` only, and
  `tools/compare_xyz.py` in the GPU repo is deliberately stdlib-only. Keep it
  that way.
* **`src/soap_turbo` is still dirty** with the `TURBOGAP_DUMP_CNK`
  instrumentation from the GPU debugging session. It is env-var gated and has
  zero effect when unset, so it cannot affect the bit-exactness comparisons,
  and it was left alone deliberately. Still item §7.3 of the GPU handoff:
  commit or revert it.

---

## 7. What is left

In the order the measurements support.

### 7.1 Bundle the contribution families — do this next

The driver holds ten additive contributions (soap, vdw, lp, pdf, sf, xrd, nd,
2b, core_pot, 3b), each with `energies_X`, `forces_X`, `virial_X` and a
`this_*` partner. Two reasons this is the next move:

* The MPI reduce block (`turbogap.f90:1971-2221`, 251 lines in one
  `#ifdef _MPIF90`)
  walks that ten-item list **three times** — to count, to pack, to unpack —
  re-evaluating the same ten predicates independently. If any pair ever
  disagrees, `counter2` shifts and energies are silently attributed to the
  wrong term. That is the same "several copies of one condition that must
  agree" shape as the ts+mbd bug in §5.1, and it should be killed the same way.
* **It is the enabler for the next extraction.** 30 of the exp-spectra block's
  73 boundary-crossing variables are exactly these families.

Sketch: a `contribution_ref` type holding pointers to the existing arrays, the
list built once, then short pack/unpack loops. Requires adding `target` to the
~30 driver arrays, which the compiler enforces.

### 7.2 Extract the remaining phases

Coupling was measured to choose the order; do not re-derive it:

| block | lines | vars used | block-only (become locals) |
|---|---|---|---|
| setup/read-input | 496 | 58 | 19 | **done** (`b7781bb`) |
| exp-spectra | 524 | 109 | 36 |
| nested-sampling | 1012 | 126 | 23 |
| md-driver | 373 | 75 | 8 |

`exp-spectra` next, after §7.1 drops its crossings from 73 to ~48. It is also
where the GPU branch diverges most (`exp_interface.f90` is 4838 changed lines),
so extracting it creates a seam exactly where the merge is hardest.

`md-driver` is a **poor standalone candidate** at only 8 block-only variables —
it is worth little until more state is bundled. `nested-sampling` is the
biggest single block but needs a regression case first; there is none.

### 7.3 The backend seam, then the merge

Once the phases exist: `gap_backend.f90` declares the interface,
`gap_backend_cpu.f90` implements it, the Makefile selects it. Device pointers,
cuBLAS handles and streams live **inside** the GPU implementation, never in the
driver. Then the GPU branch adopts the same skeleton and the merge proceeds
per-phase.

### 7.4 Open questions for whoever owns the vdW code

* The `S_xyz_inv` strain term (`forces_vdw = forces_vdw -
  local_virial_vdw_diag_corr * S_xyz_inv`) and the
  `this_local_virial_vdw_diag` correction are still commented out. They are a
  *different* correction from the one enabled in §5.2, and
  `local_virial_vdw_diag_corr` already has a live role feeding the scaling
  update, so enabling them is its own question.
* `mbd_correction_freq` defaults to 100. Given the staleness curve in §5.2
  (0.15% force error after 10 reuse steps on a near-equilibrium cell), 100 is
  probably optimistic for anything dynamic. Worth a considered default.

---

## 8. Method notes — what actually worked

**The bit-identical contract is the whole safety net.** It is far stronger and
cheaper than tolerance comparison, and it turns "did I break anything?" into a
one-command question. Protect it: when something must change, isolate it in its
own commit and convert only the affected cases.

**Extract verbatim.** Both phase extractions lifted the block with zero edits
to its body, so the diff reads as a pure move and can be reviewed against the
original. Rename and restructure *afterwards*, in separate commits.

**Check what the tests cannot.** Before moving a block, verify that every loop
counter and scratch variable it leaves behind (`i`, `j`, `k`, `time1`,
`iostatus`, …) is *written before it is next read* in the driver. A leaked
value is invisible to the suite — the vdW and exp blocks are never enabled in
the same test case, so a dependency between them would go undetected. Both
extractions were cleared this way, per variable.

**Measure coupling before choosing the next block**, rather than guessing from
size. That is what revealed `md-driver` as a poor candidate and what turned the
contribution bundle from a nice-to-have into a prerequisite.

**Add the test before the change, not after.** The `ts+mbd` correction state
had no coverage at all; adding `vdw_tsmbd_md` first (`852d41f`) is what made
the subsequent extraction and both bug fixes verifiable rather than hopeful.

**Argument counts are the honest residue.** `compute_vdw` takes 37 arguments
and `read_input_and_gap_files` 31. That is not a design failure — it is exactly
what is left once block-private state is removed, and it is the measurement
that tells you how much bundling is still owed.
