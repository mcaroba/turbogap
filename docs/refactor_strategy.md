# Strategy: de-monolithing `turbogap.f90` and merging the GPU branch

Written 2026-08-04. Companion to `docs/refactor_handoff.md` (this repo) and
`docs/gpu_fixes_handoff.md` (the GPU repo). Those two describe what was done;
this describes what to do next and why the ordering is what it is.

All measurements below were taken on `alt` against
`turbogap_master_2026` @ `95372cd` and `turbogap_gpu_commit_mahti` @ `d19603c`.

---

## 1. The merge is roughly a third of the size it appears to be

`docs/refactor_handoff.md` §2 sizes the problem as *"9591 lines over 27 hunks
spanning the whole file — hopeless as a textual merge"*. That figure counts
whitespace. Fortran free-form ignores leading whitespace, and the two branches
were edited with different indentation settings:

| | `diff` | `diff -wB` | is whitespace |
|---|---|---|---|
| `turbogap.f90`, whole file | 7697 | **3247** | 58% |
| ‣ nested-sampling phase | 1990 | 189 | 91% |
| ‣ md phase | 532 | 108 | 80% |
| ‣ snapshot loop | 535 | 122 | 77% |
| ‣ progress bar | 82 | 2 | 98% |
| ‣ multi-snapshot check | 66 | **0** | 100% |
| ‣ mode parsing | 9 | **0** | 100% |
| `read_files.f90` | 5546 | 469 | 92% |
| `exp_utils.f90` | 6705 | 1024 | 85% |
| `md.f90` | 1103 | **52** | 95% |
| `mc.f90` | 1381 | 283 | 80% |
| `neighbors.f90` | 1352 | 316 | 77% |
| `splines.f90` | 0 | 0 | — |
| `xyz.f90` | 7 | 7 | 0% |

The real merge surface is ~3250 lines in the driver, and it is concentrated in
two phases (§3). Normalising whitespace is the single cheapest thing that can
be done to this codebase, and it must be done on **both** branches with the
same formatter settings or it does not converge.

## 2. Both drivers still have the same skeleton

Extracting the top-level `!****` banners from each file gives the same phases in
the same order — the GPU branch rewrote block interiors in place and never
restructured the program. That is what makes per-phase merging viable.

| phase | master | gpu | delta | real diff (`-wB`) | gpu tokens¹ |
|---|---|---|---|---|---|
| variable definitions | 145 | 166 | +21 | 193 | 6 |
| 3b gpu locals | — | 55 | +55 | — | 22 |
| time / mpi / mode / welcome | 113 | 221 | +108 | 123 | 7 |
| read input & GAP files | 56 | 670 | +614 | 630 | 0 |
| progress bar & timers | 51 | 57 | +6 | 2 | 0 |
| multi-snapshot check | 49 | 51 | +2 | 0 | 0 |
| snapshot read / neighbours | 316 | 353 | +37 | 122 | 0 |
| **prediction** | 1614 | 2363 | +749 | **1496** | 286 |
| md driver | 373 | 381 | +8 | 108 | 0 |
| nested sampling | 1158 | 1186 | +28 | 189 | 10 |
| gpu internal subroutines | — | 315 | +315 | — | 154 |

¹ occurrences of `gpu_*`, `cuda_*`, `c_loc`, `c_ptr`, `hip*`, `cublas*`.

The sub-blocks *inside* the prediction phase also line up in order
(realloc → collect energies → vdW → [electrostatics] → core-electron BE →
exp spectra → 2b/core_pot/3b → MPI reduce → sum terms → write).

Two facts from this table drive everything:

* **475 GPU-specific tokens in the driver, zero CUDA `#ifdef`.** Confirms
  handoff §2.1: one merged `turbogap.f90` cannot be reached by sprinkling
  `#ifdef`. It needs a backend interface with build-time selection.
* **The GPU tokens are not spread out.** 286 in prediction, 154 in the three
  internal subroutines, 22 in the 3b locals, 10 in nested sampling, **0
  everywhere else**. The seam has a small, well-defined footprint.

## 3. Where the divergence actually is

**Read-input (670 lines on the GPU side) is not a GPU problem.** It contains
*zero* GPU calls, and 88% of its code lines appear verbatim in master's already
extracted `turbogap_setup.f90`. The GPU-only remainder is electrostatics
parameters, `get_time()` calls, and a few printouts. This block does not need
merging — it needs the master extraction transplanted onto it.

**Prediction (1496 real diff lines) is the hard one**, and splits into three
different kinds of divergence that must not be conflated:

| | lines | kind |
|---|---|---|
| electrostatics block | 271, GPU only | a **feature** master lacks — physics, not GPU |
| exp-spectra (Gutierrez/Johansson) | 269 master vs 757 gpu | genuine GPU port |
| 2b / core_pot / 3b | inline in master, extracted in gpu | already seamed |

## 4. Plan

### Phase 0 — shrink before refactoring

Zero-risk, verifiable, and removes more merge surface than any extraction.

**0a. Dead code in the GPU tree — done** (`4e42b78`). 21,389 lines across nine
files that were in no Makefile list, produced no object file, and were `USE`d by
nothing. Archived to `../phase0_attic/` and `../phase0_backup/` rather than
deleted, because all but one were untracked and git could not have recovered
them. Same 25 objects build before and after; `tests/gpu/run_regression.sh`
passes `CO_predict` and `CO_md`.

Still outstanding: `src/orthonormalization_kernels.cc` (541 lines) is
*commented out* of `SRC_CC` rather than absent — a deliberate disable that needs
an owner decision. `src/soap_turbo/src/*_bak.f90` (2671 lines) sit in the
submodule, which is already dirty with `TURBOGAP_DUMP_CNK`; handle both together.

**0b. Whitespace normalisation — attempted, and it is the wrong tool.**
Superseded by §6.3 and §6.4. Review with `diff -w` and merge with
`git merge -X ignore-all-space`; do not reformat.

**0c. Adopt `get_time()` on master — done** (`fef536c`). 76 call sites across
`turbogap.f90` (53), `turbogap_setup.f90` (8), `turbogap_vdw.f90` (8) and
`gap_interface.f90` (7), behind a new `src/timing.f90`. All 15 regression cases
pass, wall-clock ratios 0.97–1.00.

### Phase 1 — transplant what master already built

`turbogap_setup.f90` onto the GPU branch (§3): 670 lines become ~56, and both
branches get an identical setup phase. Then `vdw.f90` + `turbogap_vdw.f90`
wholesale — the GPU branch still carries the old 495-line pre-MBD `vdw.f90` and
per handoff §2.3 never touched it, so master's 6538-line version simply replaces
it.

### Phase 2 — the backend seam

`gap_backend.f90` declares the interface, `gap_backend_cpu.f90` and
`gap_backend_gpu.f90` implement it, the Makefile selects one. Scope it to what
§2 measured: SOAP/2b/3b/core_pot evaluation, device buffer lifetime, and the 3b
locals. Device pointers, cuBLAS handles and streams live inside the GPU
implementation and never appear in the driver.

Half of this already exists. The GPU branch's `add_2b_contribution_gpu`,
`add_core_pot_contribution_gpu` and `add_3b_contribution_gpu` (315 lines, called
from three sites) are already the right decomposition. One caveat: the GPU
handoff §7.4 calls internal procedures "the cheap route… host association means
no argument lists", and that is exactly what blocks the merge — a contained
procedure cannot move to a module. Converting them to module procedures with
explicit arguments is a necessary cost, and the resulting argument count is the
honest measure of coupling still owed (handoff §8).

### Phase 3 — `contribution_ref` bundling

Keep it where `docs/refactor_handoff.md` §7.1 put it. It now serves double duty:
the same ten-family count/pack/unpack triple-walk exists on the GPU side too
(251 lines master, ~287 GPU), so bundling on both sides converges the MPI reduce
block as a side effect.

### Phase 4 — exp-spectra

Confirmed as the worst remaining merge site: 269 lines vs 757 in the driver, and
`exp_interface.f90` is 1865 vs 3455 with a 2193-line real diff. Extract on master
first into its own module, then the GPU version becomes the backend
implementation behind the same interface.

### Phase 5 — one-sided features, decided explicitly

* **Electrostatics** (`electrostatics.f90`, 1162 lines, plus 271 in the GPU
  driver) exists only on the GPU branch and is *physics, not GPU code*. Lift it
  into its own module and merge it to master independently of the GPU work.
* **Radiation cascade stack** (`adaptive_time`, `eph_*`) is the mirror case:
  wired into master's driver (5 `use` statements), and in the GPU tree compiled
  as `SRC_STOP` but never linked — `OBJ_STOP` appears in no rule's
  prerequisites — and `use`d zero times. Master's wiring is the one to keep.
* **MC** has diverged in master's favour: master uses `params%mc_mu(mc_mu_id)`
  and `n_mc_species(mc_mu_id)` where the GPU branch passes them un-indexed, and
  master uses `random_number` where the GPU branch uses `mod(irand(), n_xyz)`.

## 5. Two corrections to the existing plan

**Do not extract nested-sampling for merge reasons.** `docs/refactor_handoff.md`
§7.2 lists it as the biggest single block (1012 lines) and notes it needs a
regression case first. It is also 91% whitespace-divergent — after Phase 0b it
is ~189 lines apart and already merges nearly free. It is not on the merge
critical path. Extract it later for maintainability, once a test exists.

**`contribution_ref` is the one step with real performance risk.** Pointer
arrays defeat the compiler's aliasing analysis, and the SOAP loop is exactly
where that matters. Measure `co_predict` immediately before and immediately
after that specific commit rather than at the end of the phase. Elsewhere,
declare extracted dummy arguments `contiguous` to preserve the assumptions host
association previously gave for free.

## 6. Verification notes

**The binary is not reproducible; do not compare binaries.** Three builds of an
unchanged GPU tree gave three different md5sums
(`1f251c64…`, `20d29991…`, `3ec89fe2…`). The bit-exact contract in
`docs/refactor_handoff.md` §4.2 is about *program output*, and the `diff`-based
comparator is the right instrument. A binary checksum would have reported a
false failure for Phase 0a, which the regression suite correctly passed.

**`make -j` is broken.** The Makefile expresses no Fortran module dependencies,
so a parallel build races and fails with `Cannot open module file
'read_files.mod'`. Builds must be serial (28 s for the GPU tree, DEBUG=0).
Adding dependency ordering is a cheap, separate improvement worth doing before
the build gets larger.

**Timing.** `co_predict` (~18 s) is the reference; ignore sub-5-second ratios
(`vdw_mbd` gave 2.01/2.01/3.34 s across repeats of one binary). Always build the
GPU tree `DEBUG=0` — the default is `DEBUG=1`, a 2.1x slowdown. Fix the negative
"Miscellaneous" bucket in the end-of-run summary early, so the in-code timers
become usable as refactor instrumentation instead of relying on wall clock alone.

### 6.3 Reformatting does not converge the branches — measured

§1 recommends normalising whitespace so the real diff becomes visible. That
recommendation was wrong. This is what the attempt showed.

`fprettify 0.3.7` with `-i 2 --disable-whitespace -l 1000` formats every module
in the tree cleanly and provably changes nothing but whitespace — `diff -wB`
between original and formatted is empty for all seventeen. But it **cannot parse
`turbogap.f90` on either branch** (§6.4), and on the files it can parse it does
not reach the goal:

| module | branches differ, raw | after formatting both | `diff -w` |
|---|---|---|---|
| `read_files.f90` | 5546 | 3337 | **469** |
| `exp_utils.f90` | 6705 | 2786 | **1024** |
| `md.f90` | 1103 | 353 | **52** |
| `mc.f90` | 1381 | 632 | **283** |
| `neighbors.f90` | 1352 | 781 | **316** |

Formatting more than halves the apparent difference but never reaches the floor,
because fprettify indents according to block structure and the two branches
genuinely have different block structure in places. `diff -w` is strictly better
on every file, costs nothing, and touches no history.

So: **do not reformat.** Review with `diff -w`, merge with
`git merge -X ignore-all-space`. That captures the whole benefit §1 was reaching
for, at zero risk, without a ~30,000-line commit that would destroy `git blame`
across the tree.

### 6.4 Sixteen statements that no Fortran tool can parse

fprettify fails on `turbogap.f90` at an argument list interrupted by the
preprocessor — the continuation is split across `#ifdef _MPIF90`, with a
different tail of arguments on each side:

```fortran
     & xyz(1:3, j_beg:j_end),&
#ifdef _MPIF90
     & this_energies_lp(i_beg:i_end), this_forces_lp, this_virial_lp, ... )
#else
     & energies_lp(i_beg:i_end), forces_lp, virial_lp, ... )
#endif
```

There are **5 such sites on master** (lines 1502, 1696, 1721, 1745, 1773) and
**11 on the GPU branch**. No Fortran-aware tool can handle these, because the
token stream depends on the preprocessor — they defeat formatters, static
analysers and IDE parsers alike.

All sixteen are the same shape, and the only thing differing across the `#ifdef`
is the `this_` prefix on an `energies_X` / `forces_X` / `virial_X` triple — i.e.
exactly the contribution families §7.1 of `docs/refactor_handoff.md` wants to
bundle. They are the most extreme instance of the "one condition written at N
sites" antipattern that produced the ts+mbd bug. **Phase 3 should eliminate all
sixteen as a by-product**, which is a further argument for doing it early; after
it, `turbogap.f90` becomes parseable by ordinary tooling for the first time.

## 7. Suggested order

1. ~~0a dead code~~ — done (`4e42b78`)
2. ~~0c `get_time()`~~ — done (`fef536c`)
3. ~~0b whitespace~~ — dropped, see §6.3; use `diff -w` and `-X ignore-all-space`
4. Phase 3 `contribution_ref` — promoted ahead of Phase 1, because §6.4 shows it
   also removes the sixteen unparseable statements. Measure `co_predict`
   immediately either side of this commit.
5. Phase 1 transplant — turns the GPU read-input block from 670 lines to 56
6. Phase 4 exp-spectra extraction on master
7. Phase 2 backend seam, then merge per phase
8. Phase 5 one-sided features, each on its own commit
