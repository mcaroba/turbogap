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

**Both remaining items closed 2026-08-05.**

`src/soap_turbo/src/*_bak.f90` — removed (submodule `856b860`, parent `6a2f772`).
2710 lines, not 2671. Archived to `../phase0_attic/soap_turbo_bak_files_20260805.tar.gz`,
though unlike the phase 0a files these were *tracked*, so git alone could recover
them. The reason they had to go first: each defines the **same module name** as
its live counterpart (`soap_turbo_desc`, `soap_turbo_radial`,
`soap_turbo_angular`), so any module-name-to-file mapping is ambiguous while they
exist — which blocks the `make -j` fix below.

`src/orthonormalization_kernels.cc` — **kept, by owner decision.** It stays
commented out of `SRC_CC` and stays in `src/` for possible reintegration. The
Makefile now says so, so it does not get re-litigated as dead code.

**0b. Whitespace normalisation — attempted, and it is the wrong tool.**
Superseded by §6.3 and §6.4. Review with `diff -w` and merge with
`git merge -X ignore-all-space`; do not reformat.
**Superseded 2026-08-05 — see §6.3.** Both trees are now fprettify-normalised
with identical settings, the raw diff has fallen to the `diff -w` floor
(26786 → 7828 lines across the shared modules), and `-X ignore-all-space` is no
longer needed for the merge to behave.

**0c. Adopt `get_time()` on master — done** (`fef536c`). 76 call sites across
`turbogap.f90` (53), `turbogap_setup.f90` (8), `turbogap_vdw.f90` (8) and
`gap_interface.f90` (7), behind a new `src/timing.f90`. All 15 regression cases
pass, wall-clock ratios 0.97–1.00.

### Phase 1 — transplant what master already built — done (setup), deferred (vdw)

On the GPU branch: `timing.f90` adopted (`9a1a2b4`), then the start-up phase
extracted into its own `turbogap_setup.f90` (`08809dd`). `turbogap.f90` there
went 5889 → 5271.

The measurement that matters: the start-up phase was **56 lines on master
against 670 on the GPU branch** — an apples-to-oranges comparison, because one
side had been extracted and the other had not, which hid the shared 88%. It is
now **16 lines against 12**, and the two modules differ by **235 lines ignoring
whitespace**. That 235 is the real divergence, and it is almost entirely two
features: electrostatics parameters on the GPU side, the eph and
electronic-stopping initialisation on master's.

**`vdw.f90` deliberately not transplanted.** It is not a one-file move: master's
version `use`s `misc`, `constants` and `nonneg_leastsq`, none of which exist in
the GPU tree. And nothing in `tests/gpu` exercises vdW at all, which is thin for
a 495 → 6538 line replacement. Same reasoning for wiring the cascade stack: the
GPU branch compiles `SRC_STOP` but `OBJ_STOP` is in no rule's prerequisites, so
it is never linked. Both are genuine merge work and want master's 18-case suite,
not the GPU branch's two.

### Phase 2 — the backend seam — done, both sides

| | |
|---|---|
| `85c1fdb` (cpu) | `gap_backend_cpu.f90`: the three contribution procedures |
| `33c377d` (cpu) | `gap_backend_begin` / `gap_backend_end` added to the interface |
| `5a42ddf`, `5e06143` (gpu) | `gap_backend_gpu.f90`: the GPU implementation |

Both branches now define `module gap_backend` with the same five public
procedures and the same argument lists; the Makefile compiles one. The driver
calls the names and never learns which implementation it got.

`turbogap.f90`: master 3552 → 3487, GPU 5129 → 4791. **Neither driver names a
device pointer in the contribution path any more.** The ten buffers the GPU
procedures used to reach by host association are module state inside
`gap_backend_gpu`; `gap_backend_begin` uploads the neighbour data once for all
three calls and `gap_backend_end` releases it, and on the CPU both are empty.

`alphas_d` is the backend's own rather than shared with the driver's SOAP path,
which removes one instance of the name reuse that hid bug 5.

**Three defects had to be fixed to get here, and only one was in the
extraction.** Two were latent in the original GPU code:

* the three procedures had an **undocumented ordering dependency** —
  `add_2b_contribution_gpu` computed `st_n_sites_double` and `st_virial`, and
  the other two read them through host association without ever setting them.
  Calling 3b without 2b first would have copied back the wrong number of bytes.
* `virial_2b_d` and `virial_core_pot_d` were released with the plain
  `gpu_free` (`hipFree`) despite being `hipMallocAsync` allocations — the same
  defect as `gpu_fixes_handoff.md` §6e, at two further sites, with the intended
  fix half-written in a trailing comment.

The third was mine, and it cost the most: the lifted setup block stopped one
line short and dropped the `cpy_htod` that uploads `neighbors_list_d`, so the
buffer was allocated and never filled. The kernel then indexed
`neighbors_list` with uninitialised device memory, producing a negative index
and out-of-bounds atomics. **Every argument was byte-identical because the bug
was a missing statement, not a wrong value** — which is why an
argument-by-argument comparison could never have found it, and why the
structural bisect found it in one step. Record that: when inputs all match and
the answer does not, stop comparing inputs.

### Phase 3 — `contribution_ref` bundling — done

Landed in three commits, plus the coverage it turned out to require first:

| | |
|---|---|
| `3cdc306` | the ten predicates evaluated once into `contrib_on` |
| `bb8e16f` | `xrd_mad`, `xrd_mad_mpi2`, `xrd_predict` — coverage for pdf/sf/xrd |
| `7390225` | `contribution_ref`, pack and unpack collapse to one loop each |
| `<fix>` | KNOWN_ISSUES #4, the uninitialised force slots |

Sequencing note worth keeping. The plan was to bundle in one step. Reading the
block first showed that four of the ten families — `pdf`, `sf`, `xrd`, `nd` —
had **no regression coverage at all**, so bundling would have restructured them
blind. Writing the exp-spectra cases had to come first, and it immediately paid
for itself by exposing KNOWN_ISSUES #4. Splitting predicate-unification from
pointer-bundling also meant the `target` attribute could be measured on its own.

Net line count barely moved (3933 → 3928): the ~110-line setup replaces ~215
lines that stated the same per-family facts three times. The win is that a
family's identity is now declared in exactly one place, so the walks cannot
disagree — and that `contribution_ref` is the aggregate Phase 4 needs.

Still outstanding on the GPU side: the same triple-walk exists there (~287
lines), and bundling it the same way converges the reduce block as a side
effect.

### Phase 4 — exp-spectra — done

Both blocks are out, into `src/turbogap_exp.f90`.

| | |
|---|---|
| `afd0c19` | scattering observables (pdf/sf/xrd/nd), 272 lines, 28 arrays become locals, 37 args |
| `03877d5` | the XPS spectrum, 161 lines, 5 arrays plus i/j/j2/k become locals, 30 args |
| `533ccf9` | the 39 declarations those two orphaned |

`turbogap.f90`: **4992 at the branch point → 3538, down 29%.**

**The driver parses.** `fprettify` accepts `turbogap.f90`, `turbogap_exp.f90`,
`turbogap_setup.f90` and `turbogap_vdw.f90` — every file in `src/`. Until this
phase the driver defeated formatters, static analysers and IDE parsers alike,
because an argument list interrupted by `#ifdef` has no token stream until the
preprocessor has run. That was the §6.4 blocker, and extracting these two blocks
is what removed it: the five statements each had their `this_`-prefixed and
un-prefixed variants split across the preprocessor, and inside a procedure the
dummy has one name either way, so the caller chooses once and the conditional
disappears from the moved code.

Both blocks are now out on the GPU branch too (`a891432` XPS, `a84be03`
spectra), so `turbogap_exp.f90` holds the same two procedures on both sides.

**`exp_interface.f90` is not the largest divergence left — it is one of the
smallest.** The 2193-line figure is a file-level diff, and a file-level diff
cannot tell a *conflict* from an *addition*. Compared procedure by procedure
with whitespace and comments stripped (`tools/proc_diff.py`), 1867 vs 3457
decomposes as:

| | lines | |
|---|---|---|
| 16 GPU-only procedures | **1898** | pure additions — the batched-device paths. No conflict; they merge by appearing. |
| 10 shared procedures | **55** real diff | the entire merge surface |
| 1 master-only procedure | 283 | `get_pdf_sf_xrd_explicitly_kde` — owed to the GPU branch |
| remainder | — | formatting |

**Eight of the ten shared procedures are identical** ignoring whitespace and
comments — including `calculate_pair_distribution`, all 479 lines of it. The
whole conflict is 55 lines in two procedures:

* `calculate_xrd` (31) and `calculate_structure_factor` (24), and most of that
  is inherent device plumbing: the `nk_d`/`k_index_d`/`j2_index_d` argument
  block, `cublas_handle`/`gpu_stream`/`gpu_host_storage`/`gpu_low_memory`, and
  their declarations. Those are what a GPU implementation *is*.

Three items in there are genuine and worth reconciling on their own, none of
them GPU-specific:

1. `preprocess_exp_data`'s i(q)-from-q\*i(q) conversion — **done** (`451735e`),
   and that procedure is now identical.
2. `virial_sf = 0.d0` in `calculate_structure_factor`, present only on the GPU
   side. Master zeroes the exp-spectra virials elsewhere (`4f99e92`); check the
   two are not both needed, or both missing a case.
3. `get_structure_factor_forces_matrix` takes `i_beg, i_end` and
   `species(i_beg:i_end)` on the GPU branch and neither on master. That is an
   MPI-decomposition difference in a routine every test exercises single-rank,
   which is exactly the shape of the unexercised multi-rank indexing bug noted
   at the end of §7.

**The measurement lesson generalises.** A whitespace-insensitive comparison
must strip whitespace *entirely*, not collapse runs of it: the GPU copy of this
file has been through fprettify, so it writes `(:, :)` for `(:,:)`,
`deallocate (x)` for `deallocate(x)`, and `a*b/c` for `a * b / c`. Collapsing
runs to a single space still counts every one of those as a difference, and
reports 757 lines of divergence where the truth is 67. Re-check the other
files' figures the same way before trusting them.

### Phase 5 — one-sided features, decided explicitly

* ~~**Electrostatics**~~ — **done (`d64d7f4`).** The "physics, not GPU code"
  reading was nearly right: the module is 24 procedures, of which exactly one,
  `calculate_batched_electrostatics` (362 lines), carries all 120 of its device
  tokens. The other 23 — including all three the driver calls — have none. So
  master's copy is the GPU copy minus that one procedure and the two `USE`
  statements it needed, with zero residual device tokens, and `gsf` takes the
  `compute_coulomb_lamichhane` path the GPU branch itself falls back to.

  Came with it: the `atomic_charge` local property threaded through
  `read_files`, `turbogap_setup` and the driver; `N_CONTRIB` 10 → 11 with
  `C_ESTAT` inserted at 3 and `C_LP..C_3B` renumbered **to match the GPU
  branch's numbering rather than appending**, so the two files agree; eleven
  keywords; `options_estat`.

  **A latent defect was found and fixed on both branches.** The block's guard
  was `estat_method /= "none"` alone, so a deck asking for electrostatics
  against a GAP with no `atomic_charge` property indexed `local_properties` with
  an uninitialised `charge_lp_index`. It segfaults; reproduced. It now warns and
  skips. That is the **fourth** defect here of the form *one condition, two
  predicates free to disagree* — after ts+mbd (§5.1), the MPI reduce triple walk
  (§7.1), and `has_vdw` against `has_local_properties`.

  **Not ported, deliberately:** the unconditional `energy_estat=` field
  `xyz.f90` writes on the GPU branch. It changes every trajectory header, which
  takes all ten baseline-compared cases red at once — a result-changing commit
  that belongs on its own with the cases blessed deliberately, per §4.2.

  **Not verified:** the electrostatics numbers. No GAP in the test data declares
  an `atomic_charge` local property, so the path cannot be exercised on *either*
  branch. The GPU branch's electrostatics has never been tested either. A
  charge-carrying GAP is the next thing this needs.
* **Radiation cascade stack** (`adaptive_time`, `eph_*`) is the mirror case:
  wired into master's driver (5 `use` statements), and in the GPU tree compiled
  as `SRC_STOP` but never linked — `OBJ_STOP` appears in no rule's
  prerequisites — and `use`d zero times. Master's wiring is the one to keep.
* **MC** has diverged in master's favour: master uses `params%mc_mu(mc_mu_id)`
  and `n_mc_species(mc_mu_id)` where the GPU branch passes them un-indexed, and
  master uses `random_number` where the GPU branch uses `mod(irand(), n_xyz)`.

### Phase 6 — the driver has reached its extraction limit, measured

Re-measured both remaining candidate blocks in master's driver with
`tools/classify.py`, after Phases 3 and 4 had done their bundling:

| block | lines | arguments | block-private locals |
|---|---|---|---|
| nested sampling | 1099 | **145** | 23 |
| md driver | 368 | **67** | 9 |

Both confirm the §5 and §7.2 judgements rather than overturning them, and now
show *why*. Splitting the arguments by kind:

| | contribution arrays | timers | everything else |
|---|---|---|---|
| nested sampling | 38 | 22 | 85 |
| md driver | 14 | 3 | 50 |

So bundling the driver's ~38 contribution arrays into the `contrib` aggregate it
already builds, and its 22 timing buckets into one derived type, would take
nested sampling from 145 arguments to **87** and the md driver from 67 to
**52**. That is a 40% cut and still not a good extraction — those blocks
genuinely read and write state that spans the whole program.

**The conclusion to record: de-monolithing `turbogap.f90` by lifting phases is
finished.** The blocks that were worth lifting have been lifted — setup, vdW,
exp-spectra, XPS, the backend seam. What remains is 3682 lines on master and
4247 on the GPU branch of code whose coupling is to the driver's *state*, not to
a phase boundary. Cutting it further means bundling that state first — the
contribution arrays and the timers are the two obvious aggregates — and that is
a different kind of change from an extraction, with its own risk and no merge
payoff on its own.

Meanwhile the merge surface is where the remaining value is. Current real
divergence, whitespace-insensitive, after the vdW and electrostatics work:

| file | master | gpu | real diff |
|---|---|---|---|
| `exp_interface.f90` | 1947 | 3760 | 2416 (but only ~55 is conflict; see §7.2) |
| `turbogap.f90` | 3682 | 4247 | 1519 |
| `turbogap_exp.f90` | 581 | 1046 | 1098 |
| `exp_utils.f90` | 3413 | 4125 | 1086 |
| `read_files.f90` | 3291 | 3285 | 458 |
| `electrostatics.f90` | 840 | 1197 | 373 |
| `vdw.f90` | 6938 | 6938 | **0** |

## 5. Two corrections to the existing plan

**Do not extract nested-sampling for merge reasons.** `docs/refactor_handoff.md`
§7.2 lists it as the biggest single block (1012 lines) and notes it needs a
regression case first. It is also 91% whitespace-divergent — after Phase 0b it
is ~189 lines apart and already merges nearly free. It is not on the merge
critical path. Extract it later for maintainability, once a test exists.

**~~`contribution_ref` is the one step with real performance risk~~ — measured,
and it is not.** The claim was that pointer arrays defeat aliasing analysis and
the SOAP loop is where that matters. Measured on `co_predict`, interleaved, four
repetitions with and without the `target` attribute on the ten families:

| | reps (s) | mean |
|---|---|---|
| without `target` | 18.17, 20.51, 21.15, 19.65 | 19.87 |
| with `target` | 19.41, 19.69, 20.98, 19.40 | 19.87 |

No regression. The reason is that the hot compute lives in callees operating on
dummy arguments, where aliasing is governed by Fortran's argument rules rather
than by `target` on the actual argument — `target` on a driver array does not
reach into `soap_turbo`. Keeping it off `energies` and `forces`, the
accumulators, costs nothing and keeps the assumption narrow.

The spread above is the more useful number: **18.2 to 21.2 s on one unchanged
binary**, a 16% range. Never trust a single wall-clock ratio here; the suite's
one-shot `x1.07` for this case was noise. `contiguous` on extracted dummy
arguments is still worth doing when blocks move out of the driver, for the
assumptions host association gave for free.

## 6. Verification notes

**The binary is not reproducible; do not compare binaries.** Three builds of an
unchanged GPU tree gave three different md5sums
(`1f251c64…`, `20d29991…`, `3ec89fe2…`). The bit-exact contract in
`docs/refactor_handoff.md` §4.2 is about *program output*, and the `diff`-based
comparator is the right instrument. A binary checksum would have reported a
false failure for Phase 0a, which the regression suite correctly passed.

**~~`make -j` is broken.~~ — fixed on the GPU tree 2026-08-05 (`6a2f772`).**
The Makefile expressed no Fortran module dependencies, so a parallel build raced
and failed with `Cannot open module file`. `tools/gen_fortran_deps.py` now
generates the 23 ordering rules into `makefiles/Makefile.deps`, included by the
Makefile. Serial 29 s unchanged; `-j12` 15 s, verified over six builds at
`-j12/-j8/-j4` plus a `deepclean` rebuild, 30 objects every time. Regression
suite after a parallel build is identical to the serial baseline.

Two traps worth keeping, both of which report success while being wrong:

* **The include must come last.** Placed before the first rule, one of the
  generated dependency lines becomes make's *default goal*: `make` builds a
  single object and **exits 0**. An exit-code check passes; only counting the
  objects catches it.
* **Duplicate module definitions must be gone first.** The generator fails loudly
  if two files define the same module, which is exactly what the `*_bak.f90`
  files did. Ordering here was forced, not chosen.

Not covered: `SRC_STOP`. `OBJ_STOP` is still in no rule's prerequisites, so it is
never built or linked, and generating dependencies for it would be meaningless
until it is wired up. Master's tree still needs the same treatment.

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

~~So: **do not reformat.**~~ **Reversed 2026-08-05** (`ca390bc` cpu,
`415f881` gpu). Both trees are now fprettify-normalised with identical settings
and a `.pre-commit-config.yaml` keeps them that way.

The measurement above was right; the conclusion drawn from it was too narrow. It
tested whether formatting **alone** converges the branches. It does not — but
that was never the choice. The choice was between formatting *neither* tree and
formatting *both with the same settings*, and the table above formats each in
isolation without asking what the pair then looks like. Measured on the pair:

| | before | after |
|---|---|---|
| `exp_utils.f90` | 6039 | **1106** |
| `read_files.f90` | 5399 | **478** |
| `turbogap.f90` | 5535 | 2372 |
| `exp_interface.f90` | 4957 | 2465 |
| `md.f90` | 983 | **55** |
| `mc.f90` | 1228 | 333 |
| `neighbors.f90` | 1328 | 339 |
| `types.f90` | 743 | 179 |

26786 → 7828, **−71%**, and those land on the `diff -w` floor in the table
above — `read_files` 478 against 469, `md` 55 against 52. The raw diff and the
semantic diff have converged, which is the thing §1 actually wanted: `git merge`
now behaves without `-X ignore-all-space`, and a file-level diff between the
branches means something again.

What made this worth the `git blame` cost was §7.2's finding: the GPU copy of
`exp_interface.f90` had already been through fprettify and master's had not, so
the pair read as 757 diverged lines where the truth was 67. Formatting drift
between the branches was not cosmetic — it was inflating the measured merge
surface by an order of magnitude and hiding what actually conflicts.

The whole change is whitespace, verified per tree with
`git diff -w --ignore-blank-lines` over `src`: one line each, in `xyz.f90`,
where fprettify drops a redundant leading `&`. Both suites pass unchanged.

Three things had to be fixed first, all of which would otherwise fail *quietly
and forever*:

* **The hygiene hooks rewrote test data** — `silicon.gap`, `silicon.xml`, the
  `.beta` tables, a regression case's input deck. Parsed inputs where a stripped
  trailing space changes what the code reads. The hooks are now allow-listed to
  `src/ tools/ docs/ makefiles/`; an exclude list fails open.
* **fprettify and `trailing-whitespace` never converge** on `eph_fdm.f90`:
  fprettify emits a trailing space after a statement-separating `;`, the fixer
  strips it, repeat. Fortran is excluded from the fixer — when a formatter and a
  fixer disagree, exactly one may own the file.
* **fprettify is not idempotent on the GPU `turbogap.f90`** as it stood, for two
  reasons in one `get_gap_soap` call: eight commented-out argument lines *inside*
  the continued argument list (dead host-side arguments the GPU replaced with
  device pointers), and the call being right-aligned to exactly column 132, where
  the alignment drifts a few spaces per pass. Both fixed; the statement's 75
  arguments are byte-identical across the rewrap.

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
sites" antipattern that produced the ts+mbd bug.

**Correction after doing Phase 3:** it did *not* eliminate them. Phase 3 changed
the MPI reduce block, whereas the sixteen are in the exp-spectra *call sites* a
few hundred lines earlier — the families are the same, the code is not.

**Resolved in Phase 4.** All five on master are gone and `turbogap.f90` parses.
The eleven on the GPU branch remain, and the same extraction removes them.

A hazard for whoever does the GPU side. These blocks contain two kinds of
`#ifdef _MPIF90` that look alike. One has an `#else` and splits an argument
list — that is the defect, and it collapses. The other has **no** `#else` and
guards an MPI-only *allocation* of the same arrays; it is an ordinary
conditional over whole statements, it parses fine, and it must survive the move.
A transformation that scans for the first `#else` will swallow it and silently
make the allocation unconditional. Match directive nesting, and assert how many
groups you collapsed against how many you kept.

### 6.5 The GPU build is not bit-reproducible, and what follows from it

Three runs of the **same unchanged GPU binary** on the `XRD_mad` case produce
three different trajectories; `compare_xyz.py` reports two differing quantities
between two such runs, appearing once MD has accumulated a step, at ~1e-12
relative. The CPU build on the identical input is reproducible. Recorded as §6c
of `docs/gpu_fixes_handoff.md` on that branch.

**This changes the method.** The bit-exact contract — output identical to a
pre-refactor baseline, compared with `diff` — cannot be used on the GPU branch,
because there is no baseline to be identical to. A change smaller than the
run-to-run spread is invisible there.

So: **do each extraction on master first, where the contract holds, and port the
verified result.** That is not a preference, it is the only place the work can be
checked. It is why Phase 2's CPU implementation was written before the GPU one,
and it is the sequencing rule for everything that remains.

### 6.6 Coverage on the GPU branch, and what it found

The GPU suite had no coverage of the exp-spectra paths at all. Adding one case
(`XRD_mad`, `515ad35`) immediately found that those paths **did not work**: a
leftover debugging `exit(0)` in `gpu_exp.cu` between the pair-distribution
evaluation and reduction kernels terminated the process with status 0 and no
trajectory (`c1d96f0`). Because the status was 0, nothing checking a return code
would ever have caught it — the case tests for the absence of a trajectory.

With that fixed, GPU and CPU agree exactly on energy, every energy component and
forces, and to 3e-08 on `local_energy`. The ~~**virial is wrong by ~90x**~~ —
**fixed (owner, 2026-08-05).** `XRD_mad` is still xfail, but against a different
and much smaller defect: a `local_energy` drift of ~1e-5 by frame 5, with
everything else agreeing. Open, deferred deliberately.

The general lesson is the one Phase 3 already taught: **write the case before
the refactor, not after.** Both times, the coverage found a real defect within
minutes of existing.

## 7. Suggested order

Done, in the order it happened:

| | |
|---|---|
| 0a dead code | `4e42b78` (gpu) |
| 0c `get_time()` | `fef536c` (cpu), `9a1a2b4` (gpu) |
| 0b whitespace | dropped, §6.3 |
| Phase 3 `contribution_ref` | `3cdc306` `7390225` `b51d832` (cpu), `4357188` (gpu) |
| exp-spectra coverage | `bb8e16f` (cpu), `515ad35` (gpu) |
| Phase 4 exp-spectra | `afd0c19` `03877d5` `533ccf9` (cpu), `a891432` (gpu, XPS only) |
| Phase 1 setup transplant | `08809dd` (gpu) |
| Phase 2 seam, CPU side | `85c1fdb` (cpu) |
| 0a remainder, `*_bak.f90` | `856b860` (gpu submodule), `6a2f772` (gpu) |
| `make -j` module deps | `6a2f772` (gpu) |
| exp-spectra virial | fixed by owner, 2026-08-05 |

Left, in the order the measurements support:

1. **The `local_energy` drift on the GPU branch**, ~1e-5 by frame 5 on
   `XRD_mad`, everything else agreeing. This is what replaced the virial bug as
   the reason that case is xfail. Deferred deliberately — not blocking.
2. **The diffraction block on the GPU branch — blocked, and on what.**
   Re-measured after Phase 2: still 701 lines and **still 71 crossing
   variables**. Phase 2 did not help, because the 2b/3b device buffers it took
   over are not the ones this block uses.

   The blocker is that the device context here is **genuinely shared**, unlike
   `gap_backend`'s. Uses inside the block against uses elsewhere in the driver:

   | | in block | outside |
   |---|---|---|
   | `gpu_stream` | 22 | **54** |
   | `n_omp` | 3 | 13 |
   | `i_beg_list` | 11 | 12 |
   | `gpu_streams` | 7 | 6 |
   | `gpu_neigh` | 6 | 6 |
   | `omp_task` | 11 | 9 |

   The batch lists are *built* at lines 1283-1364 and 1490 for the
   electrostatics and SOAP paths and consumed at 1746-1849 as well as here.
   `gap_backend_gpu` worked because its ten buffers were private to the
   2b/3b region, so they could become module state without the names changing
   and the lift stayed verbatim. That trick does not apply here.

   **What this needs is a shared GPU context module** — owning the stream, the
   cuBLAS handles, the OpenMP batch count and the batch decomposition lists —
   used by the electrostatics, SOAP and exp paths alike. That is a design
   decision spanning three subsystems rather than an extraction, and it should
   be agreed before it is started. Doing the extraction first, with 71
   arguments including ten device handles, would freeze exactly the interface
   that context module then has to undo.

   **Unblocked 2026-08-05 (`19aa6a9`). The context is eight symbols, not 71.**
   The count above conflates three different things, and separating them is what
   made the module writable:

   * **28 of the crossings are declaration-only.** Their sole appearance outside
     the block is the `type(c_ptr) :: ...` line itself — `st_x_d`, `xpdf_d`,
     `dv_d`, `pair_distribution_partial_d`, the whole `_d` family. They are
     block-private and move with the code.
   * **Four are procedures, not state** — `gpu_malloc_all`, `gpu_free_async`,
     `gpu_malloc_neighbors`, `gpu_free_neighbors` arrive via `use` at zero
     interface cost.
   * **Three are `params%` fields** — `gpu_batched`, `gpu_n_batches`,
     `gpu_max_batch_size` — already carried by `params`.

   What genuinely crosses is `gpu_stream`, `cublas_handle`, `gpu_streams`,
   `cublas_handles`, `gpu_exp`, `gpu_neigh`, `gpu_batch_storage` and
   `gpu_memory_usage`. Those are now `src/gpu_context.f90`, reached by `use`, so
   the names are unchanged and the lift stays verbatim — the `gap_backend_gpu`
   trick does apply here after all, once the private buffers are told apart from
   the shared ones.

   **`omp_task` must not go in the module**, and this is the trap. It reads like
   context and it crosses like context, but `turbogap.f90` intends it to be
   OpenMP-*private* — see the `!$omp parallel private(omp_task)` and
   `!$OMP PRIVATE(i, omp_task, ...)` directives. Shared module state there is a
   correctness bug the moment `-fopenmp` is enabled, which it currently is not,
   so nothing would catch it. The arrays it indexes (`gpu_streams`,
   `cublas_handles`) are shared and do belong in the module; the index does not.
   The batch decomposition lists are excluded for a weaker reason: they are
   recomputed per snapshot and handed to `neighbors.f90` as `intent(out)`, so
   they read as data flowing through the driver.

   **Done (`a84be03`, `ca6ee53`).** Block was lines 1979-2680, 702 lines,
   confirming the 701 above. `turbogap.f90`: **4793 → 4101, −14%**.
   `compute_exp_spectra` now sits beside `compute_exp_xps` in
   `src/turbogap_exp.f90`, same module and procedure names as the CPU branch, so
   the two are finally comparable procedure against procedure.

   **51 arguments against master's 37**, 57 driver declarations become locals,
   52 of which were then removed as orphaned (`ca6ee53`). The gap to 37 is
   entirely the batched-device machinery: the four batch-decomposition lists,
   the four per-batch bounds, `n_omp`/`omp_task`, `n_sites_temp`,
   `n_pairs_temp`, `time_exp_batched`.

   Two reductions were applied *before* lifting, not after, which is what got it
   from 64 to 51:

   * `i`, `j`, `k`, `f` became locals under §8's leaked-value rule, checked per
     variable. `j`'s first use in the block is the *labelled* `outerchk: do j =`
     loop, which a naive "is this a do-loop index" test misses; `f` turned out
     to be used nowhere else in the file at all.
   * The nine `this_`/plain contribution pairs collapse to nine dummies rather
     than eighteen, because inside a procedure the dummy has one name and the
     caller chooses once.

   That second one removed **four of the nine argument lists no Fortran tool can
   parse** (§6.4), leaving five. The new call site's own `#ifdef` is an ordinary
   whole-statement one and parses.

   Three traps, all of which report success while being wrong:

   * **Comparing the two `#ifdef` branches line by line says they differ.** They
     break their continuations in different places. Join each branch into one
     statement and normalise whitespace before comparing, or the collapse looks
     unsafe when it is not.
   * **Three `this_forces_xrd`/`this_virial_xrd` references sit outside any
     `#ifdef`** — the block allocates and zeroes them unconditionally before
     choosing which pair to hand to `collect_batched_forces`. A transformation
     that only rewrites inside directive groups leaves them behind, and they
     then have no IMPLICIT type. Both trees build with `-D _MPIF90`
     unconditionally, so the `#ifdef` branch is the one that compiles.
   * **A crossing count that does not discount declaration *continuation* lines
     overstates the argument list badly** — 89 against the real 64 on this same
     block. Discounting only the first line of a declaration is not enough.

   Tooling left behind: `tools/extract_spectra.py` (the move, with every count
   asserted), `tools/find_orphans.py` and `tools/drop_orphans.py` (the
   declaration cleanup; the former counts a name as live if it appears in *any*
   other declaration, which is the `N_CONTRIB` trap from `533ccf9`).
3. **Phase 5 one-sided features**, each on its own commit: electrostatics to
   master, the cascade stack wired on the GPU branch, master's MC fixes across.
4. `vdw.f90` + `misc`/`constants`/`nonneg_leastsq` to the GPU branch, once
   there is vdW coverage there.
5. **The merge itself.** Both drivers now share `turbogap_setup`, `turbogap_exp`
   and `gap_backend`, so the per-phase comparison the whole plan was built
   around is finally available.

Cheap, any time: an `nd` regression case (shares `calculate_xrd` with `xrd`);
the 65 remaining unreferenced declarations in the driver (mind `N_CONTRIB`,
§533ccf9); `KNOWN_ISSUES` #5, the unguarded non-root unpack; the multi-rank 3b
indexing bug (`kappas_array[i]` with the `i_beg`-relative form commented out,
never exercised because every test is single-rank).

Not on the critical path: nested-sampling extraction (§5), and the six remaining
preprocessor-interrupted statements on the GPU branch, which sit outside the exp
blocks in the earlier prediction code.
