# Session handoff — 2026-08-05/06

Work done across both trees on `alt` (`d22-0052`). Nothing has been pushed to any
remote. Companion documents: `docs/refactor_strategy.md` (what to do and why),
`docs/refactor_handoff.md` (the CPU-side de-monolithing), and
`docs/gpu_fixes_handoff.md` (the GPU correctness work).

| | |
|---|---|
| CPU branch | `refactor/modularize`, head `425f34f`, 9 commits this session |
| GPU branch | `gpu_fixes`, head `9bbebb1`, 11 commits this session |
| GPU submodule | `soap_turbo` @ `9db51a0`, parent records it |
| CPU suite | **18 passed, 0 failed**, bit-exact |
| GPU suite | **3 passed, 0 failed, 1 xfail** (`XRD_mad`, known `local_energy` drift) |

---

## 1. Read this first: five traps that report success while being wrong

Every one of these was hit during the session. They cost the most time and they
are the least likely to be re-derived.

### 1.1 "Block-private" is not "safe to make a local" when the block is in a loop

Extracting `compute_md`, nine driver locals were used *only* inside the block, so
the coupling analysis called them block-private and the first version made them
subroutine locals. **The build was clean and five regression cases failed** —
`co_md`, `vdw_tsmbd_md`, `vdw_tsmbd_md_mpi2` and the two XRD MD cases. Every
non-MD case passed.

The block runs inside the MD loop. A driver local that nothing outside touches is
still *per-run* state that persists from one step to the next; a subroutine local
resets on every call. Six of the nine were like that, and the test that finds them
is **read before written within a single call**:

```
cum_eel  gd_box_do_pos  gd_istep  restart_box_optim  target_temp  time_step_prev
```

`refactor_handoff.md` §8's leaked-value check asks whether the *driver* still
needs the value afterwards. This is the mirror question — whether the *block*
needs its own value next time round. Ask both.

### 1.2 A bare `stop` is exit status 0

Three separate places reported failure and exited successfully: `turbogap_abort`
as taken from turbogap_2.0, the unrecognised-keyword branch of `read_input_file`,
and (previously recorded) the debugging `exit(0)` in `gpu_exp.cu`. Nothing
checking a return code — a script, a queue system, the suite — could tell. All
three now exit non-zero. `turbogap_abort` uses `MPI_abort`, not `MPI_finalize`:
finalize is the orderly shutdown and expects every rank to call it, so one rank
finalizing alone leaves the others in the next collective.

### 1.3 A single system size cannot see a whole class of bug

`gap_interface.f90` counted neighbour pairs with `<` and filled them with `<=`,
so a pair exactly on `rcut_max` was written one element past six arrays. On the
7176-atom CO system the overrun lands in memory the allocator already handed out:
it corrupts the heap **silently**. Every test on the GPU branch used that one
system. On the 8-atom P4 dimer it aborts instantly with `malloc(): corrupted top
size`. The suite now has a small cell (`vdw_ts`) for exactly this reason.

### 1.4 One condition, two predicates that are free to disagree

Four defects this session, on top of the two already recorded (`ts+mbd` §5.1, the
MPI reduce triple walk §7.1):

* `gap_interface.f90` counting `<` against filling `<=` (§1.3 above).
* `has_vdw` against `has_local_properties`: a GAP declaring the deprecated inline
  `has_vdw = .true.` form never got a `local_property_models` entry, so `has_vdw`
  was true while `has_local_properties` was false, and the TS block read an array
  the allocation guard had skipped.
* The electrostatics block was entered on `estat_method` alone, without
  `valid_estat_charges`, so a deck asking for electrostatics against a GAP with no
  `atomic_charge` property indexed `local_properties` with an uninitialised
  `charge_lp_index`. Segfault, reproduced.
* `accessible_volume` and `print_vdw_forces` each had two handlers in the keyword
  chain; the second was unreachable.

When you find yourself writing a condition twice, that is the bug.

### 1.5 Measurement artefacts that look like divergence

* **A whitespace-insensitive comparison must strip whitespace entirely, not
  collapse runs of it.** `exp_interface.f90` read as 757 diverged lines where the
  truth was 67, because one copy had been through fprettify and writes `(:, :)`
  for `(:,:)`.
* **`git diff -B` is `--break-rewrites`, not "ignore blank lines".** It reported
  ~19,000 non-whitespace changes in a pure reformat. Use
  `git diff -w --ignore-blank-lines`.
* **A crossing count must discount declaration *continuation* lines**, not just
  the first line of a declaration. 89 arguments against a real 64 on the same
  block.
* **A file-level diff cannot tell a conflict from an addition.**
  `exp_interface.f90`'s 2193 was 1898 lines of GPU-only procedures that merge by
  appearing, plus 55 of actual conflict.

---

## 2. What changed

### 2.1 Build and tooling

* **`make -j` works on both trees** (`204fbb8` cpu, `6a2f772` gpu). It never did:
  the Makefiles expressed no Fortran module dependencies. `tools/gen_fortran_deps.py`
  reads the Makefile's own `SRC_*` lists and pattern rules, so one tool serves both
  trees — which matters, because master links `OBJ_STOP` and the GPU tree does not.
  Serial timing unchanged; `-j12` takes master 43 s → 26 s and the GPU tree 29 s → 15 s.

  Two traps in that change, both of which report success: the `include` must come
  **last** or a generated dependency line becomes make's default goal and `make`
  builds one object and exits 0; and duplicate module definitions must be gone
  first, which is why the dead `*_bak.f90` files had to be removed before it.

* **Phase 0a closed** (`856b860` submodule, `6a2f772` gpu). The three
  `soap_turbo/src/*_bak.f90` are gone — 2710 lines, in no Makefile, `USE`d by
  nothing, and each defining the *same module name* as its live counterpart.
  `orthonormalization_kernels.cc` is kept by owner decision, and the Makefile now
  says so.

* **pre-commit with fprettify on both trees** (`ca390bc` cpu, `415f881` gpu),
  `--indent 3 --line-length 132`, plus the venv at `~/.venvs/turbogap-tools`
  (alt is PEP-668 externally managed; do not use `--break-system-packages`).

  This **reverses** `refactor_strategy.md` §6.3's "do not reformat". That
  measurement asked whether formatting *alone* converges the branches — it does
  not, and that still stands — but the real choice was formatting *neither* tree
  or *both with identical settings*. Done on the pair, the raw cross-branch diff
  over the shared modules fell **26786 → 7828 lines, −71%**, landing on the
  `diff -w` floor. `git merge` no longer needs `-X ignore-all-space`.

  Three things had to be fixed first: the hygiene hooks were rewriting **test
  data** (`silicon.gap`, `.beta` tables, a case's input deck), so they are
  allow-listed to `src/ tools/ docs/ makefiles/` rather than excluded case by
  case; fprettify and `trailing-whitespace` never converged on `eph_fdm.f90`, so
  Fortran is excluded from the fixer and fprettify owns `.f90`; and fprettify was
  not idempotent on the GPU `turbogap.f90` until eight commented-out argument
  lines inside a continued argument list were removed and one call right-aligned
  to exactly column 132 was rewrapped.

* **Test data moved to its own repository.** `tests/fetch_test_data.sh` on both
  trees clones `https://github.com/TiganyZ/turbogap_tests` into
  `<repo>/../turbogap_tests` on first use, and both run scripts call it, so a
  fresh checkout needs no environment variable.

### 2.2 Merging the branches

* **The GPU diffraction block is out** (`a84be03`, `ca6ee53`). It was recorded as
  blocked on "71 crossing variables". That figure conflated 28 declaration-only
  crossings, four procedures that arrive via `use`, and three `params%` fields.
  The genuinely shared device context is **eight symbols**, now
  `src/gpu_context.f90`, reached by `use` so the lift stayed verbatim.
  `turbogap.f90` there went 4793 → 4101.

  `omp_task` deliberately stays out of that module: it reads like context and
  crosses like context, but the code intends it OpenMP-*private*, and shared
  module state would be a bug the moment `-fopenmp` is enabled — which it is not,
  so nothing would catch it.

* **`vdw.f90` adopted from master by the GPU branch** (`46be238`): divergence
  6398 → **0**. Also brought `misc.f90` and `src/third_party/nnls/`.

* **Electrostatics merged to master** (`d64d7f4`). 23 of the module's 24
  procedures have no device token; the one that does,
  `calculate_batched_electrostatics`, stays on the GPU branch. `N_CONTRIB` 10 → 11
  with `C_ESTAT` inserted at 3 and the rest renumbered **to match** the GPU
  branch rather than appending.

* **Input-keyword parity.** Thirteen vdW/MBD keywords master had and the GPU
  branch did not (`39a0b6d`), and the eleven `estat_*` the other way. The two
  trees now recognise the same keyword set.

### 2.3 Input handling

* **Parameters are echoed and validated** (`0833120` cpu, `d961b11` gpu), from
  turbogap_2.0's `src/read/read_utils.f90` and `src/utils/printing.f90`. New
  `src/{printing,error,read_utils}.f90`, byte-identical on both trees. ~190
  handlers per tree call `check_iostatus` and echo through `print_parameter`;
  `species`, `masses` and `radii` go through `read_parameters`, which counts the
  fields on the line first.

  Two corrections to what was taken: the field counter incremented on *every*
  blank, so `species = C  O` over-counted by one per repeated separator; and
  `read_parameters_char8` stays out of the generic interface because its dummy is
  not distinguishable from the `character(len=*)` one.

  **Behaviour change:** a malformed value used to be tolerated, because
  `iostatus` was set and never inspected, and the parameter kept its default. It
  is now fatal.

* **Keywords split into families** (`650e28d` cpu, `a711685` gpu) and
  **alphabetised** (`01b577f`, `7fd6ce5`). The ~1500-line `if/else-if` chain
  became `read_options_{general,control,md,nested,mc,vdw,estat,exp,stopping,
  local_properties,output}`, plus `read_options_gpu` on the GPU branch only, so
  that branch's five GPU-only keywords sit in one place. **Six families are now
  character-identical between the trees.**

### 2.4 The driver

* **`compute_md`** (`425f34f` cpu, `9bbebb1` gpu) into `src/turbogap_md.f90`.
  Drivers: master 3682 → 3329, GPU 4247 → 3885.

---

## 3. Where the branches stand

Whitespace-insensitive, largest first:

| file | master | gpu | real diff |
|---|---|---|---|
| `exp_interface.f90` | 1947 | 3760 | 2416 — but only ~55 is conflict (§1.5) |
| `turbogap.f90` | 3329 | 3885 | 1443 |
| `turbogap_exp.f90` | 581 | 1046 | 1098 |
| `exp_utils.f90` | 3413 | 4125 | 1086 |
| `gap.f90` | 948 | 1153 | 477 |
| `read_files.f90` | 4006 | 4022 | 394 |
| `electrostatics.f90` | 840 | 1197 | 373 |
| `turbogap_md.f90` | 524 | 521 | 157 |
| `vdw.f90` | 6938 | 6938 | **0** |

Master-only files: `gap_backend_cpu.f90`, `turbogap_vdw.f90`.
GPU-only: `fortran_cuda_interfaces.f90`, `gap_backend_gpu.f90`, `gpu_context.f90`.

---

## 4. What is left, in the order the measurements support

1. **A charge-carrying GAP.** No test data declares an `atomic_charge` local
   property, so **electrostatics has never been exercised on either branch**.
   Master's implementation is correct by construction (it is the GPU branch's,
   minus the device procedure) but untested. This is the cheapest way to turn a
   whole subsystem from assumed to checked.

2. **The `energy_estat=` output field.** The GPU branch writes it unconditionally
   in `xyz.f90`; master does not. Adding it changes every trajectory header and
   takes all ten baseline-compared cases red at once, so it is a deliberate
   result-changing commit with its own blessing — not something to slip into
   another change.

3. **`exp_interface.f90`**, nominally the largest divergence and actually the
   smallest: 1898 lines of GPU-only procedures that merge by appearing, 306 lines
   of a master-only procedure the GPU branch is owed, and ~55 lines of genuine
   conflict in two procedures — mostly the inherent device plumbing. Three
   non-GPU items in there are worth their own commits: the GPU-only
   `virial_sf = 0.d0` against master's `4f99e92`, and
   `get_structure_factor_forces_matrix` taking `i_beg, i_end` and
   `species(i_beg:i_end)` on one side only, which is an MPI-decomposition
   difference in a routine every test exercises single-rank.

4. **`turbogap_exp.f90`** at 1098. The GPU copy carries the batched device paths;
   the CPU copy does not. Same shape as `exp_interface`.

5. **Nested sampling** remains the last big block in the driver — 1099 lines, but
   **145 arguments** against 23 block-private locals. Bundling the contribution
   arrays (38 of those arguments) and the timing buckets (22) would bring it to
   ~87, which is still a poor extraction. That bundling is a different kind of
   change from a lift, and it has no merge payoff on its own. It also has no
   regression coverage.

6. **`local_energy` drift on the GPU branch**, ~1e-5 by frame 5, the only thing
   keeping `XRD_mad` xfail. Deferred deliberately by the owner.

Cheap, any time: master's `Makefile` still needs the `SRC_STOP` dependency
treatment if `OBJ_STOP` is ever linked on the GPU side; both trees have untracked
`compile.sh` variants that have been left alone throughout.

---

## 5. How to check any of this

```sh
# CPU branch — 18 cases, bit-exact against a frozen baseline binary
cd turbogap_master_2026 && make -j12 && tests/regression/run.sh

# GPU branch — compares against the CPU build
cd turbogap_gpu_commit_mahti
export HOP_ROOT=/u/74/zarrout1/unix/work/hop
export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda
make -j12 DEBUG=0                     # DEBUG defaults to 1: a 2.1x slowdown
TURBOGAP_REF_BIN=../turbogap_master_2026/bin/turbogap tests/gpu/run_regression.sh
```

Test data is fetched automatically on first run. `pre-commit` needs
`~/.venvs/turbogap-tools/bin` on `PATH`.

The GPU build is **not bit-reproducible run to run**, so the bit-exact contract
cannot be used there — do each extraction on master first, where it holds, and
port the verified result. That is not a preference; it is the only place the work
can be checked.

Debugging on `alt`: `ptrace_scope = 2`, so gdb cannot attach. `-fcheck=bounds` is
the tool of choice — patch the `DEBUG=1` branch of `makefiles/Makefile.Aalto_*`
and restore it afterwards. It named the `get_gap_soap` overrun exactly.

---

## 6. Tools left behind

In `tools/` on both trees unless noted. Each exists because the check it performs
is easy to get wrong by hand.

| | |
|---|---|
| `gen_fortran_deps.py` | module dependencies from the Makefile's own lists |
| `classify.py` | a block's arguments / locals / free names, continuation-aware |
| `proc_diff.py`, `proc_show.py` | compare two Fortran files procedure by procedure |
| `split_input_keywords.py` | the keyword-family split |
| `align_input_keywords.py` | alphabetise families, drop duplicate handlers |
| `extract_md.py` | the `compute_md` lift (per-branch argument sets) |
| `extract_spectra.py`, `find_orphans.py`, `drop_orphans.py` | GPU branch, the diffraction lift and its orphaned declarations |

A note on writing more of these: collecting `re.finditer` matches up front and
then reassigning the source string inside the loop invalidates every later
offset. It looks like it worked — three keyword families came out unsorted and
three that had been identical between the branches were left diverged. Collect
edits against the original text and apply them **last-to-first**.
