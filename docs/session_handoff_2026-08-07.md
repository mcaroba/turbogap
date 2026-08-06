# Session handoff — 2026-08-06/07

Work on both trees on `alt`. Nothing pushed to any remote. Companions:
`docs/refactor_strategy.md` (what to do and why), `docs/session_handoff_2026-08.md`
(the previous session), `docs/from_turbogap_2.0.md` (what to take from the 2.0
rewrite and what not to), `docs/TODO.md` (`state_t`/`simulation_t` and the SOAP
seam), `docs/BUILD_AND_TEST.md` (every one-liner).

| | |
|---|---|
| CPU branch | `refactor/modularize`, head `06f6a18`, 17 commits this session |
| GPU branch | `gpu_fixes`, head `c1165f8`, 20 commits this session |
| CPU suite | **21 passed, 0 failed**, bit-exact |
| GPU suite | **3 passed, 0 failed, 2 xfail** |
| `turbogap.f90` divergence | 1443 → **1081** |
| byte-identical files | 5 → **9** |

---

## 1. Read this first: how to read the divergence numbers

**A file-level diff cannot tell a conflict from an addition, and the gap is
often an order of magnitude.** Use `tools/proc_diff.py`, which attributes
divergence to individual procedures:

```sh
python3 tools/proc_diff.py ../turbogap_master_2026/src/X.f90 ../turbogap_gpu_commit_mahti/src/X.f90
```

| file | `diff -w` | real conflict |
|---|---|---|
| `exp_interface.f90` | 2414 | **81** |
| `turbogap_exp.f90` | 1121 | 422 |
| `turbogap.f90` | 1081 | — (see below) |
| `exp_utils.f90` | 1034 | 349 |
| `gap.f90` | 536 | 384 |
| `electrostatics.f90` | 374 | **0** |
| `read_files.f90` | 364 | 190 |
| `neighbors.f90` | 338 | **3** |
| `local_properties.f90` | 305 | **0** |
| `mc.f90` | 295 | 161 |
| `gap_interface.f90` | 279 | 213 |
| `turbogap_setup.f90` | 220 | 56 |
| `types.f90` | 161 | **0** |
| `turbogap_md.f90` | 157 | 74 |
| `turbogap_estat.f90` | 138 | 87 |
| `gpu_context.f90` | 132 | 29 |
| `md.f90` | 51 | 5 |
| `mpi.f90` | 43 | 9 |
| `xyz.f90` | 6 | 5 |

Byte-identical across the trees: `error`, `kinds`, `misc`, `printing`,
`read_utils`, `splines`, `timing`, `turbogap_vdw`, `vdw`.
Master-only: `gap_backend_cpu.f90`. GPU-only: `fortran_cuda_interfaces.f90`,
`gap_backend_gpu.f90`.

**Four files show a large diff and zero real conflict** — `electrostatics.f90`,
`local_properties.f90`, `types.f90` and (nearly) `neighbors.f90`. Their whole
figure is GPU-only procedures that merge by appearing. **They can merge today.**

Two caveats on the metric:

* **`turbogap.f90` reports 0 because `proc_diff` finds no shared *procedures*
  in a program file.** Ignore that cell; use the file diff there.
* **Once a file contains a seam, the total conflates two things.** `gap.f90`'s
  384 is 371 for `get_soap_energy_and_forces` (the body, which relocates at the
  end) plus 13 for the `soap_backend_begin`/`end` pair, which are the interface
  boundary and are *supposed* to differ. Read per procedure, not the total.

---

## 2. What changed

### 2.1 Structure — four seams now exist

The pattern established by `gap_backend` (one module name, two implementation
files, the Makefile picks one) is now used four times:

| seam | CPU | GPU |
|---|---|---|
| `gap_backend` | `gap_backend_cpu.f90` | `gap_backend_gpu.f90` |
| `gpu_context` | stub, empty bodies | streams, handles, batch storage |
| `turbogap_estat` | `compute_estat` | same + batched device path |
| `soap_backend_begin/end` | points at the sparse set | uploads it |

`turbogap_estat.f90` took 92 lines out of master's driver and 201 out of the
GPU's. `gpu_context_init`/`finalize` took ~50 lines of device set-up and
teardown out of the divergence.

**The `#ifdef`-split argument lists are gone.** `refactor_strategy.md` §6.4
called these the gate on any Fortran tooling working at all — 5 on master, 11
on the GPU branch, "no Fortran-aware tool can handle these". **Both drivers are
now at zero and both parse under fprettify.**

### 2.2 Bundling, from turbogap_2.0

* `times_t` — 22 (CPU) / 29 (GPU) timer declarations → 1; 13 module arguments
  → 5. `src/timing.f90` is byte-identical on both trees.
* `perform_t` — the exp-observable decisions evaluated once.
* Adjudication of everything else in `docs/from_turbogap_2.0.md`, including why
  `tg_memory` is **not** worth adopting (it collides with `c_loc`, used at ~30
  sites on the branch where the bit-exact contract does not hold).

### 2.3 Defects found and fixed

Nine, of which these are the ones that were silently wrong in shipped code:

1. **Electrostatics segfaulted on first use, on master.** `energies_estat` and
   `forces_estat` were never allocated — only the `this_` copies — so the MPI
   unpack wrote through an unallocated pointer. Unconditional for any deck with
   `estat_method /= none`.
2. **…and once that was fixed, the whole contribution was discarded.** The
   energy sum omitted `energies_estat`; the force/virial sums omitted theirs.
3. **Ten allocator mismatches on the GPU branch** — every remaining plain
   `gpu_free`. `gpu_malloc_all_blocking` is called *zero* times, so every
   `hipFree` was releasing a pool pointer. Two of them freed the source of an
   async device→host copy, valid only because `hipFree` synchronises the whole
   device.
4. **The exp-observable guards disagreed three ways.** Allocation asked
   `do_X .and. valid_X`, zeroing asked `do_X .and. exp_forces .and. valid_X`,
   accumulation asked `exp_forces .and. valid_X` — dropping `do_X`. Those are
   independent inputs, so a deck supplying exp data for an observable it had
   not switched on accumulated an unallocated array.
5. **The negative "Miscellaneous" timer** — `time_gap` accumulated *inside* the
   SOAP loop from a stamp taken before it, over a window that also enclosed an
   `mpi_reduce` already charged elsewhere.
6. **Four keyword-parity gaps** — `mbd_correction_freq`, `print_estat_forces`,
   the four XPS spectrum parameters, and the five `gpu_*` on master. In each
   case a deck setting them was silently ignored.
7. **The GPU branch was missing the standalone XPS spectrum entirely** — not
   broken, absent. Found because `proc_diff` reported two master-only
   procedures in `exp_utils.f90`.

### 2.4 Coverage

* `estat_gsf` (both suites) — the first time electrostatics has *ever* run on
  either branch. Uses the CCLi potential, the only one in the test data
  declaring an `atomic_charge` local property.
* `mad_sf_matrix` and `mad_sf_matrix_mpi2` — the structure-factor matrix path
  with forces, and the two-rank companion.
* CPU suite 18 → **21 cases**.

### 2.5 Tooling

* `tools/setup_dev_env.sh` (both) — rebuilds the hook environment; `--check`
  verifies version, pre-commit, installed hook and PATH.
* `tools/gpu_check_alloc_pairs.py` (GPU) — static allocator/deallocator
  pairing. **This is the check that matters when moving an allocation and its
  free between modules.**
* `tools/gpu_memcheck.sh` (GPU) — compute-sanitizer runner. memcheck and
  initcheck are both clean on `vdw_ts`.
* `tools/check_reformat_only.py` (GPU) — token-stream comparison, because
  `diff -w` cannot verify a C reformat.
* `clang-format` pre-commit hook for the CUDA/C/C++ sources.

---

## 3. What to do next, in order

### 3.1 Merge what is already converged — do this first, it is free

`electrostatics.f90`, `local_properties.f90`, `types.f90` and `neighbors.f90`
have **zero or three** lines of real conflict. Nothing blocks them.

### 3.2 Finish the SOAP seam

Four increments are done (`ef58814`, `a0dbf1c`, `3bc43b0` on GPU, `d005810` on
CPU). `get_soap_energy_and_forces` is down from eleven extra arguments to six,
and all six are now GPU-side *additions* rather than differences in what the
CPU needs.

**The order matters and the obvious next step is wrong.** See `docs/TODO.md`:
`n_neigh_d` is allocated inside the **`soap_turbo` submodule** and handed back
out, so giving the three remaining buffers an owner in `gap.f90` needs either
the vendored submodule writing into `gap`'s variables or a setter that is a
wash. **Seam `get_gap_soap` first** — one procedure, 410 vs 485 lines, 213 of
conflict, the last structural difference in the SOAP path. Then:

1. `soap_d` / `soap_cart_der_d` / `n_neigh_d` get an owner for free.
2. `n_sparse`, `n_pairs`, `solo_time_soap` are bookkeeping.
3. Move `get_soap_energy_and_forces` into `gap_backend_cpu`/`gap_backend_gpu`.
   **`gap.f90` converges at that point, not before** — its 371 relocates in one
   move.

### 3.3 The two remaining exp files

`turbogap_exp.f90` (422, all `compute_exp_spectra`'s batched path) and
`exp_utils.f90` (349, all `get_structure_factor_forces_matrix`). Both are one
procedure each; both are the same shape as the SOAP seam.

### 3.4 Deliberately not done

* **The `energy_estat=` xyz field.** The GPU branch writes it unconditionally;
  master does not. Adding it changes every trajectory header and takes all ten
  baseline-compared cases red at once, converting them to golden — a real
  reduction in the safety net. It is a deliberate result-changing commit with
  its own blessing, not something to slip in.
* **The GPU/CPU electrostatics disagreement** (below). Owner deferred it.
* **`simulation_t` and overlapping domain decomposition** — designed in
  `docs/TODO.md`, not to be started before the merge lands.

---

## 4. Open, and owed a decision

**The batched device electrostatics disagrees with the CPU implementation.**
`estat_gsf` is xfail on the GPU suite: forces to 1.1 eV/Å against |F|max 20.5,
virial 0.7%, `local_energy` 0.16 eV — while `energy_soap`, `2b`, `3b`,
`core_pot` and `vdw` agree to the last digit. That exactness elsewhere is what
makes the electrostatics the suspect rather than the comparison.

Which side is wrong is open. The CPU path had four defects of its own until
this session, and **neither implementation has ever had a reference**. The
`estat_gsf` case is the instrument for settling it.

**The virial debug block is guarded by `params%print_vdw_forces`** but dumps
the estat, soap, 2b, 3b, core_pot and xrd virials — one flag with six meanings.
Ported verbatim so the files converge; giving it its own flag changes what an
existing keyword does.

---

## 5. Traps that cost time here

Each of these was hit. None is obvious in advance.

**`open(p, 'w')` truncates before its argument is evaluated.** A `NameError`
while building the replacement string left `exp_utils.f90` at **zero bytes**.
Recovered with `git checkout --` because it was committed. Build the whole text
first, write to a temp file, `os.replace`. And use single-quoted heredocs or
`scp` the script — a double-quoted `ssh "..."` heredoc eats `\n`.

**A verbatim transplant only converges anything if it lands where its
counterpart already is.** Twice: `read_options_gpu` dropped at the first
plausible spot took `read_files.f90` from 389 to **444**; in the right place,
366. The XPS helpers appended at the end of the module took `exp_utils.f90`
from 1086 to **1145**; placed where master has them, 1034.

**§6.4's `#ifdef` hazard is real.** "A transformation that scans for the first
`#else` will swallow it and silently make the allocation unconditional." A
non-greedy `.*?` from `#ifdef` to `#else` runs past an intervening `#endif`.
The regex must refuse to cross a directive line, and the transform must assert
that groups removed equals groups collapsed.

**Write the check so its failure mode is a false alarm, not a false pass.**
Both new GPU tools first reported a clean run that was not — one captured
`Qs_d)` with the trailing paren and split one buffer into two entries; the
other grepped six `=` where compute-sanitizer's records use nine.

**`/usr/bin/compute-sanitizer` on alt is a stub** that dies with "Unable to
find injection library" while still producing plausible output. Use
`/usr/local/cuda*/bin/compute-sanitizer`.

**`ERROR SUMMARY` counts cuBLAS instrumentation notices**, which fire on every
run touching `cublasDgemv` — i.e. every SOAP prediction. Counted with real
faults they make "1 error" the normal state.

**A depth counter that stops at the first `end if` reaching depth 0 lands one
`end if` early** on an `if/else/end if/end if` block, on both branches. The
symptom is a stray `end if` and "Unexpected ELSE statement" from a *later*
construct.

**`gap.f90` has no default private**, so `use gpu_context, only: ...`
re-exports and collides with `gap_backend_gpu`'s own `gpu_stream`. And the
`private ::` that fixes it is an attribute declaration, so it must follow
**every** USE — just before `contains`, not beside the import.

**`gap_interface.f90` says things twice.** The `get_gap_soap` argument list
appears again in a commented-out call; `alphas(:)`/`Qs(:,:)` are declared in
two procedures. Skip lines whose first non-space char is `!`; slice to the
procedure for declarations.

**`makefiles/Makefile.deps` must be regenerated whenever a source file is added
OR removed, and whenever a `use` is added.** Otherwise `make -j` races on
`kinds.mod`. `python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps`.

**`TURBOGAP_REF_BIN` must be an absolute path** — it is resolved from the
staging directory, and a relative one fails every case with exit 127 while
looking like a harness problem.

**Do not let a rebuild race a running suite.** Stage the binary:
`cp bin/turbogap /tmp/tg_under_test && TURBOGAP_BIN=/tmp/tg_under_test ...`.

**fprettify rewrites files during `git commit` and the commit then aborts.**
That is the normal flow. Confirm with `git diff -w --ignore-blank-lines` that
the rewrite was whitespace-only, **rebuild and re-run the suite**, then commit
again.

**Commit messages containing parentheses break `ssh` heredocs.** Write to a
file and use `git commit -F`.

**Do not carve small cells out of `CLi/atoms.xyz`.** The 125k-atom parent has
sub-1.25 Å pairs and Li in pockets: a sub-box needs ~30% of its atoms pruned
and a cluster inherits the contacts (total energy came out at 5.8e7 eV). Use
`C/atoms_897.xyz` with the CLi potential instead — that is what `estat_gsf`
does.

---

## 6. How to check any of this

See `docs/BUILD_AND_TEST.md` for the full set. The short version:

```sh
export M=.../turbogap_master_2026 G=.../turbogap_gpu_commit_mahti
export PATH="$HOME/.venvs/turbogap-tools/bin:$PATH"     # or tools/setup_dev_env.sh

cd $M && make -j12 && tests/regression/run.sh           # 21, bit-exact

cd $G && export HOP_ROOT=/u/74/zarrout1/unix/work/hop \
      && export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda \
      && make -j12 DEBUG=0 \
      && TURBOGAP_REF_BIN=$M/bin/turbogap tests/gpu/run_regression.sh

cd $G && python3 tools/gpu_check_alloc_pairs.py         # 0 mismatched
cd $G && tools/gpu_memcheck.sh --tool initcheck         # clean
```

The GPU build is **not** bit-reproducible run to run, so the bit-exact contract
cannot be used there. Do each change on master first where it holds, and port
the verified result. That is not a preference; it is the only place the work
can be checked.
