# Developing TurboGAP

Where code goes and how to prove it works. Commands for building, testing and
committing are in `docs/BUILD_AND_TEST.md`.

## Layout

| layer | files | holds |
|---|---|---|
| driver | `src/turbogap.f90` | the step loop, one call per phase; no MPI, no `#ifdef` |
| features | `src/turbogap_<feature>.f90` | one module per feature, owning its state; called at fixed hooks |
| potential | `src/turbogap_evaluate.f90` | every energy and force term, in order |
| decomposition | `src/turbogap_domain.f90` | which rank owns which atoms; every broadcast and reduction of state and results |
| transport | `src/turbogap_comm.f90` | `comm_t`, `comm_bcast`, `comm_sum_to_root`, `comm_sum_all`, `comm_allgather`; serial stubs |
| kernels | `gap.f90`, `vdw.f90`, `electrostatics.f90`, `exp_utils.f90`, ... | the physics, on plain arrays |

One step: `structure_acquire` → `domain_sync_state` → `domain_build` →
`evaluate` → `md_step` or i-PI → `nested_step` → `mc_step` → `domain_sync_state`
→ `loop_end_step`.

## Adding an energy term

1. Kernel in its own module, taking plain arrays, no MPI, no TurboGAP types.
2. Its result arrays in `results_t` (`turbogap_results.f90`): allocated and
   zeroed in `results_prepare`, freed in `results_free`.
3. Call it from `evaluate`, in its place in the order.
4. Partial sums across ranks: add the family to the table in
   `domain_complete_contributions`. A global scalar the term needs mid-way goes
   through `comm_sum_all` (see `compute_estat`).
5. Keyword: default in `types.f90`, parse branch with its `!>` block in
   `read_files.f90` (`docs/keywords-howto.md`), `make docs`.
6. `SRC` in the `Makefile`; regenerate `makefiles/Makefile.deps` and
   `Makefile.deps.gpu`.

A feature that is not a term (a sampler, a bias, an output) gets its own module
with init/step/finish hooks, called from the driver in the phase it belongs to.

## Rules

- New code reaches MPI only through `turbogap_comm`. Older code still calls it
  directly: `turbogap_setup`, `turbogap_vdw`, `exp_interface`, and one call
  each in `error`, `gap_interface` and `gpu_context_cpu`.
- A rank holds centres `i_beg:i_end` and their pairs `j_beg:j_end`.
  `neighbors_list` is global (fold with `modulo(j - 1, n_sites) + 1`); arrays
  passed as a slice take the local index. Mixing the two is correct on rank 0
  and wrong everywhere else, so test at three ranks.
- A refactor is bit-exact at every rank count. A change in results is its own
  commit, with a test that shows the new result is right.

## Proving it

- Forces and virial: a finite difference of the term's energy, with
  `tests/xrd_debye/fd_gradient.py`. Isolate the term by differencing it on and
  off (`--family estat`, `vdw`) or by a potential holding only it
  (`tests/fd_gap`). Run at one rank and at three (`--ranks 3`).
- Regression: `tests/regression/run.sh`, bit-exact on the host. A new path gets
  a case; `REFERENCE=golden` when the baseline cannot run it.
- Device: `run.sh --gpu`, within a tolerance. A case may loosen it only with a
  host control -- the input shifted by 1e-9 A, host against host -- showing the
  drift is amplification (KNOWN_ISSUES 13d). Device builds are not
  bit-reproducible, so a device change is judged by the suite, not by `cmp`.
- Known defects go in `tests/regression/KNOWN_ISSUES.md`, with how to reproduce
  them.
