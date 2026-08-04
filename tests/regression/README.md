# Refactor regression suite

Guards the `turbogap.f90` modularization. The contract is **bit-identical
output** against the pre-refactor baseline: moving a block into a procedure
must not reorder floating-point operations, so not one digit should move.
That is why the comparator is `diff` and not a tolerance check — it is exact,
needs no tuning, and cannot be quietly loosened as the refactor proceeds.

A change that legitimately alters results is not a refactor and does not
belong on a commit checked by this suite.

## Running

```sh
make DEBUG=0 all                       # or whatever your arch build is
TURBOGAP_DATA_ROOT=$HOME/work/cpu_vs_gpu_tests/input \
  tests/regression/run.sh              # all cases
tests/regression/run.sh vdw_mbd        # one case
tests/regression/run.sh --list         # case names
```

Non-zero exit means at least one case differs from the baseline.

| variable | meaning |
|---|---|
| `TURBOGAP_BIN` | binary under test (default `bin/turbogap`) |
| `TURBOGAP_REF_BIN` | reference (default `tests/regression/baseline/turbogap.e6eb1aa`) |
| `TURBOGAP_DATA_ROOT` | directory holding the test systems |
| `TURBOGAP_TIME_TOL` | fail if test/ref wall-clock exceeds this ratio |
| `TURBOGAP_KEEP` | keep staging dirs for inspection |
| `TURBOGAP_BLESS` | regenerate `expected/` for golden cases |

## Two kinds of case

Most cases compare against the frozen baseline binary, and must be
bit-identical to it.

A case whose `case.conf` sets `REFERENCE=golden` compares against a checked-in
`expected/` directory instead. That is for cases the baseline cannot serve as a reference for: either it
crashes on the input, or the branch has deliberately changed the result. All
five ts+mbd cases are golden for one of those two reasons; everything else is
still compared bit-for-bit against the baseline. These are
**characterization** tests: they pin behaviour so a later refactor cannot drift
it silently, and assert nothing about whether the physics is right. Regenerate
them deliberately with `TURBOGAP_BLESS=1`, never to make a red suite go green.

Expected outputs are stored per rank count, because MPI reductions are not
associative and rank counts legitimately differ in the last digit or two.

## The baseline

`baseline/turbogap.e6eb1aa` is master HEAD at the point the refactor started;
`baseline/COMMIT` records the parent and submodule SHAs it was built from.
**Do not regenerate it from a refactored tree** — it is the only thing that
knows what the code did before, and rebuilding it from the branch under test
would make every case pass vacuously.

## Cases

| case | ranks | covers |
|---|---|---|
| `vdw_ts` | 1 | TS pairwise vdW, Hirshfeld local properties |
| `vdw_mbd` | 1 | many-body dispersion |
| `vdw_mbd_cell_mpi2` | 2 | MBD + the cross-rank Hirshfeld gradient exchange |
| `vdw_tsmbd_md` | 1 | ts+mbd: both the recompute and the reuse branch (golden) |
| `vdw_tsmbd_md_mpi2` | 2 | the same under MPI (golden) |
| `vdw_tsmbd_predict` | 1 | ts+mbd outside MD (golden) |
| `vdw_tsmbd_predict_mpi2` | 2 | the same under MPI |
| `vdw_tsmbd_mc` | 1 | ts+mbd under mc, which also never advances md_istep |
| `xps_predict` | 1 | SOAP loop, local properties, 2b/3b/core_pot, XPS spectrum |
| `xps_predict_mpi4` | 4 | the energy/force/virial MPI reductions |
| `gcmc_xps` | 1 | MC insertion/removal: `n_sites` changes, reallocation, neighbor rebuilds |
| `gcmc_xps_mpi2` | 2 | the same under MPI |
| `co_predict` | 1 | 7176-atom single point; timing reference |
| `co_predict_mpi4` | 4 | the same under MPI; timing reference |
| `co_md` | 1 | velocity Verlet, Berendsen thermostat, per-step rebuild decision |

Every case fixes `random_seed`, and `co_md` starts from explicit velocities,
so all of them are reproducible run to run. This was verified before the suite
was trusted: baseline vs. an identical rebuild passes every case.

Test data lives outside the repo (it is ~35 MB of GAP files) under
`$TURBOGAP_DATA_ROOT`:

- `vdw_P/` — P4 dimer and 108-atom cell, from `TurboGAP_School/vdw_interactions`
- `xps_opt/` — 512-atom CO with core-electron binding energies, from
  `tutorials/xps_optimization`
- `CO/` — the 7176-atom CO system also used by the GPU comparison, so timing
  numbers stay comparable with `docs/gpu_fixes_handoff.md`

A case whose data is absent is skipped, not failed.

## Coverage gaps

Read `KNOWN_ISSUES.md` before trusting a green run. The five ts+mbd cases are
characterization tests, so for those a green run means "unchanged", not
"correct". Still dormant in the vdW path: the `S_xyz_inv` strain term and the
`this_local_virial_vdw_diag` correction.
