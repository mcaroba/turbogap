# Refactor regression suite

## Running

```sh
make DEBUG=0 all                       # or whatever your arch build is
TURBOGAP_DATA_ROOT=$HOME/turbogap_tests/input \
  tests/regression/run.sh              # all cases
tests/regression/run.sh vdw_mbd        # one case
tests/regression/run.sh --list         # case names
```

Non-zero exit means at least one case differs from the baseline.

| variable             | meaning                                                          |
| -------------------- | ---------------------------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)                       |
| `TURBOGAP_REF_BIN`   | reference (default `tests/regression/baseline/turbogap.e6eb1aa`) |
| `TURBOGAP_DATA_ROOT` | directory holding the test systems                               |
| `TURBOGAP_TIME_TOL`  | fail if test/ref wall-clock exceeds this ratio                   |
| `TURBOGAP_KEEP`      | keep staging dirs for inspection                                 |
| `TURBOGAP_BLESS`     | regenerate `expected/` for golden cases                          |

## Two kinds of case

Most cases compare against the frozen baseline binary, and must be
bit-identical to it.

A case whose `case.conf` sets `REFERENCE=golden` compares against a checked-in
`expected/` directory instead. That is for inputs the baseline cannot run at
all because it crashes on them — there is no baseline output to diff. These are
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

| case                     | ranks | covers                                                                   |
| ------------------------ | ----- | ------------------------------------------------------------------------ |
| `vdw_ts`                 | 1     | TS pairwise vdW, Hirshfeld local properties                              |
| `vdw_mbd`                | 1     | many-body dispersion                                                     |
| `vdw_mbd_cell_mpi2`      | 2     | MBD + the cross-rank Hirshfeld gradient exchange                         |
| `vdw_tsmbd_md`           | 1     | ts+mbd correction state: both the recompute and the reuse branch         |
| `vdw_tsmbd_md_mpi2`      | 2     | the same under MPI                                                       |
| `vdw_tsmbd_predict`      | 1     | ts+mbd outside MD (golden; baseline crashes here)                        |
| `vdw_tsmbd_predict_mpi2` | 2     | the same under MPI                                                       |
| `vdw_tsmbd_mc`           | 1     | ts+mbd under mc, which also never advances md_istep                      |
| `xps_predict`            | 1     | SOAP loop, local properties, 2b/3b/core_pot, XPS spectrum                |
| `xps_predict_mpi4`       | 4     | the energy/force/virial MPI reductions                                   |
| `gcmc_xps`               | 1     | MC insertion/removal: `n_sites` changes, reallocation, neighbor rebuilds |
| `gcmc_xps_mpi2`          | 2     | the same under MPI                                                       |
| `co_predict`             | 1     | 7176-atom single point; timing reference                                 |
| `co_predict_mpi4`        | 4     | the same under MPI; timing reference                                     |
| `co_md`                  | 1     | velocity Verlet, Berendsen thermostat, per-step rebuild decision         |
| `xrd_mad`                | 1     | MAD: pdf, sf and xrd contributions **carrying forces**                   |
| `xrd_mad_mpi2`           | 2     | the same through pack/`mpi_reduce`/unpack                                |
| `xrd_predict`            | 1     | pdf/sf/xrd forward prediction, `exp_forces` off (KNOWN_ISSUES #4)        |

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

Read `KNOWN_ISSUES.md` before trusting a green run.

`ts+mbd` applies its correction to energies but not to forces. That is
deliberate and correct (KNOWN_ISSUES #2), but it does mean a green suite pins
the behaviour rather than endorsing any particular reading of it.

The `xrd_*` cases close what was the largest gap: `pdf`, `sf` and `xrd` had no
coverage at all, which is the reason KNOWN_ISSUES #4 went unnoticed and the
reason the contribution bundling in `docs/refactor_strategy.md` Phase 3 could
not be finished. `nd` is still uncovered — it shares the `calculate_xrd` path
with `xrd`, so it is the cheapest remaining case to add.

No coverage of nested sampling, box scaling or NPT.
