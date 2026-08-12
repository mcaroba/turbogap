# Refactor regression suite

## Running

```sh
make DEBUG=0 all                       # or whatever your arch build is
tests/regression/run.sh                # all cases
tests/regression/run.sh vdw_mbd        # one case
tests/regression/run.sh --list         # case names
```

Non-zero exit means at least one case differs from the baseline.

The test systems come from the `turbogap_tests` repository beside this one.
`run.sh` clones it through `tests/fetch_test_data.sh` on the first run, so
there is nothing to set up by hand.

| variable             | meaning                                                          |
| -------------------- | ---------------------------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)                       |
| `TURBOGAP_REF_BIN`   | reference (default `tests/regression/baseline/turbogap.e6eb1aa`) |
| `TURBOGAP_DATA_ROOT` | take the systems from here and fetch nothing                     |
| `TURBOGAP_TESTS_DIR` | where the data repository is cloned                              |
| `TURBOGAP_TIME_TOL`  | fail if test/ref wall-clock exceeds this ratio                   |
| `TURBOGAP_CASE_TIMEOUT` | seconds one run may take before it is killed (default 1800)   |
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
| `xrd_debye_mad`          | 1     | the Debye XRD route carrying forces (golden; baseline rejects keyword)  |
| `xrd_debye_mad_mpi2`     | 2     | the Debye double sum split and reduced across ranks                     |
| `mc_moves_relax`         | 1     | `swap` + `mc_relax`: gradient descent re-entered at a new `n_sites`      |
| `mc_moves_relax_mpi2`    | 2     | the same, with the per-step broadcast resizing behind it                 |
| `mc_hybrid_md`           | 1     | the `md` move under `mc_hamiltonian`: `E_kin` in the Metropolis test     |
| `mc_volume`              | 1     | the `volume` move: `modify_box` and the NPT acceptance ratio             |
| `mc_planes`              | 1     | insertions restricted to a slab by `mc_planes`                           |
| `relax_gd`               | 1     | `optimize = gd`: backtracking line search, Barzilai-Borwein step         |
| `relax_gd_box`           | 1     | variable-cell relaxation (golden; the baseline never moves the cell)     |
| `relax_gd_box_ortho`     | 1     | the same with the cell held orthorhombic                                 |
| `relax_gd_box_sheared`   | 1     | the same from a triclinic cell, for the off-diagonal lattice gradient    |
| `mc_hybrid_md_maxwell`   | 1     | the same walk with Maxwell-Boltzmann momenta rather than the old draw    |
| `mc_reverse_xps`         | 1     | reverse MC: the experimental mismatch inside the Metropolis test         |
| `mc_molecule`            | 1     | grand-canonical exchange of a whole rigid molecule, and `mc_mu_reference` |

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

`xps_opt/atoms_sheared.xyz` is `atoms.xyz` with a 4%/3%/2% affine shear applied
to the cell and the positions together, generated by `tools/shear_xyz.py`. The
fractional coordinates are unchanged, so it is the same structure under a
strain the variable-cell relaxation has to undo.

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

The `mc_*` and `relax_*` cases close the second-largest gap. Before them the
only Monte Carlo coverage was `gcmc_xps`, which uses `move`, `insertion` and
`removal` and nothing else: `swap`, `volume`, `md`, `relax`, `mc_hamiltonian`,
`mc_relax_after` and the `mc_planes` geometry were all uncovered, and so were
both relaxation drivers. That is how KNOWN_ISSUES #6 through #10 survived. In
particular `mc_hybrid_md` is the only case that reads `E_kinetic` in an
acceptance ratio, and `relax_gd_box_sheared` the only one whose cell is not
nearly cubic.

They are also the reason most of the new cases are `REFERENCE=golden`: the
baseline binary does not merely differ on them, it hangs (KNOWN_ISSUES #6).

Two things these cases cannot check, which have suites of their own:

- **`tests/mc_mpi_neighbors`** — that the neighbour lists the walk carries
  through insertions and removals are the ones a cold start would build, on
  every rank. A wrong neighbour list is reproducible, so a stored reference
  would defend it rather than catch it; that test compares the walk against a
  fresh prediction instead. It is what found KNOWN_ISSUES #11.
- **`tests/velocity_draw`** — that `velocity_distribution = "maxwell"` really
  is Maxwell-Boltzmann. Both draws are built to have the right mean kinetic
  energy, so only the *shape* separates them, and a diff against stored output
  says nothing about shape.

Still no coverage of nested sampling, box scaling, the Berendsen barostat, or
`mc_relax_opt = gd-box` (relaxing the cell inside a Monte Carlo walk, which
runs but is not pinned — and where KNOWN_ISSUES #11 is still open for the
cell).
