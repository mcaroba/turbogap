# GAP forces and virial against finite differences

Does each contribution's force equal −dE/dr, and its virial −dE/dε, for its own
energy?

```sh
make
tests/fd_gap/run.sh
tests/fd_gap/run.sh 3b core_pot      # or any subset of the legs
```

The test systems come from the `turbogap_tests` repository beside this one;
`tests/fetch_test_data.sh` clones it, and `run.sh` calls that itself if it is
not there yet. Non-zero exit means at least one leg failed. Missing test data
is a skip, not a failure.

| variable             | meaning                                           |
| -------------------- | ------------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)        |
| `TURBOGAP_DATA_ROOT` | take the systems from here and fetch nothing      |
| `TURBOGAP_TESTS_DIR` | where the data repository is cloned               |
| `TURBOGAP_PYTHON`    | interpreter for the reference (default `python3`) |
| `TURBOGAP_KEEP`      | keep the staging directory                        |

The reference needs numpy. On `alt` the system python3 has none:

```sh
TURBOGAP_PYTHON=$HOME/.venvs/turbogap-profiling/bin/python tests/fd_gap/run.sh
```

## One leg per contribution, and why

`soap_turbo`, `distance_2b`, `angle_3b` and `core_pot` are four separate
derivative implementations that happen to be summed into the same `forces` and
`virial` arrays. Every regression case in the tree checks the total, and a
factor or a sign wrong in one of them is invisible in the total on any system
where another dominates — on this cell the 3b energy is 0.55 eV against SOAP's
4115, and the core potential is identically zero.

Each leg therefore builds a potential file holding **only** that contribution's
blocks (`reduce_gap.py`). The run's total energy, forces and virial are then
that contribution's own, which is what lets the check be a plain
finite-difference of the total with nothing extra printed out of the Fortran.

`fd_gradient.py` does the differencing, in `--absolute` mode. For the
experimental families it isolates a term by running twice at different energy
scales and subtracting; a GAP family has no such knob, so isolation is by
reduced potential instead.

Forces and virial are independent claims and both are checked. Forces alone
would pass on a potential whose energy depended on the cell in a way the virial
did not capture, which is exactly the defect `tests/pdf_virial` exists for on a
different route.

## The legs

| leg | h | tol | notes |
| --- | --- | --- | --- |
| `soap` | 1e-3 | 1e-4 | excludes atom 326, see below |
| `2b` | 1e-3 | 1e-4 | |
| `3b` | 1e-2 | 1e-3 | forces are O(0.02) eV/Å here, so the difference needs a bigger step to clear the print resolution |
| `core_pot` | 2e-4 | 1e-4 | cell compressed to 0.75, and `--pick largest` |

The numbers are not uniform, and the reasons are worth knowing before changing
them.

**The floor.** The energy is printed to 8 decimals, so a central difference over
`h` is good to about 1e-8/h eV/Å. Divided by the size of the forces that gives
the relative floor a tolerance has to clear: ~1e-6 for `soap` and `2b`, but 5e-4
for `3b` at h = 1e-3 — above the tolerance, which is why that leg uses a larger
step. Its energy surface is shallow enough that the extra truncation error costs
less than the signal gained.

**The ceiling.** Truncation error grows as h². `core_pot` is a spline through a
steep repulsive wall, and at h = 1e-3 the central difference's own error is
2e-4 relative — the tolerance, not the code. h = 2e-4 puts it at 7e-6.

**Sparse contributions need `--pick largest`.** Most atoms feel no core
potential even in the compressed cell, so a randomly chosen atom has an analytic
force of exactly zero and a finite difference of exactly zero: the leg passes
having checked nothing. Picking the largest forces is the only useful choice
there, and it is also how the leg shows action-reaction pairs — atoms 9 and 24
come out with exactly opposite forces.

**`core_pot` needs a compressed cell.** The core potential is the short-range
repulsion the GAP is not fitted to describe, so on an equilibrated structure it
is identically zero. `scale_xyz.py` applies a uniform strain to cell and
positions until it switches on; 0.75 gives 0.62 eV and forces up to 19 eV/Å.

## The excluded atom

The `soap` leg passes `--exclude-atoms 326`. That atom's x-component force
disagrees with a central difference by 0.5%, and the reason is a kink in the
energy at that geometry rather than a wrong force: the energy is smooth either
side, the analytic force sits between the two one-sided slopes, and 24
components on 8 other atoms agree to 5e-6.

It is an open question with a written-up reproduction — `KNOWN_ISSUES.md` #12 —
not a tolerance problem. It is excluded by name so this leg stays a working gate
on the other 511 atoms instead of a permanent red. Remove the exclusion when
picking it up.

## What it does not cover

- `vdw`, and the local-property models. Both have their own derivative code.
- The experimental families — `tests/xrd_debye` and `tests/pdf_virial` cover
  those, with the same tool in its differencing mode.
- Any potential but `xps_opt`'s CO GAP. A contribution's derivative can be right
  for one set of hyperparameters and wrong for another; `compress_soap`,
  `radial_enhancement` and the non-polynomial `scaling_mode` are all untested
  here.
