# Neighbour lists along the Monte Carlo path

Does a Monte Carlo walk that inserts, removes, swaps and displaces atoms end
each step with the neighbour lists a cold start would have built — on every
rank?

```sh
make
tests/mc_mpi_neighbors/run.sh
tests/mc_mpi_neighbors/run.sh r1        # or any subset of the legs
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
| `TURBOGAP_PYTHON`    | interpreter for the comparison (default `python3`) |
| `TURBOGAP_KEEP`      | keep the staging directory                        |

The comparison needs no numpy.

## Why this is not a regression case

`tests/regression` compares output against stored output, so it catches a
change. A wrong neighbour list is not a change: it is perfectly reproducible,
and it would have been baked into the stored reference on the day it was
recorded. The suite would then defend it.

What this test compares instead is one computation against a *different*
computation of the same thing.

1. `turbogap mc` runs the walk on N ranks with `mc_write_xyz`, so every
   accepted configuration reaches `mc_all.xyz` carrying the forces and local
   energies the MC path computed for it.
2. `turbogap predict` reads `mc_all.xyz` back on one rank. Prediction handles
   concatenated frames and builds every neighbour list from scratch for each,
   which is the reference.
3. `compare_frames.py` compares per-atom forces and local energies frame by
   frame.

Round-off agreement means the lists the walk carried were the lists a cold
start would have built. Disagreement localises to the frame and the atom, which
is usually enough to say which move type broke it.

## The legs

| leg | ranks | what an isolated failure means |
| --- | ----- | ------------------------------- |
| `r1` | 1 | the rebuild itself is wrong; not a communication problem |
| `r2` | 2 | `r1` passing and this failing puts it in the per-step broadcast or the slice split |
| `r4` | 4 | a split fine at 2 ranks and wrong at 4 — an empty or boundary slice |

## Tolerances

`1e-5` eV/Å on force components, `1e-6` eV on local energies. One neighbour
missing from one atom moves that atom's force by O(0.1) eV/Å, five orders of
magnitude above this, so the tolerance is not what decides the verdict. It is
there to absorb two things that are legitimate: MPI reductions are not
associative, so the rank count changes the last digit or two, and the
coordinates make the round trip through `mc_all.xyz` at `F16.8`.

## What it does not cover

- The walk's *statistics*. This says each configuration's energy was computed
  correctly, not that the configurations were sampled from the right
  distribution.
- `mc_relax` and the `md` move. Both are integrator bursts, and they used to
  leave the stored configuration one step ahead of the energy stored with it
  (`KNOWN_ISSUES.md` #11) — which this test is what found. They are left out of
  the deck so that a future regression there fails *that* issue's cases rather
  than confusing this one; add a leg if you want them here too.
- Supercells. `indices > 1` forces a rebuild every step through a different
  branch.
