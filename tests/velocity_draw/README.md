# The velocity draw

Does `randomize_velocities` draw from the distribution its keyword names?

```sh
make
tests/velocity_draw/run.sh
tests/velocity_draw/run.sh maxwell      # or any subset of the legs
```

The test systems come from the `turbogap_tests` repository beside this one;
`tests/fetch_test_data.sh` clones it, and `run.sh` calls that itself if it is
not there yet. Non-zero exit means at least one leg failed. Missing test data
is a skip, not a failure.

| variable             | meaning                                            |
| -------------------- | -------------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)         |
| `TURBOGAP_DATA_ROOT` | take the systems from here and fetch nothing       |
| `TURBOGAP_TESTS_DIR` | where the data repository is cloned                |
| `TURBOGAP_PYTHON`    | interpreter for the check (default `python3`)      |
| `TURBOGAP_KEEP`      | keep the staging directory                         |

The check needs no numpy.

## Why a temperature check would not do

`velocity_distribution` picks between two draws.

- `"uniform"` — each component uniform on [0,1), the centre-of-mass velocity
  removed, then everything scaled so the kinetic energy is exactly
  (3/2)(N−1)kT. Historical, and the default so that existing trajectories
  reproduce.
- `"maxwell"` — each component from N(0, kT/m), centre-of-mass velocity
  removed, **no rescaling**, so the kinetic energy fluctuates as it should.

Both are constructed to have the right *mean* kinetic energy. A temperature
check therefore passes on either, which is exactly why the uniform draw sat
under `mc_hamiltonian` unnoticed: it looks right by every summary statistic
anyone was printing.

What separates them is the shape. The cheapest statistic that sees it is the
excess kurtosis of the mass-scaled components √mᵢ·vᵢⱼ, which is 0 for a
Gaussian and −1.2 for a box, and which scaling, centre-of-mass removal and the
choice of temperature all leave alone. So the test needs no reference data and
no tolerance tuning:

```
== maxwell ==   T = 624.7 K   excess kurtosis = +0.160
== uniform ==   T = 600.0 K   excess kurtosis = -1.162
```

Note the temperatures as well as the kurtoses. `"uniform"` hits 600 K to the
digit because it is scaled to; `"maxwell"` does not, because it is not, and
should not.

## Why it matters

Under `mc_hamiltonian` the momenta are part of the state the Metropolis test
accepts or rejects. Detailed balance holds only if they are drawn from the
Maxwell-Boltzmann distribution that test assumes. With `"uniform"` the walk
still runs and still moves; what it samples is not the canonical ensemble.

## The legs

| leg | what it asserts |
| --- | --------------- |
| `maxwell` | Gaussian components, and a temperature near but not equal to `t_beg` |
| `uniform` | box-shaped components, and a temperature equal to `t_beg` |
| `reject` | an unknown `velocity_distribution` aborts, *and does so naming that keyword* — a non-zero exit alone would also be produced by falling over for some unrelated reason |
