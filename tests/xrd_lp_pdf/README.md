# The Lorentz-polarization factor on the pdf route

Does `xrd_lorentz_polarization = .true.` weight the pattern the pdf → partial
structure factor → `I(q)` chain produces, and does the MAD gradient carry the
same weight?

```sh
make DEBUG=0 all
TURBOGAP_DATA_ROOT=$HOME/work/cpu_vs_gpu_tests/input tests/xrd_lp_pdf/run.sh
tests/xrd_lp_pdf/run.sh forward     # just the pattern
tests/xrd_lp_pdf/run.sh gradient    # just the forces
```

Non-zero exit means at least one leg failed. Missing test data is a skip, not
a failure.

| variable             | meaning                                      |
| -------------------- | -------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)   |
| `TURBOGAP_DATA_ROOT` | directory holding the test systems           |
| `TURBOGAP_RANKS`     | ranks for the MPI leg (default 2)            |
| `TURBOGAP_PYTHON`    | interpreter for the reference (default `python3`) |
| `TURBOGAP_KEEP`      | keep the staging directory                   |

On `alt` the system python3 has no numpy:

```sh
TURBOGAP_PYTHON=$HOME/.venvs/turbogap-profiling/bin/python tests/xrd_lp_pdf/run.sh
```

## What this tests that `tests/xrd_debye` does not

`tests/xrd_debye` checks the same factor on the *other* route, against a
pattern rebuilt from scratch in numpy. The two routes apply it in different
places: the Debye one inside `calculate_xrd_debye`, this one in
`calculate_xrd`, after the MPI reduce and on top of whatever `xrd_output`
asked for. The pdf route also has a separate gradient — `dy/dr` comes out of
`get_structure_factor_forces{,_matrix}` rather than out of the Debye
derivative — and that is the part a wrong sign or a missing factor would
break silently, because the pattern would still look right.

## The two legs

**forward** — the same deck run twice, with the factor off and on, and
`lp_reference.py` asserting `y_lp == y_raw * LP(theta)` with `LP` recomputed
from the q column of the file and nothing else. Differencing two runs is what
makes this independent of the pdf chain: whatever that chain produces cancels,
so the leg says nothing about it and everything about the weight applied on
top. Run over the raw intensity, `q*F(q)`, neutron diffraction, and 2 ranks —
the MPI leg is the one that would catch a factor applied per rank rather than
once, since each rank owns only its slice of the q grid.

Agreement is checked to `1e-6` relative to the peak of the weighted pattern,
against an observed `6e-8`. The floor is the q column being written `F20.8`:
`LP` goes as `1/sin^2(theta)`, so at `LP ~ 100` a 5e-9 error in `q` moves it
by more than the 5e-9 in `y` does.

**gradient** — `fd_gradient.py`, shared with `tests/xrd_debye`, on an
LP-weighted MAD deck: displace one atom at a time and compare
`-(E(+h) - E(-h))/2h` against the reported force. Observed worst deviation
`2e-7` relative, tolerance `1e-4`.

`exp_energy_scales` is 1.0 here where the `xrd_mad` case uses 100. `LP` peaks
around 117 over this q range and the residual enters the energy squared, so at
scale 100 the virial overflows the field width of the xyz writer and comes
back with two numbers run together. The gradient check does not depend on the
scale — it divides out.

## The virial is not checked, and why

`--skip-virial`. The pdf route's virial does not agree with a strain
derivative of its own energy, and that has nothing to do with this factor:

| deck                       | analytic `W_00` | finite difference |
| -------------------------- | --------------- | ----------------- |
| pdf route, LP off          | 34.46           | 16.75             |
| pdf route, LP on           | 28475.12        | 3621.60           |

The LP-off row is reproduced bit-for-bit by `bin/turbogap` built from the
commit *before* the factor was carried to this route, so it predates the
change. The Debye route passes the same check to `3e-6`, which is what the
difference between the two points at: `get_structure_factor_forces{,_matrix}`
accumulates `sum_pairs F (x) r` over both orderings of each pair, and the pdf
it differentiates is normalised by the number density `N/V`, so the energy
depends on the cell explicitly and not only through the interatomic distances
that expression telescopes into. Neither term is the Lorentz-polarization
factor's to fix.

Consequence: MAD forces on the pdf route are right, MAD *stress* on it is not,
with or without this factor. `tests/regression/cases/xrd_mad` pins the number
that comes out; it does not check it against anything.

## Test data

`$TURBOGAP_DATA_ROOT/xrd_mad`, the 91-atom melted graphite-O system and the
glassy-carbon structure factor the `xrd_mad` regression case uses. Fetched by
`tests/fetch_test_data.sh`.
