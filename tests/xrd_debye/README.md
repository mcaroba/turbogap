# Debye scattering acceptance test

Does `xrd_debye = .true.` compute the Debye scattering equation, and is its
gradient the gradient of its energy?

```sh
make DEBUG=0 all
tests/xrd_debye/run.sh
tests/xrd_debye/run.sh forward     # just the pattern
tests/xrd_debye/run.sh gradient    # just the forces and virial
```

Non-zero exit means at least one leg failed. Missing test data is a skip, not
a failure.

| variable             | meaning                                      |
| -------------------- | -------------------------------------------- |
| `TURBOGAP_BIN`       | binary under test (default `bin/turbogap`)   |
| `TURBOGAP_DATA_ROOT` | take the systems from here and fetch nothing |
| `TURBOGAP_TESTS_DIR` | where the data repository is cloned          |
| `TURBOGAP_RANKS`     | ranks for the MPI leg (default 2)            |
| `TURBOGAP_PYTHON`    | interpreter for the reference (default `python3`) |
| `TURBOGAP_KEEP`      | keep the staging directory                   |

It needs a `python3` with `numpy`, and `mpirun` for the MPI leg only. On `alt`
the system python3 has no numpy, so run it as

```sh
TURBOGAP_PYTHON=$HOME/.venvs/turbogap-profiling/bin/python tests/xrd_debye/run.sh
```

## Why this is not in tests/regression

Same reason `tests/dipole` is not. That suite compares byte-for-byte against
`baseline/turbogap.e6eb1aa`, a binary frozen before this keyword existed — it
exits with *"I do not recognize the input file keyword xrd_debye"* and cannot
produce a reference at all. Its comparator is plain `diff`, which has no
tolerance, and the question here needs one.

`tests/regression/cases/xrd_debye_mad` and `..._mpi2` do live in that suite, as
golden characterization cases. They pin the output so a later refactor cannot
drift it silently. They say nothing about whether the physics is right. This
test is what says that.

## The two legs

**forward** — `debye_reference.py` recomputes the pattern with an explicit
double loop over pairs in numpy, reusing nothing from the Fortran but the
tabulated form factors, and diffs it against `xrd_prediction.dat`. Run over:
the raw intensity, `q*F(q)`, `F(q)` with the Lorch window, the
Lorentz-polarization factor at `P = 0.8` (neither 0 nor 1, so a dropped
parameter shows up), neutron diffraction, and the same thing under MPI where
the double sum is split across ranks and reduced.

Agreement is checked to `1e-7` absolute, which is the floor set by
`xrd_prediction.dat` being written `F20.8`. Getting a tolerance that tight
requires one deliberate quirk in the reference: it rounds the form-factor
constants through float32, because the Fortran table is declared `real(dp)`
but written with literals like `2.31000` that carry no `d0` suffix, so each is
rounded before it is widened. Without that the two disagree in the 8th
significant digit for that reason alone.

**gradient** — `fd_gradient.py` displaces one atom at a time and compares
`-(E(+h) - E(-h)) / 2h` with the reported force, then applies a symmetric
homogeneous strain and compares against the virial. The experimental term is
isolated by running each geometry twice, once at the deck's `exp_energy_scales`
and once at zero: everything else is identical between the two and cancels, so
the test does not depend on the GAP forces being right and needs no extra
diagnostic printed from the Fortran. It also checks that the forces sum to
zero, which they must, the energy being a function of interatomic distances.

The default tolerance is `1e-4` relative, against an observed deviation around
`3e-6`. The gap is deliberate: the finite difference is limited by
`energy_xrd` being printed to 8 decimals, so at `E ~ 1e3` and `h = 1e-3` a few
parts in `1e6` is the noise floor and a tighter tolerance would make the test
flaky rather than sharper. A real sign or factor error is off by a factor, not
by parts per million.

The gradient leg runs turbogap 62 times and takes a couple of minutes.

## Test data

`$TURBOGAP_DATA_ROOT/xrd_mad`, the same 91-atom melted graphite-O system and
glassy-carbon structure factor the `xrd_mad` regression case uses. Fetched by
`tests/fetch_test_data.sh`.
