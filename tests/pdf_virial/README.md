# The MAD virial on the pdf route

Does the virial of a MAD run driven through the pair distribution agree with a
strain derivative of that run's own energy?

This is the GPU branch's copy of the CPU branch's suite. The two differ in two
places: there is no Lorentz-polarization factor here, so the `xrd` leg runs
without it, and the matrix force route accumulates its pair virial in
`kernel_exp_force_virial_collection` (`src/gpu/mad_xrd.cu`) rather than on the
host, which makes that leg a check of the device kernel.

```sh
make
tests/pdf_virial/run.sh
tests/pdf_virial/run.sh xrd sf      # or any subset of the legs
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
| `TURBOGAP_RANKS`     | ranks for the MPI leg (default 2)                 |
| `TURBOGAP_PYTHON`    | interpreter for the reference (default `python3`) |
| `TURBOGAP_KEEP`      | keep the staging directory                        |

On `alt` the system python3 has no numpy:

```sh
TURBOGAP_PYTHON=$HOME/.venvs/turbogap-profiling/bin/python tests/pdf_virial/run.sh
```

## What was wrong

The forces on this route agreed with finite differences of the energy all
along. The virial did not, by a factor of two on the shear components and by
more than that on the diagonal. Two independent causes:

**Every pair was counted twice.** `sum_a F_a (x) r_a` telescopes into one term
per *unordered* pair, `f_ij (x) r_ij`. The neighbour list is a full list, so
`get_structure_factor_forces{,_matrix}` and `get_pair_distribution_forces` walk
each pair twice, once from either end, and both visits contribute the same
`f (x) r` — the force and the separation change sign together. The accumulation
now carries `0.25`: a half to symmetrise and a half for the second visit. The
matrix route's accumulation is the device kernel's `loc_viri`, which is built
over the same `k_list` of ordered pairs and carries the same `0.25`.

**The pattern depends on the cell directly, not only through the distances.**
Each partial is normalised by the number density `N/V`, so a homogeneous strain
moves the prediction through the volume as well as through the interatomic
distances, and `sum_a F_a (x) r_a` cannot see that. Where the volume enters
depends on the observable:

- `S_ab(q) - delta_ab = c_factor * sum_r M(q,r) ( g_ab(r) - 1 )`, with the
  `N/V` in `c_factor` and a compensating `V` inside `g_ab`. The two cancel in
  the `g_ab` term, which is why the missing piece is not the whole thing — only
  the subtracted `-1` survives, going as `1/V`. That is what
  `add_structure_factor_volume_virial` puts back, on the diagonal, for the xrd,
  nd and structure-factor observables.
- `g(r)` itself is proportional to `V` outright, so `V dy/dV` is the pattern.
- `D(r) = 4 pi rho r ( g(r) - 1 )` has the volume cancel in the `g` term again,
  leaving `4 pi rho r`.

Both pdf cases are added by `calculate_pair_distribution` and its GPU twin
`gpu_calculate_pair_distribution`, which hold the assembled pattern; the
structure-factor case is added inside the force routine, which is where the
per-channel weights already live. In the matrix route that means after the
device virial has been copied back, since the two halves land in the same
array.

The volume term belongs to the cell rather than to the pairs a rank owns, so
exactly one rank may add it — hence the `rank == 0` guards, and hence the `mpi`
leg, which is the only thing that would catch every rank adding it.

## The legs

Each runs `fd_gradient.py` over a deck: it applies a homogeneous strain to the
cell and every position, differences the energy, and compares against the
reported virial. Forces are checked in the same run.

| leg | what it reaches |
| --- | --- |
| `xrd` | the batched device route, which is what a default deck runs |
| `xrd-plain` | `get_structure_factor_forces`, the host route, `gpu_batched` off |
| `sf` | both routines with `do_xrd` off |
| `pdf` | `get_pair_distribution_forces`, `g(r)` |
| `pdf-dr` | the same with `pair_distribution_output = 'D(r)'` |
| `mpi` | the virial at 1 rank against `TURBOGAP_RANKS` |

All six run. Four of them could not until `gpu_batched = .false.` was given a
route to take again — see below, and `src/turbogap_exp.f90`.

`get_structure_factor_forces_matrix` has no leg. It reads the device pair lists
that only the batched pdf builds, so it is reachable from neither route: the
batched one calls `get_structure_factor_forces_matrix_original` instead, and the
unbatched one has no device lists to give it. It carries the fix, unexercised.

The `mpi` leg compares with a tolerance rather than byte for byte. The pair half
is a device reduction whose summation order is not fixed, and splitting the
pairs across ranks reorders it again, so the last digit moves between runs.
What the leg is looking for — the volume term added once per rank instead of
once — is a whole term, percent-sized, not 1e-8.

Tolerance is `1e-4` throughout. What each leg actually lands at is
finite-difference resolution, not physics: the frame carries 8 decimals, so a
leg whose energy is large and whose virial is small resolves the derivative
less well, and the device reductions do not fix their summation order.

## Four of these legs used to be impossible, and why they are not

Until `src/turbogap_exp.f90` grew an unbatched branch, `compute_exp_spectra`
assembled the *total* pattern for exactly one configuration: `gpu_batched` on,
with `do_xrd` or `do_nd`. Everything else fell through to the similarity block,
which assigned an unallocated array into `y_pred` and segfaulted. That took out
`sf`, `pdf` and `pdf-dr`, whose labels had no assembled total at all, and
`xrd-plain`, which reaches the `exp_utils` routines by turning `gpu_batched`
off — the very thing that removed the only assembly of `y_xrd`.

The batched route collects the partials and stops; nothing turned them into
`y_pair_distribution`. The unbatched branch calls `calculate_pair_distribution`,
the host implementation, which does. The pattern assembly — the
`calculate_structure_factor` / `calculate_xrd` / `calculate_nd` calls — is now
shared by both routes, and still runs before the batched forces, which build
their residual from `y_xrd`.

`gpu_batched = .false.` forces `structure_factor_matrix_forces = .false.` and
says so. The matrix route reads the device pair lists the batched pdf builds,
and the host pdf builds none.

The `xrd` and `xrd-plain` legs agree to the last printed digit, which is worth
more than either alone: the batched device route and the host route are
independent implementations of the same virial, and neither was checked against
the other before.

`run.sh` still probes for the old failure and skips with the reason rather than
failing, so the suite stays honest on a branch where it comes back.

`--strain 1e-5`, against the harness default of `1e-4`. The pair distribution
is cut hard at `pair_distribution_rcut`, so a pair crossing that radius takes a
finite bite out of the pattern; at `1e-4` a few pairs on this system cross and
the finite difference picks up about a percent of noise.

## One case is still inconsistent, deliberately

Ask for `pair_distribution_output = 'D(r)'` and hand it a target that is *not*
declared `exp_input_type = 'D(r)'`, and `preprocess_exp_data` builds the target
by applying `4 pi rho r ( g - 1 )` to the file — using the density of the step
it is on, every step. The fit target then breathes with the cell, the energy
depends on the volume through `y_exp` as well as through `y`, and the virial
does not carry that term.

The `pdf-dr` leg declares its input, which is the well-posed usage, and passes.
The other case looks more like a question about whether an experimental target
should move at all than one about the virial, so it is left alone rather than
papered over.

## Test data

`$DATA_ROOT/xrd_mad`: the 91-atom melted graphite-O system and the
glassy-carbon structure factor the `xrd_mad` regression case uses. The targets
for the `sf`, `pdf` and `pdf-dr` legs are generated by `make_targets.py` from a
prediction run, since no file in the test set is on those grids.
