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
| `xrd-matrix` | `get_structure_factor_forces_matrix`, `gpu_batched` off |
| `xrd-plain` | `get_structure_factor_forces`, the host implementation |
| `sf` | both routines with `do_xrd` off |
| `pdf` | `get_pair_distribution_forces`, `g(r)` |
| `pdf-dr` | the same with `pair_distribution_output = 'D(r)'` |
| `mpi` | the virial at 1 rank against `TURBOGAP_RANKS` |

Only `xrd` and `mpi` actually run here; the rest report why they cannot, below.
`xrd` is not a narrow leg for that — the batched route is the only route a deck
on this branch can take, so it is the whole of what a user sees.

The `mpi` leg compares with a tolerance rather than byte for byte. The pair half
is a device reduction whose summation order is not fixed, and splitting the
pairs across ranks reorders it again, so the last digit moves between runs.
What the leg is looking for — the volume term added once per rank instead of
once — is a whole term, percent-sized, not 1e-8.

Tolerance is `1e-4` throughout. What each leg actually lands at is
finite-difference resolution, not physics: the frame carries 8 decimals, so a
leg whose energy is large and whose virial is small resolves the derivative
less well, and the device reductions do not fix their summation order.

## Four legs cannot run here, and it is not the virial

`compute_exp_spectra` assembles the *total* pattern for exactly one
configuration: `gpu_batched` on, with `do_xrd` or `do_nd`. Everything else
falls through to the similarity block, which assigns an unallocated array into
`y_pred` and segfaults — `src/turbogap_exp.f90:1024`, `:1026` and `:1028` for
`y_pair_distribution`, `y_structure_factor` and `y_xrd` in turn. So:

- `sf`, `pdf` and `pdf-dr` cannot run: those labels have no assembled total at
  all. The batched path collects the partials and stops there, and there is no
  other path.
- `xrd-matrix` and `xrd-plain` cannot run: they exist to reach the routines in
  `exp_utils`, which means turning `gpu_batched` off, which is what removes the
  only thing that assembles `y_xrd`.

All of it predates the virial work — the binary from before dies at the same
lines on the same decks. Adding `do_xrd = .true.` to a pdf deck does not help
either: what is missing is the assembly, not the trigger. Guarding the
dereference would stop the crash without making the observable work, so
nothing here does that.

What it costs: `get_structure_factor_forces`, `get_structure_factor_forces_matrix`
and `get_pair_distribution_forces` in `src/exp_utils.f90`, and the volume term
in `calculate_pair_distribution` and `gpu_calculate_pair_distribution`, are
unreachable on this branch and therefore unexercised. They carry the fix anyway,
as a line-for-line port of the CPU branch's change, where `tests/pdf_virial`
does run those legs and they pass. The legs are kept rather than deleted so
they come back on their own the day those routes are reachable.

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
