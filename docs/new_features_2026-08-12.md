# New in TurboGAP — August 2026

Everything here is new or changed behaviour a user can reach from an input
file. Defects fixed alongside it are in `tests/regression/KNOWN_ISSUES.md`
#6–#11; the working notes are in `docs/session_handoff_2026-08-12.md`.

| | |
| --- | --- |
| [Variable-cell relaxation](#variable-cell-relaxation) | `optimize = gd-box` rewritten; new `gd_box_weight` |
| [Maxwell-Boltzmann momenta](#maxwell-boltzmann-momenta) | new `velocity_distribution` |
| [Chemical potential reference](#chemical-potential-reference) | new `mc_mu_reference` |
| [Molecular grand-canonical moves](#molecular-grand-canonical-moves) | new `mc_molecule_files` |
| [Readable input echo](#readable-input-echo) | output format change |
| [New test suites](#new-test-suites) | `fd_gap`, `mc_mpi_neighbors`, `velocity_draw` |

---

## Variable-cell relaxation

`optimize = "gd-box"` and `"gd-box-ortho"` relax the positions and the lattice
**together**, in the preconditioned variables of Gubler, Finkler, Schaefer and
Goedecker, *Efficient variable cell shape geometry optimization* (2023).

```
optimize      = 'gd-box'      ! or 'gd-box-ortho' to keep the cell orthorhombic
gd_box_weight = 2.0           ! how heavily the lattice counts against the atoms
max_opt_step  = 0.1           ! largest atomic displacement in one step, A
e_tol         = 1.d-8
f_tol         = 0.005
p_tol         = 0.5           ! bar
```

The optimizer works in

```
q_i = A_0 A^-1 x_i                        positions in the reference cell
A~  = w sqrt(N) A diag(1/|a0|,1/|b0|,1/|c0|)   lattice on the atoms' scale
```

so that a lattice step moves the structure affinely and leaves the atoms alone,
and a single Barzilai-Borwein step length under an Armijo-Goldstein line search
serves both blocks. `w` is `gd_box_weight`; raise it to make the cell move more
slowly relative to the atoms.

**What changed for you.** The old scheme alternated — positions to convergence,
then the box to convergence, then back. On a 512-atom amorphous carbon cell the
position half never reached its tolerance, so the box half never ran and the
cell never moved at all: `gd-box` produced a trajectory bit-identical to `gd`.
Same system, same 80 steps:

| | old | new |
| --- | --- | --- |
| `gd-box` energy | −4321.801 eV | **−4321.967 eV** |
| `gd-box` pressure | −950 bar | **−23 bar** |
| `gd-box-ortho` energy | −4321.600 eV | **−4321.641 eV** |
| `gd-box-ortho` pressure | −2142 bar | **−32 bar** |

`optimize = "gd"` is untouched and bit-identical to before.

**Two things to know.**

- `max_opt_step_eps` is now read and ignored. There is one step length and it
  comes from `max_opt_step`; use `gd_box_weight` to change the balance.
- `fix_atom` constrains fractional coordinates, not Cartesian ones. Under a
  moving cell "this atom does not move" has no frame-independent meaning, so a
  fixed atom rides the cell. Pinning an atom in the laboratory frame *and*
  relaxing the cell is not expressible.

---

## Maxwell-Boltzmann momenta

```
randomize_velocities  = .true.
velocity_distribution = 'maxwell'    ! or 'uniform' (default)
t_beg = 300
```

- `"maxwell"` draws each velocity component from N(0, kT/m), removes the
  centre-of-mass velocity, and **does not rescale** — the kinetic energy
  fluctuates, as it should.
- `"uniform"` is the historical draw: components uniform on [0,1), then scaled
  so the kinetic energy is exactly (3/2)(N−1)kT. It remains the default so that
  existing trajectories reproduce.

**Use `"maxwell"` with `mc_hamiltonian`.** There the momenta are part of the
state the Metropolis test accepts or rejects, and detailed balance holds only
if they come from the distribution that test assumes. Under `"uniform"` the
walk runs and moves, but what it samples is not the canonical ensemble.

Both draws give the same *mean* kinetic energy, which is why the difference is
easy to miss. `tests/velocity_draw` separates them by the kurtosis of the
components — 0 for a Gaussian, −1.2 for a box:

```
maxwell   T = 624.7 K   excess kurtosis = +0.160
uniform   T = 600.0 K   excess kurtosis = -1.162
```

Note the temperatures too: `"uniform"` hits 600 K to the digit because it is
scaled to; `"maxwell"` does not, and should not.

---

## Chemical potential reference

```
mc_mu_reference = 'e0'       ! or 'absolute' (default)
```

An insertion's energy change carries the e0 of the atom that arrived, so with
`"absolute"` the `mc_mu` you quote is measured from zero. With `"e0"` that
reference is added back into the acceptance ratio and `mc_mu` is measured from
the isolated species instead:

```
insertion:  p ∝ exp(-(dE - mu - e0_s)/kT)
removal:    p ∝ exp(-(dE + mu + e0_s)/kT)
```

For a molecule, `e0_s` is summed over its atoms.

---

## Molecular grand-canonical moves

A chemical potential can now exchange a whole rigid molecule instead of a
single atom.

```
n_mc_mu           = 1
mc_species        = 'CO'                ! a label; the atoms come from the file
mc_molecule_files = 'co_molecule.xyz'   ! 'none' for an entry that is an atom
mc_mu             = -1.0
mc_mu_reference   = 'e0'
mc_min_dist       = 1.2
mc_max_dist       = 3.0
```

`co_molecule.xyz` is a plain xyz file:

```
2
CO molecule; the centre of mass is computed on read
C   0.00000000   0.00000000   0.00000000
O   0.00000000   0.00000000   1.12800000
```

- Every species in the file must appear in `species`, or the run stops.
- **Insertion** places a uniformly random orientation (Shoemake's quaternion
  method — sampling Euler angles would bias orientations) with the centre of
  mass drawn the way an atomic insertion draws its site. `mc_min_dist` and
  `mc_max_dist` are checked for every atom of the molecule against the atoms
  already present, never against each other.
- **Removal** takes all of one molecule's atoms together.
- The acceptance ratio treats the molecule as a point particle of the total
  mass, so `mc_mu` is the chemical potential **of the molecule**. The
  rotational partition function is a temperature-dependent constant absorbed
  into it.
- `mc.log`'s `N_mc_species` column counts molecules, not atoms.
- An atom belonging to a molecule is not a candidate for an *atomic* removal
  under some other chemical potential.

**Do not combine a molecular chemical potential with `swap`.** A swap changes
one atom's species, and a `move` displaces one atom; either applied inside a
molecule leaves something that is no longer what `mc_mu` prices. Nothing stops
you, and nothing warns.

`wrap_pbc` will split a molecule across a cell face in the written xyz. That is
cosmetic — every consumer uses the minimum image, and
`tests/mc_mpi_neighbors` confirms the forces are unaffected.

---

## Readable input echo

The "Checking input file..." block now prints one `name = value` line per
keyword, aligned, with no `|` end cap and nothing truncated:

```
 atoms_file                       = atoms.xyz
 species                          = C O
 masses                           = 12.01000000 15.99000000
 mc_types                         = insertion removal move
 mc_acceptance                    = 0.5000000000 0.2500000000 0.2500000000
 mc_mu                            = -1.000000000   ! eV
 mc_molecule_files                = co_molecule.xyz
```

The point is that you can paste any of it back into an input file. Three things
changed to make that true:

- **Nothing is truncated.** Names longer than 20 characters used to lose their
  tails (`write_pair_distribut`, `mc_planes_restrict_t`), and so did values
  longer than the 42-character line — silently, including long file paths.
- **List keywords print on one line.** `species = C O` rather than a header
  followed by one value per line, which read back as nothing.
- **Units are comments.** `mc_mu = -1.0   ! eV`, which the input parser ignores.

Roughly thirty keywords that were parsed but never echoed now are, most of them
list-valued: `e0`, `mc_types`, `mc_acceptance`, `mc_mu`, `mc_species`,
`mc_swaps`, `mc_relax_after`, `mc_planes`, `exp_labels`, `exp_data_files`,
`exp_energy_scales`, `vdw_c6_ref` and the rest. The `.gap` potential-file
parser still does not echo its keywords; that is a separate, much larger,
stream.

---

## New test suites

Three suites outside `tests/regression`, because each checks something a diff
against stored output cannot.

```sh
tests/fd_gap/run.sh              # forces and virial vs finite differences
tests/mc_mpi_neighbors/run.sh    # neighbour lists along the MC path
tests/velocity_draw/run.sh       # the shape of the velocity distribution
```

**`fd_gap`** checks each GAP contribution's forces (F = −dE/dr) and virial
(W = −dE/dε) against a finite difference of its own energy, one leg per
contribution — `soap`, `2b`, `3b`, `core_pot`. They are four separate
derivative implementations summed into the same arrays, so an error in one is
invisible in the total wherever another dominates: on this cell the 3b energy is
0.55 eV against SOAP's 4115, and the core potential is identically zero. Each
leg builds a potential file holding only that contribution's blocks, which makes
the total energy, forces and virial that contribution's own.

```
soap       max rel dev  3.6e-06 (forces)  7.6e-05 (virial)
2b                      6.3e-06           8.3e-08
3b                      1.8e-04           2.7e-06
core_pot                7.1e-06           3.9e-06
```

Two things came out of building it.

*A potential file with no `soap_turbo` block segfaulted before the first
energy*, because five sites read `soap_turbo_hypers(:)` and the array is only
allocated when there is at least one block. A pure 2b, 3b or core_pot potential
now runs — which is what made the per-contribution legs possible in the first
place.

*One atom has a kinked SOAP energy.* Atom 326 of `xps_opt/atoms.xyz` has an
x-force 0.5% away from the central difference, while 24 components on 8 other
atoms agree to 5e-6. The energy is smooth either side of that point and the
analytic force sits between the two one-sided slopes, so it is a discontinuity
rather than a wrong force. It is open, written up as `KNOWN_ISSUES.md` #12, and
excluded by name from the `soap` leg so that leg stays a working gate.

**`mc_mpi_neighbors`** checks that the neighbour lists a Monte Carlo walk
carries through insertions and removals — atomic and molecular, at 1, 2 and 4
ranks — are the lists a cold start would build. A wrong neighbour list is
perfectly reproducible, so a stored reference would defend it rather than catch
it; the test runs the walk, feeds `mc_all.xyz` back through `turbogap predict`,
and compares per-atom forces and local energies frame by frame.

**`velocity_draw`** checks the shape of the velocity distribution, as above.

All three take `TURBOGAP_BIN`, `TURBOGAP_DATA_ROOT` and `TURBOGAP_KEEP` like
the rest of `tests/`. `fd_gap` needs numpy; on `alt` the system python has
none, so pass `TURBOGAP_PYTHON=$HOME/.venvs/turbogap-profiling/bin/python`.

`tests/regression` gained twelve cases and a `TURBOGAP_CASE_TIMEOUT`, so a
defect that hangs a run fails that case instead of stalling the whole suite
with nothing on the terminal to say which one is stuck.

---

## Committing

The tree is a fixed point of `.pre-commit-config.yaml`: `pre-commit run` makes
no changes, so a commit passes first time rather than needing a second one to
pick up what the hooks rewrote. Two things had to be settled to get there.

`tools/keyword_docs.py` wrote `docs/keywords.md` with a blank last line and
pre-commit's `end-of-file-fixer` stripped it, so the two rewrote each other on
every commit, for ever. The generator now emits exactly one trailing newline —
the same class of conflict the config already warns about for fprettify against
`trailing-whitespace`, one hook over.

`src/read_files.f90` was not a fprettify fixed point: the formatter wanted the
`!> @kw` doc blocks between `else if` branches indented to the body level, and
there are ~340 of them. That normalisation is in, which is why the diff carries
about 2600 lines of pure re-indentation in that file. Nothing else moved: the
suite is unchanged by it.
