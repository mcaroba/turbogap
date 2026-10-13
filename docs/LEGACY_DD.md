# Legacy domain decomposition

`dd_type = legacy` splits the cell into a grid of spatial domains for the
**evaluation** of energies, forces and the virial. It is the algorithm of the
2024–25 `domain-decomposition` branch (M. Louhivuori, T. Zarrouk), moved behind
the `turbogap_domain` seam and with that branch's defects fixed. It is called
legacy because a different decomposition is being developed separately; the
two are meant to coexist under `dd_type`.

## What it does

Rank 0 stays the owner of the structure, as in the replicated default. At
every neighbour-list rebuild, the decomposition:

1. bins the atoms on rank 0 into a `dd_grid` of domains along the three cell
   vectors, by fractional coordinate, with the cuts placed so every slab holds
   the same number of atoms (`dd_balance = count`, the default) or the same
   volume (`volume`);
2. sends each domain's atoms to that domain's root rank;
3. exchanges halos between neighbouring domains: the atoms within `rcut_max`
   of each face. This is the legacy branch's six-sweep exchange. Atoms received
   along one axis are forwarded along the later ones, which covers edges and
   corners;
4. has each domain's ranks build their neighbour lists and evaluate the owned
   atoms, split between them as `i_beg:i_end`;
5. sends the forces and energies that landed on halo atoms back to their owners,
   in reverse sweep order;
6. gathers the owned results to rank 0 in the original atom order.

Steps 4–6 happen every step. Between rebuilds, which rank 0's skin test
decides exactly as in the replicated mode, steps 1–3 are replaced by a
refresh: each domain keeps its atoms and receives their new positions, and the
neighbour lists are only re-measured. The halo is `rcut_max` wide, which
includes `neighbors_buffer`, so it stays complete for as long as the lists do.

From then on rank 0 has whole arrays again, and MD, the optimisers and the
output run unchanged.

`dd_grid` left at its default `0 0 0` gives one domain per rank, laid out at
the first build as the factorisation with the smallest halo for the cell. An
axis that is not cut has no halo, so 8 domains of a cubic cell become `2 4 1`,
not `2 2 2`.

Ranks are grouped into domains by `dd_grid_affinity`: `block` puts consecutive
ranks in one domain, `cyclic` deals them out round-robin. With more ranks than
domains, a domain's ranks share its atoms the way the replicated mode shares all
of them.

## Why it is shaped like this

The seam's contract is that rank 0 holds the structure and the domain layer
decides what each rank sees and how results come back. Keeping to that contract
means nothing outside `turbogap_domain` learns about domains:

- `domain_sync_state` does nothing; the structure is distributed in
  `domain_build`;
- `domain_build` stashes rank 0's global structure and puts the domain's local
  one (owned atoms first, then halo) into `state`, for `evaluate` to read;
- `domain_complete_contributions` sends the halo forces home, gathers the owned
  results, resizes rank 0's result arrays to the whole system, and restores the
  global structure.

The legacy branch also integrated MD on each domain and migrated atoms across
faces (`migrate` in its `decompose.f90`). That is not ported; the refresh and
the rebalancing rebuilds take its place. It means every MD feature (thermostats, GLE,
barostats, electronic stopping, the optimisers, output) has to learn about
domains, and that is the scope of the other decomposition. Here the
per-step cost is one scatter and one gather of O(N) on rank 0. The replicated
mode pays a broadcast of O(N) to every rank and a reduction of
O(N × families) every step, so this is never worse on communication.

## What runs how

Nothing is refused. Each feature either runs on the domains, runs on a
whole-structure view while the potential stays on the domains, or makes the
run evaluate replicated, with a line in the output saying which and why.

| feature | how |
| --- | --- |
| SOAP, `distance_2b`, `angle_3b`, `core_pot`, local properties | domains |
| TS van der Waals | domains; Hirshfeld volumes are copied from owners to halo copies (`legacy_forward`) |
| electrostatics | domains; charges copied to halos, and the neutralising mean charge and `phi_sum` taken over all atoms |
| XPS | domains for the forces; the spectrum from every atom's binding energy, gathered |
| dipoles | domains; the total dipole summed over owned atoms |
| MAD-IR (ACF, aux, XL) | domains; dmu/dr and the XL force are returned from halos and gathered into global order, the XL weights reordered into local order, and after the evaluation every rank is given the whole structure and results (`replicate_after`), as the IR code expects |
| pair distribution, structure factor, XRD, ND, Debye XRD | the whole-structure view: the structure on every rank with a replicated neighbour list, as the replicated mode does. They are global by nature, so only the potential is decomposed |
| MC, nested sampling, i-PI | domains; the moves and the exchange run on rank 0, which keeps the whole structure. Nested sampling uses `replicate_after` |
| `write_soap`, `write_derivatives` | domains; serial-only in every mode, as before |
| MBD, `ts+mbd` | replicated: its clusters reach beyond `rcut_max` and need polarisabilities and volume derivatives that other domains own |
| a cell that needs a supercell, a grid too fine for the cell | replicated, from the step it is found on |
| a build without MPI, a run without `do_prediction` | replicated |

A domain may be empty, owned and halo alike; it waits.

## Defects fixed from the legacy branch

1. Only the SOAP forces on halo atoms were sent back; the `distance_2b`,
   `angle_3b`, `core_pot` and vdW forces on them were dropped.
2. The reduction inside a domain covered owned atoms only, so with more than one
   rank per domain the halo forces of every rank but the root were lost.
3. The reverse exchange ran in forward sweep order. A contribution on a corner
   atom reaches the intermediate domain after that domain has already sent its
   share back, so it was lost.
4. When the halo arrays grew mid-exchange, the new rows of the send mask were
   left uninitialised.
5. `prepare_grid_distribute` declared `integer :: j = 1`, which is implicitly
   `save`, so a second distribution counted from a stale index. It also
   misassigned counts when a domain was empty.
6. Duplicate removal and the force lookup were linear searches, quadratic in the
   number of atoms; both are now an O(1) id-to-slot map.
7. Errors used a bare `stop`, which exits with status 0; they now go through
   `turbogap_abort`.
8. `dd_grid` had no default.
9. Atoms outside the cell were binned by Cartesian projection against borders
   that only cover the cell, so their domain was undefined. They are now binned
   by wrapped fractional coordinate.
10. The barostat rescaled the domain borders by mixing vector lengths with
    projected widths, which is wrong for a triclinic cell. The borders are now
    recomputed from the cell every step.
11. vdW, MC and the experimental terms ran and gave wrong numbers; they are now
    refused.
12. `migrate` (not ported) could grow its arrays to less than it needed, and its
    debug check had `intent(out)` on the mask it was checking.

## Testing

`tests/legacy_dd/run.sh` runs each input decomposed and replicated and
compares `trajectory_out.xyz` number by number (rtol = atol = 1e-7, since
summation order differs). A leg passes only if its run stayed decomposed to
the end. The 22 legs cover:
- grids of 2, 3, 4, 6 and 8 domains, more than one rank per domain, both
  affinities, a sheared cell, the automatic grid, count and volume cuts, a
  slab with vacuum and two empty domains;
- MD that reuses the domains between rebuilds, and that rebuilds every step;
- TS vdW, electrostatics, both XRD routes on the whole-structure view,
  dipoles, and 40 MD steps of the MAD-IR bias;
- a grid too fine for the cell, which must fall back and still match.

It passes with `TURBOGAP_BIN=bin-gpu/turbogap` at the same tolerance. Three
mutations are caught by it: the reverse exchange in forward order, the
Hirshfeld volumes not forwarded to halos, and the mean charge summed over
local sites.

The regression suite also runs decomposed: `TURBOGAP_INPUT_EXTRA="dd_type =
legacy"` appends that line to the input under test only, and
`TURBOGAP_MIN_RANKS` raises every case's rank count, so `run.sh --gpu`, with a
second copy of the host binary as the reference, compares decomposed against
replicated within each case's tolerance. At 4 ranks, 43 of 44 pass; relax_gd
amplifies the reordered sums to 3e-5, inside the 1e-4 its device tolerance and
its own 1e-9 A position-shift control already need. Many of its cells are
too small to split and fall back, which is why the legs above use larger
ones. With the default `dd_type = none` the suite is bit-identical to
`origin/master`.

## Profiling

The run ends with a per-rank table: owned and local sites, this rank's slice
and pairs, and its time in GAP, the neighbour build and each decomposition
phase (`scatter`, `halo`, `share` inside the neighbour bucket; `reverse`,
`gather` inside E & F). A rank that waits for slower neighbours shows the wait
in `reverse`, so that column is the load-imbalance meter. A `PROFILE=1` build
also marks each phase as an NVTX range (`dd_scatter`, `dd_halo`, ...).

## Measured on one workstation

On alt (Core Ultra 5 225: 6 performance and 4 efficiency cores; RTX 5070), CO
GAP, two MD evaluations per run, one domain per rank. Step-0 energies agree
across every run to 5e-8 eV.

Host build, total seconds (replicated / legacy):

| atoms | 1 rank | 2 | 4 | 8 | 10 |
| --- | --- | --- | --- | --- | --- |
| 57k | 71.2 / 72.4 | 38.7 / 39.4 | 20.6 / 21.1 | 12.8 / 12.7 | 11.7 / 11.9 |
| 194k | 237 / 241 | 129 / 130 | 67.5 / 69.2 | 42.2 / 41.2 | 40.5 / 39.4 |
| 308k | — | 201 / 205 | 108 / 108 | 68.2 / 65.5 | 65.5 / 63.0 |

Device build, 194k atoms, median of three (ranks share the one GPU): 65.2 /
65.0 s on 1 rank, 56.9 / 58.5 on 2, 36.1 / 38.1 on 4.

What the numbers say:

- The two modes scale the same on one node, and that is expected. The
  replicated mode already splits the evaluation over ranks; a decomposition
  changes what each rank holds and sends, not how much it computes.
- With 8–10 ranks on the larger cells legacy is 3–4% faster. Its GAP time is
  lower (53.6 against 60.8 s at 308k, 8 ranks) because a domain's atoms are
  close together in memory.
- Aggregate memory is 10–19% lower at 8–10 ranks (6.0 against 7.4 GB at 308k).
  Per rank the saving is small: the SOAP descriptors dominate, and both modes
  hold them only for their own slice.
- The decomposition's own phases are small: scatter, halo, reverse and gather
  take under 0.2 s per step at 308k atoms, not counting the time a rank waits
  in `reverse` for a slower neighbour.
- Efficiency falls to 0.6 at 10 ranks in both modes because four of the cores
  are efficiency cores. Their ranks take 8% longer in GAP and twice as long in
  the neighbour build, and the other six wait for them, in `reverse` here and
  in the reduction in the replicated mode.
- On the device, kernel and copy time are identical in the two modes
  (`nsys`), and the GPU is busy less than half the wall time: runs are bound by
  host-side work. Legacy's deficit there is all neighbour building.

## Improvements made after that study

1. **Neighbour lists reused across steps** (follows rank 0's skin test). The
   device deficit is gone: 194k atoms, five MD steps, medians, replicated
   against legacy, 99.5 / 98.9 s on 1 rank, 81.7 / 81.4 on 2, 76.9 / 78.2 on 4.
   On the host, 308k atoms and ten MD steps, legacy is now 7.5% faster:
   356.4 / 329.5 s on 8 ranks, 347.7 / 321.7 on 10, with half the neighbour
   time (8.3 / 4.8 s and 7.1 / 3.9).
2. **Cuts balanced by atom count** (`dd_balance = count`, the default). On
   the 57k cell with half of it vacuum, 8 ranks along c, ten MD steps:
   replicated 63.6 s, count 60.7, automatic grid 61.5. Equal volume crashed
   there, on its four empty domains; that was a latent bug in
   `results_prepare`, fixed, and empty domains now have a test leg.
3. **Grid chosen automatically** (`dd_grid` left at `0 0 0`): the
   factorisation with the smallest estimated halo; an uncut axis has none.
4. **Host work on the device's critical path.** The two fixes from
   `perf/gpu-host-overhead` are cherry-picked here: no pair distribution in a
   GAP-only device run, and no zeroing of an unused host array. Giving the
   idle cores to OpenMP (`OMP_NUM_THREADS` = cores / ranks) takes a 1-rank
   device run from 99.5 to 90.7 s, mostly the neighbour build (8.7 to 1.3 s),
   and matters little at 2 and 4 ranks.
5. **Rank 0's serial share** turned out to be small. The 1.2 s per step
   reported first was trajectory output at the first and last steps; the
   integration itself is about 0.2 µs per atom per step, 65 ms at 308k atoms.

## Not yet done

- MBD on the domains: a halo as wide as its cluster radius, a forward exchange
  of the SCS polarisabilities, and the volume-derivative transfer moved onto
  the domain map.
- `dd_integrate = domain`: integrating on the domains, so that rank 0 no
  longer holds the structure between output steps.
