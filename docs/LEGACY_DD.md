# Legacy domain decomposition

`dd_type = legacy` splits the cell into a grid of spatial domains for the
**evaluation** of energies, forces and the virial. It is the algorithm of the
2024–25 `domain-decomposition` branch (M. Louhivuori, T. Zarrouk), moved behind
the `turbogap_domain` seam and with that branch's defects fixed. It is called
legacy because a different decomposition is being developed separately; the
two are meant to coexist under `dd_type`.

## What it does

Rank 0 stays the owner of the structure, as in the replicated default. Every
step, the decomposition:

1. bins the atoms on rank 0 into a `dd_grid` of domains along the three cell
   vectors, by fractional coordinate;
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

From then on rank 0 has whole arrays again, and MD, the optimisers and the
output run unchanged.

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
faces. That is not ported. It means every MD feature (thermostats, GLE,
barostats, electronic stopping, the optimisers, output) has to learn about
domains, and that is the scope of the other decomposition. Here the
per-step cost is one scatter and one gather of O(N) on rank 0. The replicated
mode pays a broadcast of O(N) to every rank and a reduction of
O(N × families) every step, so this is never worse on communication.

## Supported, and refused

Supported: SOAP GAP with `distance_2b`, `angle_3b` and `core_pot`;
energies, forces and the virial; local properties as output; static
prediction, multi-frame prediction, MD and the optimisers. It needs an MPI
build, host or device.

Each of these is refused at start-up with a message, not computed wrongly:

| refused | why |
| --- | --- |
| van der Waals (TS, MBD) | the energy of a halo atom's pair needs that atom's Hirshfeld volume, which comes from its own neighbourhood; that needs a second halo layer or a forward exchange |
| electrostatics, `exp_*`, XPS, pair distribution / structure factor / XRD / ND | global sums, not local ones |
| dipoles, MAD-IR | the dipole is reduced in `evaluate` before the gather |
| MC, nested sampling | `n_sites` changes on rank 0; this probably works but is untested |
| i-PI | untested |
| `write_soap`, `write_derivatives` | per-atom files written by rank 0 from its local slice |
| a cell that needs a supercell (`indices > 1`) | halos would need images |
| a domain narrower than `rcut_max` | a halo would need atoms from beyond the adjacent domain |

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

`tests/legacy_dd/run.sh` compares each decomposition with the replicated run of
the same input, per atom: energies, forces and the virial for a single point,
and the trajectory of a short MD run. It covers grids of 2, 3, 4 and 8
domains, more than one rank per domain, both affinities and a triclinic cell.
Summation order differs, so the check is a relative tolerance, not equality.
Run it with `TURBOGAP_BIN=bin-gpu/turbogap` for the device build; it passes there
at the same tolerance. The regression suite runs in the default
`dd_type = none` and must not move.

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

## What to improve, in order of measured return

1. **Reuse neighbour lists across steps.** Legacy rebuilds every step;
   the replicated mode rebuilds only when an atom has moved half the
   `neighbors_buffer`. That is legacy's whole deficit on the device (+1.6 to
   +3.3 s at 194k) and 2–2.5x the host neighbour time. Keep the domain and halo
   membership between rebuilds, and send only the halo positions each step.
2. **Balance the domains by work, not volume.** Equal-volume slabs are fine
   for this homogeneous cell (owned counts within 1.3%) but cannot absorb
   slower cores, and would be badly unbalanced for a slab with vacuum or a
   cluster. Place the borders at atom-count quantiles, optionally weighted by
   each rank's measured GAP time.
3. **Choose the grid automatically.** The halo costs local/owned = 1.6–2.1 at
   8–10 domains here; flatter grids (`5 2 1`) cost more than near-cubic ones.
4. **Move host work off the device's critical path.** Not specific to the
   decomposition: one GPU shared by 4 ranks reaches 1.6–1.8x in either mode.
   This branch also still computes a pair distribution in a GAP-only device
   run; `perf/gpu-host-overhead` removes it.
5. **Rank 0's serial share.** Integration, binning and unpacking are O(N) on
   rank 0 every step (about 1.2 s per step at 308k atoms), the same in both
   modes. Small at 10 ranks; it caps scaling at a few hundred.

## Not yet done

In rough order of value:

- a forward exchange of per-site properties, which unlocks TS vdW (the two-layer
  halo comes with it);
- dipoles and MAD-IR: sum the owned dipoles over domains before `evaluate`
  reads them;
- MC and nested sampling: test, then allow.
