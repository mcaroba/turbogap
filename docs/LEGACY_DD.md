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
prediction, multi-frame prediction, MD and the optimisers. It needs an MPI host
build.

Each of these is refused at start-up with a message, not computed wrongly:

| refused | why |
| --- | --- |
| van der Waals (TS, MBD) | the energy of a halo atom's pair needs that atom's Hirshfeld volume, which comes from its own neighbourhood; that needs a second halo layer or a forward exchange |
| electrostatics, `exp_*`, XPS, pair distribution / structure factor / XRD / ND | global sums, not local ones |
| dipoles, MAD-IR | the dipole is reduced in `evaluate` before the gather |
| MC, nested sampling | `n_sites` changes on rank 0; this probably works but is untested |
| i-PI | untested |
| `write_soap`, `write_derivatives` | per-atom files written by rank 0 from its local slice |
| the device build | untested |
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
and the thermodynamic log of a short MD run. It covers grids of 2, 4 and 8
domains, more than one rank per domain, both affinities and a triclinic cell.
Summation order differs, so the check is a relative tolerance, not equality.
The regression suite runs in the default `dd_type = none` and must not move.

## Not yet done

In rough order of value:

- a forward exchange of per-site properties, which unlocks TS vdW (the two-layer
  halo comes with it);
- dipoles and MAD-IR: sum the owned dipoles over domains before `evaluate`
  reads them;
- MC and nested sampling: test, then allow;
- the device build;
- a test at scale for the scatter/gather cost against the replicated broadcast.
