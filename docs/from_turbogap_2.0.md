# What to take from turbogap_2.0, and what not to

Written 2026-08-06. Companions: `docs/refactor_strategy.md` (what to do and
why), `docs/session_handoff_2026-08.md` (what was done), `docs/refactor_handoff.md`
(the CPU-side de-monolithing).

`~/software/turbogap_2.0` is a from-scratch modular rewrite. It is
experimental, its numerics are not trusted, and it does not carry electronic
stopping or the MBD rewrite. None of that matters here: what is worth taking is
its **type discipline**, and that transfers without any of its physics.

This document adjudicates each candidate, answers the graph report's
`control_t` / `state_t` betweenness question, and sets out what `simulation_t`
has to look like if overlapping domain decomposition is coming.

---

## 1. Adjudication

| feature | verdict | why |
|---|---|---|
| `printing.f90`, `error.f90`, `read_utils.f90`, parameter echo + validation | **taken** (`0833120`/`d961b11`) | |
| keyword families, one reader per subsystem | **taken** (`650e28d`/`a711685`) | six families now character-identical across the branches |
| `get_time` wrapper | **taken** (`fef536c`/`9a1a2b4`) | |
| **`times_t`** — one type for every wall-clock bucket | **TAKE — done, §4** | perf-neutral, removes 13 arguments, converges signatures, fixes the negative Miscellaneous |
| **`perform_t`** — per-iteration decisions evaluated once | **TAKE — next** | the systematic cure for the defect family that has cost six bugs; half of it already exists as `contrib_on` |
| **pure `decide_*` / `check_exit` functions** | **TAKE**, with `perform_t` | |
| **`check_options_*`** — the third leg of the reader triad | **TAKE** | where the electrostatics segfault and `has_vdw` vs `has_local_properties` would have been caught |
| `initialize_options_*` — mode to defaults | TAKE, cheap | master decides mode consequences at scattered sites |
| `split_t` — `{i_beg, i_end, j_beg, j_end}` | TAKE, cheap | four arguments to one in nearly every signature |
| `energy_t` — the ~11 energy scalars | TAKE, **after the merge** | it changes what `xyz.f90` writes, so it is a result-changing commit (§4.2 rules) |
| `calculation_t` — value type per contribution | **SKIP** | `contribution_ref` (`7390225`) already banked the payoff; converting to a value type is a whole-driver rewrite for no additional merge benefit |
| `state_t` — the full bundle | **DEFER into `simulation_t`** | §3 |
| **`tg_memory` / `tg_alloc`** tracked allocation | **SKIP wholesale**, §2 | |
| neighbour-list overallocation (the `used_dims` idea alone) | TAKE, standalone, later | a real MD/MC win that needs no wrapper type |
| `src/<subsystem>/` directory layout | **SKIP until the merge lands** | churn that breaks the per-file cross-branch diff the whole plan is built on |
| per-subsystem `_types.f90` + options types | TAKE gradually | this is the actual fix for `control_t`'s betweenness, §3 |
| `fypp` code generation | SKIP | a build dependency for one generated file |

### The three that are worth the most, and why they are the same thing

`perform_t`, `check_options_*` and `contribution_ref` all attack one defect
shape: **a condition that exists at more than one site and is free to
disagree with itself.** The session handoff lists six instances, and this
document's own work found a seventh (§4). turbogap_2.0's answer is
structural — a decision is evaluated exactly once, into a named field, and
everything downstream reads the field. That is the single most transferable
idea in that tree and it costs nothing at runtime, because every one of these
decisions is made outside the hot loops.

---

## 2. Why `tg_memory` is not worth it here

It is the most impressive thing in turbogap_2.0 and the wrong thing to adopt.
Four reasons, in descending order of weight:

1. **It collides with `c_loc`.** The GPU branch hands driver arrays to kernels
   by `c_loc(array)` at ~30 sites. Under `tg_memory` the array is a component
   of a wrapper type and its logical extent (`%used_dims`) is not its physical
   extent (`%dims`), so every one of those sites needs re-auditing for
   contiguity, `target` attribution and the right length. That audit lands
   entirely on the branch where **the bit-exact contract does not hold** (§6.5)
   — the one place the work cannot be checked.
2. **The overallocation is the valuable half and it does not need the
   wrapper.** Keeping a buffer larger than its logical extent to avoid
   reallocation churn is a genuine MD/MC win. It can be done to the neighbour
   list alone, with a plain `n_used` integer beside it, in one file.
3. **It is a whole-tree change with no merge payoff.** Every other item in §1
   converges the two branches. This one touches every array on both.
4. **turbogap_2.0 has not finished it either** — the commit history says
   "added allocation now make it pervasive", and the report's 15 near-identical
   `tg_alloc`/`tg_memory` communities are the shape of a rollout in progress.

**Take the accounting without the wrapper.** A `memory_t` counter updated by a
thin `tg_allocate` that still yields a plain `allocatable` gives the peak/total
tracking with zero change to how arrays are indexed or passed. That is the part
worth having, and it is a day's work rather than a month's.

---

## 3. The betweenness question, and the alternative

The graph report flags `control_t` (52 edges, betweenness 0.053) and `state_t`
(42 edges, 0.035) as cross-community bridges and asks whether they should be
split. **They should not, and they are not the same problem.** High betweenness
is not the defect; the *direction* of the edges is.

### 3.1 `control_t` — the betweenness is a real defect

`control_t` reaches into the physics: `calculate_pdf`, `calculate_xrd`,
`calculate_sf`, `calculate_gap_soap`, `check_local_properties`, even
`print_times`. Each of those reads two to five flags out of a ~60-field type.
Three consequences:

* every physics module recompiles when any unrelated option is added;
* a kernel is *able* to re-derive a decision the main loop already made, which
  is exactly the defect family in §1;
* the type cannot be split, because splitting it means touching every kernel —
  which is why the report's suggestion is unactionable as posed.

### 3.2 The alternative: do not split it, stop passing it

Three rules, none of which requires touching the type:

1. **`control_t` is read-only below the loop head.** Declare the dummy
   `intent(in)` everywhere after `check_options_*` has run. turbogap_2.0 does
   this by convention; make the compiler enforce it.
2. **Every per-iteration decision is evaluated once, at the loop head, into
   `perform_t`.** Kernels take `perform_t`, or a single `logical`, never
   `control_t`.
3. **Numeric options a kernel genuinely needs live in that subsystem's own
   options type** — `options_vdw_t`, `options_exp_t`, … declared in the
   subsystem's own module. This is turbogap_2.0's own stated convention
   ("feature options-type lives in `<feature>.f90` only", "no module
   variables"), and its `control_t` betweenness is high precisely because it
   has not finished applying it.

`perform_t` then inherits the betweenness, and **that is fine**, because:

> A hub type that `use`s nothing, holds no behaviour, and is written by exactly
> one procedure is a **data bus**. A hub type that is written from many places
> and read from many places is a **god object**. Betweenness cannot tell them
> apart; the two counts below can.

**The measurement to use instead of betweenness**, for any candidate hub:

| | `control_t` today | `perform_t` |
|---|---|---|
| modules it `use`s | 0 | 0 |
| distinct procedures that **write** it | reader + `check_options` + every kernel taking it non-`intent(in)` | 1 |

Drive the second row to 1 and the god object becomes a bus, without splitting
anything. That is a far cheaper refactor than the one the report implies, and
it is the one to do.

### 3.3 `state_t` — the betweenness is inherent and should be left alone

`state_t` *is* the physical system; every physics module reads positions. That
is not accidental coupling and reducing it is not a goal. Splitting it by
subject — positions here, velocities there, species elsewhere — makes every
call site longer and buys nothing. That is the 145-argument nested-sampling
problem run in reverse, and §Phase 6 has already measured where it ends up.

The right move is the opposite of splitting: make `state_t` **nestable**, so
the widest node in the graph is one type *at every scale* rather than a
different bundle per scope. Which is the next section.

---

## 4. `simulation_t` and overlapping domain decomposition

The motivation is worth stating so the design can be judged against it. The
current parallel model is **replicated data**: every rank holds full positions,
work is split by atom index (`split_tasks` giving `i_beg:i_end`), and results
come back through an `mpi_reduce` over `n_sites`. Memory per rank is O(N) and
communication per step is O(N) *regardless of rank count*, so the model caps
both the system size and the useful number of ranks. Spatial decomposition with
a halo fixes both. That, not tidiness, is what `simulation_t` is for.

### 4.1 The shape

```fortran
!  A domain is a spatial region owned by one group of ranks, plus a halo of
!  width rcut_max + buffer.  It is a simulation for the short-ranged additive
!  terms, and only for those -- see 4.3.
type domain_t
   integer  :: id   = 0
   integer  :: comm = MPI_COMM_NULL      ! the ranks sharing this domain
   integer  :: root = 0                  ! rank within comm that aggregates
   integer  :: n_owned = 0, n_halo = 0
   real(dp) :: lo(3) = 0.d0, hi(3) = 0.d0
   real(dp) :: halo_width = 0.d0         ! rcut_max + buffer
   integer,  allocatable :: global_index(:)  ! owned+halo -> global site id
   logical,  allocatable :: is_owned(:)      ! contributes to the global sum
end type

type simulation_t
   type(domain_t)         :: domain
   type(state_t)          :: state     ! owned + halo sites only
   type(neighbors_t)      :: neighbors
   type(split_t)          :: split     ! intra-domain rank split
   type(contribution_ref) :: contrib(1:N_CONTRIB)
   type(times_t)          :: time
end type
```

Four of those six fields either exist today or are the near-term work in §1.
That is deliberate: **the bundling owed now is the down payment on
`simulation_t`**, and doing it first turns the eventual change from a redesign
into a struct that names types already in the tree.

### 4.2 Rule one: `control_t` / `params` is NOT a field

One immutable copy, broadcast once, passed alongside. N domains each carrying
their own option set is the "one condition, N predicates free to disagree"
antipattern promoted to type scale, and it would stay invisible until two
domains disagreed — on a large run, in a way no small test reproduces. This is
the constraint most likely to be violated by the instinct that everything
belongs in one type, and it is the one that must not be.

### 4.3 Rule two: declare which contributions are domain-local

They are not all local, and conflating them is how a decomposition silently
produces wrong physics rather than a crash:

* **domain-local** — soap, 2b, 3b, core_pot, local properties, TS pairwise vdW.
  For these a domain genuinely is a simulation.
* **global** — pdf / sf / xrd / nd are structure-wide sums over all pairs and
  all q-points; MBD is a global eigenproblem; electrostatics is long-ranged;
  XPS depends on a global spectrum normalisation; nested sampling and MC
  insert/remove are whole-configuration moves that change `n_sites`.

`contrib_on` and the `C_*` indices already enumerate exactly these eleven
families, so the encoding is a companion `logical, parameter ::
C_IS_LOCAL(1:N_CONTRIB)` beside them, and the reduce block — which since
`7390225` already walks the family list **once** — is what routes each family
to an intra-domain or a global reduction. The Phase 3 bundling is what makes
this a table lookup rather than a redesign.

### 4.4 Rule three: the communicators are the cheap part

`MPI_Comm_split(MPI_COMM_WORLD, domain%id, rank, domain%comm)`, then a
communicator over the domain roots, is standard and small. The cost is the halo
exchange and the ownership bookkeeping — `is_owned` masks so a halo atom's
energy is not double-counted. That is the same shape as the `gap_interface`
overrun that counted pairs with `<` and filled them with `<=`, at much larger
scale, and it will not announce itself: a double-counted halo atom is a wrong
number, not a crash.

So, per §1.3 of the session handoff: **write the multi-domain regression case
before the code**, on a *small* cell, with at least one atom whose halo image
lives on another domain. A single system size cannot see this class of bug, and
the CO system is far too large to see it at all.

### 4.5 Rule four: on the GPU branch a batch and a domain are the same thing

`gpu_context.f90` already holds `gpu_streams(:)` and `cublas_handles(:)`, and
the driver builds `i_beg_list` and the per-batch bounds at two sites and
consumes them at three. That is a proto-decomposition: a list of index ranges,
each with its own stream. If `simulation_t` is designed so that a GPU batch is
a `simulation_t` with `comm = MPI_COMM_SELF`, the batching machinery and the
domain machinery become one thing, and the diffraction block's 51-against-37
argument gap — which is *entirely* batch bookkeeping — closes as a side effect.
If it is not designed that way, the GPU branch grows a second, parallel
decomposition that has to be kept in agreement with the first, forever.

**This is the decision to make before any of it is written**, and it is the
only part of `simulation_t` that the merge sequencing has to accommodate.

### 4.6 Sequencing

Do not start `simulation_t` before the merge lands. It changes every argument
list on both branches at once, which is the one thing the whole plan has been
avoiding. Build its fields as free-standing types now — `times_t` (done),
`split_t`, `perform_t`, alongside the existing `contrib` — each verified
bit-exact on master and ported, per §6.5.

---

## 5. What was implemented: `times_t`

`src/timing.f90` now declares one `times_t` carrying the union of both
branches' buckets, plus `time_start`, `time_end` and `sum_times`. The file is
byte-identical on the two trees, as `printing.f90`, `error.f90` and
`read_utils.f90` already are.

| | before | after |
|---|---|---|
| timer declarations in `turbogap.f90` | 22 (CPU) / 29 (GPU) | 1 |
| timer arguments across the extracted modules | 13 | 5 (one `times_t` each) |
| hand-written `x(3) = x(3) + x(2) - x(1)` sites | 28 (CPU) | 0 |

**CPU suite: 18 passed, 0 failed, bit-exact**, wall-clock ratios 0.88–1.06 —
which is inside the 16% single-binary spread §5 of `refactor_strategy.md`
measured, so read it as "unchanged", not as a speed-up.

Twelve of the GPU tree's accumulate sites do **not** collapse to `time_end`,
because the idiom there is broken up by interleaved comments and blank lines
that the tool's pattern does not match. They were renamed and are correct;
they are simply not tidied. Collapsing them is worth a follow-up pass and
nothing depends on it.

`tools/bundle_timers.py` did the mechanical part with every count asserted, and
is kept because three of its checks are things that are easy to get wrong by
hand and quiet when wrong:

* **`time_step` and `time_step_prev` are not timers.** They are the MD
  integration step in fs. A blanket `s/time_/time%/` destroys both, and the
  code still compiles.
* **A bucket may be promoted from scalar to `(1:3)` and left unsubscripted.**
  `time_neigh` and `time_gap` were plain scalars. After promotion, a `write`
  with an `F13.3` edit descriptor accepts the whole rank-1 array without
  complaint and the format reverts — the symptom is garbled output several
  lines further down, nowhere near the change. Hit during this work; the tool
  now asserts against it.
* **An accumulate line missing its own `(3)` term assigns instead of
  accumulating.** `turbogap_vdw.f90:720` did, so vdW time reported only the
  *last* call rather than the total, on both branches. `time_end` always
  accumulates, so collapsing the idiom fixes it — the tool reports rather than
  hides that.

### 5.1 The negative "Miscellaneous", diagnosed

`docs/refactor_handoff.md` §6 records this as "a double-counting bug in the
timing buckets" that "makes the printed totals useless". It is two bugs:

* `time_gap` was accumulated **inside** the SOAP descriptor loop from a stamp
  taken *before* the loop, so an N-descriptor potential charged nearly the
  whole prediction phase N times.
* that window enclosed an `mpi_reduce` already charged to `time_mpi_ef`, and
  Miscellaneous subtracted both.

`time%gap` is now charged to the SOAP loop alone, once. The reallocation and
the reduce fall to Miscellaneous, which is honest.

The accounting rule is now stated in `timing.f90` and enforced in one place:
**buckets nest, `sum_times` adds only the parents, and a bucket whose nesting
is not established is not added.** An over-large Miscellaneous is a
correct-signed statement of ignorance; a negative one is a lie, and that lie is
why the in-code timers could not be used as refactor instrumentation.

Two holes are left deliberately, both now visible in Miscellaneous rather than
hidden in a wrong subtraction: the prediction-phase reallocation is untimed,
and the end-of-loop neighbour-array teardown opens an interval it never closes.
