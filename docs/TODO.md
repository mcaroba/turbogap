# TODO: `state_t`, `simulation_t`, and overlapping domain decomposition

Written 2026-08-06. Companion to `docs/from_turbogap_2.0.md` §3–§4, which argues
*why*; this says *what to declare and how to use it*.

**Nothing here is to be started before the merge lands.** It changes argument
lists on both branches at once, which is the one thing the whole plan has been
avoiding, and the GPU branch cannot check it bit-exactly (§6.5 of
`refactor_strategy.md`). Declaring the types changes no behaviour and is free;
migrating consumers is not.

---

## 0. You already have `state_t`. It is called `image`.

`src/types.f90:250` defines `type image`; `turbogap.f90` mentions `images(` **97
times**; `from_properties_to_image` / `from_image_to_properties` marshal the
driver's loose variables in and out of it. Above that pair, in the original
author's hand:

> `! In time I should make the image data type the default way to store these`
> `! properties!!!!!!!`

So this is not a new abstraction. It is finishing one that is half-built and
currently second-class: the truth lives in the driver's loose variables, and an
`image` is a *snapshot* copied out when MC or nested sampling needs one and
copied back afterwards. Every round trip is a chance for the copy and the
original to disagree — this codebase's characteristic defect wearing a hat.

**`image` already nests along the *ensemble* axis** (nested-sampling walkers, MC
trial/current). Domain decomposition adds a second, *spatial* axis to the same
type. That is the argument for one nestable type rather than three parallel
ones.

## 1. Why `image` is not nestable as written

It mixes three things that nest three different ways:

| kind | fields in `image` | how it nests spatially |
|---|---|---|
| per-site extensive | `positions`, `positions_prev`, `velocities`, `forces`, `forces_prev`, `masses`, `species`, `energies`, `local_properties`, `fix_atom`, `xyz_species` | decomposable — a domain holds a slice **plus halo** |
| the cell | `a_box`, `b_box`, `c_box`, `indices` | **replicated**, one writer, read-only to domains |
| reduced observables | `energy`, `e_kin`, `energy_exp` | meaningful **only at the level that did the reduce** |

Along the ensemble axis all three are per-walker, so a flat type works — which
is why nobody has noticed. Along the spatial axis a flat type forces one meaning
on all three. Concretely: if a domain's `state` has `%instant_temp`, is that the
global temperature replicated, or this domain's partial? Both are defensible.
Silently having both is not.

## 2. The declarations

Split **by how each part nests**, not by subject. Note this is the opposite of
what the graph report's betweenness question implies: betweenness stays high and
should — see `from_turbogap_2.0.md` §3.3.

```fortran
!**************************************************************************
!   The decomposable part: sites this rank holds, owned and halo alike.
!
!   n_local counts everything in the arrays; n_owned counts the prefix that
!   contributes to global sums.  Keeping both, plus is_owned, is what stops a
!   halo atom's energy being counted by two domains -- a wrong number, not a
!   crash.  This is the same shape as gap_interface counting pairs with < and
!   filling with <=, at much larger scale, so the count and the fill must come
!   from ONE place: n_owned.
   type sites_t
      integer :: n_local = 0                    !! owned + halo
      integer :: n_owned = 0                    !! contributes to global sums
      integer, allocatable :: global_index(:)   !! local -> global site id
      logical, allocatable :: is_owned(:)       !! redundant with n_owned when
                                                !! owned sites are the prefix;
                                                !! keep ONE as authoritative

      real(dp), allocatable :: positions(:, :)
      real(dp), allocatable :: positions_prev(:, :)
      real(dp), allocatable :: velocities(:, :)
      real(dp), allocatable :: forces(:, :)
      real(dp), allocatable :: forces_prev(:, :)
      real(dp), allocatable :: masses(:)
      real(dp), allocatable :: energies(:)
      real(dp), allocatable :: local_properties(:, :)
      integer,  allocatable :: species(:)
      logical,  allocatable :: fix_atom(:, :)
      character*8, allocatable :: xyz_species(:)
   end type sites_t

!**************************************************************************
!   Replicated, one writer.  Under NPT the box changes and that change is
!   BROADCAST, never recomputed per domain -- two domains deriving the same
!   volume independently is the antipattern at type scale.
   type cell_t
      real(dp) :: a_box(1:3) = [1.d0, 0.d0, 0.d0]
      real(dp) :: b_box(1:3) = [0.d0, 1.d0, 0.d0]
      real(dp) :: c_box(1:3) = [0.d0, 0.d0, 1.d0]
      real(dp) :: volume = 1.d0
      integer  :: indices(1:3) = 1
   end type cell_t

!**************************************************************************
!   Outputs of a reduction.  Filled ONLY at the level whose reduce produced
!   them; left at zero elsewhere, so reading one at the wrong level gives an
!   obviously wrong answer rather than a plausible one.
   type observables_t
      integer  :: level = LEVEL_NONE            !! see the parameters below
      real(dp) :: energy = 0.d0
      real(dp) :: e_kin = 0.d0
      real(dp) :: energy_exp = 0.d0
      real(dp) :: instant_temp = 0.d0
      real(dp) :: instant_pressure = 0.d0
      real(dp) :: virial(1:3, 1:3) = 0.d0
   end type observables_t

   integer, parameter :: LEVEL_NONE = 0, LEVEL_DOMAIN = 1, LEVEL_GLOBAL = 2

!**************************************************************************
!   What `image` should become.  Migrating `image` to this and deleting the
!   two marshalling routines is the FIRST consumer to convert -- it removes
!   97 call sites' worth of copying and is testable on its own.
   type state_t
      type(sites_t)       :: sites
      type(cell_t)        :: cell
      type(observables_t) :: obs
   end type state_t
```

### What must stay out

The mirror of the `control_t` rule (`from_turbogap_2.0.md` §3.2):

* **`params` / `control_t`** — policy, one immutable copy, passed alongside.
* **the GAP hypers** (`soap_turbo_hypers`, `distance_2b_hypers`, …) — model, not
  state; read-only and shared.
* **the neighbour list** — it *is* per-domain, but it is *derived* from state
  plus a cutoff, so it belongs in `simulation_t` **beside** `state_t`, not
  inside it. `image` already gets this right: it carries no neighbour list.

### The `_prev` trap

`n_sites_prev`, `volume_prev`, `indices_prev`, `E_kinetic_prev` exist for change
detection. Under nesting, *"did `n_sites` change"* is a **global** question (an
MC insert anywhere changes it) while *"did my domain's site count change"* is a
**local** one — same name, different question, both true at different levels.

This is already visible: `n_sites /= n_sites_prev .or. params%do_mc` appears at
`turbogap.f90` 1039, 1132 and 1176. Exactly one place must decide "rebuild", and
that place must say which level it is asking about.

---

## 3. `domain_t` and `simulation_t`

```fortran
!**************************************************************************
!   A spatial region owned by one group of ranks, plus a halo of width
!   rcut_max + buffer.  It is a simulation for the short-ranged additive
!   terms and ONLY for those -- see §4.
   type domain_t
      integer  :: id   = 0
      integer  :: comm = MPI_COMM_NULL      !! ranks sharing this domain
      integer  :: root = 0                  !! rank in comm that aggregates
      integer  :: n_domains = 1
      real(dp) :: lo(1:3) = 0.d0            !! this domain's spatial extent
      real(dp) :: hi(1:3) = 0.d0
      real(dp) :: halo_width = 0.d0         !! rcut_max + buffer
      integer  :: coords(1:3) = 0           !! position in the 3-D grid
      integer  :: neighbour_rank(1:26) = MPI_PROC_NULL   !! 3^3-1 directions
   end type domain_t

   type simulation_t
      type(domain_t)         :: domain
      type(state_t)          :: state       !! owned + halo sites only
      type(neighbors_t)      :: neighbors
      type(split_t)          :: split       !! intra-domain rank split
      type(contribution_ref) :: contrib(1:N_CONTRIB)
      type(times_t)          :: time
   end type simulation_t
```

Four of those six fields already exist or are the near-term work: `times_t`
landed (`d72580f`), `contrib` landed (`7390225`), `split_t` is a trivial bundle
of `i_beg/i_end/j_beg/j_end`, and `neighbors_t` is a bundle of what
`neighbors.f90` already produces. **The bundling owed now is the down payment on
this**, which is why it is worth doing on its own terms.

---

## 4. Which contributions are domain-local — declare it, do not assume it

Conflating these is how a decomposition produces wrong physics rather than a
crash.

| | families | why |
|---|---|---|
| **domain-local** | `C_SOAP`, `C_2B`, `C_3B`, `C_CP`, `C_LP`, TS pairwise vdW | short-ranged and additive; a halo of `rcut_max` is sufficient |
| **global** | `C_PDF`, `C_SF`, `C_XRD`, `C_ND` | structure-wide sums over *all* pairs and all q-points |
| **global** | MBD | a global eigenproblem, not a sum |
| **global** | `C_ESTAT` | long-ranged |
| **global** | XPS | depends on a global spectrum normalisation |
| **not decomposable at all** | nested sampling, MC insert/remove | whole-configuration moves that change `n_sites` |

`contrib_on` and the `C_*` indices already enumerate exactly these eleven
families, so the encoding is a companion table beside them:

```fortran
   logical, parameter :: C_IS_LOCAL(1:N_CONTRIB) = [ &
      .true.,  &   ! C_SOAP
      .false., &   ! C_VDW    -- TS is local, MBD is not; see below
      .false., &   ! C_ESTAT
      .true.,  &   ! C_LP
      .false., &   ! C_PDF
      .false., &   ! C_SF
      .false., &   ! C_XRD
      .false., &   ! C_ND
      .true.,  &   ! C_2B
      .true.,  &   ! C_CP
      .true.   ]   ! C_3B
```

`C_VDW` is the one that cannot be answered by the table alone — TS is pairwise
and local, MBD is not, and they share a family index. Either split the family or
make the entry a function of `params%vdw_type`. **Decide this before writing the
reduce**, not during.

The MPI reduce block since `7390225` already walks the family list **once**, so
routing each family to an intra-domain or a global reduction is a table lookup
rather than a redesign. That is the concrete return on Phase 3.

---

## 5. Overlapping domain decomposition: the actual work

### 5.1 Communicators

Two levels, built once after `read_input_and_gap_files` and rebuilt only if the
decomposition changes:

```fortran
!  1. Split the world into domains.  Ranks within a domain share its sites and
!     split work over them exactly as split_tasks does today.
   call MPI_Comm_split(MPI_COMM_WORLD, dom%id, rank, dom%comm, ierr)
   call MPI_Comm_rank(dom%comm, rank_in_domain, ierr)
   dom%root = 0

!  2. A communicator over the domain roots, for the global reductions in §4.
!     Non-root ranks pass MPI_UNDEFINED and get MPI_COMM_NULL back.
   colour = merge(0, MPI_UNDEFINED, rank_in_domain == dom%root)
   call MPI_Comm_split(MPI_COMM_WORLD, colour, dom%id, comm_roots, ierr)

!  3. A Cartesian communicator for the halo exchange, which gives the 26
!     neighbour ranks and handles periodicity for free.  periods = .true. is
!     what makes the wrap-around a communicator property rather than something
!     the exchange code has to special-case.
   call MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, .true., cart_comm, ierr)
   call MPI_Cart_coords(cart_comm, rank, 3, dom%coords, ierr)
```

The two-level reduce then replaces today's single `mpi_reduce` over
`MPI_COMM_WORLD`:

```
   local terms   -> MPI_Reduce over dom%comm     to dom%root
   global terms  -> MPI_Allreduce over comm_roots, then broadcast down dom%comm
```

**This is the cheap part.** Budget almost nothing for it.

### 5.2 The halo exchange, which is where the bugs are

Per iteration, before the neighbour-list build:

1. Each domain sends the positions/species of its owned sites within
   `halo_width` of each face to the 26 neighbours, and receives theirs.
2. Received sites append to the local arrays **after** the owned prefix, so
   `1:n_owned` stays the owned set and `n_owned+1:n_local` is halo. Keeping
   owned as a contiguous prefix is what lets every sum be
   `sum(x(1:n_owned))` rather than a masked reduction, and it removes the
   `is_owned` array as a second source of truth. **Keep the prefix invariant;
   do not keep both.**
3. Forces computed on halo sites must be sent **back** and accumulated onto
   their owner (Newton's third law), or dropped, depending on whether the
   kernel computes each pair once or twice. `gap_interface` currently does the
   latter. Getting this wrong is silent.

**Write the test before the code.** A small cell — the P4 dimer or the 108-atom
vdW cell — split over 2 and 4 domains, with at least one atom whose halo image
lives on another domain, compared against the single-domain result. §1.3 of
`session_handoff_2026-08.md`: a single system size cannot see this class of bug,
and the 7176-atom CO system is far too large to see it at all.

### 5.3 What breaks and needs a decision

* **MC insert/remove changes `n_sites` globally.** A domain cannot decide alone.
  Either serialise the move through `comm_roots`, or keep MC replicated and
  decompose only the energy evaluation.
* **Nested sampling** is the ensemble axis, not the spatial one. A walker is a
  `simulation_t`; walkers may themselves be decomposed. Two nestings, and the
  communicator hierarchy becomes three levels.
* **The box under NPT** is replicated state that changes. One writer, broadcast.
* **`rcut_max` sets `halo_width`.** If the box shrinks below `2*rcut_max` the
  decomposition is invalid and must fall back to replicated. The driver already
  has the analogous check for the supercell; reuse that logic, do not re-derive
  it.

### 5.4 On the GPU branch: a batch and a domain are the same abstraction

`src/gpu_context.f90` already holds `gpu_streams(:)` and `cublas_handles(:)`,
and the driver builds `i_beg_list`/`i_end_list`/`j_beg_list`/`j_end_list` at two
sites and consumes them at three. **That is a proto-decomposition**: a list of
index ranges, each with its own stream.

Design `simulation_t` so a GPU batch is a `simulation_t` with
`comm = MPI_COMM_SELF`, and the batching machinery and the domain machinery
become one thing — the diffraction block's 51-against-37 argument gap, which is
entirely batch bookkeeping, closes as a side effect. Do not, and the GPU branch
grows a second decomposition to keep in agreement with the first forever.

**This is the only part of the design the merge sequencing has to accommodate,
and it should be agreed before anything is written.**

---

## 6. Order of work

1. Declare the types. No behaviour change, no risk. Can happen any time.
2. Migrate `image` to `state_t` and delete `from_properties_to_image` /
   `from_image_to_properties`. One consumer, 97 call sites, testable on the
   18-case suite, and it pays for itself immediately.
3. `split_t`, then `neighbors_t`. Cheap bundles; both are `simulation_t` fields.
4. Measure nested sampling again with `tools/classify.py`. Phase 6 got it to
   87 arguments with contributions and timers bundled and called that still not
   worth extracting; per-site state is the third bundle and *may* tip it.
   **That is a projection, not a measurement — measure before relying on it.**
5. `simulation_t` and the decomposition, after the merge, with §5.2's test
   written first.

## 7. Performance rule, so this does not get relitigated

Array components of a derived type are as contiguous as bare allocatables, and
passing `state%sites%positions` to a `real(dp), intent(in) :: x(:,:)` dummy is
not a copy. This was measured on `contribution_ref`: 19.87 s either way, because
the hot compute lives in callees operating on **dummy arguments**, where
aliasing is governed by Fortran's argument rules rather than by anything the
driver does.

**Bundling stops at the phase boundary. Kernels keep taking plain arrays.**
Nothing in `soap_turbo`, `gap.f90` or `vdw.f90` learns about these types.
