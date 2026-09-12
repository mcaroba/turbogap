# TurboGAP through LAMMPS, and the Kokkos backend that gets it there

## Why there is a second device backend

TurboGAP's device code is CUDA, written once in the HIP API and compiled for
CUDA through the `src/hop` headers. That is enough to run TurboGAP's own
driver on a GPU, and not enough to run TurboGAP *inside* LAMMPS.

LAMMPS's accelerator package for GPUs is **KOKKOS**. A `pair_style` in that
package is handed its neighbour list, positions and types as `Kokkos::View`s
living in LAMMPS's memory space, and is expected to launch its kernels in
LAMMPS's execution space. A potential that can only be evaluated by launching
raw `<<<>>>` kernels against buffers it allocated itself cannot be that pair
style — not because the arithmetic is wrong, but because it cannot be handed
LAMMPS's data.

So the kernels have to be expressible in Kokkos. `-D_KOKKOS` is that step, and
it is deliberately *not* a second copy of the sources: one kernel body, two
dispatches, chosen in `src/gpu/gpu_backend.h`. A kernel cannot be fixed on one
path and left broken on the other, which is the failure mode a parallel
`src/kokkos/` tree would have had.

## The architecture this is aiming at

The model is [`symmetrix`](https://github.com/wcwitt/symmetrix), which does
exactly this for MACE, and whose layout is worth copying because it is already
known to work against LAMMPS `develop`:

```
libturbogap_kokkos/          a Kokkos library that knows nothing about LAMMPS
  CMakeLists.txt
  source/turbogap_kokkos.hpp

pair_turbogap/
  install.sh                 symlinks into LAMMPS, appends to its CMakeLists
  pair_turbogap.cpp/.h       the plain pair style        -> lammps/src/
  pair_turbogap_kokkos.cpp/.h  the Kokkos one            -> lammps/src/KOKKOS/
```

`install.sh` in `symmetrix` is nine lines: four symlinks and an
`add_subdirectory` appended to `lammps/cmake/CMakeLists.txt`. LAMMPS is then
built with `-D PKG_KOKKOS=ON -D TURBOGAP_KOKKOS=ON`. Nothing is forked.

## The part that is already true

The interesting thing about this port is how little translation the *data*
needs. `symmetrix` hands MACE flat one-dimensional views:

```cpp
Kokkos::View<int*>    num_neigh, first_neigh, neigh_indices, neigh_types;
Kokkos::View<double*> xyz, r;
```

That is TurboGAP's device layout already. Its kernels take
`n_neigh`, `n_neigh_index`, `neighbors_list`, `rjs` and `xyz` as flat device
arrays with exactly that meaning — a CSR neighbour list with the separations
precomputed. The `extern "C"` entry points in `gap_gpu.h` and `mad_gpu.h` take
raw pointers to them.

So the bridge is `Kokkos::View<double*, MemoryTraits<Unmanaged>>` wrapped
around a pointer LAMMPS owns, rather than a conversion. That is why the port
is worth doing in this direction: it is mostly re-expressing launches, not
re-expressing physics.

## What is done, and what is not

Done:

- `src/gpu/gpu_backend.h` — the dispatch layer. `TG_LAMBDA`,
  `tg_parallel_for`, `tg_inclusive_scan_int`, `tg_reduce_sum_int`, the
  execution-space instance adopted from the stream `gpu_context` already owns,
  and Kokkos initialisation tied to the card `cuda_set_device` picked.
- `src/gpu/gpu_scan.cu` — the reduction and the inclusive scan, which under
  `-D_KOKKOS` are `Kokkos::parallel_reduce` and `Kokkos::parallel_scan`
  instead of a hand-written block reduction and Blelloch scan. Both are
  integer, so the two backends agree exactly.
- `KOKKOS=1`, its own object tree, `tools/install_kokkos.sh`,
  `tools/verify_kokkos.sh`.
- `src/gpu/gpu_scatter.cu` — all nine passes of the deterministic scatter.
  They are flat maps and per-chunk serial sums, so the conversion does not
  touch the ordering the file exists to guarantee.

## About `src/soap_turbo_gpu`

Worth stating plainly, because the name misleads: **the GPU `soap_turbo`
submodule contains no device code at all.** It is 4100 lines of Fortran, and
every kernel it drives lives in this repository:

| what `soap_turbo_gpu` calls | where the kernel is |
| --- | --- |
| `gpu_radial_poly3`, `gpu_radial_poly3gauss` | `gap_soap_radial.cu` |
| `gpu_radial_poly3operator` | `gap_soap_radial_operator.cu` |
| `gpu_get_plm_array_global`, `gpu_get_exp_coeff_array`, `gpu_get_cnk` | `gap_soap_angular.cu` |
| `gpu_get_sqrt_dot_p`, `gpu_soap_normalize`, `gpu_get_derivatives` | `gap_soap_descriptor.cu` |
| `gpu_get_soap_der` | `gap_soap_forces.cu` |

So porting `soap_turbo_gpu` needs **no change to the submodule**: its Fortran is
already backend-agnostic, and the work is converting those five files here.

## What is converted, and what is not

| file | state |
| --- | --- |
| `gpu_scan.cu` | Kokkos `parallel_reduce` and `parallel_scan`; the hand-written primitives stay as the CUDA fallback |
| `gpu_scatter.cu` | all 9 passes |
| `mad_xrd.cu` | all 5 |
| `gap_predict.cu` | all 4 |
| `gap_2b.cu` | both |
| `gap_soap_radial.cu` | all 3 |
| `gap_soap_angular.cu` | all 6, one through `tg_parallel_for_3d` |
| `mad_electrostatics.cu` | 2 of 3 |
| `mad_pdf.cu` | 4 of 10 |
| `gap_soap_descriptor.cu` | 2 of 9 |
| `gap_soap_radial_operator.cu` | none |
| `gap_soap_forces.cu` | none |
| `gap_3b.cc` | none |

### Why the rest are not converted

Four reasons, and they are not equally serious. The first is the one that
matters.

**It would change the numbers.** These kernels finish with a `__shared__` tree
reduction over doubles -- `sh[tid] += sh[tid + s]`, halving until one value
remains:

- `kernel_electrostatics_gsf` (`mad_electrostatics.cu`)
- `kernel_reduce_pair_distribution` (`mad_pdf.cu`)
- all three kernels in `gap_soap_forces.cu`
- `cuda_get_soap_der_two_one` (`gap_soap_descriptor.cu`)

Floating-point addition is not associative, so a Kokkos team reduction summing
the same values in a different order gives a different answer. Converting these
is a physics change wearing a refactor's clothes, and it would land as exactly
the last-digit drift the device already has too much of. Each needs a
`TeamPolicy` that reproduces *this* tree order, and a check of its own. The
reason is written beside each kernel, so nobody converts one by reflex.

**A block cooperating without reducing.** `gap_soap_radial_operator.cu` stages
`A_d` into dynamic shared memory and gives each thread scratch;
`naive_transpose_soap_rad_azi_pol` is a tiled transpose. Both use `__shared__`
purely to share data, with one barrier and no accumulation, so there is no
numerical hazard at all -- they need `TeamPolicy` with `team_scratch` and
nothing more.

**A grid index that is part of the addressing.** The rest of
`gap_soap_descriptor.cu` -- `_der_two_two`, `_thr_one`, `_thr_two`,
`cuda_soap_normalize`, `cuda_get_derivatives_new_new` -- index by block as well
as thread. `MDRangePolicy` or an explicit flattening; mechanical, just not
done.

**No Kokkos equivalent.** `gap_3b.cc`'s main kernel is a template carrying
`__maxnreg__` (through `TG_3B_REGCAP`) to cap its register count. Kokkos offers
`LaunchBounds`, which is `__launch_bounds__` -- a different mechanism, bounding
occupancy rather than registers. Converting would silently drop the cap on
CUDA >= 12.4, where it is live. Worth doing only alongside a measurement
showing the register count did not run away.

### About `tg_parallel_for_3d`

`cuda_get_cnk_one_new_new` indexes site, radial `n` and angular `k`, carried on
the grid's x, y and z. Kokkos tiles an `MDRangePolicy` differently, so the
order the triples are visited is not the same. That is safe here because each
triple runs a serial neighbour loop and writes one `cnk` entry, reading nothing
another triple writes -- and would be wrong the moment that stopped holding.
The CUDA path keeps the original y/z grid exactly.


## Then the part that is not kernels

1. **A compute entry point that does not come from Fortran.** Every device
   entry point today is reached through `bind(C)` from a Fortran driver that
   also owns the neighbour list, the GAP hyperparameters and the SOAP
   compression. A pair style has none of that. This needs a C++ object that
   holds the model and exposes `compute(positions, neigh, types) -> (energy,
   forces)`.
3. **The `.gap` file, read from C++.** Presently parsed by
   `src/read_files.f90`. `symmetrix` sidesteps the equivalent problem by
   exporting the model to JSON from Python first, which is worth copying.
4. **The pair style itself**, and `install.sh`.

Step 2 is the real work, not step 1. The kernels are mechanical; the model
object is a genuine interface design, and it is the thing to think hardest
about before writing any of it.

## Checking a Kokkos build

Not with a plain diff. The device binary is not reproducible against itself —
run `estat_gsf` twice with one unchanged CUDA binary and `energy_soap` moves in
the tenth digit. A Kokkos-against-CUDA comparison would charge the port for all
of that.

`tools/verify_kokkos.sh` runs the suite twice instead: the CUDA binary against
a copy of itself, which measures the noise, and the Kokkos binary against the
CUDA one. What matters is the difference of the two failure sets.

```sh
tools/verify_kokkos.sh                    # every case
tools/verify_kokkos.sh estat_gsf xrd_predict
```
