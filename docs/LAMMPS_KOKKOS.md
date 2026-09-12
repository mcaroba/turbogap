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

## The remaining kernels, by what each one needs

51 kernels are still raw launches. They are not equally hard; cheapest first:

| file | kernels | needs |
| --- | --- | --- |
| `gap_predict.cu` | 4 | flat maps — `tg_parallel_for` as it stands |
| `gap_2b.cu` | 2 | flat maps |
| `mad_xrd.cu` | 6 | flat maps |
| `gap_soap_radial.cu` | 4 | flat maps. The first of the SOAP group |
| `mad_pdf.cu` | 10 | flat, except one shared-memory reduction |
| `mad_electrostatics.cu` | 3 | two shared-memory reductions |
| `gap_soap_angular.cu` | 8 | one 3-D grid — needs an MDRange policy |
| `gap_soap_radial_operator.cu` | 1 | one shared-memory kernel |
| `gap_soap_descriptor.cu` | 10 | 4 shared-memory, one 3-D grid, tiled transposes |
| `gap_soap_forces.cu` | 3 | 9 shared-memory arrays across 3 kernels. The hardest |

Two things the dispatch layer does not have yet, which the bottom half of that
table needs:

- **`tg_parallel_for_2d` / `_3d`.** Two kernels launch
  `dim3((n + tpb - 1) / tpb, n_max, k_max)`. Flattening the index by hand would
  work and would read badly; Kokkos has `MDRangePolicy` for exactly this.
- **A team policy with scratch.** `<<<n_sites, tpb>>>` with a block's threads
  cooperating over one atom is a `TeamPolicy`, and `__shared__` becomes
  `team.team_scratch(0)`. This is the real work in `gap_soap_forces.cu`.

Neither is hard. Both should be added when the first kernel actually needs one,
and checked against the control the same way, rather than written speculatively.

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
