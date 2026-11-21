# Developing TurboGAP

Where code goes, how to build, test and commit it, and how to document a
keyword.

## Once per clone

```sh
git submodule update --init --recursive
tools/setup_dev_env.sh                 # venv, pinned formatters, git hook
export PATH="$HOME/.venvs/turbogap-tools/bin:$PATH"
tools/setup_dev_env.sh --check         # exits non-zero if anything is missing
```

`src/soap_turbo` and `src/soap_turbo_gpu` are one repository on two branches,
not an old and a new version. Their `get_soap` signatures differ, so each build
takes its own and the architecture makefile says which.

## Build

```sh
./compile_cpu.sh                       # bin/turbogap
./compile_gpu.sh                       # bin-gpu/turbogap
KOKKOS=1 ./compile_gpu.sh              # bin-kokkos-gpu/turbogap; needs KOKKOS_ROOT
```

Both scripts take `TURBOGAP_ARCH`, `JOBS` and `DEBUG`. Every variant builds
into its own object tree, so all of them coexist and the suite can run one
against another: `OPENMP=0` (`bin-noomp/`), `DEBUG=1` (`-G` on the device,
about 2.1x slower; never time it), `PROFILE=1` (line info and NVTX markers for
Nsight), `KOKKOS=1`. Extra flags go through `F90_EXTRA` with a
`BUILD_TAG_EXTRA` naming their tree:

```sh
make F90_EXTRA=-fcheck=bounds BUILD_TAG_EXTRA=-bounds
```

A device architecture makefile is one with `GPU = 1`: it selects the `_gpu`
backends, pulls in `src/gpu/` and puts `-D _GPU` in `PP`. Kokkos is built by
`tools/install_kokkos.sh` and must be linked with the `nvcc` that built it;
`docs/LAMMPS_KOKKOS.md` has why it exists and what is missing.

**After adding or removing a source file**, regenerate both dependency files,
or `make -j` races and fails with `Cannot open module file`:

```sh
python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps
python3 tools/gen_fortran_deps.py . --gpu ST_DIR=src/soap_turbo_gpu/src \
    > makefiles/Makefile.deps.gpu
```

The `include` of them must stay last in the `Makefile`: placed before the first
rule, a dependency line becomes the default goal and `make` builds one object
and exits 0. CI fails on a stale file.

`make clean` removes objects, modules and binaries; `make deepclean` removes
the build directories. If `-j` fails after a clean, rebuild once serially: that
is a stale dependency file, not a flaky build.

## Layout

| layer | files | holds |
|---|---|---|
| driver | `src/turbogap.f90` | the step loop, one call per phase; no MPI, no `#ifdef` |
| features | `src/turbogap_<feature>.f90` | one module per feature, owning its state; called at fixed hooks |
| potential | `src/turbogap_evaluate.f90` | every energy and force term, in order |
| decomposition | `src/turbogap_domain.f90` | which rank owns which atoms; every broadcast and reduction of state and results |
| transport | `src/turbogap_comm.f90` | `comm_t`, `comm_bcast`, `comm_sum_to_root`, `comm_sum_all`, `comm_allgather`; serial stubs |
| kernels | `gap.f90`, `vdw.f90`, `electrostatics.f90`, `exp_utils.f90`, ... | the physics, on plain arrays |
| device | `src/gpu/` | CUDA/HIP kernels, written once for both backends; `src/gpu/README.md` |

One step: `structure_acquire` → `domain_sync_state` → `domain_build` →
`evaluate` → `md_step` or i-PI → `nested_step` → `mc_step` → `domain_sync_state`
→ `loop_end_step`.

## Adding an energy term

1. Kernel in its own module, taking plain arrays, no MPI, no TurboGAP types.
2. Its result arrays in `results_t` (`turbogap_results.f90`): allocated and
   zeroed in `results_prepare`, freed in `results_free`.
3. Call it from `evaluate`, in its place in the order.
4. Partial sums across ranks: add the family to the table in
   `domain_complete_contributions`. A global scalar the term needs mid-way goes
   through `comm_sum_all` (see `compute_estat`).
5. Keyword: default in `types.f90`, parse branch with its `!>` block in
   `read_files.f90` (below), `make docs`.
6. `SRC` in the `Makefile`; regenerate both dependency files.

A feature that is not a term (a sampler, a bias, an output) gets its own module
with init/step/finish hooks, called from the driver in the phase it belongs to.

## Rules

- New code reaches MPI only through `turbogap_comm`. Older code still calls it
  directly: `turbogap_setup`, `turbogap_vdw`, `exp_interface`, and one call
  each in `error`, `gap_interface` and `gpu_context_cpu`.
- A rank holds centres `i_beg:i_end` and their pairs `j_beg:j_end`.
  `neighbors_list` is global (fold with `modulo(j - 1, n_sites) + 1`); arrays
  passed as a slice take the local index. Mixing the two is correct on rank 0
  and wrong everywhere else, so test at three ranks.
- A refactor is bit-exact at every rank count. A change in results is its own
  commit, with a test that shows the new result is right.
- Fatal errors go through `turbogap_abort()`. A bare `stop` exits 0, and has
  hidden failures before.

## Documenting a keyword

`docs/keywords.md`, `docs/keywords.html` and `turbogap --help` are generated
from `src/read_files.f90` by `tools/keyword_docs.py`. Put a `!>` block directly
above the branch that parses the keyword:

```fortran
      !> @kw mc_move_max
      !> Largest displacement of a single-atom "move". The trial displacement
      !> is drawn uniformly up to this.
      !> @units A
      !> @modes mc
      !> @needs mc_types
      else if (keyword == 'mc_move_max') then
```

then `make docs`. `make check-docs` fails if a keyword has no block or an
output is stale; the pre-commit hook runs it on any change to `read_files.f90`
or `types.f90`.

| tag | meaning |
| --- | --- |
| `@kw` | required, first: every spelling the branch tests, in the order it tests them. The generator aborts if they disagree |
| `@units` | free text; omit for dimensionless and logical keywords |
| `@modes` | which of `predict`, `md`, `mc`, `soap` it acts in; omit for all. Drives `turbogap --help <mode>`. Not for `.gap` keywords |
| `@needs`, `@see` | keywords this one needs to have an effect; related keywords |
| `@status` | `ignored` (parsed, warned about, discarded) or `deprecated` |
| `@type`, `@default` | only where the generator cannot find the field itself |

Do not write the group, the fields, the type, the default or the side effects:
the generator reads them from the enclosing `read_options_*` subroutine, the
`params%<field>` assignments and `types.f90`. Potential-file keywords are
documented the same way inside `read_gap_hypers`, keyed by their `.gap` block.

An unknown keyword in `input` aborts the run; an unknown keyword in a `.gap`
block is skipped without a word, so a misspelling there shows up as a wrong
answer.

## Testing

```sh
tests/regression/run.sh                    # host, bit for bit against the baseline
tests/regression/run.sh --gpu              # device, against the host build, within a tolerance
tests/regression/run.sh --both             # both; fails if either does
tests/regression/run.sh vdw_ts co_predict  # named cases
tests/regression/run.sh --list
tests/<name>/run.sh                        # the standalone tests
```

The first run fetches the test data (`tests/fetch_test_data.sh`) to
`../turbogap_tests`; `tests/regression/make_baseline.sh` builds the frozen
baseline binary. `tests/regression/README.md` explains the cases and their
`case.conf`.

| variable | |
|---|---|
| `TURBOGAP_BIN` | binary under test |
| `TURBOGAP_REF_BIN` | reference binary, absolute path |
| `TURBOGAP_KEEP=1` | keep the staging directory |
| `TURBOGAP_BLESS=1` | regenerate a golden case's `expected/`, deliberately; never to turn a red suite green |
| `TURBOGAP_TIME_TOL` | fail if test over reference wall time exceeds this |

`run.sh` reads the binary as it goes, so do not rebuild while it runs; copy the
binary and point `TURBOGAP_BIN` at the copy. Read the suite's exit code, not
its last lines.

What a new feature needs:

- Forces and virial: a finite difference of the term's energy, with
  `tests/xrd_debye/fd_gradient.py`. Isolate the term by differencing it on and
  off (`--family estat`, `vdw`) or by a potential holding only it
  (`tests/fd_gap`). Run at one rank and at three (`--ranks 3`).
- Regression: a new path gets a case; `REFERENCE=golden` when the baseline
  cannot run it.
- Device: `run.sh --gpu`. The device sums in a different order, so it is not
  reproducible even against itself; judge a device change by the suite, never
  by `cmp`. A case may loosen `GPU_RTOL`/`GPU_ATOL` only with a host control --
  the input shifted by 1e-9 A, host against host -- showing the drift is
  amplification (KNOWN_ISSUES 13d). A Kokkos build is checked with
  `tools/verify_kokkos.sh`, which measures that noise first and subtracts it.
- Known defects go in `tests/regression/KNOWN_ISSUES.md`, with how to reproduce
  them.

Device checks outside the suite:

```sh
python3 tools/gpu_check_alloc_pairs.py     # every device buffer freed by its matching call
tools/gpu_memcheck.sh                      # a case under compute-sanitizer
tools/gpu_memcheck.sh --tool initcheck     # uninitialised device reads
```

## Committing

The pre-commit hook formats (`fprettify` for Fortran, `clang-format` for
`src/gpu`) and regenerates the keyword docs. When it changes a file, the commit
aborts; that is the normal flow:

```sh
git add <paths>                            # name them; never git add -A
git commit -F msg.txt                      # may abort, having reformatted
git diff -w --ignore-blank-lines           # Fortran: must be empty
tools/check_reformat_only.py HEAD src/gpu/*.cu   # CUDA/C++: layout only
git add <paths> && git commit -F msg.txt   # now it lands
```

`-w` settles it for Fortran, where `fprettify` only moves whitespace within a
line; `clang-format` also breaks lines, so CUDA needs the token comparison.
After an aborted commit, `git checkout HEAD -- <file>` discards a change;
`git checkout -- <file>` restores the hook's rewrite.
