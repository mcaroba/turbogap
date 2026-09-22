# The TurboGAP code

**TurboGAP** (c) 2018-2026 by **Miguel A. Caro** and others (see "contributors" below
for detailed authorship info).

[www.turbogap.fi](http://www.turbogap.fi)

## Contributors, copyright and license

The following people, listed chronologically, have contributed code or ideas to the **TurboGAP**
code. Those whose names are in bold have contributed code _to the master branch_ (relevant for
purpose of copyright; each file in the **TurboGAP** repo that contains code has a copyright
statement at the beginning).

- **Miguel A. Caro** (Aalto University)
- Patricia Hernández-León (Aalto University)
- Suresh Kondati Natarajan (formerly @ Aalto University)
- **Albert P. Bartók-Pártay** (Warwick University)
- Eelis V. Mielonen (formerly @ Aalto University)
- **Heikki Muhli** (Aalto University)
- **Mikhail Kuklin** (formerly @ Aalto University)
- Gábor Csányi (University of Cambridge)
- **Jan Kloppenburg** (Aalto University)
- **Richard Jana** (Aalto University)
- **Tigany Zarrouk** (Aalto University)
- **Uttiyoarnab Saha** (Aalto University)

**TurboGAP** is licensed under the Academic Software License (ASL), an "available source"
non-commercial license. This means that you are free to use and distribute the code for
non-commercial academic research (or teaching) under the terms of the license. See
`LICENSE.md` for details. If you want to obtain a commercial license for **TurboGAP**, please
contact Miguel Caro (mcaroba@gmail.com).

Some third-party code is included with **TurboGAP** for convenience, under the
`src/third_party` directory. These codes are licensed independently from **TurboGAP** and their
respective licenses have been verified to be compatible for redistribution with **TurboGAP**.
They may be redistributed separately from **TurboGAP** under their respective licenses.
Refer to each piece of software in that subdirectory for further information.

The **soap_turbo** submodule, under `src/soap_turbo`, is a separate distribution from
**TurboGAP**, but it is required
for running **TurboGAP**, since it contains the `soap_turbo` routines. These routines are
copyright (c) of Miguel A. Caro and they are also distributed under the ASL. Therefore, you
can freely use this code for non-commercial academic research or teaching. If you want to
obtain a commercial license for **soap_turbo** please contact Miguel Caro (mcaroba@gmail.com).

`git clone --recursive` pulls three submodules, none of which needs any
further setup:

| submodule            | from                          | why                                          |
| -------------------- | ----------------------------- | -------------------------------------------- |
| `src/soap_turbo`     | `TiganyZ/soap_turbo` `master` | the CPU SOAP routines                        |
| `src/soap_turbo_gpu` | `TiganyZ/soap_turbo` `gpu`    | the GPU port of soap_turbo                   |
| `src/hop`            | `cschpc/hop` `master`         | HIP/CUDA portability headers for a gpu build |

Which `soap_turbo` is compiled is set by `ST_DIR`, defaulting to
`src/soap_turbo/src` and overridden to `src/soap_turbo_gpu/src` by a GPU arch
makefile. `HOP_ROOT` defaults to `src/hop`, so a device build needs nothing in
the environment; setting it still overrides.

## Overview of the code

The **TurboGAP** code consists of a series of Fortran routines written by Dr. Miguel A. Caro
and others. These routines are designed to efficiently and
accurately build many-body atomic descriptors and carry out other related computations
employed in machine-learning approaches to atomistic modeling, such as to run molecular
dynamics simulations. The main current functionality is the computation of SOAP-type
descriptors [1], which were originally introduced by Bartók, Csányi et al. [2] in the
context of the Gaussian Approximation Potential framework [3].

**TurboGAP** is a primitive but efficient interface to the **soap_turbo** library.
This native interface is currently restricted to limited functionality but reasonably
fast (can outperform QUIP+LAMMPS in most situations); however it may be buggy, is
undocumented, and can be "temperamental". If you want to use
**soap_turbo** routines to run molecular dynamics or to carry out other simulations involving
heavy use of CPU power, without the worries of using the native interface, you are advised to
use QUIP. However, some new or experimental features (e.g.,
full support for van der Waals corrections) may only be available via the native interface.
If you're feeling adventurous, and know what you're doing, you are more than welcome to
use the **TurboGAP** interface, and feedback can be sent to Miguel Caro (mcaroba@gmail.com)
or left on the Issues section of the Github page.

[1] M.A. Caro. [Phys. Rev. B 100, 024112
(2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.024112).  
[2] A.P. Bartók, R. Kondor, G. Csányi. [Phys. Rev. B 87, 184115
(2013)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.184115).  
[3] A.P. Bartók, M.C. Payne, R. Kondor, G. Csányi. [Phys. Rev. Lett. 104, 136403
(2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.136403).

## Installation

_**tl;dr, CPU**_ (standard linux, prerequisites installed):

```sh
git clone --recursive https://github.com/mcaroba/turbogap.git
cd turbogap
./compile_cpu.sh
export PATH="$(realpath bin):$PATH"
```

_**tl;dr, GPU**_ (CUDA or ROCm, on top of the above):

```sh
./compile_gpu.sh
export PATH="$(realpath bin-gpu):$PATH"
```

Remember you must specify the architecture of your GPU: so you must change the `Makefile.<TURBOGAP_ARCH>` to reflect your machine.

Both binaries are built into separate directories: `bin/turbogap`
and `bin-gpu/turbogap`.

_**tl;dr, GPU through Kokkos**_ (optional; needed to run TurboGAP from LAMMPS):

```sh
tools/install_kokkos.sh              # builds Kokkos, ~15 min, once
export KOKKOS_ROOT=$HOME/.local/kokkos-4.7.04
KOKKOS=1 ./compile_gpu.sh
export PATH="$(realpath bin-kokkos-gpu):$PATH"
```

This is the same device code, dispatching its kernels through Kokkos instead of
raw CUDA launches, into its own `bin-kokkos-gpu/`. It exists because LAMMPS's
accelerator package is Kokkos, so a `pair_style` calling TurboGAP has to run in
LAMMPS's execution space on LAMMPS's memory -- see `docs/LAMMPS_KOKKOS.md` for
what that needs and how far it has got. `tools/install_kokkos.sh` also installs
`cmake` from PyPI if the machine has none.

`compile_cpu.sh` and `compile_gpu.sh` default to `Ubuntu_gfortran_mpi` and
`Aalto_gfortran_openblas_hip_cuda`; pick another from `makefiles/` with
`TURBOGAP_ARCH=<name> ./compile_cpu.sh`. Each script fetches the submodules,
checks for the compilers it needs, and refuses an architecture meant for the
other script rather than producing a binary that is not what you asked for.

_**tl;dr, developer tooling**_ (only needed to commit, or to run the python
test drivers):

```sh
tools/setup_dev_env.sh
export PATH="$HOME/.venvs/turbogap-tools/bin:$PATH"
tools/setup_dev_env.sh --check
```

### What gets installed

A `--recursive` clone brings three submodules:

| submodule            | what it is                                                                                                                                                                                                                                |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/soap_turbo`     | the SOAP routines the **host** build uses                                                                                                                                                                                                 |
| `src/soap_turbo_gpu` | the SOAP routines the **device** build uses. A second implementation, not a later version: its `get_soap` takes device pointers and it represents the descriptor compression differently, so both are needed and each build takes its own |
| `src/hop`            | the HIP-on-CUDA headers the device build compiles against, so one source tree targets both vendors. No `HOP_ROOT` to set                                                                                                                  |

`tools/setup_dev_env.sh` creates a virtual environment at
`~/.venvs/turbogap-tools` holding, all pinned:

| tool                        | why                                                                                                                                 |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `fprettify`, `clang-format` | the formatting the pre-commit hook applies. Pinned because an unpinned formatter rewrites files nobody touched                      |
| `pre-commit`                | the hook itself: formatting, the generated keyword reference, whitespace                                                            |
| `numpy`                     | what the python test drivers compare against. Without it `tests/ir_fft`, `tests/mad_ir` and `tests/ipi_socket` skip rather than run |
| `i-PI`                      | drives TurboGAP over a socket in `turbogap ipi` mode, for path-integral and other advanced dynamics.                                |
| `cmake` | only installed by `tools/install_kokkos.sh`, and only when the machine has none of its own. Kokkos has no other build system |

Kokkos is **not** part of that environment: it is a from-source CUDA build of
about fifteen minutes, and only a `KOKKOS=1` device build needs it. Install it
separately with `tools/install_kokkos.sh`, or ask for it at the same time with
`tools/setup_dev_env.sh --with-kokkos`. It lands in `~/.local/kokkos-<version>`
and is found through `KOKKOS_ROOT`.

### Prerequisites

To use TurboGAP you must have the libraries `openmpi`, `lapack` + `blas` (or `openblas`)
and the Fortran compilers necessary for your system (e.g. `gcc`).

_For more detailed info on download, installation, etc., you can visit the
[TurboGAP wiki](https://turbogap.fi/wiki/index.php/Installation)._

### Getting the code

To get the **TurboGAP** code and the necessary **soap_turbo** routines, do a recursive
`git clone` with `--depth=1` to speed up things:

```sh
git clone --recursive --depth=1 http://github.com/mcaroba/turbogap.git
```

### Building

To build the **TurboGAP** binary and library, you need to select the options
that best match your architecture by exporting the environment variable
`TURBOGAP_ARCH=<your_architecture>` which is used in the `Makefile` e.g.:

```sh
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
```

This will include the architecture specific examples found in `makefiles`, e.g.
the above will include the specification for that architecture,

```sh
include makefiles/Makefile.Ubuntu_gfortran_mpi
```

A list of example makefiles is provided under the `makefiles/` directory for different systems.

Once you are happy with your `Makefile`, to build the code just type

```sh
make j<N_processes>
```

where `<N_processes>` is the number of processes you wish to build the code with.

Then add `turbogap/bin` to your path.

```sh
turbogap_dir=$(realpath bin)
export PATH="$turbogap_dir:$PATH"
echo "export PATH=\"${turbogap_dir}:\$PATH\"" >> ~/.bashrc
```

If you need to rebuild the code,
you can `make clean; make` or `make deepclean; make`.

## Running TurboGAP

**TurboGAP** can be used to run static (single-point) calculations (`turbogap
predict`) or molecular dynamics (`turbogap md`) or (`turbogap mc`) for
Monte-Carlo. For details, documentation and up-to-date information refer to the
[TurboGAP wiki](http://turbogap.fi).

If you need help with TurboGAP modes or its keywords, you can either consult
the html document (which is searchable)

```sh
# On Linux
xdg-open docs/keywords.html
```

```sh
# On Mac
open docs/keywords.html
```

or the markdown `docs/keywords.md`, or one can use the `--help`:

```sh
turbogap --help
```

or

```sh
turbogap --help <topic>
```

where topic is a simulation mode (`predict`, `md`, `mc`, `soap`) or `gap`.

Given a mode, the keywords for the `input` file are listed with descriptions
and their dependencies, filtered to those that do something in that mode.
`turbogap --help gap` lists the keywords of the potential (`.gap`) file
instead, grouped by the block they belong to (`soap_turbo`, `distance_2b`,
`angle_3b`, `core_pot`). With no topic at all, both files are listed.

## TurboGAP Tutorials

For simple test cases of various simulation modes, one can consult `tests/regression/cases`.

For practical use cases (i.e. useful simulations for understanding a particular
research question) one can follow the tutorials which were ran for the TurboGAP
School, and adapt them for their machine: [TurboGAP
School](https://github.com/mcaroba/TurboGAP_School).

Simple tutorials can be found in the [TurboGAP Tutorials](https://github.com/TiganyZ/turbogap_tutorials) repository `git clone https://github.com/TiganyZ/turbogap_tutorials.git`. Further tutorials and all documentation can be found on the [Turbogap Website](https://turbogap.fi/wiki/index.php/Tutorials).

## Testing TurboGAP

The test systems -- trajectories and fitted potentials -- are not in this
repository. They live in [TurboGAP
Tests](https://github.com/TiganyZ/turbogap_tests), and
`tests/fetch_test_data.sh` clones them beside the source tree. Every suite
calls it when the clone is missing, so a first run on a fresh checkout works;
`TURBOGAP_DATA_ROOT` points at a checkout kept somewhere else.

One source builds four ways, and what you can test is decided by what the
machine has. Each build has its own object tree and its own `bin`, so they
coexist and can be run against one another:

| build | build it with | needs | binary |
| --- | --- | --- | --- |
| host | `./compile_cpu.sh` | gfortran, MPI, LAPACK/BLAS | `bin/turbogap` |
| serial | `make TURBOGAP_ARCH=Ubuntu_gfortran BUILD_TAG_EXTRA=-serial` | the same, without MPI | `bin-serial/turbogap` |
| device | `./compile_gpu.sh` | a CUDA or HIP toolchain, and a GPU to run on | `bin-gpu/turbogap` |
| Kokkos | `KOKKOS=1 ./compile_gpu.sh` | the above, and Kokkos at `KOKKOS_ROOT` built by the same nvcc release | `bin-kokkos-gpu/turbogap` |

With no GPU on the machine the last two are out of reach, and the host build
is what you test -- which is also all that CI can do; see below.

### Testing a CPU build

Once per clone: the test data, and the frozen baseline binary that the
regression suite compares against.

```sh
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
make -j4
tests/fetch_test_data.sh
tests/regression/make_baseline.sh
```

Then the two kinds of test. The regression suite asks whether any output
moved, byte for byte, against that baseline:

```sh
tests/regression/run.sh --list                       # the cases
tests/regression/run.sh                              # all of them
TURBOGAP_KEEP=1 tests/regression/run.sh estat_gsf    # one, keeping its run directory
```

The physics suites ask whether the answer is right, each against something
that is not TurboGAP -- an independent implementation, an analytic result, or
a finite difference of the code's own energy:

```sh
for suite in tests/*/run.sh; do
   case $suite in tests/gpu_zero_trunc/*|tests/regression/*) continue ;; esac
   echo "== $suite"
   "$suite" || echo "FAILED: $suite"
done
```

`tests/gpu_zero_trunc` is the device's own and is left out there; the rest run
on the host build. All of them honour:

| variable | |
| --- | --- |
| `TURBOGAP_BIN` | the binary under test (default `bin/turbogap`), so the same suite can be pointed at any of the four |
| `TURBOGAP_PYTHON` | the interpreter for the references; it needs `numpy`, and `tests/ir_fft` also wants `scipy` and says so loudly when it is absent |
| `TURBOGAP_DATA_ROOT` | where the test systems are |
| `TURBOGAP_KEEP` | keep the staging directory to look at afterwards |

```sh
TURBOGAP_BIN=$(realpath bin-gpu/turbogap) tests/ipi_pimd/run.sh
```

Read the exit status rather than the last lines: piping a suite through `tail`
replaces its status with `tail`'s, and every run then looks like a success.

### Testing a GPU build

Needs a CUDA (or HIP) toolchain and a device to run on. `HOP_ROOT` needs no
setting -- `src/hop` is a submodule. The same case list runs against the device
binary; build both, then:

```sh
tests/regression/run.sh --gpu      # every case, on bin-gpu/turbogap
tests/regression/run.sh --both     # the host pass, then the device pass
tests/gpu_zero_trunc/run.sh        # the device's own suites
tests/gpu/run_regression.sh        # whole trajectories, rather than the decks
```

`--gpu` compares the device binary against the **host build of the same
source**, not against the frozen baseline: the baseline is a snapshot of an
older commit, and comparing against it would fold every intended change since
into the same number as the host/device difference.

The comparison is numerical rather than bit-exact, through
`tests/regression/compare_tol.py`. The device sums the SOAP batches and the
cuBLAS reductions in a different order, so the last digits move on a run that is
entirely correct; a bit-exact diff would call every device run a failure.
Non-numeric tokens -- species columns, `Properties=` strings -- must still match
exactly, so a structural change is still caught. The default tolerance is
`rtol = atol = 1e-6`, overridable per case with `GPU_RTOL` and `GPU_ATOL` in its
`case.conf`, or globally with `TURBOGAP_GPU_RTOL` and `TURBOGAP_GPU_ATOL`.

### Testing a Kokkos build

Kokkos comes from `tools/install_kokkos.sh` (or `tools/setup_dev_env.sh
--with-kokkos`) and is found through `KOKKOS_ROOT`. It must be built by the
same CUDA release that compiles TurboGAP -- another release links and then
fails on a runtime symbol -- and the Makefile refuses a mismatch it can see.

Not with a plain diff against the CUDA binary. The device binary is not
reproducible run to run -- run `estat_gsf` twice with one unchanged binary and
`energy_soap` moves in the tenth digit -- so a straight comparison would charge
the Kokkos backend for differences that were there without it.

```sh
tools/verify_kokkos.sh             # every case; add case names to narrow it
```

That runs the suite twice: the CUDA binary against a copy of itself, which
measures how much the device moves on its own, and the Kokkos binary against
the CUDA one. What it reports is the difference of the two failure sets, so a
case listed at the end differs *because of* the Kokkos backend and not because
of the device.

### What CI covers, and what it leaves to you

GitHub's hosted runners have no GPU and no CUDA, so they test the host portion
and check the rest only as far as a dry run reaches. The device and Kokkos
binaries are yours to build and test on a machine that has them, before the
pull request.

- `build`: the host and serial builds, the generated dependency and keyword
  files being current, and the device and Kokkos builds *resolving* --
  `make -n` against a stub prefix, which catches the wiring rotting while
  nobody has a GPU to notice, without compiling anything.
- `tests`: the physics suites, on the host build.
- `regression`, on a pull request: the case list against the merge base, built
  in the same job. Red there means this branch changes some output, which a
  deliberate fix to the physics does too -- say so in the pull request rather
  than making it green.
- `gpu-regression`: the device suites, on a self-hosted runner labelled `gpu`.
  Without such a runner the job stays skipped.

## Developing TurboGAP

To develop for TurboGAP, one can install the development tools python environment (which has fprettify and pre-commit and so on which is installed through uv) such that formatting is preserved. Make a new branch or fork and then once ready submit a pull request. Please add tests for your new feature in the ` tests/` folder such that the CI interface can test upon pushing.

If you're including new source files, make sure to run the script which
regenerates the `makefiles/Makefile.deps`.

```sh
python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps
```

This may require you to first install a python environment if you don't have a
reasonable version of python. This can be done using the `source
tools/setup_dev_env.sh` which will install one for you using `uv` which will be
placed in `$HOME/.venvs`

To add new keywords in to TurboGAP, please follow the reference for adding them
in `docs/keywords-howto.md`. Following this format allows for automated
documentation and help information for usage.

### Debugging

Debug flags can be enabled by exporting

```sh
DEBUG=1
```

This will create a `bin-dbg` folder which contains the (now slower) binary.

## Attribution

When using **TurboGAP**, you should give attribution to the
**TurboGAP** author(s). The appropriate way to do that is to provide a link to the
[**TurboGAP** website](https://www.turbogap.fi) and, if you publish results obtained
using **TurboGAP** or the **soap_turbo** library,
even if it is through one of its external interfaces, you should cite:

> **Miguel A. Caro**. _Optimizing many-body atomic descriptors for enhanced computational
> performance of machine learning based interatomic potentials_. [Phys. Rev. B 100, 024112
> (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.024112).

In addition, you should cite any other relevant literature and code websites (e.g., the
original SOAP/GAP papers) as appropriate.
