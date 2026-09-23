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

## Overview of the code

The **TurboGAP** code is a specialized sampling/prediction engine which is made
for use with Gaussian Approximation Potentials that use many-body `soap_turbo`
descriptors.

Unique features are included in the code (e.g., Molecular Augmented Dynamics and experimental spectra prediction and
full support for van der Waals corrections).
Feedback can be sent to Tigany Zarrouk (tiganyzarrouk@gmail.com) or Miguel Caro (mcaroba@gmail.com)
or left on the Issues section of the Github page.

## Installation Quick Start

More detailed information on installation is found further on.

_**tl;dr, CPU**_ (standard linux, prerequisites installed):

```sh
git clone --recursive https://github.com/mcaroba/turbogap.git
cd turbogap
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
./compile_cpu.sh
export PATH="$(realpath bin):$PATH"
```

_**tl;dr, GPU**_ (CUDA or ROCm, on top of the above):

```sh
git clone --recursive https://github.com/mcaroba/turbogap.git
cd turbogap
export TURBOGAP_ARCH=Aalto_gfortran_openblas_hip_cuda
./compile_gpu.sh
export PATH="$(realpath bin-gpu):$PATH"
```

_**tl;dr, GPU through Kokkos**_ (optional; needed to run TurboGAP from LAMMPS):

```sh
# If you haven't installed Kokkos you can use the script in turbogap/tools
tools/install_kokkos.sh              # builds Kokkos, ~15 min, once
export KOKKOS_ROOT=$HOME/.local/kokkos-4.7.04 # Change this to where your kokkos installation is
KOKKOS=1 ./compile_gpu.sh
export PATH="$(realpath bin-kokkos-gpu):$PATH"
```

_**tl;dr, developer tooling**_ (only needed to commit, or to run the python
test drivers):

```sh
tools/setup_dev_env.sh
export PATH="$HOME/.venvs/turbogap-tools/bin:$PATH"
tools/setup_dev_env.sh --check
```

### Prerequisites

To use TurboGAP for CPU architectures, you must have the libraries `openmpi`, `lapack` + `blas` (or `openblas`)
and the Fortran compilers necessary for your system (e.g. `gcc`).

For GPU installations, you can use the CUDA/ROCm which is done via the `hop`
library that is included, or via Kokkos which is used for LAMMPS integration.

_For more detailed info on download, installation, etc., you can visit the
[TurboGAP wiki](https://turbogap.fi/wiki/index.php/Installation)._

### Getting the code

To get the **TurboGAP** code and the necessary **soap_turbo** routines, do a recursive
`git clone` with `--depth=1` to speed up things:

```sh
git clone --recursive --depth=1 http://github.com/mcaroba/turbogap.git
```

### What gets installed

A `--recursive` clone brings three submodules:

| submodule                                                 | what it is                                                          |
| --------------------------------------------------------- | ------------------------------------------------------------------- |
| `src/soap_turbo`                                          | CPU SOAP turbo routines uses                                        |
| `src/soap_turbo_gpu`                                      | GPU SOAP routines                                                   |
| `src/hop`                                                 | Header Only Porting for interoperability on CUDA/ROCm architectures |
| `tools/setup_dev_env.sh` creates a virtual environment at |
| `~/.venvs/turbogap-tools` holding, all pinned:            |

| tool                        | why                                                                                                                                 |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `fprettify`, `clang-format` | the formatting the pre-commit hook applies. Pinned because an unpinned formatter rewrites files nobody touched                      |
| `pre-commit`                | the hook itself: formatting, the generated keyword reference, whitespace                                                            |
| `numpy`                     | what the python test drivers compare against. Without it `tests/ir_fft`, `tests/mad_ir` and `tests/ipi_socket` skip rather than run |
| `i-PI`                      | drives TurboGAP over a socket in `turbogap ipi` mode, for path-integral and other advanced dynamics.                                |
| `cmake`                     | only installed by `tools/install_kokkos.sh`, and only when the machine has none of its own. Kokkos has no other build system        |

### Building Manually

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
predict`) or molecular dynamics (`turbogap md`) or (`turbogap mc`) or (`turbogap ipi`) for
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

The data for tests can be by cloning the [TurboGAP Tests](https://github.com/TiganyZ/turbogap_tests). This can be done using the script in `tests/regression/fetch_test_data.sh`.

Instead of doing a shallow clone (with depth=1) you can do the full clone so you have access to the baseline commit for regression tests.

```sh
git clone --recursive http://github.com/mcaroba/turbogap.git
cd turbogap
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
make -j4
turbogap_dir=$(realpath bin)
export PATH="$turbogap_dir:$PATH"
cd tests/regression
./make_baseline.sh
./fetch_test_data.sh
TURBOGAP_KEEP=1 ./run.sh
```

and the tests will be found in `$TMPDIR/turbogap_regression.xxxxx`.

### Testing a GPU build

The same case list runs against the device binary. Build both, then:

```sh
tests/regression/run.sh --gpu      # every case, on bin-gpu/turbogap
tests/regression/run.sh --both     # the host pass, then the device pass
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

```sh
tools/verify_kokkos.sh             # every case; add case names to narrow it
```

Other tests can be done by running the scripts in the `tests/<test_name>/run.sh` directories respectively.
Each honours `TURBOGAP_BIN`, so the same suite can be pointed at either build:

```sh
TURBOGAP_BIN=$(realpath bin-gpu/turbogap) tests/ipi_pimd/run.sh
```

## Developing TurboGAP

Follow the details of how the code is structured in `docs/DEVELOPMENT.md`. Install the developer tooling as stated above.

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
