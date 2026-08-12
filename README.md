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

_**tl;dr**_ (for standard linux builds with the prerequisites installed):

```sh
git clone --recursive --depth=1 http://github.com/mcaroba/turbogap.git
cd turbogap
export TURBOGAP_ARCH=Ubuntu_gfortran_mpi
make -j4
turbogap_dir=$(realpath bin)
export PATH="$turbogap_dir:$PATH"
echo "export PATH=\"${turbogap_dir}:\$PATH\"" >> ~/.bashrc
```

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

If you need help with TurboGAP modes or its keywords, you can either consult the html document

```sh
open docs/keywords.html
```

or the markdown `docs/keywords.md`, or one can use the `--help`:

:

```sh
turbogap --help
```

or

```sh
turbogap --help <mode>
```

where mode is any of the above (`predict`, `md`, `mc`).

The keywords for the `input` file will be listed with descriptions and their
dependencies and be filtered according to the simulation mode.

## TurboGAP Tutorials

For the latest features in TurboGAP, one can follow the tutorials which were
ran for the TurboGAP School, and adapt them for their machine: [TurboGAP
School](https://github.com/mcaroba/TurboGAP_School).

Simple tutorials can be found in the [TurboGAP Tutorials](https://github.com/TiganyZ/turbogap_tutorials) repository `git clone https://github.com/TiganyZ/turbogap_tutorials.git`. Further tutorials and all documentation can be found on the [Turbogap Website](https://turbogap.fi/wiki/index.php/Tutorials).

## Testing TurboGAP

The data for tests can be by cloning the [TurboGAP Tests](https://github.com/TiganyZ/turbogap_tests). One can clone the repository `git clone https://github.com/TiganyZ/turbogap_tests.git` in the directory where the TurboGAP repo is stored and execute the regression tests.

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
