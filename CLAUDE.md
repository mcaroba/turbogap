# TurboGAP Codebase Guide for AI Assistants

**Last Updated:** 2025-11-18
**Repository:** TurboGAP - Machine Learning Interatomic Potentials
**Primary Language:** Fortran 90 (18,221 lines), Python utilities (14 scripts)
**License:** Academic Software License (ASL) - Non-commercial use only

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Codebase Architecture](#codebase-architecture)
3. [Key Conventions](#key-conventions)
4. [Development Workflow](#development-workflow)
5. [Build System](#build-system)
6. [Testing Infrastructure](#testing-infrastructure)
7. [Common Tasks](#common-tasks)
8. [Important Notes for AI Assistants](#important-notes-for-ai-assistants)

---

## Project Overview

**TurboGAP** is a scientific computing package for machine learning-based interatomic potentials, developed primarily by Dr. Miguel A. Caro and contributors at Aalto University. It provides:

- **SOAP Descriptors**: Smooth Overlap of Atomic Positions for many-body atomic environments
- **GAP Framework**: Gaussian Approximation Potentials for molecular modeling
- **MD/MC Simulations**: Molecular dynamics and Monte Carlo sampling
- **Experimental Observables**: Structure factors, XRD, XPS, pair distribution functions
- **Radiation Damage**: Electronic stopping power for radiation cascade simulations

### Key References

When working with this code, be aware of these foundational papers:

- **Main Citation**: Miguel A. Caro, *Phys. Rev. B* 100, 024112 (2019)
- **SOAP Descriptors**: Bartók et al., *Phys. Rev. B* 87, 184115 (2013)
- **GAP Framework**: Bartók et al., *Phys. Rev. Lett.* 104, 136403 (2010)

### Website and Documentation

- **Website**: http://www.turbogap.fi
- **Wiki**: https://turbogap.fi/wiki/
- **GitHub**: http://github.com/mcaroba/turbogap

---

## Codebase Architecture

### Directory Structure

```
TurboReax/
├── src/                          # Fortran source code
│   ├── turbogap.f90              # Main program (4,322 lines)
│   ├── read_files.f90            # I/O handling (3,032 lines)
│   ├── exp_utils.f90             # Experimental properties (2,982 lines)
│   ├── neighbors.f90             # Neighbor lists, PBC (719 lines)
│   ├── gap.f90                   # GAP energy/forces (826 lines)
│   ├── gap_interface.f90         # GAP file reading (784 lines)
│   ├── exp_interface.f90         # Observable gradients (1,865 lines)
│   ├── md.f90                    # Molecular dynamics (800 lines)
│   ├── mc.f90                    # Monte Carlo (855 lines)
│   ├── vdw.f90                   # Van der Waals (495 lines)
│   ├── xyz.f90                   # XYZ I/O (477 lines)
│   ├── types.f90                 # Data structures (418 lines)
│   ├── local_properties.f90      # Hirshfeld volumes (286 lines)
│   ├── mpi.f90                   # MPI helpers (197 lines)
│   ├── splines.f90               # Spline interpolation (163 lines)
│   ├── stopping/                 # Electronic stopping (5 modules, 2,086 lines)
│   │   ├── electronic_stopping.f90      # SRIM stopping tables
│   │   ├── adaptive_time.f90            # Adaptive timestep
│   │   ├── eph_beta.f90                 # EPH beta parameters
│   │   ├── eph_fdm.f90                  # Finite difference method
│   │   └── eph_electronic_stopping.f90  # Temperature-dependent stopping
│   ├── soap_turbo/               # SOAP library (Git submodule)
│   │   └── src/                  # https://github.com/libAtoms/soap_turbo.git
│   └── third_party/
│       └── bussi_thermostat/     # Bussi velocity rescaling
├── tools/                        # Python utilities
│   ├── add_tags/                 # Prepare training databases
│   ├── compress_indices/         # SOAP descriptor compression
│   └── quip_xml_to_gap/          # Convert QUIP XML to TurboGAP format
├── tutorials/                    # Complete example workflows
│   ├── PtAu_NPs_MDS/            # Pt-Au nanoparticle MD (16 MB)
│   ├── creating_a-COx/          # Amorphous carbon-oxygen (23 MB)
│   └── xps_optimization/         # XPS-guided structure optimization (25 MB)
├── tests/                        # Test cases with reference data
│   └── stopping/                 # 5 electronic stopping test cases (64 MB)
├── makefiles/                    # Platform-specific build configs
├── docs/                         # Documentation files
├── Makefile                      # Main build file
├── README.md                     # Project overview
└── LICENSE.md                    # Academic Software License
```

### Module Dependency Graph

```
turbogap.f90 (main program)
├── types.f90                    # Core data structures
├── neighbors.f90                # Neighbor lists, distances
│   └── types.f90
├── gap.f90                      # GAP energy/force calculation
│   └── splines.f90             # Potential interpolation
├── gap_interface.f90            # Read GAP files
├── soap_turbo_* (submodule)    # SOAP descriptor computation
├── md.f90                       # Molecular dynamics integrator
├── mc.f90                       # Monte Carlo sampler
├── vdw.f90                      # Van der Waals corrections
├── exp_utils.f90                # Experimental observables
│   ├── types.f90
│   ├── splines.f90
│   └── xyz_module.f90
├── exp_interface.f90            # Observable gradients
├── local_properties.f90         # Hirshfeld volumes, charges
├── read_files.f90               # All input/output
│   ├── neighbors.f90
│   ├── vdw.f90
│   ├── xyz_module.f90
│   └── soap_turbo_compress_module
├── stopping/* (5 modules)       # Electronic stopping power
├── mpi.f90                      # Parallel computing
└── bussi.f90                    # Thermostat
```

### Key Data Structures (types.f90)

```fortran
type input_parameters              ! Main simulation configuration
  - do_md, do_mc                  ! Simulation type flags
  - n_steps, time_step            ! Temporal parameters
  - T_target, P_target            ! Thermodynamic targets
  - exp_data(:)                   ! Array of experimental observables
  - xps, pair_distribution, structure_factor, xrd_params
  - soap_turbo_hypers             ! SOAP descriptor settings

type soap_turbo                    ! SOAP descriptor parameters
  - l_max, n_max                  ! Angular/radial basis size
  - rcut_hard, rcut_soft          ! Cutoff radii
  - atom_sigma_r, atom_sigma_t    ! Smearing widths
  - alphas(:), alphas_compress(:) ! GAP coefficients
  - basis(:)                      ! Descriptor basis

type exp_data_container            ! Experimental observable data
  - structure_factor, pair_distribution, xrd, xps, neutron_diffraction

type image                         ! Atomic configuration snapshot
  - n_atoms                       ! Number of atoms
  - positions(:,:)                ! Cartesian coordinates
  - velocities(:,:)               ! Atomic velocities
  - forces(:,:)                   ! Forces on atoms
  - species(:)                    ! Atomic types
  - lattice(1:3,1:3)             ! Periodic cell vectors
```

---

## Key Conventions

### 1. File Header Format

All Fortran source files begin with a standard copyright header:

```fortran
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, <filename>, is copyright (c) <years>, <authors>
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

**When creating new files**: Always include this header with appropriate author/year information.

### 2. Fortran Coding Style

- **Language Standard**: Fortran 90
- **Free-form source**: `.f90` extension
- **Indentation**: 2 spaces per level
- **Line length**: Generally <132 characters
- **Variable declarations**: Explicit `implicit none` in all modules/programs
- **Real precision**: `real*8` for double precision (8 bytes)
- **Array indexing**: 1-based (Fortran convention)
- **Intent specifications**: Always declare intent (in/out/inout) for subroutine arguments

Example:
```fortran
subroutine calculate_energy(positions, n_atoms, species, energy, forces)
  implicit none
  real*8, intent(in) :: positions(:,:)
  integer, intent(in) :: n_atoms
  integer, intent(in) :: species(:)
  real*8, intent(out) :: energy
  real*8, intent(out) :: forces(:,:)
  ! ... implementation
end subroutine
```

### 3. Naming Conventions

- **Modules**: Lowercase, descriptive (e.g., `neighbors`, `exp_utils`)
- **Subroutines**: Lowercase with underscores (e.g., `get_soap_energy_and_forces`)
- **Variables**: Lowercase with underscores (e.g., `n_atoms`, `rcut_max`)
- **Type names**: Lowercase with underscores (e.g., `input_parameters`, `soap_turbo`)
- **Constants**: Lowercase with underscores (e.g., `kB = 8.6173303d-5`)

### 4. Comment Style

- **Intent comments**: Describe what subroutine does, not how
- **Algorithm notes**: Use `!` for inline comments
- **TODO markers**: Use `! TODO:` or `! WE NEED TO...`
- **Section separators**: Use multiple `!` for visual separation

Example from code:
```fortran
! This subroutine reads in the XYZ file
! WE NEED TO WRITE A PROPER EXTXYZ READER...
subroutine read_xyz(filename, ase_format, all_atoms, ...)
```

### 5. Input File Format

Input files use a keyword-based format with case-insensitive matching:

```
# Comment lines start with #

# Simulation type
do_md
n_steps = 10000
time_step = 0.5

# Temperature control
T_target = 300.0

# Structure file
structure_file = silicon.xyz

# GAP model
gap_file = model.gap
```

**Important**: Keywords are matched with `index(line, 'keyword')`, so partial matches may occur. Be specific.

### 6. File Formats

- **Atomic structures**: Extended XYZ format (ASE-compatible)
- **GAP models**: Binary `.gap` files (use `tools/quip_xml_to_gap/` to convert from QUIP XML)
- **Stopping tables**: Text `.dat` files (SRIM format)
- **Trajectories**: XYZ format with velocities/forces
- **Logs**: `thermo.log` for thermodynamic data, `ElectronicEnergyLoss.txt` for stopping

---

## Development Workflow

### Git Workflow

**Current Repository State:**
- **Main branch**: `master` (not explicitly shown, but standard)
- **Development branch**: `claude/claude-md-<session-id>`
- **Submodule**: `src/soap_turbo` tracking `https://github.com/libAtoms/soap_turbo.git`

**Working with Git:**

1. **Always develop on the designated Claude branch**:
   ```bash
   git checkout claude/claude-md-<session-id>
   ```

2. **Initialize submodules** (required for building):
   ```bash
   git submodule update --init --recursive
   ```

3. **Commit messages** should be descriptive:
   ```bash
   git commit -m "Fix neighbor list bug for triclinic cells"
   ```

4. **Push to remote** (use `-u` for new branches):
   ```bash
   git push -u origin claude/claude-md-<session-id>
   ```

5. **Retry logic for network failures**: If push/pull fails, retry up to 4 times with exponential backoff (2s, 4s, 8s, 16s)

**Recent Development Focus** (from git log):
- SOAP descriptor improvements (poly3operator radial basis)
- XPS functionality enhancements
- Bug fixes for parallel execution
- Radiation cascade simulation improvements

### Making Changes

#### When Modifying Fortran Code:

1. **Understand dependencies**: Check module `use` statements
2. **Preserve interfaces**: Don't break existing subroutine signatures unless necessary
3. **Test compilation**: Always compile after changes
4. **Check for side effects**: Fortran modules have shared state

#### When Modifying Python Tools:

1. **Dependencies**: Tools require `ase` (Atomic Simulation Environment), `numpy`, `BeautifulSoup4`
2. **Input/output**: Tools read/write ASE-compatible XYZ files
3. **Test on small data**: Use tutorial data for validation

---

## Build System

### Compilation Process

**Standard Build:**
```bash
# Clean previous build
make clean

# Build library and executables
make

# Deep clean (removes all build directories)
make deepclean
```

**Output Locations:**
- Executable: `bin/turbogap`
- Library: `lib/libturbogap.a`
- Modules: `include/*.mod`
- Objects: `build/*.o`

### Platform-Specific Configuration

The build uses platform-specific makefiles in `makefiles/`. Current default:
```makefile
include makefiles/Makefile.Ubuntu_gfortran_mpi
```

**Available Configurations:**
1. `Makefile.Ubuntu_gfortran` - Single-processor, gfortran
2. `Makefile.Ubuntu_gfortran_mpi` - MPI parallel, gfortran (DEFAULT)
3. `Makefile.CSC-Puhti_gfortran_mkl_mpi` - CSC Puhti HPC, Intel MKL
4. `Makefile.CSC-Mahti_gfortran_openblas_mpi` - CSC Mahti HPC, OpenBLAS

**Key Makefile Variables:**
```makefile
F90 = mpif90                    # Fortran compiler
PP = -cpp -D _MPIF90           # Preprocessor flags
F90_OPTS = -fPIC -O3           # Optimization flags
LIBS = -llapack -lblas         # Required libraries
```

### Dependencies

**Required:**
- Fortran 90 compiler (gfortran, ifort)
- BLAS (Basic Linear Algebra Subprograms)
- LAPACK (Linear Algebra Package)

**Optional:**
- MPI library (mpich, openmpi) - for parallel execution
- Python 3 with ASE - for utility scripts

**Compilation Order** (managed by Makefile):
1. Third-party code (`bussi_thermostat`)
2. SOAP_turbo modules (from submodule)
3. Electronic stopping modules
4. Core TurboGAP modules
5. Main program

### Debugging Build

**Enable debugging**:
Edit the makefile to use debug flags:
```makefile
F90_OPTS=-fPIC -O0 -fcheck=bounds -g -fcheck=all -Wall -fbacktrace
```

**Common build errors**:
- **Missing soap_turbo**: Run `git submodule update --init --recursive`
- **LAPACK not found**: Install `liblapack-dev` or specify path in makefile
- **MPI errors**: Ensure `mpif90` wrapper is in PATH, or use non-MPI makefile

---

## Testing Infrastructure

### Test Organization

**Location**: `tests/stopping/` - 5 complete test cases for electronic stopping

Each test directory contains:
- `input` - Simulation configuration
- `silicon.xyz` - Initial atomic structure
- `*.gap` - GAP model files
- `*.beta`, `*.fdm`, `*.dat` - Stopping power data
- `readme` - Test description and expected outputs

### Running Tests

**Single test**:
```bash
cd tests/stopping/1/
../../bin/turbogap md
```

**Expected outputs**:
- `thermo.log` - Thermodynamic properties (time, temp, energy, pressure)
- `ElectronicEnergyLoss.txt` - Electronic energy loss over time
- `trajectory_out.xyz` - Atomic trajectory
- Console output with timing information

**Validation**:
Compare outputs with reference data mentioned in `readme` files. Key metrics:
- Energy conservation (for NVE)
- Temperature fluctuations (for NVT)
- Electronic stopping magnitude
- Simulation performance (time per step)

### Test Cases Overview

1. **Test 1**: Binary collision with adaptive timestep
2. **Test 2**: Temperature-dependent electronic stopping
3. **Test 3**: Constant beta stopping model
4. **Test 4**: Quadratic beta model
5. **Test 5**: Finite difference method implementation

---

## Common Tasks

### 1. Running a Molecular Dynamics Simulation

```bash
# Prepare input file (see examples in tests/ or tutorials/)
cat > input << EOF
do_md
n_steps = 1000
time_step = 1.0
T_target = 300.0
structure_file = my_structure.xyz
gap_file = my_model.gap
EOF

# Run simulation
/path/to/bin/turbogap md
```

### 2. Running a Monte Carlo Simulation

```bash
# Input file for MC
cat > input << EOF
do_mc
n_steps = 10000
T_target = 300.0
P_target = 0.0
structure_file = initial.xyz
gap_file = model.gap
EOF

# Run
turbogap mc
```

### 3. Converting QUIP XML to TurboGAP Format

```bash
cd tools/quip_xml_to_gap/
python make_gap_files.py \
  --xml model.xml \
  --output_dir ../../my_model/
```

This creates `.gap` files readable by TurboGAP.

### 4. Preparing Training Data

```bash
cd tools/add_tags/
python add_tags.py \
  --input training_data.xyz \
  --output tagged_data.xyz \
  --energy_sigma 0.001 \
  --force_sigma 0.01
```

Adds regularization tags for GAP fitting.

### 5. Adding a New Observable Type

**Steps**:

1. **Define data structure** in `src/types.f90`:
   ```fortran
   type my_new_observable
     real*8 :: parameter1
     real*8, allocatable :: data(:)
   end type
   ```

2. **Add reading routine** in `src/read_files.f90`:
   ```fortran
   subroutine read_my_observable(file, obs)
     ! Read from file
   end subroutine
   ```

3. **Implement calculation** in `src/exp_utils.f90`:
   ```fortran
   subroutine get_my_observable(positions, species, obs, result)
     ! Calculate observable
   end subroutine
   ```

4. **Add gradient calculation** in `src/exp_interface.f90`:
   ```fortran
   subroutine get_my_observable_forces(positions, species, obs, forces, virial)
     ! Calculate derivatives for force matching
   end subroutine
   ```

5. **Integrate into main program** in `src/turbogap.f90`:
   ```fortran
   if (params%do_my_observable) then
     call get_my_observable(...)
   end if
   ```

6. **Update Makefile** if new files created:
   ```makefile
   SRC := ... my_observable.f90
   ```

### 6. Analyzing Simulation Output

**Python analysis** (using ASE):
```python
from ase.io import read

# Read trajectory
atoms_list = read('trajectory_out.xyz', index=':')

# Analyze
for atoms in atoms_list:
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    # ... analyze
```

**Extract thermodynamic data**:
```bash
# From thermo.log
awk '{print $1, $3}' thermo.log > energy_vs_time.dat
```

---

## Important Notes for AI Assistants

### 1. Code Maturity and Stability

- **Core MD/MC**: Stable and well-tested
- **GAP calculations**: Production-ready
- **SOAP descriptors**: Well-established (via soap_turbo submodule)
- **Experimental observables**: Active development area
- **Native interface**: Described as "temperamental" in README
- **Documentation**: Limited inline documentation; rely on wiki and papers

### 2. Common Pitfalls

#### Fortran-Specific:

- **Array bounds**: Fortran uses 1-based indexing, not 0-based
- **Column-major ordering**: Arrays are stored column-major (opposite of C/Python)
- **Intent violations**: Don't modify `intent(in)` parameters
- **Module dependencies**: Order matters in compilation
- **Implicit typing**: Always use `implicit none`

#### TurboGAP-Specific:

- **Submodule required**: Code won't compile without `soap_turbo` submodule
- **Input file parsing**: Uses simple string matching, can have false positives
- **File paths**: Many routines assume files are in current directory
- **MPI configuration**: MPI code is conditionally compiled with `-D _MPIF90`
- **Periodic boundaries**: Always active; non-periodic systems need large cell

### 3. Performance Considerations

- **Neighbor lists**: Rebuilt every step; major bottleneck for large systems
- **SOAP calculation**: Most expensive per-atom operation
- **MPI scaling**: Good for large systems (>1000 atoms), overhead for small
- **Memory**: SOAP descriptors can require significant memory for large basis sets
- **I/O**: Trajectory writing can be slow for large systems

### 4. Scientific Accuracy Requirements

When modifying force/energy calculations:

- **Energy conservation**: NVE simulations should conserve energy to ~1e-6 eV/atom
- **Force consistency**: Forces must be negative gradient of energy
- **Virial correctness**: Pressure calculation requires correct virial implementation
- **Symmetry preservation**: Don't break rotational/translational invariance
- **Units**: TurboGAP uses eV, Angstrom, femtoseconds internally

### 5. Interaction with External Tools

**QUIP Integration**: TurboGAP can use GAP models trained with QUIP, but requires conversion via `tools/quip_xml_to_gap/`

**ASE Integration**: Trajectory files are ASE-compatible XYZ format

**LAMMPS**: TurboGAP is NOT a LAMMPS interface; it's a standalone code

**SRIM**: Electronic stopping tables can be imported from SRIM

### 6. License Compliance

**Academic Software License (ASL)**:
- Non-commercial use only
- Free for academic research and teaching
- Commercial use requires separate license
- Proper attribution required (cite Phys. Rev. B paper)
- Third-party code in `src/third_party/` has separate compatible licenses

**When adding new code**:
- Include proper copyright header
- Ensure any external code is ASL-compatible
- Document sources and licenses

### 7. Debugging Strategies

**Compilation errors**:
```bash
# Use debug flags
make clean
# Edit makefile to enable debugging
make 2>&1 | tee build.log
```

**Runtime errors**:
- Check `implicit none` is used
- Verify array bounds with `-fcheck=bounds`
- Use `-fbacktrace` to get stack traces
- Check MPI configuration if parallel run fails

**Scientific errors**:
- Verify energy conservation for NVE
- Check forces are smooth (no discontinuities)
- Validate against known test cases
- Compare with QUIP results if available

### 8. Documentation and Help

**When stuck**:
1. Check the [TurboGAP Wiki](https://turbogap.fi/wiki/)
2. Look at tutorial examples in `tutorials/`
3. Search for similar functionality in existing code
4. Contact: Miguel Caro (mcaroba@gmail.com)
5. GitHub Issues: https://github.com/mcaroba/turbogap/issues

**Citation template for reference**:
```
@article{Caro2019,
  author = {Caro, Miguel A.},
  title = {Optimizing many-body atomic descriptors for enhanced computational
           performance of machine learning based interatomic potentials},
  journal = {Phys. Rev. B},
  volume = {100},
  pages = {024112},
  year = {2019},
  doi = {10.1103/PhysRevB.100.024112}
}
```

### 9. Code Review Checklist

Before committing changes, verify:

- [ ] Compilation succeeds with no warnings
- [ ] At least one test case runs successfully
- [ ] Energy conservation checked (if modifying forces/energy)
- [ ] Copyright header present in new files
- [ ] Variable intent declared for all subroutine arguments
- [ ] No breaking changes to existing interfaces (or documented if necessary)
- [ ] Memory allocated arrays are deallocated
- [ ] MPI-related code wrapped in `#ifdef _MPIF90`
- [ ] Comments explain non-obvious algorithms
- [ ] Commit message is descriptive

### 10. Quick Reference Commands

```bash
# Build from scratch
git submodule update --init --recursive
make clean
make

# Run test
cd tests/stopping/1/
../../bin/turbogap md

# Convert QUIP model
cd tools/quip_xml_to_gap/
python make_gap_files.py --xml model.xml --output_dir output/

# Prepare training data
cd tools/add_tags/
python add_tags.py --input data.xyz --output tagged.xyz

# Check code statistics
find src/ -name "*.f90" -exec wc -l {} + | sort -n

# Search for function/subroutine
grep -r "subroutine get_soap" src/
```

---

## Appendix: File Size Reference

**Source Code Statistics** (as of 2025-11-18):

| Category | Files | Lines |
|----------|-------|-------|
| Main modules | 14 | ~12,000 |
| Main program | 1 | 4,322 |
| Stopping modules | 5 | 2,086 |
| SOAP_turbo (submodule) | 5 | ~8,000 |
| Third-party | 1 | ~200 |
| **Total Fortran** | **26** | **~18,221** |
| Python utilities | 14 | ~1,500 |

**Test Data**: ~64 MB (5 test cases)
**Tutorial Data**: ~64 MB (3 complete workflows)

---

## Version History

- **2025-11-18**: Initial CLAUDE.md created based on codebase analysis
- Recent commits focus on: SOAP poly3operator, XPS functionality, parallel execution fixes

---

**End of CLAUDE.md**

For the most current information, always check:
- Website: http://www.turbogap.fi
- Wiki: https://turbogap.fi/wiki/
- GitHub: http://github.com/mcaroba/turbogap
