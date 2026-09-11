# Copyright (c) 2020-2023 by Albert Bartók and Miguel Caro

SHELL = /bin/sh

# Include user-modifiable variables from a customizable file.
# Check the makefiles/ directory for a list of tested architectures

TURBOGAP_ARCH ?= Ubuntu_gfortran_mpi

# Where the HOP headers live. hop is a submodule, so a recursive clone needs
# no HOP_ROOT in the environment; setting one still overrides this.
HOP_ROOT ?= $(CURDIR)/src/hop
export HOP_ROOT

include makefiles/Makefile.$(TURBOGAP_ARCH)

# Which soap_turbo checkout the SOAP routines come from. The CPU and GPU lines
# are separate submodules: their get_soap signatures differ (the GPU one takes
# device pointers, a cuBLAS handle and a stream) and so do their module
# dependencies. A GPU arch makefile sets this to src/soap_turbo_gpu/src.
ST_DIR ?= src/soap_turbo/src

# The module graph differs with ST_DIR, so the generated deps do too.
DEPS_FILE ?= makefiles/Makefile.deps

# ------------------------------------------------------------------- DEBUG
#
# A bounds-checked, symbol-carrying build, in its own object tree.
#
#     make            ->  build/       bin/turbogap
#     make DEBUG=1    ->  build-dbg/   bin-dbg/turbogap
#
# The separate tree is the point, not tidiness. make decides whether to
# rebuild from TIMESTAMPS, not from flags, so building DEBUG=1 over an
# existing build/ recompiles only what changed and links the rest of the
# optimised objects -- while reporting success. The reverse is worse: the
# optimised binary everything else names would then contain debug objects, and
# nothing in the build would say so. The GPU branch measured that: a stale
# DEBUG=1 object tree ran a kernel 4.7x slower than the same source rebuilt
# clean and shifted an energy by ~1e-10, which reads as a numerical regression
# from an unrelated change.
#
# DEBUG=0 is left untagged deliberately. bin/turbogap is what the regression
# suite, the extra suites and the docs all name.
#
# These flags are appended after the architecture's, so the last -O wins and
# anything arch-specific (-march, and so on) survives. -fcheck=all names the
# file and line of an out-of-range array access instead of leaving a glibc
# "free(): invalid next size" and a backtrace of addresses.
DEBUG ?= 0
ifeq ($(DEBUG),1)
  F90_OPTS += -O0 -g -fbacktrace -fcheck=all
  BUILD_TAG := $(BUILD_TAG)-dbg
endif

# ---------------------------------------------------------------- PROFILE=1
#
# A build the NVIDIA profilers can read. Kept here rather than in the
# architecture makefiles so every device architecture gets it from one place,
# and so it is applied AFTER the include -- these are additions to whatever
# optimisation level the architecture chose, never a replacement for it.
#
#   -lineinfo   maps device instructions back to source lines, so an ncu report
#               names the line that stalled rather than the SASS address. It
#               does not change codegen.
#   -g (host)   host symbols, so the nsys sampler shows function names.
#   -D_NVTX     compiles in the phase markers in src/nvtx.f90, driven from the
#               timing buckets. Without them a timeline is kernels with no
#               statement of which phase launched them.
#
# DEBUG=1 is the wrong thing to profile and the mistake is easy to make, so it
# is refused rather than warned about: -G disables essentially all device
# optimisation, and docs/BUILD_AND_TEST.md measures that at a 2.1x slowdown.
# Numbers from such a build do not rank the kernels in the order the real build
# ranks them, which is the only thing a profile is for.
PROFILE ?= 0
ifeq ($(PROFILE),1)
  ifeq ($(DEBUG),1)
    $(error PROFILE=1 needs DEBUG=0. -G disables device optimisation (~2.1x \
      slower) and reorders the kernel ranking, so a DEBUG profile does not \
      describe the build you run. Use: make PROFILE=1 DEBUG=0)
  endif
  F90_OPTS += -g -fno-omit-frame-pointer
  PP += -D_NVTX
  CU += -lineinfo -g
  CC += -lineinfo -g
  BUILD_TAG := $(BUILD_TAG)-profile
endif

# ----------------------------------------------------------------- OPENMP=1
#
# Host threading. No architecture makefile passes -fopenmp, so _OPENMP has been
# undefined on every build that ever shipped and the !$OMP directives were
# comments. The neighbour build and the per-step pair geometry in
# src/neighbors.f90 carry them now, and they are the reason to turn this on:
# with one MPI rank per GPU a single core built the whole list while the device
# idled, and a node has tens of cores per rank going spare.
#
# Off by default all the same. It is a build-time choice that changes which
# code runs, so it gets its own object tree and its own bin, and the regression
# suite keeps testing the untagged one.
#
# gfortran needs -fopenmp at BOTH compile and link: it selects the runtime
# library as well as enabling the directives, and the !$ sentinel lines are
# conditionally compiled only when it is present.
#
#     make OPENMP=1               build-omp/   bin-omp/turbogap
OPENMP ?= 0
ifeq ($(OPENMP),1)
  F90_OPTS += -fopenmp
  LIBS += -fopenmp
  BUILD_TAG := $(BUILD_TAG)-omp
endif

# ------------------------------------------------------- an ad-hoc variant
#
# Extra compiler flags, and a tag so the variant gets its own object tree.
# Both are needed together: without the tag, make sees unchanged sources and
# relinks the previous variant's objects, so the "comparison" is a binary
# against itself.
#
#     make F90_EXTRA='-ffpe-trap=invalid,zero,overflow' BUILD_TAG_EXTRA=-fpe
#     make F90_EXTRA='-finit-real=snan -finit-integer=-99999999' BUILD_TAG_EXTRA=-poison
#     make F90_EXTRA=-pg BUILD_TAG_EXTRA=-prof
#
# Use this rather than overriding F90_OPTS on the command line. A command-line
# assignment replaces the variable outright, including the `-J $(INC_DIR)`
# appended below, so gfortran writes every .mod into the working directory
# instead of the build's include dir. Those stray .mod files are then found
# ahead of the real ones by every later build in that tree, which fails with
# "is not a member of the structure" against a field that plainly is -- and
# deleting the build dir does not fix it, because the stale modules are not in
# it. That happened here on 2026-08-12 and cost half an hour.
CC += $(CC_EXTRA)
CU += $(CU_EXTRA)
F90_OPTS += $(F90_EXTRA)
BUILD_TAG := $(BUILD_TAG)$(BUILD_TAG_EXTRA)

# ---------------------------------------------------------------- the trees
#
# One object tree per flag combination. Every rule, plus clean and deepclean,
# is written in terms of these four variables, and the generated
# makefiles/Makefile.deps in terms of $(BUILD_DIR), so this is the only place
# the split has to be stated. Overriding any of them on the command line still
# works, for a one-off tree with no flag change behind it.
BUILD_DIR=build$(BUILD_TAG)
BIN_DIR=bin$(BUILD_TAG)
INC_DIR=include$(BUILD_TAG)
LIB_DIR=lib$(BUILD_TAG)

# Do not change anything below this line
##########################################################

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR)

PROGRAMS := turbogap



# ------------------------------------------------------------------- GPU=1
#
# A device build. Set by the device architecture makefiles, which also put
# -D _GPU in PP so the guarded Fortran compiles, and point ST_DIR at the GPU
# soap_turbo.
#
# Two modules have one name and two implementation files, chosen here rather
# than by #ifdef, because the two bodies share almost nothing: gpu_context owns
# the streams and cuBLAS handles on one side and nothing at all on the other,
# and gap_backend either uploads the neighbour data once for the three
# contribution calls or does nothing. The driver calls the same names on both.
GPU ?= 0
ifeq ($(GPU),1)
  GPU_CONTEXT := gpu_context_gpu.f90
  GAP_BACKEND := gap_backend_gpu.f90
# Device sources, all under src/gpu/ and named after what they compute. See
# src/gpu/README.md: gpu_* is infrastructure both sides use, gap_* is the
# interatomic potential and mad_* the experimental-data path.
  SRC_CUDA := gpu_memory.cu gpu_blas.cu gpu_scan.cu \
              gap_predict.cu gap_soap_radial.cu gap_soap_radial_operator.cu gap_soap_angular.cu \
              gap_soap_descriptor.cu gap_soap_forces.cu gap_2b.cu \
              mad_pdf.cu mad_xrd.cu mad_electrostatics.cu gpu_scatter.cu
# gap_3b stays .cc, and so is compiled by $(CC) rather than $(CU), because it
# needs -std=c++20 for <numbers> and <bit>; only $(CC) passes it.
#
# orthonormalization_kernels.cc is deliberately disabled, not dead: it is kept
# in src/gpu/ for possible reintegration. Re-enable by moving it before the '#'.
  SRC_CC := gap_3b.cc # orthonormalization_kernels.cc
else
  GPU_CONTEXT := gpu_context_cpu.f90
  GAP_BACKEND := gap_backend_cpu.f90
  SRC_CUDA :=
  SRC_CC :=
endif

SRC := printing.f90 error.f90 read_utils.f90 nvtx.f90 timing.f90 misc.f90 electrostatics.f90 constants.f90 gle.f90 ipi_socket.f90 ipi_driver.f90 mad_ir.f90 mad_ir_xl.f90 ir_fft.f90 ir_fft_io.f90 nonneg_leastsq.f90 splines.f90 types.f90 $(GPU_CONTEXT) neighbors.f90 neighbors_skin.f90 gap.f90 vdw.f90		\
	local_properties.f90 exp_utils.f90  xyz.f90 md.f90 ir_auxiliary_dynamics.f90 mc.f90 read_files.f90	\
	$(GAP_BACKEND) gap_interface.f90 mpi.f90 exp_interface.f90 turbogap_exp.f90 turbogap_md.f90 turbogap_vdw.f90 turbogap_estat.f90 turbogap_setup.f90

# kinds must build before everything, including SRC_STOP, so it gets its own
# group placed first in every prerequisite list.
SRC_BASE := kinds.f90 keyword_help.f90
SRC_TP_BT := resamplekin.f90
ifeq ($(GPU),1)
  SRC_TP_BT += fortran_cuda_interfaces.f90
endif
SRC_ST := soap_turbo_functions.f90 soap_turbo_radial.f90 soap_turbo_angular.f90 \
          soap_turbo.f90 soap_turbo_compress.f90
SRC_STOP := adaptive_time.f90 electronic_stopping.f90 eph_beta.f90 eph_fdm.f90 \
            eph_electronic_stopping.f90

OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC)))
OBJ_BASE := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_BASE)))
OBJ_TP_BT := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_TP_BT)))
OBJ_ST := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_ST)))
OBJ_STOP := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_STOP)))
OBJ_CUDA := $(addprefix $(BUILD_DIR)/,$(patsubst %.cu,%.o,$(SRC_CUDA)))
OBJ_CC := $(addprefix $(BUILD_DIR)/,$(patsubst %.cc,%.o,$(SRC_CC)))

PROG := $(addprefix $(BIN_DIR)/,$(PROGRAMS))

.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: default all programs clean deepclean libturbogap docs check-docs

default: libturbogap programs

all: default

# Regenerate the input-keyword reference from the !> blocks above each keyword
# in src/read_files.f90: docs/keywords.md, docs/keywords.html and
# src/keyword_help.f90, the module `turbogap --help` prints from. Run after
# adding, renaming or documenting a keyword.
docs:
	python3 tools/keyword_docs.py

# Fail if a keyword has no !> block or one of those three outputs is stale.
# Also runs from .pre-commit-config.yaml.
check-docs:
	python3 tools/keyword_docs.py --check

clean:
	rm -rf $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(OBJ_CUDA) $(OBJ_CC) $(INC_DIR)/*.mod $(PROG)

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

.SECONDEXPANSION:
.SECONDARY: $(OBJS)

programs: $(PROG)

libturbogap: $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(OBJ_CUDA) $(OBJ_CC) ${LIB_DIR}
	ar scr $(LIB_DIR)/libturbogap.a $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(OBJ_CUDA) $(OBJ_CC)

$(BIN_DIR)/%: src/%.f90 $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(OBJ_CUDA) $(OBJ_CC) | $$(@D)
	$(F90) $(PP) $(F90_OPTS) $< -o $@ $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(OBJ_CUDA) $(OBJ_CC) $(LIBS)

# Device sources. The headers are listed as prerequisites explicitly: they are
# not in the generated deps (which is Fortran modules only) and nvcc is not
# asked to emit depfiles, so without this an edit to gpu_common.h would rebuild
# nothing -- the same silent staleness the BUILD_TAG note above describes.
GPU_HEADERS := $(wildcard src/gpu/*.h)

$(BUILD_DIR)/%.o: src/gpu/%.cu $(GPU_HEADERS) | $$(@D)
	$(CU) $(CUDA_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/gpu/%.cc $(GPU_HEADERS) | $$(@D)
	$(CC) $(CC_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/stopping/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/third_party/bussi_thermostat/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/third_party/nnls/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: $(ST_DIR)/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR): ${INC_DIR}
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(INC_DIR):
	mkdir -p $@

$(LIB_DIR):
	mkdir -p $@

# Fortran module dependencies. Compiling a file that USEs a module requires that
# module's .mod file to exist already, and make cannot infer that from the source
# tree. Without this include the build works only because the SRC lists happen to
# be in a workable order for a *serial* build; -j races and fails with
# "Cannot open module file".
#
# Included last on purpose: an include placed before the first rule would make one
# of these dependency lines the default goal, and `make` would silently build a
# single object and exit 0.
#
# Regenerate after adding, removing or moving a USE / module. DEPS_FILE is
# per-arch because ST_DIR selects a soap_turbo whose module graph differs:
#     python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps
#     make TURBOGAP_ARCH=<gpu arch> ... then regenerate into Makefile.deps.gpu
include $(DEPS_FILE)
