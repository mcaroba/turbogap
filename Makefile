# Copyright (c) 2020-2023 by Albert Bartók and Miguel Caro

SHELL = /bin/sh

# Include user-modifiable variables from a customizable file.
# Check the makefiles/ directory for a list of tested architectures

ifndef TURBOGAP_ARCH
$(error TURBOGAP_ARCH is not defined. Please export it before running make to pick a makefile from makefiles/)
endif

ifndef HOP_ROOT
$(error HOP_ROOT is not defined. Please export it before running make!)
endif
include makefiles/Makefile.${TURBOGAP_ARCH}



# ----------------------------------------------------------------- OPENMP=1
#
# Host threading, and with it the per-thread GPU streams.
#
# No architecture makefile has ever passed -fopenmp, so _OPENMP has always been
# undefined: gpu_context_init's thread-count block was dead, n_omp came out 1,
# exactly one stream was created, and the !$OMP PARALLEL DO over the batched
# pdf and xrd loops was a comment. Everything the code says about "n_omp
# streams for the batched gpu calculation" described an intent, not a build --
# nsys measures peak device concurrency of 1 on every test case.
#
# Off by default, deliberately. Turning it on makes two loops in
# turbogap_exp.f90 run concurrently, and the only test case that reaches them
# (XRD_mad) is an xfail, so the suite cannot tell a regression there from the
# drift it already tolerates. Compare trajectory_out.xyz between an OPENMP=0
# and an OPENMP=1 build before trusting it on a new path.
#
# gfortran needs -fopenmp at BOTH compile and link: it selects the runtime
# library as well as enabling the directives, and the !$ sentinel lines are
# conditionally compiled only when it is present.
OPENMP ?= 0
ifeq ($(OPENMP),1)
  F90_OPTS += -fopenmp
  LIBS += -fopenmp
  BUILD_TAG := $(BUILD_TAG)-omp
endif

# ------------------------------------------------------------------ DEBUG
#
# DEBUG gets its own object tree for the same reason PROFILE and OPENMP do, and
# this one has already cost real time.
#
# The architecture makefiles default DEBUG to 1, which compiles device code with
# -g -G and Fortran with -g3 and no -O. `make DEBUG=0` afterwards recompiles
# NOTHING, because make decides from timestamps and no source changed -- so the
# "optimised" build links whatever debug objects were lying in build/ and
# reports success.
#
# Measured, on 2026-08-09: a build/ tree whose 3b_final.o was three days old and
# DEBUG=1 ran the three-body kernel 4.7x slower than the same source rebuilt
# clean (173 s against 37 s on 125k atoms), and its stale soap_turbo objects
# shifted the SOAP energy by ~1e-10 relative -- enough to look like a numerical
# regression from an unrelated change. Every timing measured against that tree
# was wrong, and nothing in the build said so.
#
# Tagging the tree makes the two builds physically separate, so the flag cannot
# lie about what is in the binary:
#     make DEBUG=0   ->  build/       bin/turbogap
#     make DEBUG=1   ->  build-dbg/   bin-dbg/turbogap
#
# DEBUG=0 is left untagged deliberately: it is the build everything else refers
# to -- bin/turbogap is what the regression suite, the cron job and the docs all
# name.
ifeq ($(DEBUG),1)
  BUILD_TAG := $(BUILD_TAG)-dbg
endif

# ---------------------------------------------------------------- PROFILE=1
#
# A build the NVIDIA profilers can read. Kept here rather than in the
# architecture makefiles so every CUDA architecture gets it from one place, and
# so it is applied AFTER the include -- these are additions to whatever
# optimisation level the architecture chose, never a replacement for it.
#
# What each flag buys, and why none of them is the debug build:
#
#   -lineinfo   maps device instructions back to source lines. This is what
#               makes an ncu report say which line stalled, instead of which
#               SASS address. It does NOT change codegen.
#   -g (host)   host symbols, so the nsys CPU sampler and the backtraces in an
#               ncu report show TurboGAP function names rather than addresses.
#               Host-side only; it does not touch device code.
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

# ---------------------------------------------------------------- the trees
#
# One object tree per flag combination, and this is not tidiness. make decides
# whether to rebuild from TIMESTAMPS, not from flags, so `make PROFILE=1` over
# an existing build/ relinks the old objects and produces a binary with no line
# info and no NVTX ranges -- while reporting success. The same trap applies to
# OPENMP=1, where the stale objects would be the ones with the directives
# compiled out, i.e. exactly the thing being tested.
#
# Splitting the trees makes the builds independent, and in particular means
# neither can invalidate bin/turbogap while the regression suite is reading it.
# Every rule, plus clean and deepclean, is written in terms of these four
# variables, and the generated Makefile.deps in terms of $(BUILD_DIR), so this
# is the only place the split has to be stated.
#
#     make                        build/       bin/turbogap
#     make OPENMP=1               build-omp/   bin-omp/turbogap
#     make PROFILE=1 DEBUG=0      build-profile/ ...
#     make OPENMP=1 PROFILE=1 DEBUG=0          build-omp-profile/ ...
BUILD_DIR=build$(BUILD_TAG)
BIN_DIR=bin$(BUILD_TAG)
INC_DIR=include$(BUILD_TAG)
LIB_DIR=lib$(BUILD_TAG)

# ------------------------------------------------------- tuning a device build
#
# Extra flags for the device compilers, and a tag so a tuned build gets its own
# object tree. Both are needed together: without the tag, make sees unchanged
# sources and relinks the previous variant's objects, so the "comparison" is a
# binary against itself.
#
#     make CC_EXTRA=-DTG_3B_WARPS_PER_BLOCK=1 BUILD_TAG_EXTRA=-3bw1
#     make CC_EXTRA=-DTG_3B_WARPS_PER_BLOCK=8 BUILD_TAG_EXTRA=-3bw8
#
# TG_3B_WARPS_PER_BLOCK is the one knob this exists for today: how many atoms a
# block of the three-body kernel covers, one per warp. 1 reproduces the geometry
# that kernel had before 2026-08-08 and is the reference to diff against.
CC += $(CC_EXTRA)
CU += $(CU_EXTRA)
F90_OPTS += $(F90_EXTRA)
BUILD_TAG := $(BUILD_TAG)$(BUILD_TAG_EXTRA)

# F90_EXTRA is how you get a bounds-checked build without disturbing bin/:
#
#     make F90_EXTRA='-fcheck=all -g -fbacktrace' BUILD_TAG_EXTRA=-chk
#
# which names the file and line of an out-of-range Fortran array access instead
# of leaving a glibc "free(): invalid next size" and a backtrace of addresses.
# Note what it CANNOT see: a cpy_dtoh writing more bytes into a host array than
# it holds is a memcpy through c_loc, invisible to Fortran bounds checking. For
# those, valgrind or an inspection of the byte counts is the tool.

# Do not change anything below this line
##########################################################

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR)

PROGRAMS := turbogap


# Device sources, all under src/gpu/ and named after what they compute. See
# src/gpu/README.md for the partition; in short, gpu_* is infrastructure both
# sides use, gap_* is the interatomic potential (2b, 3b, soap_turbo) and mad_*
# is the experimental-data path (pdf, structure factor, electrostatics).
#
# Every one of these was carved out of the three monolithic files this list used
# to name -- cuda_wrappers.cu (2684 lines), gpu_exp.cu (1709) and 3b_final.cc --
# so a change to one kernel family no longer recompiles the other twelve.
SRC_CUDA := gpu_memory.cu gpu_blas.cu gpu_scan.cu \
            gap_predict.cu gap_soap_radial.cu gap_soap_angular.cu \
            gap_soap_descriptor.cu gap_soap_forces.cu gap_2b.cu \
            mad_pdf.cu mad_xrd.cu mad_electrostatics.cu
# gap_3b stays .cc, and so is compiled by $(CC) rather than $(CU), because it
# needs -std=c++20 for <numbers> and <bit>; only $(CC) passes it.
#
# orthonormalization_kernels.cc is deliberately disabled, not dead: it is kept in
# src/gpu/ for possible reintegration. Re-enable by moving it before the '#'.
SRC_CC :=  gap_3b.cc # orthonormalization_kernels.cc
SRC := printing.f90 error.f90 read_utils.f90 nvtx.f90 timing.f90 misc.f90 constants.f90 nonneg_leastsq.f90 splines.f90 types.f90 gpu_context.f90 neighbors.f90 gap.f90 vdw.f90 local_properties.f90 exp_utils.f90 \
       xyz.f90 md.f90 mc.f90 read_files.f90 \
       gap_backend_gpu.f90 gap_interface.f90 mpi.f90 exp_interface.f90 turbogap_exp.f90 turbogap_md.f90 turbogap_vdw.f90 turbogap_estat.f90 turbogap_setup.f90 \
       electrostatics.f90
# kinds must build before every other group.
SRC_BASE := kinds.f90
SRC_TP_BT := resamplekin.f90 fortran_cuda_interfaces.f90
SRC_ST := soap_turbo_functions.f90 soap_turbo_radial.f90 soap_turbo_angular.f90 \
          soap_turbo.f90 soap_turbo_compress.f90
SRC_STOP := adaptive_time.f90 electronic_stopping.f90 eph_beta.f90 eph_fdm.f90 \
            eph_electronic_stopping.f90

OBJ_CUDA := $(addprefix $(BUILD_DIR)/,$(patsubst %.cu,%.o,$(SRC_CUDA)))
OBJ_CC := $(addprefix $(BUILD_DIR)/,$(patsubst %.cc,%.o,$(SRC_CC)))
OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC)))
OBJ_BASE := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_BASE)))
OBJ_TP_BT := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_TP_BT)))
OBJ_ST := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_ST)))
OBJ_STOP := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_STOP)))

PROG := $(addprefix $(BIN_DIR)/,$(PROGRAMS))

.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: default all programs clean deepclean libturbogap

default: libturbogap programs

all: default

clean:
	rm -rf $(OBJ_BASE) $(OBJ_CUDA) $(OBJ_CC) $(OBJ_TP_BT) $(OBJ_ST)  $(OBJ) $(INC_DIR)/*.mod $(PROG)

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

.SECONDEXPANSION:
.SECONDARY: $(OBJS)

programs: $(PROG)


libturbogap: $(OBJ_BASE) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(OBJ_CUDA) ${LIB_DIR}
	ar scr $(LIB_DIR)/libturbogap.a $(OBJ_BASE) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ)  $(OBJ_CUDA)

# Device sources. The headers are listed as prerequisites explicitly: they are
# not in Makefile.deps (which is Fortran modules only) and nvcc is not asked to
# emit depfiles, so without this an edit to gpu_common.h would rebuild nothing
# -- the same class of silent-staleness trap the BUILD_TAG note above describes.
GPU_HEADERS := $(wildcard src/gpu/*.h)

$(BUILD_DIR)/%.o: src/gpu/%.cu $(GPU_HEADERS) | $$(@D)
	$(CU) $(CUDA_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/gpu/%.cc $(GPU_HEADERS) | $$(@D)
	$(CC) $(CC_OPTS) -c $< -o $@
$(BIN_DIR)/%: src/%.f90 $(OBJ_BASE) $(OBJ_TP_BT) $(OBJ_ST)  $(OBJ) $(OBJ_CUDA) $(OBJ_CC) | $$(@D)
	$(F90) $(PP) $(F90_OPTS) $< -o $@ $(OBJ_BASE) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(OBJ_CUDA) $(OBJ_CC) $(LIBS)

$(BUILD_DIR)/%.o: src/stopping/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/third_party/nnls/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/third_party/bussi_thermostat/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/soap_turbo/src/%.f90 | $$(@D)
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
# Regenerate after adding, removing or moving a USE / module:
#     python3 tools/gen_fortran_deps.py . > makefiles/Makefile.deps
include makefiles/Makefile.deps
