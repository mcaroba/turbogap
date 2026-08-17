# Copyright (c) 2020-2023 by Albert Bartók and Miguel Caro

SHELL = /bin/sh

# Include user-modifiable variables from a customizable file.
# Check the makefiles/ directory for a list of tested architectures

TURBOGAP_ARCH ?= Ubuntu_gfortran_mpi

include makefiles/Makefile.$(TURBOGAP_ARCH)

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


SRC := printing.f90 error.f90 read_utils.f90 timing.f90 misc.f90 electrostatics.f90 constants.f90 mad_ir.f90 nonneg_leastsq.f90 splines.f90 types.f90 gpu_context.f90 neighbors.f90 gap.f90 vdw.f90		\
	local_properties.f90 exp_utils.f90  xyz.f90 md.f90 mc.f90 read_files.f90	\
	gap_backend_cpu.f90 gap_interface.f90 mpi.f90 exp_interface.f90 turbogap_exp.f90 turbogap_md.f90 turbogap_vdw.f90 turbogap_estat.f90 turbogap_setup.f90

# kinds must build before everything, including SRC_STOP, so it gets its own
# group placed first in every prerequisite list.
SRC_BASE := kinds.f90 keyword_help.f90
SRC_TP_BT := resamplekin.f90
SRC_ST := soap_turbo_functions.f90 soap_turbo_radial.f90 soap_turbo_angular.f90 \
          soap_turbo.f90 soap_turbo_compress.f90
SRC_STOP := adaptive_time.f90 electronic_stopping.f90 eph_beta.f90 eph_fdm.f90 \
            eph_electronic_stopping.f90

OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC)))
OBJ_BASE := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_BASE)))
OBJ_TP_BT := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_TP_BT)))
OBJ_ST := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_ST)))
OBJ_STOP := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_STOP)))

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
	rm -rf $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(INC_DIR)/*.mod $(PROG)

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

.SECONDEXPANSION:
.SECONDARY: $(OBJS)

programs: $(PROG)

libturbogap: $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) ${LIB_DIR}
	ar scr $(LIB_DIR)/libturbogap.a $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ)

$(BIN_DIR)/%: src/%.f90 $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) | $$(@D)
	$(F90) $(PP) $(F90_OPTS) $< -o $@ $(OBJ_BASE) $(OBJ_STOP) $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o: src/stopping/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/third_party/bussi_thermostat/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/third_party/nnls/%.f90 | $$(@D)
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
