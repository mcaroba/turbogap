# Copyright (c) 2020-2023 by Albert Bartók and Miguel Caro

SHELL = /bin/sh

# Include user-modifiable variables from a customizable file.
# Check the makefiles/ directory for a list of tested architectures

TURBOGAP_ARCH ?= Ubuntu_gfortran_mpi

include makefiles/Makefile.$(TURBOGAP_ARCH)


# Default locations for various files
BUILD_DIR=build
BIN_DIR=bin
INC_DIR=include
LIB_DIR=lib

# Do not change anything below this line
##########################################################

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR)

PROGRAMS := turbogap


SRC := printing.f90 error.f90 read_utils.f90 timing.f90 misc.f90 electrostatics.f90 constants.f90 nonneg_leastsq.f90 splines.f90 types.f90 gpu_context.f90 neighbors.f90 gap.f90 vdw.f90		\
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
