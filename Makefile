# Copyright (c) 2020-2021 by Albert Bartók and Miguel Caro

SHELL = /bin/sh

# User-modifiable variables
F90=gfortran
MPIF90=mpif90
#F90=ifort
LAPACK_LIB_DIR=

# Default locations for various files
BUILD_DIR=build
BIN_DIR=bin
INC_DIR=include
LIB_DIR=lib

# Edit F90_OPTS, F90_MOD_DIR_OPT and LIBS as required
#########################################################
ifeq ($(MPIF90),mpif90)
F90=mpif90
F90_OPTS=-fPIC -O3
PP=-cpp -D _MPIF90
#F90_OPTS=-fPIC -O3 -fcheck=bounds -g -fcheck=all -Wall
F90_MOD_DIR_OPT=-J
#LIBS=-L$(LAPACK_LIB_DIR) -llapack -lopenblas
LIBS=-llapack -lblas
endif

ifeq ($(F90),gfortran)
F90_OPTS=-fPIC -O3
#F90_OPTS=-fPIC -O3 -fcheck=bounds -g -fcheck=all -Wall
PP=-cpp
F90_MOD_DIR_OPT=-J
#LIBS=-L$(LAPACK_LIB_DIR) -llapack -lopenblas
LIBS=-llapack -lblas
endif

ifeq ($(F90),ifort)
F90_OPTS=-fPIC -O3 -fp-model precise
PP=-fpp
F90_MOD_DIR_OPT=-module
LIBS=-L$(LAPACK_LIB_DIR) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
endif

# Do not change anything below this line
##########################################################

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR)

PROGRAMS := turbogap

SRC := splines.f90 types.f90 neighbors.f90 gap.f90 vdw.f90 read_files.f90 md.f90 \
       gap_interface.f90 mpi.f90
SRC_TP_BT := resamplekin.f90
SRC_GT := functions.f90 radial.f90 angular.f90 soap.f90

OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC)))
OBJ_TP_BT := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_TP_BT)))
OBJ_GT := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_GT)))

PROG := $(addprefix $(BIN_DIR)/,$(PROGRAMS))



.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: default all programs clean deepclean libturbogap

default: libturbogap programs

all: default

clean:
	rm -rf $(OBJ_TP_BT) $(OBJ_GT) $(OBJ) $(INC_DIR)/*.mod $(PROG)

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

.SECONDEXPANSION:
.SECONDARY: $(OBJS)

programs: $(PROG)

libturbogap: $(OBJ_TP_BT) $(OBJ_GT) $(OBJ) ${LIB_DIR}
	ar scr $(LIB_DIR)/libturbogap.a $(OBJ_TP_BT) $(OBJ_GT) $(OBJ)

$(BIN_DIR)/%: src/%.f90 $(OBJ_TP_BT) $(OBJ_GT) $(OBJ) | $$(@D)
	$(F90) $(PP) $(F90_OPTS) $< -o $@ $(OBJ_TP_BT) $(OBJ_GT) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o: src/third_party/bussi_thermostat/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/GAP_turbo/src/%.f90 | $$(@D)
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
