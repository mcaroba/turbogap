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
F90_OPTS=-fPIC -O3 -fcheck=bounds -g -fcheck=all -Wall
F90_MOD_DIR_OPT=-J
#LIBS=-L$(LAPACK_LIB_DIR) -llapack -lopenblas
LIBS=-llapack -lblas #-L/u/27/muhlih1/unix/psblas3/lib -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base  -llapack -lblas
endif

ifeq ($(F90),gfortran)
F90_OPTS=-fPIC -O3
F90_OPTS=-fPIC -O3 -fcheck=bounds -g -fcheck=all -Wall
PP=-cpp
F90_MOD_DIR_OPT=-J
#LIBS=-L$(LAPACK_LIB_DIR) -llapack -lopenblas
LIBS=-llapack -lblas #-L/u/27/muhlih1/unix/psblas3/lib -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base  -llapack -lblas
endif

ifeq ($(F90),ifort)
F90_OPTS=-fPIC -O3 -fp-model precise
PP=-fpp
F90_MOD_DIR_OPT=-module
LIBS=-L$(LAPACK_LIB_DIR) -llapack -lblas
endif

# Do not change anything below this line
##########################################################

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR)

PROGRAMS := turbogap

SRC := misc.f90 constants.f90 nonneg_leastsq.f90 splines.f90 types.f90 neighbors.f90 gap.f90 vdw.f90 read_files.f90 md.f90 \
       gap_interface.f90 mpi.f90 xyz.f90
SRC_TP_BT := resamplekin.f90
SRC_ST := soap_turbo_functions.f90 soap_turbo_radial.f90 soap_turbo_angular.f90 \
          soap_turbo.f90 soap_turbo_compress.f90

OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC)))
OBJ_TP_BT := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_TP_BT)))
OBJ_ST := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC_ST)))

PROG := $(addprefix $(BIN_DIR)/,$(PROGRAMS))



.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: default all programs clean deepclean libturbogap

default: libturbogap programs

all: default

clean:
	rm -rf $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(INC_DIR)/*.mod $(PROG)

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

.SECONDEXPANSION:
.SECONDARY: $(OBJS)

programs: $(PROG)

libturbogap: $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) ${LIB_DIR}
	ar scr $(LIB_DIR)/libturbogap.a $(OBJ_TP_BT) $(OBJ_ST) $(OBJ)

$(BIN_DIR)/%: src/%.f90 $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) | $$(@D)
	$(F90) $(PP) $(F90_OPTS) $< -o $@ $(OBJ_TP_BT) $(OBJ_ST) $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o: src/third_party/bussi_thermostat/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/soap_turbo/src/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@
$(BUILD_DIR)/%.o: src/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@ #-I/u/27/muhlih1/unix/psblas3/modules

$(BUILD_DIR): ${INC_DIR}
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(INC_DIR):
	mkdir -p $@

$(LIB_DIR):
	mkdir -p $@
