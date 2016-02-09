########################################################################
# Parallel Computing - Makefile
# Erika Fabris <fabriser@dei.unipd.it>
# Thomas Gagliardi <gagliard@dei.unipd.it>
# Marco Zanella <marco.zanella.9@studenti.unipd.it>


########################################################################
# Configuration
CC      = mpicc
CCOPT   = -Wall -Wextra -pedantic -ansi -Wno-long-long -c
INCPATH =
OPIMIZE = -O2 -funroll-loops -fomit-frame-pointer
DEBUG   =
LIBPATH =
LIBS    = -lm -llapack

BIN     = ../bin/

# Uses cc as compiler, if specified
ifdef cc
	CC = $(cc)
endif

# Includes debug symbols, if asked to
ifdef debug
	DEBUG = -g
endif

# Uses different options for mpcc (assumes xlc-like)
ifeq ($(CC), mpcc)
	CCOPT   = -qlanglvl=stdc89 -qlibansi -qlibmpi -c -DIBM -DHPC
	INCPATH = -I/usr/lpp/ppe.hpct/include
	LIBPATH = -L/usr/linux/lib -L/usr/lpp/xlf/lib -L/usr/lpp/ppe.hpct/lib
	LIBS    = -lm -llapack -lxlopt -lxl -lxlf90 -lblas -lhpc -lpmapi
endif

CCFLAGS = $(INCPATH) $(CCOPT) $(OPTIMIZE) $(DEBUG)
LDFLAGS = $(LIBPATH) $(LIBS)



########################################################################
# Dependencies
PROJ   = generator lapack qr_iterative qr_householder divide_et_impera \
         divide_et_impera_mpi
OBJ    = st_matrix.o stopwatch.o

all: $(PROJ)

generator: $(OBJ) generator.o

lapack: $(OBJ) solver_lapack.o
qr_iterative: $(OBJ) algorithms/basic_la.o algorithms/qr_iterative.o solver_qr.o
qr_householder: $(OBJ) algorithms/basic_la.o algorithms/qr_householder.o solver_householder.o
divide_et_impera: $(OBJ) algorithms/root_finding.o algorithms/basic_la.o \
                  algorithms/divide_et_impera.o solver_dei.o

divide_et_impera_mpi: $(OBJ) algorithms/root_finding.o algorithms/basic_la.o \
                      algorithms/divide_et_impera_mpi.o solver_dei_mpi.o

install: $(PROJ)

.PHONY: clean doc



########################################################################
# Recipes
%.o: %.c
	@echo "Compiling $@..."
	@$(CC) $(CCFLAGS) $^ -o $@

$(PROJ):
	@echo "Linking $@..."
	@$(CC) -o $@ $^ $(LDFLAGS)
	@echo "done."

install:
	@echo "Copying under $(BIN)"
	@cp $(PROJ) $(BIN)

clean:
	@echo "Cleaning..."
	@rm -fR *.dat *.o algorithms/*.o $(PROJ)

doc:
	@echo "Generating documentation..."
	@doxygen ../doc/Doxyfile
