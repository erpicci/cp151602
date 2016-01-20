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
OPIMIZE = -O2
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
	CCOPT   = -qlanglvl=stdc89 -qlibansi -qlibmpi -c -DIBM
	LIBPATH = -L/usr/linux/lib -L/usr/lpp/xlf/lib
	LIBS    = -lm -llapack -lxlopt -lxl -lxlf90 -lblas
endif

CCFLAGS = $(CCOPT) $(OPTIMIZE) $(DEBUG)
LDFLAGS = $(LIBPATH) $(LIBS)



########################################################################
# Dependencies
PROJ   = generator parallel lapack
SERIAL = st_matrix.o serial_chronometer.o
MPI    = st_matrix.o mpi_chronometer.o

all: $(PROJ)

generator: st_matrix.o generator.o

lapack: $(SERIAL) lapack.o

# This is here for demo purposes
parallel: $(MPI) mpi_chronometer.o parallel.o

install: $(PROJ)

.PHONY: clean doc



########################################################################
# Recipes
%.o: %.c
	@echo "Compiling $@..."
	@$(CC) $(CCFLAGS) $^ -o $@

$(PROJ):
	@echo "Linking $@..."
	@$(CC) $(LDFLAGS) -o $@ $^
	@echo "done."

install:
	@echo "Copying under $(BIN)"
	@cp $(PROJ) $(BIN)

clean:
	@echo "Cleaning..."
	@rm -fR *.o $(PROJ)

doc:
	@echo "Generating documentation..."
	@doxygen ../doc/Doxyfile
