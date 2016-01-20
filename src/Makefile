########################################################################
# Parallel Computing - Makefile
# Erika Fabris <fabriser@dei.unipd.it>
# Thomas Gagliardi <gagliard@dei.unipd.it>
# Marco Zanella <marco.zanella.9@studenti.unipd.it>


########################################################################
# Configuration
CC      = mpcc
CCFLAGS = -O2 -qlanglvl=stdc89 -qlibansi -qlibmpi
LDFLAGS = -L/usr/linux/lib -lm -llapack

BIN     = ../bin/



########################################################################
# Dependencies
PROJ   = generator parallel #lapack
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
	@$(CC) $(CCFLAGS) -c $^ -o $@

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
