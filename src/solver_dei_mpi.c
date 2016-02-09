/*
 * Copyright 2016 Erika Fabris, Thomas Gagliardi, Marco Zanella
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this file. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Computes eigenvalues of a symmetric tridiagonal matrix.
 * Uses the divide et impera approach with MPI.
 * @file solver_dei_mpi.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include <stdio.h>
#include <string.h>

#include "st_matrix.h"
#include "algorithms/divide_et_impera_mpi.h"
#include "stopwatch.h"
#include "utils.h"


/**
 * Divide et impera-based solver.
 * @param[in] argc ARGument Counter
 * @param[in] argv ARGument Vector
 * @retval EXIT_SUCCESS Normal termination of the program
 * @retval EXIT_FAILURE Some error occurred
 */
int main(const int argc, char * const argv[]) {
    int mpi_rank;

    MPI_Init((int *) &argc, (char ***) &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    
    /* Root node */
    if (ROOT == mpi_rank) {
        st_matrix_t M           = st_matrix_load(stdin);
        const unsigned int size = st_matrix_size(M);
        stopwatch_t stopwatch   = stopwatch_create("Divide_et_Impera_solver_MPI");
        double *eigenvalues;
        
        unsigned int i;

        /* Allocates resources */
        SAFE_MALLOC(eigenvalues, double *, size * sizeof(double));


        /* Computes eigenvalues */
        stopwatch_start(stopwatch, 0, "Compute eigenvalues");
        divide_et_impera_mpi(M, eigenvalues);
        stopwatch_stop(stopwatch, 0);


        /* Prints results */
        printf("Eigenvalues:\n[");
        for (i = 0; i < size - 1; ++i) {
            printf("%g, ", eigenvalues[i]);
        }
        printf("%g]\n", eigenvalues[i]);


        /* Frees memory */
        st_matrix_delete(&M);
        free(eigenvalues);
        stopwatch_delete(&stopwatch);
    }
    
    /* Slave nodes */
    else {
        divide_et_impera_mpi(NULL, NULL);
    }
    
    MPI_Finalize();

    return EXIT_SUCCESS;
}
