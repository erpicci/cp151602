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
 * Uses LAPACK library to obtain eigenvalues.
 * @note LAPACK routines are used as a black box tool
 * @file lapack.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include <stdio.h>
#include <string.h>

#include "st_matrix.h"
#include "lapack.h"
#include "eigenvalues.h"
#include "chronometer.h"
#include "utils.h"


/**
 * Computes eigenvalues.
 * @param[in] M            Symmetric tridiagonal matrix
 * @param[out] eigenvalues Array of eigenvalues
 */
void compute_eigenvalues(st_matrix_t M, double *eigenvalues) {
    char jobz = 'N';
    int info, n = st_matrix_size(M);
    double *d, *e, *work;

    d = eigenvalues;
    SAFE_MALLOC(e, double *, n * sizeof(double));
    SAFE_MALLOC(work, double *, 2 * n * sizeof(double));

    memcpy(d, st_matrix_diag(M), n * sizeof(double));
    memcpy(e, st_matrix_subdiag(M), (n - 1) * sizeof(double));

    dstev(&jobz, &n, d, e, NULL, &n, work, &info);

    free(e);
    free(work);
}



/**
 * LAPACK-based solver.
 * @param[in] argc ARGument Counter
 * @param[in] argv ARGument Vector
 * @retval EXIT_SUCCESS Normal termination of the program
 * @retval EXIT_FAILURE Some error occurred
 */
int main(const int argc, char * const argv[]) {
    st_matrix_t M = st_matrix_load(stdin);
    const unsigned int size = st_matrix_size(M);
    double *eigenvalues, time;
    unsigned int i;
    chronometer_t chronometer = chronometer_create();

    (void) argc;
    (void) argv;

    /* Allocates resources */
    SAFE_MALLOC(eigenvalues, double *, size * sizeof(double));


    /* Computes eigenvalues */
    chronometer_start(chronometer);
    compute_eigenvalues(M, eigenvalues);
    time = chronometer_stop(chronometer);


    /* Prints results */
    printf("Eigenvalues:\n[");
    for (i = 0; i < size - 1; ++i) {
        printf("%g, ", eigenvalues[i]);
    }
    printf("%g]\n", eigenvalues[i]);
    printf("Time: %g ms\n", time);


    /* Frees memory */
    st_matrix_delete(&M);
    free(eigenvalues);
    chronometer_delete(&chronometer);

    return EXIT_SUCCESS;
}
