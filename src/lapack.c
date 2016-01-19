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
#include "eigenvalues.h"
#include "chronometer.h"
#include "utils.h"


/**
 * LAPACK dstev_.
 * Computes all eigenvalues and, optionally, eigenvectors of a real
 * symmetric tridiagonal matrix A.
 * @param[in]      jobz 'N' for eigenvalues, 'V' for eigenvalues and eigenvectors
 * @param[in]      n    Order of the matrix, N >= 0
 * @param[in, out] d    Array of size n; on entry, the n diagonal
 *                      elements of the tridiagonal matrix A. On exit,
 *                      if INFO = 0, the eigenvaules in ascending order
 * @param[in, out]  e   On entry, the (n-1) subdiagonal elements of the
 *                      tridiagonal matrix A, stored in elements 1 to
 *                      n-1 of e; e(n) need not be set, but is used by
 *                      the routine. On exit, the contents of e are
 *                      destroyed
 * @param[out]     z    If jobz = 'V', then if INFO = 0, Z contains the
 *                      orthonormal eigenvectors of the matrix A, with
 *                      the i-th column of z holding the eigenvector
 *                      associated with D(i). If jobz = 'N', then z is
 *                      not referenced
 * @param[in]      ldz  The leading dimension of the array z. ldz >= 1
 *                      and if jobz = 'N', ldz >= max(1,n)
 * @param[in, out] work If jobz = 'N', work is not referenced
 * @param[out]     info = 0: successful exit,
 *                      < 0: if info = -i, the i-th argument had an
 *                           illegal value,
 *                      > 0: if info = i, the algorithm failed to
 *                           converge; i off-diagonal elements of e did
 *                           not converge to zero
 */
extern void
dstev_(char *jobz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);



void compute_eigenvalues(st_matrix_t M, double *eigenvalues) {
    char jobz = 'N';
    int info, n = st_matrix_size(M);
    double *d, *e, *work;

    d = eigenvalues;
    SAFE_MALLOC(e, double *, n * sizeof(double));
    SAFE_MALLOC(work, double *, 2 * n * sizeof(double));

    memcpy(d, st_matrix_diag(M), n * sizeof(double));
    memcpy(e, st_matrix_subdiag(M), (n - 1) * sizeof(double));

    dstev_(&jobz, &n, d, e, NULL, &n, work, &info);

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
