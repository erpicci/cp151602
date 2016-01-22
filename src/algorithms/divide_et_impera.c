/*
 * Copyright 2016 Erika Fabris, Thomas Gagliardi, Marco Zanella
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this file. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Divide et impera eigenvalue algorithm.
 * Implements the basic divide et impera eigenvalue algorithm
 * (Applied Numerical Linear Algebra @cite demmel1997).
 * @file divide_et_impera.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "divide_et_impera.h"
#include "lapack.h" /* Mock */
#include <string.h>
#include "../utils.h"

void divide_et_impera(const st_matrix_t M, double *eigenvalues) {
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
