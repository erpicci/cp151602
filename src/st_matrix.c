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
 * Symmetric tridiagonal matrix.
 * @file st_matrix.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include <stdlib.h>
#include <string.h>

#include "st_matrix.h"
#include "utils.h"


/**
 * Reads a single value from a file.
 * Calls fscanf and checks its return value.
 * @param[in]  fp   File pointer
 * @param[in]  type Type of the data (as the format string of fscanf)
 * @param[out] ptr  Pointer to which store data
 */
#define READ(fp, type, ptr)                                 \
{                                                           \
    int n = fscanf(fp, type, ptr);                          \
    if (1 != n) {                                           \
        fprintf(stderr, "Cannot read from file [%s:%d]\n",  \
                __FILE__, __LINE__);                        \
        exit(EXIT_FAILURE);                                 \
    }                                                       \
}


/**
 * Structure of a symmetric triangular matrix.
 * Elements of diagonal and subdiagonal are stored in contigous memory
 */
struct st_matrix_s {
    unsigned int size;  /**< Size of the matrix */
    double *elems;      /**< Elements on diagonals */
};



st_matrix_t st_matrix_create(const unsigned int size) {
    st_matrix_t M;

    SAFE_MALLOC(M, st_matrix_t, sizeof(struct st_matrix_s));
    SAFE_MALLOC(M->elems, double *, (2 * size - 1) * sizeof(double));
    M->size = size;

    return M;
}



void st_matrix_delete(st_matrix_t *M) {
    if (NULL == M) {
        return;
    }

    if (NULL != (*M)->elems) {
        free((*M)->elems);
    }

    free(*M);
    *M = NULL;
}



unsigned int st_matrix_size(const st_matrix_t M) {
    return (NULL != M) ? M->size : 0;
}



double
st_matrix_get(const st_matrix_t M, const unsigned int i, const unsigned int j) {
    const unsigned int min = (i < j) ? i : j,
                       max = (i < j) ? j : i,
                       d   = max - min;

    return (d <= 1 && !(NULL == M || i >= M->size || j >= M->size))
         ? M->elems[M->size * d + min]
         : 0.0;
}



double *st_matrix_diag(const st_matrix_t M) {
    return (NULL != M) ? M->elems : NULL;
}



double *st_matrix_subdiag(const st_matrix_t M) {
    return (NULL != M) ? (M->elems + M->size) : NULL;
}



double *st_matrix_superdiag(const st_matrix_t M) {
    return (NULL != M) ? (M->elems + M->size) : NULL;
}



st_matrix_t st_matrix_load(FILE *fp) {
    st_matrix_t M;
    unsigned int size, i = 0;
    double *elems;

    /* Reads size of the matrix */
    READ(fp, "%u", &size);
    M = st_matrix_create(size);
    elems = M->elems;

    /* Reads elements in main diagonal and subdiagonal */
    for (i = 0; i < 2 * size - 1; ++i) {
        READ(fp, "%lf", elems + i);
    }

    return M;
}



st_matrix_t st_matrix_save(const st_matrix_t M, FILE *fp) {
    unsigned int i, size = st_matrix_size(M);
    double *diag    = st_matrix_diag(M),
           *subdiag = st_matrix_subdiag(M);

    /* Writes size of the matrix */
    fprintf(fp, "%u\n", size);

    /* Writes diagonal of the matrix */
    for (i = 0; i < size; ++i) {
        fprintf(fp, "%g ", diag[i]);
    }

    /* Writes subdiagonal of the matrix */
    fprintf(fp, "\n");
    for (i = 0; i < size - 1; ++i) {
        fprintf(fp, "%g ", subdiag[i]);
    }

    return M;
}



st_matrix_t st_matrix_to_dense(const st_matrix_t M, double dst[]) {
    unsigned int i;
    const unsigned int size = st_matrix_size(M);
    double *diag    = st_matrix_diag(M),
           *subdiag = st_matrix_subdiag(M);

    memset(dst, 0, size * size * sizeof(double));
    for (i = 0; i < size; ++i) {
        dst[i * size + i] = diag[i];
    }
    for (i = 0; i < size - 1; ++i) {
        dst[(i + 1) * size + i] = subdiag[i];
        dst[i * size + i + 1]   = subdiag[i];
    }

    return M;
}
