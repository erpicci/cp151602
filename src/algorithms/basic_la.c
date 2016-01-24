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
 * Basic linear algebra utilities.
 * Basic utilities to handle vector and matrix multiplication.
 * @file basic_la.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "basic_la.h"

#include <stdio.h>
#include <string.h>

/** Size of a block. */
#define BLK_SIZE 32


/**
 * Matrix block multiplication.
 * @param[out] C   Result block
 * @param[in]  A   First block
 * @param[in]  B   Second block
 * @param[in]  ldc Leading dimension of C
 * @param[in]  lda Leading dimension of A
 * @param[in]  ldb Leading dimension of B
 */
void matrix_multiply_block(
    double C[], const double A[], const double B[],
    const unsigned int ldc, const unsigned int lda, const unsigned int ldb) {
    unsigned int i, j, k;

    const double *B_k = B;
    for (k = 0; k < BLK_SIZE; ++k) {
        double *C_i = C;
        const double *Ap = A + k;

        for (i = 0; i < BLK_SIZE; ++i) {
            const double A_i_k = *Ap;
            for (j = 0; j < BLK_SIZE; ++j) {
                C_i[j] += A_i_k * B_k[j];
            }
            C_i += ldc;
            Ap += lda;
        }
        B_k += ldb;
    }
}



double
scalar_product(const double v[], const double w[], const unsigned int n) {
    unsigned int i;
    double sum = 0.0;

    for (i = 0; i < n; ++i) {
        sum += v[i] * w[i];
    }

    return sum;
}



void matrix_multiply(
    double C[], const double A[], const double B[],
    const unsigned int m, const unsigned int n, const unsigned int q) {
    unsigned int i, j, k;

    memset(C, 0, m * n * sizeof(double));

    /* Standard multiplication if matrix cannot be divided in blocks */
    if (m % BLK_SIZE || n % BLK_SIZE || q % BLK_SIZE) {
        for (k = 0; k < q; ++k) {
            for (i = 0; i < m; ++i) {
                for (j = 0; j < n; ++j) {
                    C[i * n + j] += A[i * q + k] * B[k * n + j];
                }
            }
        }
        return;
    }

    /* Otherwise, block multiplication */
    for (k = 0; k < q; k += BLK_SIZE) {
        const double *B_k = B + k * n;

        for (i = 0; i < m; i += BLK_SIZE) {
            const double *A_i_k = A + i * q + k;
            double *C_i = C + i * n;

            for (j = 0; j < n; j += BLK_SIZE) {
                matrix_multiply_block(
                    C_i + j, A_i_k, B_k + j,
                    n, q, n
                );
            }
        }
    }
}
