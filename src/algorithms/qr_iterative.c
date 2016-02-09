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
 * QR iterative eigenvalue algorithm.
 * Interface for the QR iterative eigenvalue algorithm proposed
 * by Givens as described in The algebraic eigenvalue problem
 * @cite wilkinson1965 and Applied Numerical Linear Algebra
 * @cite demmel1997.
 * @file qr_iterative.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "qr_iterative.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "../utils.h"
#include "basic_la.h"



/** Machine precision. */
#define EPSILON 1e-6



/** Returns sign of an arithmetic expression.
 * @param[in] v Value or expression to evaluate
 * @retval  1.0 If expression is positive
 * @retval  0.0 If expression is zero
 * @retval -1.0 If expression is negative
 */
#define sign(v) \
    (v > EPSILON) ? 1.0 : ((v < -EPSILON) ? -1.0 : 0.0)



/** Multiplies two squared matrices: C = A * B.
 * @param[out] C    Result matrix
 * @param[in]  A    First operand
 * @param[in]  B    Second operand
 * @param[in]  size Size of the matrices
 */
#define multiply(C, A, B, size) \
    matrix_multiply(C, A, B, size, size, size)



/** Pre-multiplies a squared matrix: A = B * A.
 * Uses a pre-allocated buffer to improve performance. Buffer must be
 * (at least) of the same size of the matrices.
 * @param[in, out] A      Matrix to be pre-multiplied
 * @param[in]      B      Operand
 * @param[in]      size   Size of the matrices
 * @param[out]     buffer Support memory buffer
 */
#define pre_multiply(A, B, size, buffer)              \
{                                                     \
    multiply(buffer, B, A, size);                     \
    memcpy(A, buffer, size * size * sizeof(double));  \
}



/**
 * Post-multiplies a squared matrix: A = A * B.
 * Uses a pre-allocated buffer to improve performance. Buffer must be
 * (at least) of the same size of the matrices.
 * @param[in, out] A      Matrix to be post-multiplied
 * @param[in]      B      Operand
 * @param[in]      size   Size of the matrices
 * @param[out]     buffer Support memory buffer
 */
#define post_multiply(A, B, size, buffer)             \
{                                                     \
    multiply(buffer, A, B, size);                     \
    memcpy(A, buffer, size * size * sizeof(double));  \
}



/**
 * Computes C = A * B, assuming C will be simmetric and tridiagonal.
 * Ony values on main diagonal and super/sub diagonals are computed,
 * resulting in a faster multiplication than regular matrix to matrix
 * one.
 * @param[out] C    Result matrix
 * @param[in]  A    First operand
 * @param[in]  B    Second operand
 * @param[in]  size Size of the matrices
 */
#define tridiagonal_multiply(C, A, B, size)                  \
{                                                            \
    unsigned int i, k;                                       \
    for (i = 0; i < size; ++i) {                             \
        double sum = 0.0;                                    \
        for (k = 0; k < size; ++k) {                         \
            sum += A[i * size + k] * B[k * size + i];        \
        }                                                    \
        C[i * size + i] = sum;                               \
    }                                                        \
    for (i = 0; i < size - 1; ++i) {                         \
        double sum = 0.0;                                    \
        for (k = 0; k < size; ++k) {                         \
            sum += A[(i + 1) * size + k] * B[k * size + i];  \
        }                                                    \
        C[(i + 1) * size + i] = sum;                         \
        C[i * size + i + 1] = sum;                           \
    }                                                        \
}



/**
 * Computes a "good sigma" for the QR iterative algorithm.
 * The closer the value of sigma to an eigenvalue of M, the faster the
 * convergence of the QR iterative algorithm.
 * @param[in] M    Matrix
 * @param[in] size Size of the matrix
 * @return A value close to an eigenvalue of M
 */
static double
compute_sigma(const double M[], const unsigned int size) {
    const double a_1 = M[(size - 2) * size + size - 2],
                 a_2 = M[(size - 1) * size + size - 1],
                 b   = M[(size - 1) * size + size - 2],
                 d   = 0.5 * (a_1 + a_2);

    return a_2 + d - sign(d) * sqrt(d * d + b * b);
}



/**
 * Computes cosine and sine of a rotation angle.
 * Angle is chosen such that a Givens rotation will set b to zero.
 * @param[out] c Cosine of the angle
 * @param[out] s Sine of the angle
 * @param[in]  a First value
 * @param[in]  b Second value
 */
static void
givens(double *c, double *s, const double a, const double b) {
    const double t = b / a;
    *c = 1 / sqrt(1 + t * t);
    *s = *c * t;
}



/**
 * Computes a Givens rotation matrix.
 * G_ij will be the Givens rotation matrix which will rotate the i-th
 * and j-th element of a vector clockwise by an angle theta such that
 * cos(theta) = c and sin(theta) = s.
 * @param[out] G_ij Givens rotation matrix
 * @param[in]  size Size of the matrix
 * @param[in]  i    First index
 * @param[in]  j    Second index
 * @param[in]  c    Cosine of the rotation angle
 * @param[in]  s    Sine of the rotation angle
 */
static void G(
    double G_ij[], const unsigned int size,
    const unsigned int i, const unsigned int j,
    const double c, const double s) {
    unsigned int k;

    /* G_ij = I */
    memset(G_ij, 0, size * size * sizeof(double));
    for (k = 0; k < size; ++k) {
        G_ij[k * size + k] = 1.0;
    }

    G_ij[i * size + i] = c;
    G_ij[i * size + j] = s;
    G_ij[j * size + i] = -s;
    G_ij[j * size + j] = c;
}



/**
 * QR decomposition: A = Q * R.
 * Uses Givens' rotation to obtain the QR decomposition of a matrix.
 * Assumes input matrix is a tridiagonal matrix.
 * @param[out] Q    Orthonormal matrix
 * @param[out] R    Upper triangular matrix
 * @param[in]  A    Input matrix
 * @param[in]  size Size of the matrix
 */
static void
QR(double Q[], double R[], const double A[], const unsigned int size) {
    unsigned int i;
    double c, s, *G_i, *buffer;


    SAFE_MALLOC(G_i, double *, size * size * sizeof(double));
    SAFE_MALLOC(buffer, double *, size * size * sizeof(double));


    /* Q' = I */
    memset(Q, 0, size * size * sizeof(double));
    for (i = 0; i < size; ++i) {
        Q[i * size + i] = 1.0;
    }

    /* R = A */
    memcpy(R, A, size * size * sizeof(double));


    /* Rotates matrix, setting to 0 every element in the subdiagonal */
    for (i = 0; i < size - 1; ++i) {
        if (fabs(R[(i + 1) * size + i]) < EPSILON) continue;

        givens(&c, &s, R[i * size + i], R[(i + 1) * size + i]);
        G(G_i, size, i, i + 1, c, s);
        
        /* Q' = G_i * Q' */
        matrix_multiply(buffer, G_i + i * size, Q, 2, size, size);
        memcpy(Q + i * size, buffer, 2 * size * sizeof(double));
        
        /* R = G_i * R */
        matrix_multiply(buffer, G_i + i * size, R, 2, size, size);
        memcpy(R + i * size, buffer, 2 * size * sizeof(double));
    }
    matrix_transpose(Q, size, size);


    /* Frees memory */
    free(G_i);
    free(buffer);
}



void
qr_iterative(const st_matrix_t M, double *eigenvalues, double *eigenvectors) {
    double *T, *Q, *R;
    unsigned int i;
    const unsigned int size = st_matrix_size(M);

    SAFE_MALLOC(T, double *, size * size * sizeof(double));
    SAFE_MALLOC(Q, double *, size * size * sizeof(double));
    SAFE_MALLOC(R, double *, size * size * sizeof(double));

    /* Computes also eigenvectors, if not NULL */
    if (NULL != eigenvectors) {
        memset(eigenvectors, 0, size * size * sizeof(double));
        for (i = 0; i < size; ++i) {
            eigenvectors[i * size + i] = 1.0;
        }
    }

    st_matrix_to_dense(M, T);
    while (!is_diagonal(T, size)) {
        const double sigma = compute_sigma(T, size);
        
        /* T -= sigma * I */
        for (i = 0; i < size; ++i) {
            T[i * size + i] -= sigma;
        }

        QR(Q, R, T, size);
        tridiagonal_multiply(T, R, Q, size);
        

        /* T += sigma * I */
        for (i = 0; i < size; ++i) {
            T[i * size + i] += sigma;
        }

        /* Computes also eigenvectors, if needed: E = E * Q */
        if (NULL != eigenvectors) {
            post_multiply(eigenvectors, Q, size, R);
        }
    }

    /* eigenvalues = diag(T) */
    for (i = 0; i < size; ++i) {
        eigenvalues[i] = T[i * size + i];
    }


    /* Frees memory */
    free(T);
    free(Q);
    free(R);
}
