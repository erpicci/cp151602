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
 * QR iterative eigenvalue algorithm (Householder, parallel matrix multiplication by row).
 * QR iterative eigenvalue algorithm which uses Householder
 * transformations as described in The algebraic eigenvalue problem
 * @cite wilkinson1965 and Applied Numerical Linear Algebra
 * @cite demmel1997 .
 * Matrix multiplications are parallelized by rows.
 * @file qr_householder_r.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "qr_householder.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../utils.h"
#include "basic_la.h"
#include "mpi_la.h"



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
    rowsMultiply(C, A, B, size, size, size)



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
 * Computes vector v.
 * v = [0, 0, ... alpha, -T(i, k), -T(i + 1, k), ... -T(size, k)]
 * @param[out] v     Vector
 * @param[in]  k     Current iteration
 * @param[in]  T     Matrix
 * @param[in]  size  Size of the matrix
 * @param[in]  alpha Alpha
 */
static void
get_v(double v[], const unsigned int k, const unsigned int size, const double T[], const double alpha);


/**
 * Returns alpha.
 * alpha = +- norm(T(k:size, k)). 
 * @param[in] k    Current iteration
 * @param[in] size Size of the matrix
 * @param[in] T    Matrix
 * @return Alpha
 */
static double
get_alpha(const unsigned int k, const unsigned int size, const double T[]);


/**
 * QR decomposition: A = Q * R.
 * Uses Householder transformations.
 * @param[out] Q    Orthonormal matrix
 * @param[out] R    Upper triangular matrix
 * @param[in]  size Size of the matrices
 * @param[in]  A    Input matrix
 */
static void
qr(double Q[], double R[], const unsigned int size, const double A[]);


/**
 * Computes Q = I - 2 * v' * v. 
 * @param[out] Q    Output matrix
 * @param[in]  size Size of the matrix
 * @param[in]  v    Vector
 */
static void
get_q(double Q[], const unsigned int size, const double v[]);


/**
 * Returns norm of a vector.
 * Considers only elements whose index is greather of equal to k.
 * @param[in] v    Vector
 * @param[in] k    Starting index
 * @param[in] size Size of the vector
 * @return Norm of the vector
 */
static double
get_norm_v(const double v[], const unsigned int k, const unsigned int size);


/**
 * Computes a "good sigma" for the QR iterative algorithm.
 * The closer the value of sigma to an eigenvalue of M, the faster the
 * convergence of the QR iterative algorithm.
 * @param[in] M    Matrix
 * @param[in] size Size of the matrix
 * @return A value close to an eigenvalue of M
 */
static double
compute_sigma(const double M[], const unsigned int size);



static void
qr(double Q[], double R[], const unsigned int size, const double A[]) {
    int mpi_rank;
    unsigned int i, k;
    double *H, *v, *temp;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Bcast((unsigned int *) &size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (ROOT == mpi_rank) {
        SAFE_MALLOC(H, double *, size * size * sizeof(double));
        SAFE_MALLOC(v, double *, size * sizeof(double));
        SAFE_MALLOC(temp, double *, size * size * sizeof(double));

        /* Q = I */
        memset(Q, 0, size * size * sizeof(double));
        for (i = 0; i < size; ++i) {
            Q[i * size + i] = 1.0;
        }
    
        /* R = A */
        memcpy(R, A, size * size * sizeof(double));

        /* Main for loop */
        for (k = 0; k < size - 1; ++k) {
            double normv;
            const double alpha = get_alpha(k, size, R);
            get_v(v,k, size, R, alpha);
            normv = get_norm_v(v, k, size);
        
            /* If norm of v is zero, no operation is needed */
            if (normv < 1e-6) {
                continue;
            }

            /* Normalizes v */
            for(i = k; i < size; i++) {
                v[i] = v[i] / normv;
            }
        
            /* H = I - 2 * v' * v */
            get_q(H, size, v); 

            /* Q = Q * H */
            post_multiply(Q, H, size, temp);
        
            /* R = H * R */
            pre_multiply(R, H, size, temp);
        }
    
        /* Frees memory */
        free(H);
        free(v);
        free(temp);
    }

    /* Slave nodes waits incoming multiplications */
    else {
        double *B;
        SAFE_MALLOC(B, double *, size * size * sizeof(double));
        for (k = 0; k < size - 1; ++k) {
            rowsMultiply(NULL, NULL, B, size, size, size);
            rowsMultiply(NULL, NULL, B, size, size, size);
        }
        free(B);
    }
}



static double
get_alpha(const unsigned int k, const unsigned int size, const double T[]) {
    unsigned int j;
    double res = 0.0;
    
    for (j = k; j < size; ++j) {
        res += T[j * size + k] * T[j * size + k];
    }
    res = sqrt(res); 
    
    return (T[k * size + k] > 1e-6) ? -res : res;
}



static void
get_v(double v[], const unsigned int k, const unsigned int size, const double T[], const double alpha) {
    unsigned int i;

    memset(v, 0, k * sizeof(double));
    
    v[k] = -T[k * size + k] + alpha;
    
    for (i = k + 1; i < size; ++i) {
        v[i] = -T[i * size + k];
    }
}



static double
get_norm_v(const double v[], const unsigned int k, const unsigned int size) {
    unsigned int i;
    double ret = 0.0;
    
    for (i = k; i < size; ++i) {
        ret += v[i] * v[i];
    }
    
    return sqrt(ret);
}



static void
get_q(double Q[], const unsigned int size, const double v[]) {
    unsigned int i;

    matrix_multiply(Q, v, v, size, size, 1);
    for(i = 0; i < size * size; ++i) {
        Q[i] *= -2;
    }
    for (i = 0; i < size; ++i) {
        Q[i * size + i] += 1;
    }
}



static double
compute_sigma(const double M[], const unsigned int size) {
    const double a_1 = M[(size - 2) * size + size - 2],
                 a_2 = M[(size - 1) * size + size - 1],
                 b   = M[(size - 1) * size + size - 2],
                 d   = 0.5 * (a_1 + a_2);

    return a_2 + d - sign(d) * sqrt(d * d + b * b);
}



void
qr_householder_r(const st_matrix_t M, double *eigenvalues, double *eigenvectors) {
    double *T, *Q, *R;
    int mpi_rank;
    unsigned int i, stop;
    const unsigned int size = st_matrix_size(M);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (ROOT == mpi_rank) {
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
            stop = 0;
        
            MPI_Bcast(&stop, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

            /* T -= sigma * I */
            for (i = 0; i < size; ++i) {
                T[i * size + i] -= sigma;
            }

            qr(Q, R, size, T);
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
        stop = 1;
        MPI_Bcast(&stop, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        /* eigenvalues = diag(T) */
        for (i = 0; i < size; ++i) {
            eigenvalues[i] = T[i * size + i];
        }


        /* Frees memory */
        free(T);
        free(Q);
        free(R);
    }

    else {
        while (1) {
            MPI_Bcast(&stop, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            if (stop) break;

            qr(NULL, NULL, 0, NULL);
        }
    }
}
