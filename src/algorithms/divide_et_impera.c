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
 * Interface for the basic divide et impera eigenvalue algorithm proposed
 * by Cuppen (A divide and conquer method for the symmetric tridiagonal
 * eigenproblem @cite cuppen1980, Applied Numerical Linear Algebra
 * @cite demmel1997).
 * @file divide_et_impera.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "divide_et_impera.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "../utils.h"
#include "root_finding.h"
#include "basic_la.h"



/** Parameters for the secular equation. */
struct params_s {
    unsigned int n;  /**< Size of the matrix */
    double rho;      /**< Element on the subdiagonal */
    double *d;       /**< Eigenvalues */
    double *usqr;    /**< Eigenvector (elements are squared) */
};



/**
 * Returns sign of a value.
 * @param[in] v Value
 * @retval +1.0 If v is positive
 * @retval  0.0 If v is zero
 * @retval -1.0 if v is negative
 */
static double sign(const double v) {
    return (v > 0.0) ? 1.0 : ((v < 0.0) ? -1.0 : 0.0);
}



/**
 * Function for the secular equation.
 * @param[in] lambda Independent variable
 * @param[in] args    Parameters
 * @return Value of the function for given lambda
 */
static double secular_function(const double lambda, const void *args) {
    register double value = 0.0;
    const struct params_s *params = (struct params_s *) args;
    register int i = params->n - 1;
    register const double * const usqr = params->usqr;
    register const double * const d    = params->d;

    value += *usqr / (*d - lambda);
    for (; i; i -= 1) {
        value += usqr[i] / (d[i] - lambda);
    }

    return 1.0 + params->rho * value;
}



/**
 * Merges two couples of arrays.
 * Arrays a1 and a2 are merged keeping the resulting array sorted, b1
 * and b2 are merged preserving the same order of the elements in a1
 * and a2.
 * Both a and b must be large enough to hold n1 + n2 elements.
 * @param[out] a  Result main vector
 * @param[out] b  Result companion vector
 * @param[in]  a1 First main vector
 * @param[in]  b1 First companion vector
 * @param[in]  n1 Size of first vector
 * @param[in]  a2 Second main vector
 * @param[in]  b2 Second companion vector
 * @param[in]  n2 Size of second vector
 */
void merge(
    double a[], double b[],
    const double a1[], const double b1[], const unsigned int n1,
    const double a2[], const double b2[], const unsigned int n2) {
    const unsigned int n = n1 + n2;
    unsigned int k, i = 0, j = 0;

    for (k = 0; k < n; ++k) {
        if (a1[i] >= a2[j]) {
            a[k] = a1[i];
            b[k] = b1[i];
            if (++i == n1) break;
        }
        else {
            a[k] = a2[j];
            b[k] = b2[j];
            if (++j == n2) break;
        }
    }

    k += 1;
    memcpy(a + k, a1 + i, (n1 - i) * sizeof(double));
    memcpy(b + k, b1 + i, (n1 - i) * sizeof(double));
    memcpy(a + k, a2 + j, (n2 - j) * sizeof(double));
    memcpy(b + k, b2 + j, (n2 - j) * sizeof(double));
}



/**
 * Computes eigenvalues and eigenvectors of a symmetric tridiagonal matrix.
 * Computation is recursive with a divide et impera approach.
 * Lambda must be large enough to hold n elements, Q must be large
 * enough to hold 2 * n elements.
 * @param[out] lambda Array of eigenvalues
 * @param[out] Q      Array of eigenvectors (only first and last)
 * @param[in]  d      Main diagonal
 * @param[in]  e      First subdiagonal
 * @param[in]  n      Size of the matrix
 */
static void
eig_rec(double lambda[], double Q[], double d[], const double e[], const unsigned int n) {
    /* Base case: n = 1 */
    if (n == 1) {
        *lambda = *d;
        *Q      = 1.0;
        return;
    }

    /* Recursive case: n > 1 */
    {
    unsigned int i, j;
    const unsigned int n1 = n / 2,
                       n2 = n - n1;
    double *alpha = lambda,
           *L1, *Q1, *u1,
           *L2, *Q2, *u2,
           *Qp,
           *D, *u, *usqr;
    const double rho = e[n1 - 1];
    struct params_s params;

    SAFE_MALLOC(L1, double *, n1 * sizeof(double));
    SAFE_MALLOC(Q1, double *, 2 * n1 * sizeof(double));
    SAFE_MALLOC(L2, double *, n2 * sizeof(double));
    SAFE_MALLOC(Q2, double *, 2 * n2 * sizeof(double));


    /* [Qi, Li] = eigen(Ti) */
    d[n1 - 1] -= fabs(rho);
    d[n1]     -= fabs(rho);

    eig_rec(L1, Q1, d, e, n1);
    eig_rec(L2, Q2, d + n1, e + n1, n2);

    d[n1 - 1] += fabs(rho);
    d[n1]     += fabs(rho);


    /* Conquer */
    SAFE_MALLOC(u1, double *, n1 * sizeof(double));
    SAFE_MALLOC(u2, double *, n2 * sizeof(double));
    SAFE_MALLOC(D, double *, n * sizeof(double));
    SAFE_MALLOC(u, double *, n * sizeof(double));
    SAFE_MALLOC(usqr, double *, n * sizeof(double));
    SAFE_MALLOC(Qp, double *, n * n * sizeof(double));

    /* u = [+- last row of Q1, first row of Q2] */
    for (i = 0; i < n1; i++) {
        u1[i] = sign(rho) * Q1[(n1 != 1) * n1 + i];
    }
    memcpy(u2, Q2, n2 * sizeof(double));


    /* Finds eigenvalues of D + rho * u * u' */
    merge(D, u, L1, u1, n1, L2, u2, n2);
    for (i = 0; i < n; ++i) {
        usqr[i] = u[i] * u[i];
    }
    params.n    = n;
    params.rho  = fabs(rho);
    params.d    = D;
    params.usqr = usqr;

    for (i = 1; i < n; ++i) {
        lambda[i] = bisection(D[i], D[i - 1], secular_function, &params);
    }
    lambda[0] = bisection(D[0], DBL_MAX, secular_function, &params);
    /* END */


    /* D = [L1 0; 0 L2] */
    memcpy(D, L1, n1 * sizeof(double));
    memcpy(D + n1, L2, n2 * sizeof(double));


    /* u = [u1 u2] */
    memcpy(u, u1, n1 * sizeof(double));
    memcpy(u + n1, u2, n2 * sizeof(double));


    /* Finds eigenvectors Qprime of D + rho * u * u' */
    for (i = 0; i < n; i++) {
        double norm = 0.0;
        for (j = 0; j < n; j++) {
            double e = u[j] / (D[j] - alpha[i]);
            Qp[j * n + i] = e;
            norm += e * e;
        }
        norm = 1.0 / sqrt(norm);
        for (j = 0; j < n; ++j) {
            Qp[j * n + i] *= norm;
        }
    }


    /* Computes eigenvectors as [Q1 0; 0 Q2] * Qprime (only first and last) */
    matrix_multiply(Q, Q1, Qp, 1, n, n1);
    matrix_multiply(Q + n, Q2 + (n2 != 1) * n2, Qp + n1 * n, 1, n, n2);


    /* Frees memory */
    free(L1);
    free(L2);
    free(Q1);
    free(Q2);
    free(u1);
    free(u2);
    free(u);
    free(usqr);
    free(D);
    free(Qp);
    }
}



void divide_et_impera(const st_matrix_t M, double *eigenvalues) {
    const unsigned int n = st_matrix_size(M);
    double *d = st_matrix_diag(M),
                 *e = st_matrix_subdiag(M);
    double *evecs;

    SAFE_MALLOC(evecs, double *, n * n * sizeof(double));
    eig_rec(eigenvalues, evecs, d, e, n);
    free(evecs);
}
