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

#include <string.h>
#include <math.h>
#include <float.h>

#include "../utils.h"
#include "root_finding.h"
#include "basic_la.h"



/** Parameters for the secular equation. */
struct params_s {
    unsigned int n;  /** Size of the matrix */
    double rho;      /** Element on the subdiagonal */
    double *d;       /** Eigenvalues */
    double *usqr;    /** Eigenvector (elements are squared) */
};



double sign(double v) {
    return (v > 0.0) ? 1.0 : ((v < 0.0) ? -1.0 : 0.0);
}



/**
 * Function for the secular equation.
 * @param[in] lambda Independent variable
 * @param[in] params Parameters
 * @return Value of the function for given lambda
 */
double secular_function(double lambda, void *args) {
    unsigned int i;
    double value = 0.0;
    struct params_s *params = (struct params_s *) args;

    for (i = 0; i < params->n; ++i) {
        value += params->usqr[i] / (params->d[i] - lambda);
    }

    return 1.0 + params->rho * value;
}

double secular_function_prime(double lambda, void *args) {
    unsigned int i;
    double value = 0.0;
    struct params_s *params = (struct params_s *) args;

    for (i = 0; i < params->n; ++i) {
        value += params->usqr[i] / ((params->d[i] - lambda) * (params->d[i] - lambda));
    }

    return params->rho * value; 
}



/**
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
 * @param[in]  d     Main diagonal
 * @param[in]  e     First subdiagonal
 * @param[in]  n     Size of the matrix
 * @param[out] evals Array of eigenvalues
 * @param[out] evecs Array of eigenvectors
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
    SAFE_MALLOC(Q1, double *, n1 * n1 * sizeof(double));
    SAFE_MALLOC(L2, double *, n2 * sizeof(double));
    SAFE_MALLOC(Q2, double *, n2 * n2 * sizeof(double));
    SAFE_MALLOC(u1, double *, n1 * sizeof(double));
    SAFE_MALLOC(u2, double *, n2 * sizeof(double));
    SAFE_MALLOC(D, double *, n * sizeof(double));
    SAFE_MALLOC(u, double *, n * sizeof(double));
    SAFE_MALLOC(usqr, double *, n * sizeof(double));
    SAFE_MALLOC(Qp, double *, n * n * sizeof(double));


    /* [Qi, Li] = eigen(Ti) */
    d[n1 - 1] -= fabs(rho);
    d[n1]     -= fabs(rho);

    eig_rec(L1, Q1, d, e, n1);
    eig_rec(L2, Q2, d + n1, e + n1, n2);

    d[n1 - 1] += fabs(rho);
    d[n1]     += fabs(rho);


    /* u = [+- last row of Q1, first row of Q2] */
    for (i = 0; i < n1; i++) {
        u1[i] = sign(rho) * Q1[(n1 - 1) * n1 + i];
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


    /* Computes eigenvectors as [Q1 0; 0 Q2] * Qprime */
/*
    matrix_multiply(Q,          Q1, Qp,          n1, n, n1);
    matrix_multiply(Q + n1 * n, Q2, Qp + n1 * n, n2, n, n2);
*/
matrix_multiply(Q,      Q1, Qp,     1, n, n1);
matrix_multiply(Q + (n - 1) * n, Q2 + (n2 - 1) * n2, Qp + (n - n2) * n, 1, n, n2);


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
