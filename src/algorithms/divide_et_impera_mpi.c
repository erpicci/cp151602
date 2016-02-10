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
 * Divide et impera parallel eigenvalue algorithm.
 * Interface for the parallel divide et impera eigenvalue algorithm proposed
 * by Cuppen (A divide and conquer method for the symmetric tridiagonal
 * eigenproblem @cite cuppen1980, Applied Numerical Linear Algebra
 * @cite demmel1997).
 * @file divide_et_impera_mpi.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "divide_et_impera_mpi.h"

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
 * Conquer phase.
 * @param[out] lambda Eigenvalues
 * @param[out] Q      First and last eigenvectors
 * @param[in]  rho    Correction factor
 * @param[in]  L1     Eigenvalues of first sub-matrix
 * @param[in]  Q1     First and last eigenvectors of first sub-matrix
 * @param[in]  n1     Size of first sub-matrix
 * @param[in]  L2     Eigenvalues of second sub-matrix
 * @param[in]  Q2     First and last eigenvectors of second sub-matrix
 * @param[in]  n2     Second submatrix
 */
static void
conquer(
    double lambda[], double Q[], const double rho,
    const double L1[], const double Q1[], const unsigned int n1,
    const double L2[], const double Q2[], const unsigned int n2) {
    const unsigned int n = n1 + n2;
    unsigned int i, j;
    double *u1, *u2, *D, *u, *usqr, *Qp, *alpha = lambda;
    struct params_s params;

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
    free(u1);
    free(u2);
    free(D);
    free(u);
    free(usqr);
    free(Qp);
}






static void
multi_conquer_master(
    double lambda[], const int count[], const int displ[], const double rho,
    const double L1[], const double Q1[], const unsigned int n1,
    const double L2[], const double Q2[], const unsigned int n2) {
    const unsigned int n = n1 + n2;
    unsigned int i;
    double *u1, *u2, *D, *u, *usqr;
    struct params_s params;

    SAFE_MALLOC(u1, double *, n1 * sizeof(double));
    SAFE_MALLOC(u2, double *, n2 * sizeof(double));
    SAFE_MALLOC(D, double *, (n + 1) * sizeof(double));
    SAFE_MALLOC(u, double *, n * sizeof(double));
    SAFE_MALLOC(usqr, double *, n * sizeof(double));



    /* u = [+- last row of Q1, first row of Q2] */
    for (i = 0; i < n1; i++) {
        u1[i] = sign(rho) * Q1[(n1 != 1) * n1 + i];
    }
    memcpy(u2, Q2, n2 * sizeof(double));


    /* Finds eigenvalues of D + rho * u * u' */
    merge(D + 1, u, L1, u1, n1, L2, u2, n2);
    for (i = 0; i < n; ++i) {
        usqr[i] = u[i] * u[i];
    }



    D[0] = DBL_MAX;
    MPI_Bcast((double *) &rho, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(D, n + 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(usqr, n, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);


    params.n    = n;
    params.rho  = fabs(rho);
    params.d    = D + 1;
    params.usqr = usqr;

    for (i = 0; i < (unsigned int) count[ROOT]; ++i) {
        lambda[i] = bisection(D[i + 1], D[i], secular_function, &params);
    }

    MPI_Gatherv(
        lambda, count[ROOT], MPI_DOUBLE,
        lambda, (int *) count, (int *) displ, MPI_DOUBLE,
        ROOT, MPI_COMM_WORLD);



    /* END */


    /* Frees memory */
    free(u1);
    free(u2);
    free(D);
    free(u);
    free(usqr);
}



static void
multi_conquer_slave(const unsigned int n, const unsigned int count, const unsigned int disp) {
    double rho, *D, *usqr, *lambda;
    struct params_s params;
    unsigned int i;

    SAFE_MALLOC(D, double *, (n + 1) * sizeof(double));
    SAFE_MALLOC(usqr, double *, n * sizeof(double));
    SAFE_MALLOC(lambda, double *, count * sizeof(double));

    MPI_Bcast(&rho, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(D, n + 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(usqr, n, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);


    params.n    = n;
    params.rho  = fabs(rho);
    params.d    = D + 1;
    params.usqr = usqr;

    for (i = 0; i < count - 1; ++i) {
        lambda[i] = bisection(D[i + disp + 1], D[i + disp], secular_function, &params);
    }

    MPI_Gatherv(
        lambda, count - 1, MPI_DOUBLE,
        NULL, NULL, NULL, MPI_DOUBLE,
        ROOT, MPI_COMM_WORLD);



    free(D);
    free(usqr);
    free(lambda);
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
    const unsigned int n1 = n / 2,
                       n2 = n - n1;
    const double rho      = e[n1 - 1];
    double *L1, *Q1, *L2, *Q2;

    SAFE_MALLOC(L1, double *, n1 * sizeof(double));
    SAFE_MALLOC(Q1, double *, 2 * n1 * sizeof(double));
    SAFE_MALLOC(L2, double *, n2 * sizeof(double));
    SAFE_MALLOC(Q2, double *, 2 * n2 * sizeof(double));


    /* Divide: [Qi, Li] = eigen(Ti) */
    d[n1 - 1] -= fabs(rho);
    d[n1]     -= fabs(rho);

    eig_rec(L1, Q1, d, e, n1);
    eig_rec(L2, Q2, d + n1, e + n1, n2);

    d[n1 - 1] += fabs(rho);
    d[n1]     += fabs(rho);


    /* Conquer: [lambda, Q] = conquer(L1, Q1, L2, Q2) */
    conquer(lambda, Q, rho, L1, Q1, n1, L2, Q2, n2);


    /* Frees memory */
    free(L1);
    free(L2);
    free(Q1);
    free(Q2);
    }
}



void divide_et_impera_mpi(const st_matrix_t M, double *eigenvalues) {
    unsigned int n = st_matrix_size(M);
    int mpi_rank, mpi_size, i, *count, *displ, *displ2;
    double *evecs;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


    /* Shares size of the problem */
    MPI_Bcast(&n, 1, MPI_UNSIGNED, ROOT, MPI_COMM_WORLD);


    /* Every node computes how many elements it will work on */
    SAFE_MALLOC(count,  int *, mpi_size * sizeof(int));
    SAFE_MALLOC(displ,  int *, mpi_size * sizeof(int));
    SAFE_MALLOC(displ2, int *, mpi_size * sizeof(int));

    {
    const int avg = n / mpi_size,
              ext = n % mpi_size;

    count[0] = avg + (0 < ext);
    displ[0] = 0;
    displ2[0] = 0;

    for (i = 1; i < mpi_size; ++i) {
        count[i] = avg + (i < ext);
        displ[i] = displ[i - 1] + count[i - 1];
        displ2[i] = 2 * displ[i];
    }
    }
    


    /* DIVIDE */
    {
    double *d, *e, *tmp;
    SAFE_MALLOC(evecs, double *, n * n * sizeof(double));
    SAFE_MALLOC(d, double *, count[mpi_rank] * sizeof(double));
    SAFE_MALLOC(e, double *, count[mpi_rank] * sizeof(double));

    /* Master node: modifies matrix with rank-1 correction */
    if (ROOT == mpi_rank) {
        double *d = st_matrix_diag(M);
        SAFE_MALLOC(tmp, double * ,n * sizeof(double));
        memcpy(tmp, st_matrix_subdiag(M), (n - 1) * sizeof(double));

        for (i = 0; i < mpi_size - 1; ++i) {
            double rho = fabs(tmp[displ[i + 1] - 1]);
            d[displ[i + 1] - 1] -= rho;
            d[displ[i + 1]]     -= rho;
        }
    }

    /* Slave nodes: prepare to receive data */
    else {
        SAFE_MALLOC(eigenvalues, double *, count[mpi_rank] * sizeof(double));
    }


    /* Shares modified matrix, every node works on its portion */
    MPI_Scatterv(
        st_matrix_diag(M), count, displ, MPI_DOUBLE,
        d, count[mpi_rank], MPI_DOUBLE,
        ROOT, MPI_COMM_WORLD);
    MPI_Scatterv(
        tmp, count, displ, MPI_DOUBLE,
        e, count[mpi_rank], MPI_DOUBLE,
        ROOT, MPI_COMM_WORLD);
    eig_rec(eigenvalues, evecs, d, e, count[mpi_rank]);


    /* Master node: restores original matrix (undoes rank-1 correction) */
    if (ROOT == mpi_rank) {
        for (i = 0; i < mpi_size - 1; ++i) {
            double *d = st_matrix_diag(M);
            double rho = fabs(e[displ[i] + count[i] - 1]);
            d[displ[i] + count[i] - 1] += rho;
            d[displ[i + 1]] += rho;
        }
        free(tmp);
    }
    free(d);
    free(e);
    }



    /* CONQUER */
    {
    int k;
    unsigned int n1 = count[mpi_rank], n2;
    double *E, *L, *Q, *L1, *L2, *Q1, *Q2;

    SAFE_MALLOC(E, double *, n * sizeof(double));
    SAFE_MALLOC(L, double *, n * sizeof(double));
    SAFE_MALLOC(Q, double *, 2 * n * sizeof(double));
    SAFE_MALLOC(L1, double *, n * sizeof(double));
    SAFE_MALLOC(Q1, double *, 2 * n * sizeof(double));
    SAFE_MALLOC(L2, double *, n * sizeof(double));
    SAFE_MALLOC(Q2, double *, 2 * n * sizeof(double));

    /* Shares subdiagonal, to compute rho locally */
    if (ROOT == mpi_rank) {
        memcpy(E, st_matrix_subdiag(M), (n - 1) * sizeof(double));
    }
    MPI_Bcast(E, n, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    /* Initializes Li and Qi */
    memcpy(L, eigenvalues, count[mpi_rank] * sizeof(double));
    memcpy(Q, evecs, 2 * count[mpi_rank] * sizeof(double));


    /* Tree-shaped gathering */
    for (k = 1; k < mpi_size; k *= 2) {
        int is_sender = (mpi_rank & k) && (mpi_rank % k == 0),
            is_recv   = !(mpi_rank & k) && (mpi_rank % k == 0) && (mpi_rank + k < mpi_size),
            is_last   = k * 2 >= mpi_size;

        /* Sends size, Li and Qi */
        if (is_sender) {
            int to  = mpi_rank - k;

            MPI_Send(&n1, 1, MPI_UNSIGNED, to, 0, MPI_COMM_WORLD);
            MPI_Send(L, n1, MPI_DOUBLE, to, 1, MPI_COMM_WORLD);
            MPI_Send(Q, 2 * n1, MPI_DOUBLE, to, 2, MPI_COMM_WORLD);
        }

        /* Receives size, Li and Qi, and merges them with its Li and Qi */
        else if (is_recv) {
            const unsigned int from = mpi_rank + k;
            const double rho = E[n1 - 1];
            MPI_Status status;

            /* Prepares L1 and Q1 */
            memcpy(L1, L, n1 * sizeof(double));
            memcpy(Q1, Q, 2 * n1 * sizeof(double));

            /* Receives n2, L2 and Q2 */
            MPI_Recv(&n2, 1, MPI_UNSIGNED, from, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(L2, n2, MPI_DOUBLE, from, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(Q2, 2 * n2, MPI_DOUBLE, from, 2, MPI_COMM_WORLD, &status);

            /* Merges results */
            if (!is_last) {
                conquer(L, Q, rho, L1, Q1, n1, L2, Q2, n2);
                n1 += n2;
            }
        }


        /* Special speedup for last iteration */
        if (is_last) {
            if (ROOT == mpi_rank) {
                const double rho = E[n1 - 1];
                multi_conquer_master(L, count, displ, rho, L1, Q1, n1, L2, Q2, n2);
            }
            else {
                multi_conquer_slave(n, count[mpi_rank] + 1, displ[mpi_rank]);
            }
        }
    }


    /* Master node: copies result to destination */
    if (0 == mpi_rank) {
        memcpy(eigenvalues, L, n * sizeof(double));
    }


    /* Frees memory */
    free(E);
    free(L);
    free(Q);
    free(L1);
    free(L2);
    free(Q1);
    free(Q2);
    }
    free(evecs);
    free(count);
    free(displ);
    free(displ2);
}
