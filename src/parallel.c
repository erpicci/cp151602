#include <stdio.h>
#include <mpi.h>
#include <math.h>

#include "st_matrix.h"
#include "stopwatch.h"
#include "utils.h"



/**
 * Computes eigenvalues and eigenvectors of a 2x2 symmetric, tridiagonal matrix.
 * @param[in]  d       Main diagonal
 * @param[in]  e       First subdiagonal
 * @param[out] eigvals Array of eigenvalues
 * @param[out] eigvect Array of eigenvectors
 */
static void
eig_2x2(const double *d, const double *e, double *eigvals, double *eigvect) {
    const double d_1  = d[0],
                 d_2  = d[1],
                 b    = e[0],
                 tr   = d_1 + d_2,
                 det  = d_1 * d_2 - b * b,
                 sqrD = sqrt(tr * tr - 4.0 * det),
                 l1   = 0.5 * (tr - sqrD),
                 l2   = 0.5 * (tr + sqrD);

    /* Computes eigenvalues */
    eigvals[0] = l1;
    eigvals[1] = l2;

    /* Computes eigenvectors */
    eigvect[0] = 1.0;
    eigvect[1] = (l1 - d_1) / b;
    eigvect[2] = 1.0;
    eigvect[3] = (l2 - d_1) / b;
}



/**
 * Computes eigenvalues and eigenvectors of a symmetric tridiagonal matrix.
 * Computation is recursive, with a divide-and-conquer approach.
 * @param[in]  d       Main diagonal
 * @param[in]  e       First subdiagonal
 * @param[in]  n       Size of the matrix
 * @param[out] eigvals Array of eigenvalues
 * @param[out] eigvect Array of eigenvectors
 */
static void
eig_rec(const double *d, const double *e, const unsigned int n, double *eigvals, double *eigvect) {
    if (n == 2) {
        eig_2x2(d, e, eigvals, eigvect);
        return;
    }
}



void compute_eigenvalues(st_matrix_t M, double *eigenvalues) {
    const unsigned int size = st_matrix_size(M);
    const double *diag = st_matrix_diag(M),
                 *sub  = st_matrix_subdiag(M);
    double *eigenvectors;
    SAFE_MALLOC(eigenvectors, double *, size * size * sizeof(double));
    eig_rec(diag, sub, size, eigenvalues, eigenvectors);
    free(eigenvectors);
}



int main(int argc, char *argv[]) {
    st_matrix_t M = st_matrix_load(stdin);
    const unsigned int size = st_matrix_size(M);
    double *eigenvalues;
    unsigned int i;
    stopwatch_t stopwatch = stopwatch_create("Sample task");

    (void) argc;
    (void) argv;

    /* Initializes environment */
    MPI_Init(&argc, &argv);
    SAFE_MALLOC(eigenvalues, double *, size * sizeof(double));


    /* Computes eigenvalues */
    stopwatch_start(stopwatch, 0, "Compute eigenvalues");
    compute_eigenvalues(M, eigenvalues);
    stopwatch_stop(stopwatch, 0);


    /* Prints results */
    printf("Eigenvalues:\n[");
    for (i = 0; i < size - 1; ++i) {
        printf("%g, ", eigenvalues[i]);
    }
    printf("%g]\n", eigenvalues[i]);
    

    /* Frees memory */
    MPI_Finalize();
    st_matrix_delete(&M);
    free(eigenvalues);
    stopwatch_delete(&stopwatch);
    
    return 0;
}
