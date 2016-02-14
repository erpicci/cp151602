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
 * Parallel linear algebra utilities.
 * Basic utilities to handle vector and matrix multiplication, exploiting
 * parallelism with MPI.
 * @file mpi_la.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "../utils.h"
#include "basic_la.h"



/** MPI tag: tells a slave to compute a product */
#define COMPUTE 999

/** MPI tag: tells the master a slave is sending a result */
#define RESULT 1000

/** MPI tag: tells a slave there are no more rows to compute */
#define STOP 1001


void rowsMultiply(
    double C[], const double A[], const double B[],
    const unsigned int m, const unsigned int n, const unsigned int q) {


    int mpi_size, mpi_rank;
    int *rows, *offsets;
    int *sendcounts, *sdispls;
    int *recvcounts, *rdispls;
    int i;
    double *Ai, *Ci;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


    SAFE_MALLOC(rows, int *, mpi_size * sizeof(int));
    SAFE_MALLOC(offsets, int *, mpi_size * sizeof(int));
    SAFE_MALLOC(sendcounts, int *, mpi_size * sizeof(int));
    SAFE_MALLOC(sdispls, int *, mpi_size * sizeof(int));
    SAFE_MALLOC(recvcounts, int *, mpi_size * sizeof(int));
    SAFE_MALLOC(rdispls, int *, mpi_size * sizeof(int));



    /* Shares size of matrices */
    MPI_Bcast((void *) &m, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *) &n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *) &q, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);



    /* Offsets and counts */
    for (i = 0; i < mpi_size; ++i) {
        rows[i] = (m / mpi_size + (i < (int) m % mpi_size));;
    }

    offsets[0] = 0;
    for (i = 1; i < mpi_size; ++i) {
        offsets[i] = offsets[i - 1] + rows[i - 1];
    }

    for (i = 0; i < mpi_size; ++i) {
        sendcounts[i] = rows[i]    * q;
        sdispls[i]    = offsets[i] * q;

        recvcounts[i] = rows[i]    * n;
        rdispls[i]    = offsets[i] * n;
    }


    /* Sends rows of A */
    SAFE_MALLOC(Ai, double *, sendcounts[mpi_rank] * sizeof(double));
    MPI_Scatterv((void *) A, sendcounts, sdispls, MPI_DOUBLE,
                 Ai, sendcounts[mpi_rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    /* Sends B */
    MPI_Bcast((void *) B, n * q, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    /* Performs multiplication */
    SAFE_MALLOC(Ci, double *, recvcounts[mpi_rank] * sizeof(double));
    memcpy(Ci, Ai, recvcounts[mpi_rank] * sizeof(double));
    matrix_multiply(Ci, Ai, B, rows[mpi_rank], n, q);


    /*  Gathers Ci */
    MPI_Gatherv((void *) Ci, recvcounts[mpi_rank], MPI_DOUBLE,
                C, recvcounts, rdispls, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
}



/**
 * Reads a row from a matrix.
 * @param[out] row Row to read
 * @param[in]  A   Matrix to read from
 * @param[in]  dim Leading dimension of the matrix
 * @param[in]  num Index of the row
 */
void get_row(double *row, const double *A, const int dim, const int num);


/**
 * Writes a row to a matrix.
 * @param[in]  row Row to write
 * @param[out] C   Matrix to write to
 * @param[in]  dim Leading dimension of the matrix
 * @param[in]  num Index of the row
 */
void set_row(const double *row, double *C, const int dim, const int num);


void workerPoolMultiply(
    double C[],  const double A[],  const double B[], const unsigned int n)
{
    int i, k, righericev;
      int mpi_rank, mpi_size; 
    MPI_Status statusRic, statusSlave;
    double *row, *temprow;


    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    /* Shares dimension n and matrix B */
    MPI_Bcast((void *) &n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *) B, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    SAFE_MALLOC(row, double *, n * sizeof(double));
    SAFE_MALLOC(temprow, double *, n * sizeof(double));


    /* Master sends rows, one by one, to slaves.
     * Master node tracks which processor is assigned to which row(s). Master node is
     * assigned to no rows, and acts as an orchestrator.
     */
    if (0 == mpi_rank) {
       int *whosework;
       SAFE_MALLOC(whosework, int *, mpi_size*sizeof(int));
       whosework[0] = -1;

        for (k = 0; k < (int) n; k += mpi_size - 1) {
            for (i = 0; i < mpi_size - 1 && k + i < (int) n; ++i) {
                whosework[i + 1] = i + k;
                get_row(row, A, n, i + k);
                MPI_Send(row, n, MPI_DOUBLE, i + 1, COMPUTE, MPI_COMM_WORLD);
            }
            righericev = 0;
            while (righericev < i) {
                MPI_Recv(temprow, n, MPI_DOUBLE, MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &statusRic);
                set_row(temprow, C, n, whosework[statusRic.MPI_SOURCE]);
                ++righericev;
            }
        }

        /* Master informs every node that all the rows have been sent */
        for(i = 1; i < mpi_size; ++i) {
            MPI_Send(0, 0, MPI_INT, i, STOP, MPI_COMM_WORLD);
        } 
        free(whosework);
    }

    /* Slaves: receive rows of A and compute Ai * B */
    else {
        while (1) {
            MPI_Recv(temprow, n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &statusSlave);

            /* Stops if there are no more rows */
            if (statusSlave.MPI_TAG == STOP)
                break;

            matrix_multiply(row, temprow, B, 1, n, n);
            MPI_Send(row, n, MPI_DOUBLE, 0, RESULT, MPI_COMM_WORLD);
        }
    }

    /* Frees memory */
    free(row);
    free(temprow);
}



void get_row(double *row, const double *A, const int dim, const int num)
{
    memcpy(row, &A[dim * num] , dim * sizeof(double));
}



void set_row(const double *row, double *C, const int dim, const int num)
{
    memcpy(&C[dim * num] , row,  dim * sizeof(double));
}
