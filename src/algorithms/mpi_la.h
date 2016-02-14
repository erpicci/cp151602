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
#ifndef _MPILA_H_
#define _MPILA_H_


/**
 * Performs matrix multiplication C = A * B.
 * Multiplication is parallelized so that every processor computes a
 * row of C.
 * @param[out] C Result matrix
 * @param[in]  A First operand
 * @param[in]  B Second operand
 * @param[in]  m Number of rows of A
 * @param[in]  n Number of columns of B
 * @param[in]  q Number of columns of A and rows of B
 */
void rowsMultiply(
    double C[], const double A[], const double B[],
    const unsigned int m, const unsigned int n, const unsigned int q);



/**
 * Performs matrix multiplication C = A * B.
 * Matices are assumed to be squared. Multiplication is parallelized in
 * worker pool fashion: master node shares B, then sends one row to each
 * slave. When a slave finishes its work, it is assigned a new row.
 * @param[out] C Result matrix
 * @param[in]  A First operand
 * @param[in]  B Second operand
 * @param[in]  n Size of the matrices
 */
void workerPoolMultiply(
    double C[],  const double A[],  const double B[], const unsigned int n);

#endif  /* _MPILA_H_ */
