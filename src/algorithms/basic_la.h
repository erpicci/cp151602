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
 * @file basic_la.h
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#ifndef _BASIC_LA_H_
#define _BASIC_LA_H_

/**
 * Computes scalar product v * w.
 * @param[in] v First vector
 * @param[in] w Second vector
 * @param[in] n Size of v and w
 * @return Scalar product between v and w
 */
double
scalar_product(const double v[], const double w[], const unsigned int n);


/**
 * Tells whether a square matrix is diagonal.
 * @param[in] M Matrix
 * @param[in] n Size of the matrix
 * @retval 1 if M is diagonal
 * @retval 0 if M is not diagonal
 */
unsigned int is_diagonal(const double M[], const unsigned int n);


/**
 * Transposes a matrix.
 * @param[in, out] M Matrix to transpose
 * @param[in]      m Number of rows of M
 * @param[in]      n Number of columns of M
 */
void
matrix_transpose(double M[], const unsigned int m, const unsigned int n);


/**
 * Matrix-matrix multiplication C = A * B.
 * A is m x q, B is q x n, hence C is m x n.
 * @param[out] C Result matrix
 * @param[in]  A First matrix
 * @param[in]  B Second matrix
 * @param[in]  m Number of rows of A and C
 * @param[in]  n Number of columns of B and C
 * @param[in]  q Number of columns of A and rows of B
 */
void matrix_multiply(
    double C[], const double A[], const double B[],
    const unsigned int m, const unsigned int n, const unsigned int q);

#endif  /* _BASIC_LA_H_ */
