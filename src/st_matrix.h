/*
 * Copyright 2016 Erika Fabris, Thomas Gagliardi, Marco Zanella
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this file. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Symmetric tridiagonal matrix.
 * @file st_matrix.h
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#ifndef _ST_MATRIX_H_
#define _ST_MATRIX_H_

#include <stdio.h>

/** Simmetric tridiagonal matrix. */
typedef struct st_matrix_s *st_matrix_t;


/**
 * Creates a new symmetric tridiagonal matrix.
 * Allocates memory for the matrix.
 * @param[in] size Size of the matrix
 * @return New matrix
 */
st_matrix_t st_matrix_create(const unsigned int size);


/**
 * Deletes a symmetric tridiagonal matrix.
 * Deallocates memory and frees resources.
 * @param[out] M Matrix to delete
 */
void st_matrix_delete(st_matrix_t *M);


/**
 * Returns size of a symmetric tridiagonal matrix.
 * @param[in] M Symmetric tridiagonal matrix
 * @return Size of the matrix
 */
unsigned int st_matrix_size(const st_matrix_t M);


/**
 * Returns element at row i, column j.
 * Returns 0.0 if position is invalid.
 * @param[in] M Symmetric tridiagonal matrix
 * @param[in] i Row index
 * @param[in] j Column index
 * @return Element at row i, column j
 */
double
st_matrix_get(const st_matrix_t M, const unsigned int i, const unsigned int j);


/**
 * Returns diagonal of a symmetric trdiagonal matrix.
 * @param[in] M Symmetric tridiagonal matrix
 * @return Main diagonal of the matrix
 */
double *st_matrix_diag(const st_matrix_t M);


/**
 * Returns subdiagonal of a symmetric tridiagonal matrix.
 * @param[in] M Matrix
 * @return First subdiagonal of the matrix
 */
double *st_matrix_subdiag(const st_matrix_t M);


/**
 * Returns superdiagonal of a symmetric tridiagonal matrix.
 * @param[in] M Symmetric tridiagonal matrix
 * @return First superdiagonal of the matrix
 */
double *st_matrix_superdiag(const st_matrix_t M);


/**
 * Loads a symmetric tridiagonal matrix from a file.
 * @param[in,out] fp Pointer to file
 * @return New matrix, load from file
 */
st_matrix_t st_matrix_load(FILE *fp);


/**
 * Saves a symmetric tridiagonal matrix to a file.
 * @param[in]  M  Symmetric tridiagonal matrix
 * @param[out] fp Pointer to file
 * @return M itself
 */
st_matrix_t st_matrix_save(const st_matrix_t M, FILE *fp);


/**
 * Converts a symmetric tridiagonal matrix to dense representation.
 * @param[in]  M   Symmetric tridiagonal matrix
 * @param[out] dst Destination array for dense matrix
 * @return M itself
 */
st_matrix_t st_matrix_to_dense(const st_matrix_t M, double dst[]);

#endif  /* _ST_MATRIX_H_ */
