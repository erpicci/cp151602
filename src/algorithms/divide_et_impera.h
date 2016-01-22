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
 * Interface for the basic divide et impera eigenvalue algorithm
 * (Applied Numerical Linear Algebra @cite demmel1997).
 * @file divide_et_impera.h
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#ifndef _DIVIDEETIMPERA_H_
#define _DIVIDEETIMPERA_H_

#include "../st_matrix.h"

/**
 * Computes eigenvalues of a symmetric tridiagonal matrix.
 * Uses the divide et impera algorithm.
 * @param[in] M            Symmetric tridiagonal matrix
 * @param[out] eigenvalues Array of eigenvalues
 */
void divide_et_impera(const st_matrix_t M, double *eigenvalues);

#endif  /* _DIVIDEETIMPERA_H_ */
