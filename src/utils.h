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
 * General utility macros.
 * Offers additional macros to interact with MPI APIs.
 * @file utils.h
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/**
 * Calls an MPI API and checks its return value.
 * If MPI API returns an error value (not equal to MPI_SUCCESS), prints
 * an error message to standard error, finalizes MPI and aborts the
 * program with EXIT_FAILURE.
 * @param[in] call MPI API call
 */
#define MPI_CALL(call)                                         \
{                                                              \
    int err = call;                                            \
    if (MPI_SUCCESS != err) {                                  \
        char msg[MPI_MAX_ERROR_STRING];                        \
        int len;                                               \
        MPI_Error_string(err, msg, &len);                      \
        fprintf(stderr, "MPI Error (%d): %s [%s in %s:%d]\n",  \
                err, msg, #call, __FILE__, __LINE__);          \
        MPI_Finalize();                                        \
        exit(EXIT_FAILURE);                                    \
    }                                                          \
}



/**
 * Allocates memory and checks for error.
 * If malloc fails, prints and error and aborts the program.
 * @param[out] ptr  Pointer to allocated memory
 * @param[in]  type Type of pointer
 * @param[in]  size of memory to allocate (in bytes)
 */
#define SAFE_MALLOC(ptr, type, size)                         \
{                                                            \
    ptr = (type) malloc(size);                               \
    if (NULL == ptr) {                                       \
        fprintf(stderr, "Cannot allocate memory [%s:%d]\n",  \
                __FILE__, __LINE__);                         \
        exit(EXIT_FAILURE);                                  \
    }                                                        \
}

#endif  /* _MPI_UTILS_H_ */
