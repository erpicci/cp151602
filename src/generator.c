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
 * Instance generator.
 * Generates symmetric tridiagonal matrices.
 * @file generator.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include "utils.h"
#include "st_matrix.h"


/**
 * Shows helper.
 * Prints helper on standard output.
 * @param[in] argc ARGument Counter
 * @param[in] argv ARGument Vector
 */
static void show_helper(const int argc, char * const argv[]) {
    (void) argc;

    printf("INSTANCE GENERATOR\n");
    printf("Use this tool to generate an instance of symmetric, tridiagonal matrix.\n");
    printf("Instance is print to standard output.\n");
    printf("Instance is generated randomly, considering an uniform distribution of random values\n");
    printf("\n");
    printf("Usage:\n");
    printf("  %s -n <int> -m <num> -M <num> -h\n", argv[0]);
    printf("\n");
    printf("Options:\n");
    printf("  -n <int>\t Size of the matrix (default: 5)\n");
    printf("  -m <num>\t Minimum allowed value (default: -1.0)\n");
    printf("  -M <num>\t Maximum allowed value (default: +1.0)\n");
    printf("  -h      \t Prints this helper and exit\n");
}



/**
 * Generates an instance of symmetric, tridiagonal matrix.
 * Generates a symmetric tridiagonal matrix and prints it to standard
 * output. Reads matrix size from standard input.
 * @param[in] argc ARGument Counter
 * @param[in] argv ARGument Vector
 * @retval EXIT_SUCCESS Normal termination of the program
 * @retval EXIT_FAILURE Some error occurred
 */
int main(const int argc, char * const argv[]) {
    int size = 5, opt, i;
    double min = -1.0, max = +1.0, *diag, *subdiag;
    st_matrix_t M;


    /*******************************************************************
     * Command line options parsing.
     ******************************************************************/
    while ((opt = getopt(argc, argv, "n:m:M:h")) != -1) {
        switch (opt) {
        case 'n': size = atoi(optarg); break;
        case 'm': min  = atof(optarg); break;
        case 'M': max  = atof(optarg); break;
        case 'h':
           show_helper(argc, argv);
           exit(EXIT_SUCCESS);
           break;
        default:
           fprintf(stderr, "Unrecognized option. Run with -h to see the helper.\n");
           exit(EXIT_FAILURE);
        }
    }


    /*******************************************************************
     * Generates a symmetric tridiagonal matrix randomly.
     ******************************************************************/
    srand(time(NULL) + getpid());

    M       = st_matrix_create(size);
    diag    = st_matrix_diag(M);
    subdiag = st_matrix_subdiag(M);

    for (i = 0; i < size - 1; ++i) {
        diag[i]    = (double) rand() / RAND_MAX * (max - min) + min;
        subdiag[i] = (double) rand() / RAND_MAX * (max - min) + min;
    }
    diag[i] = (double) rand() / RAND_MAX * (max - min) + min;


    /*******************************************************************
     * Saves matrix on standard output.
     ******************************************************************/
    st_matrix_save(M, stdout);


    /* Frees memory. */
    st_matrix_delete(&M);

    return EXIT_SUCCESS;
}
