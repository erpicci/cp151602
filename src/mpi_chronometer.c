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
 * Chronometer utility.
 * Implements a chronometer suitable for an MPI environment.
 * @file mpi_chronometer.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "chronometer.h"

#include <mpi.h>

#include "utils.h"

/** Structure of an MPI chronometer. */
struct chronometer_s {
    double start;   /**< Start time */
    double finish;  /**< Finish time */
};



chronometer_t chronometer_create() {
    chronometer_t chronometer;

    SAFE_MALLOC(chronometer, chronometer_t, sizeof(struct chronometer_s));

    return chronometer;
}



void chronometer_delete(chronometer_t *chronometer) {
    if (NULL == chronometer || NULL == *chronometer) {
        return;
    }

    free(*chronometer);
    *chronometer = NULL;
}



chronometer_t chronometer_start(chronometer_t chronometer) {
    double time;

    if (NULL == chronometer) {
        return NULL;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    time = MPI_Wtime();
    MPI_CALL(MPI_Reduce(&time, &(chronometer->start), 1, MPI_DOUBLE,
                        MPI_MIN, 0, MPI_COMM_WORLD));

    return chronometer;
}



double chronometer_stop(chronometer_t chronometer) {
    double time;

    if (NULL == chronometer) {
        return 0.0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    time = MPI_Wtime();
    MPI_CALL(MPI_Reduce(&time, &(chronometer->finish), 1, MPI_DOUBLE,
                        MPI_MAX, 0, MPI_COMM_WORLD));

    return (chronometer->finish - chronometer->start) * 1e3;
}
