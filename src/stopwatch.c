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
 * Stopwatch utility.
 * Offers basic functions to measure time.
 * If HPC macro is defined, uses HPC library instead of default
 * implementation.
 * @file stopwatch.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "stopwatch.h"

#include <stdio.h>

#ifdef HPC
#include <libhpc.h>
#else
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

#define NAMELENGTH 32      /**< Maximum length of a name */
#define FILENAMELENGTH 64  /**< Maximum length of a filename */
#define MEASURES 128       /**< Maximum number of measures */

/**
 * Structure of a measure.
 */
struct measure_s {
    char name[NAMELENGTH];  /**< Name of the measure */
    unsigned int used;      /**< Tells whether the measure was used */
    struct timeval start;   /**< Start time */
    struct timeval stop;    /**< Stop time */
};


/**
 * Structure of a stopwatch.
 */
struct stopwatch_s {
    char name[NAMELENGTH];                /**< Name of the stopwatch */
    struct measure_s measures[MEASURES];  /**< Measure in the stopwatch */
};
#endif  /* HPC */



/**
 * Creates a new stopwatch.
 * @param[in] task_name Name of the task to measure
 * @return A new stopwatch
 */
stopwatch_t stopwatch_create(const char *task_name) {
#ifdef HPC
    hpmInit(0, task_name);
    return NULL;
#else
    unsigned int i;
    stopwatch_t sw;

    SAFE_MALLOC(sw, stopwatch_t, sizeof(struct stopwatch_s));
    strncpy(sw->name, task_name, NAMELENGTH);

    /* Marks every measure as not used */
    for (i = 0; i < MEASURES; ++i) {
        sw->measures[i].used = 0;
    }

    return sw;
#endif
}



/**
 * Deletes a stopwatch.
 * Saves measures in the stopwatch to a file.
 * @param[out] S Address of the stopwatch to delete
 */
void stopwatch_delete(stopwatch_t *S) {
#ifdef HPC
    (void) S;
    hpmTerminate(0);
#else
    unsigned int i;
    FILE *fp;
    char filename[FILENAMELENGTH];
    struct measure_s *m;

    /* Returns if stopwatch is invalid */
    if (NULL == S || NULL == *S) {
        return;
    }

    /* Opens file to write */
    sprintf(filename, "%s.dat", (*S)->name);
    fp = fopen(filename, "w");
    if (NULL == fp) {
        return;
    }

    /* Writes every measure */
    for (i = 0; i < MEASURES; ++i) {
        double d_t;
        m = (*S)->measures + i;

        if (!m->used) continue;

        d_t = (m->stop.tv_sec  - m->start.tv_sec)  * 1e+3
            + (m->stop.tv_usec - m->start.tv_usec) * 1e-3;

        fprintf(fp, "[%u] %s: %g ms\n", i, m->name, d_t);
    }

    /* Memory free */
    fclose(fp);
    free(*S);
    *S = NULL;
#endif
}



/**
 * Starts a measure.
 * @param[in, out] S    Stopwatch to use
 * @param[in]      id   Identifier of the measure
 * @param[in]      name Name of the measure
 */
void stopwatch_start(stopwatch_t S, const unsigned int id, const char *name) {
#ifdef HPC
    (void) S;
    hpmStart(id, name);
#else
    struct measure_s *m;

    /* Returns if stopwatch or id are invalid */
    if (NULL == S || id >= MEASURES) {
        return;
    }

    m = S->measures + id;
    m->used = 1;
    strncpy(m->name, name, NAMELENGTH);
    gettimeofday(&(m->start), NULL);
#endif
}



/**
 * Stops a measure.
 * @param[in, out] S  Stopwatch to use
 * @param[in]      id Identifier of the measure to stop
 */
void stopwatch_stop(stopwatch_t S, const unsigned int id) {
#ifdef HPC
    (void) S;
    hpmStop(id);
#else
    struct measure_s *m;

    /* Returns if stopwatch or id are invalid */
    if (NULL == S || id >= MEASURES) {
        return;
    }

    m = S->measures + id;
    m->used = 1;
    gettimeofday(&(m->stop), NULL);
#endif
}
