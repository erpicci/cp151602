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
 * @file stopwatch.h
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

/** Type of a stopwatch. */
typedef struct stopwatch_s *stopwatch_t;


/**
 * Creates a new stopwatch.
 * @param[in] task_name Name of the task to measure
 * @return A new stopwatch
 */
stopwatch_t stopwatch_create(const char *task_name);


/**
 * Deletes a stopwatch.
 * Saves measures in the stopwatch to a file.
 * @param[out] S Address of the stopwatch to delete
 */
void stopwatch_delete(stopwatch_t *S);


/**
 * Starts a measure.
 * @param[in, out] S    Stopwatch to use
 * @param[in]      id   Identifier of the measure
 * @param[in]      name Name of the measure
 */
void stopwatch_start(stopwatch_t S, const unsigned int id, const char *name);


/**
 * Stops a measure.
 * @param[in, out] S  Stopwatch to use
 * @param[in]      id Identifier of the measure to stop
 */
void stopwatch_stop(stopwatch_t S, const unsigned int id);

#endif  /* _STOPWATCH_H_ */
