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
 * This header file provides an interface for a chronometer.
 * @file chronometer.h
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 * @deprecated Chronometer utility will be removed as soon as possible.
 * Please use the stopwatch utility in stopwatch.h instead.
 */
#ifndef _CHRONOMETER_H_
#define _CHRONOMETER_H_

/** Chronometer. */
typedef struct chronometer_s *chronometer_t;


/**
 * Creates a new chronometer.
 * @return A new chronometer
 * @deprecated Chronometer utility will be removed as soon as possible.
 * Please use the stopwatch utility in stopwatch.h instead.
 */
chronometer_t chronometer_create();


/**
 * Deletes a chronometer.
 * @param[out] chronometer Pointer to chronometer
 * @deprecated Chronometer utility will be removed as soon as possible.
 * Please use the stopwatch utility in stopwatch.h instead.
 */
void chronometer_delete(chronometer_t *chronometer);


/**
 * Starts the chronometer.
 * Resets any previous measurement.
 * @param[out] chronometer Chronometer to start
 * @return Chronometer itself
 * @deprecated Chronometer utility will be removed as soon as possible.
 * Please use the stopwatch utility in stopwatch.h instead.
 */
chronometer_t chronometer_start(chronometer_t chronometer);


/**
 * Stops the chronometer and returns elapsed time.
 * Returns time in milliseconds.
 * @param[out] chronometer Chronometer to stop
 * @return Elapsed time
 * @deprecated Chronometer utility will be removed as soon as possible.
 * Please use the stopwatch utility in stopwatch.h instead.
 */
double chronometer_stop(chronometer_t chronometer);

#endif  /* _CHRONOMETER_H_ */
