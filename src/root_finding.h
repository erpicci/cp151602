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
 * Computes roots of a function.
 * @file root_finding.h
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#ifndef _ROOTFINDING_H_
#define _ROOTFINDING_H_

/**
 * Bisection method.
 * @param[in] min  Minimum value
 * @param[in] max  Maximum value
 * @param[in] f    Function to evaluate
 * @param[in] args Additional arguments to pass to f
 * @return Root of f between min and max
 */
double bisection(double min, double max, double (*f)(double, void *), void *args);


/**
 * Regula falsi method.
 * @param[in] min  Minimum value
 * @param[in] max  Maximum value
 * @param[in] f    Function to evaluate
 * @param[in] args Additional arguments to pass to f
 * @return Root of f between min and max
 */
double regula_falsi(double min, double max, double (*f)(double, void*), void *args);


/**
 * Newton method.
 * @param[in] x    Starting point
 * @param[in] f    Function to evaluate
 * @param[in] fp   First derivative of f
 * @param[in] args Additional arguments to pass to f and fp
 * @return Root of f between min and max
 */
double newton(double x, double (*f)(double, void *), double (*fp)(double, void *), void *args);

#endif  /* _ROOTFINDING_H_ */
