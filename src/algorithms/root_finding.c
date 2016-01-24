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
 * @file root_finding.c
 * @author Erika Fabris <fabriser@dei.unipd.it>
 * @author Thomas Gagliardi <gagliard@dei.unipd.it>
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @copyright GNU GPLv3 <http://www.gnu.org/licenses/gpl-3.0.txt>
 */
#include "root_finding.h"

#include <math.h>

/** Allowed tolerance. */
static const double EPS = 1e-6;



double bisection(
    double min, double max,
    double (*f)(const double, const void *), const void *args) {
    double f_x, x;

    while (max - min > EPS) {
        x   = (max + min) * 0.5;
        f_x = (*f)(x, args);

        if (f_x > EPS) {
            max = x;
        }
        else if (f_x < -EPS) {
            min = x;
        }
        else {
            break;
        }
    }
    
    return x;
}



double regula_falsi(
    double min, double max,
    double (*f)(const double, const void *), const void *args) {
    int side = 0;
    double x, f_x,
           f_min = (*f)(min, args),
           f_max = (*f)(max, args);

    while (1) {
        x   = (f_min * max - f_max * min) / (f_min - f_max);
        f_x = (*f)(x, args);
        if (fabs(max - min) < EPS * fabs(max + min)) {
            break;
        }

        if (f_x * f_max > 0) {
            max   = x;
            f_max = f_x;
            if (side == -1) {
                f_min *= 0.5;
            }
            side = -1;
        }
        else if (f_min * f_x > 0) {
            min   = x; 
            f_min = f_x;
            if (side == +1) {
                f_max *= 0.5;
            }
            side = +1;
        }
        else {
            break;
        } 
    }
    return x;
}



double newton(
    double x,
    double (*f)(const double, const void *),
    double (*fp)(const double, const void *),
    const void *args) {
    double f_x  = (*f)(x, args),
           fp_x = (*fp)(x, args);
    
    while (fabs(f_x) > EPS && fabs(fp_x) > EPS) {
        x    = x - f_x / fp_x;
        f_x  = (*f)(x, args);
        fp_x = (*fp)(x, args);
    }
    
    return x;
}
