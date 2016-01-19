#include <math.h>

double bisection(double min, double max, double (*f)(double, void *), void *args) {
    const double TOL = 1e-6;
    double f_x, x;
    
    while (max - min > TOL) {
        x   = (max + min) * 0.5;
        f_x = (*f)(x, args);
        
        if (f_x > TOL) {
            max = x;
        }
        else if (f_x < -TOL) {
            min = x;
        }
        else {
            break;
        }
    }
    
    return x;
}



double regula_falsi(double min, double max, double (*f)(double, void*), void *args) {
    const double TOL = 1e-6;
    int side = 0;
    double x, f_x,
           f_min = (*f)(min, args),
           f_max = (*f)(max, args);

    while (1) {
        x   = (f_min * max - f_max * min) / (f_min - f_max);
        f_x = (*f)(x, args);
        if (fabs(max - min) < TOL * fabs(max + min)) {
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



double newton(double x, double (*f)(double, void *), double (*fp)(double, void *), void *args) {
    const double TOL = 1e-6;
    double f_x  = (*f)(x, args),
           fp_x = (*fp)(x, args);
    
    while (fabs(f_x) > TOL && fabs(fp_x) > TOL) {
        x    = x - f_x / fp_x;
        f_x  = (*f)(x, args);
        fp_x = (*fp)(x, args);
    }
    
    return x;
}