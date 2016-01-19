#include <stdio.h>
#include "root-finding.h"


double y_x(double x, void *p) {
    return x;
}

double d1(double x, void *p) {
    return 1.0;
}

double y_x2_3x_n2(double x, void *p) {
    return x * x + 3 * x - 2;
}

double d2(double x, void *p) {
    return 2.0 * x + 3.0;
}



int main() {
    printf("Bisection\n");
    printf("Zero of y = x: %g\n", bisection(-25, +10, y_x, NULL));
    printf("Zero of y = x^2 + 3x - 2: %g\n", bisection(0, 1, y_x2_3x_n2, NULL));

    printf("Regula falsi\n");
    printf("Zero of y = x: %g\n", regula_falsi(-20, +10, y_x, NULL));
    printf("Zero of y = x^2 + 3x - 2: %g\n", regula_falsi(0, 1, y_x2_3x_n2, NULL));

    printf("Newton\n");
    printf("Zero of y = x: %g\n", newton(-20, y_x, d1, NULL));
    printf("Zero of y = x^2 + 3x - 2: %g\n", newton(1, y_x2_3x_n2, d2, NULL));

    return 0;
}
