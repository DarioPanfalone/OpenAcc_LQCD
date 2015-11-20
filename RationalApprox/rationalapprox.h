#ifndef RATIONAL_APPROX_H_
#define RATIONAL_APPROX_H_

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define MAX_APPROX_ORDER 19


typedef struct RationalApprox_t{
    int exponent_num;
    int exponent_den;
    int approx_order;
    double lambda_min;
    double lambda_max; // usually equal to one
    int gmp_remez_precision;
    double error;
    double RA_a0;
    double RA_a[MAX_APPROX_ORDER];
    double RA_b[MAX_APPROX_ORDER];
} RationalApprox;




char* rational_approx_filename(int approx_order, int exponent_num, int exponent_den, double lambda_min);

void rationalapprox_read(RationalApprox* rational_approx);


void rationalapprox_save(const char* nomefile, RationalApprox* rational_approx);


void rescale_rational_approximation(RationalApprox *in, RationalApprox *out, double *minmax);



#endif
