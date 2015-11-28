#ifndef RATIONAL_APPROX_H_
#define RATIONAL_APPROX_H_

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
void rationalapprox_read_custom_nomefile(RationalApprox* rational_approx, char* nomefile);


void rationalapprox_save(const char* nomefile, RationalApprox* rational_approx);

// NOTE: this function assumes that in->lambda_max is 1.
// during the rescaling, only max = minmax[1] is considered, but a check 
// is performed to see if min = minmax[0] is less than the out->lambda_min,
// in which case a warning will be issued.
void rescale_rational_approximation(RationalApprox *in, RationalApprox *out, double *minmax);
// this function rescales the rational approximation 
// to obtain out->lambda_max =1
void renormalize_rational_approximation(RationalApprox *in, RationalApprox *out);



#endif
