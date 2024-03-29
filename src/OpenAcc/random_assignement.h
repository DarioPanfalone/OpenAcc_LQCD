#ifndef RANDOM_ASSIGNEMENT_H_
#define RANDOM_ASSIGNEMENT_H_


#include "../OpenAcc/struct_c_def.h"
#include "./single_types.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

double double_gauss(void);
void two_double_gauss(double *r);
d_complex d_complex_gauss(void);

void generate_vec3_soa_gauss(__restrict vec3_soa * const vect);

void generate_Momenta_gauss_comp(__restrict thmat_soa * const mom);

void generate_random_su3(single_su3 * m, double f);

void generate_Conf_cold_comp(__restrict su3_soa * const conf,double factor);


void generate_Momenta_gauss(__restrict thmat_soa * const mom);

void generate_Conf_cold(__restrict su3_soa * const conf,double factor);

void generate_vec3_soa_z2noise(__restrict vec3_soa * const vect);


#endif
