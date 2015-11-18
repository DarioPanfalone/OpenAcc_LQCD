#ifndef FERMIONIC_UTILITIES_H_
#define FERMIONIC_UTILITIES_H_

#include "../DbgTools/debug_macros_glvarcheck.h"

#pragma acc routine seq
static inline double scal_prod_loc_1double(  __restrict vec3_soa * const in_vect1,
					     __restrict vec3_soa * const in_vect2,
					     const int idx_vect
                                          ) {

  double sum  = creal(in_vect1->c0[idx_vect]) * creal(in_vect2->c0[idx_vect]) + cimag(in_vect1->c0[idx_vect]) * cimag(in_vect2->c0[idx_vect]);
         sum += creal(in_vect1->c1[idx_vect]) * creal(in_vect2->c1[idx_vect]) + cimag(in_vect1->c1[idx_vect]) * cimag(in_vect2->c1[idx_vect]);
         sum += creal(in_vect1->c2[idx_vect]) * creal(in_vect2->c2[idx_vect]) + cimag(in_vect1->c2[idx_vect]) * cimag(in_vect2->c2[idx_vect]);
  return sum;
}


#pragma acc routine seq
static inline double l2norm2_loc( __restrict vec3_soa * const in_vect1,
				  const int idx_vect
				  ) {
  double sum = creal(in_vect1->c0[idx_vect])*creal(in_vect1->c0[idx_vect])+cimag(in_vect1->c0[idx_vect])*cimag(in_vect1->c0[idx_vect]);
  sum += creal(in_vect1->c1[idx_vect])*creal(in_vect1->c1[idx_vect])+cimag(in_vect1->c1[idx_vect])*cimag(in_vect1->c1[idx_vect]);
  sum += creal(in_vect1->c2[idx_vect])*creal(in_vect1->c2[idx_vect])+cimag(in_vect1->c2[idx_vect])*cimag(in_vect1->c2[idx_vect]);
  return sum;
}



d_complex scal_prod_global(  const __restrict vec3_soa * const in_vect1, const __restrict vec3_soa * const in_vect2 );
double real_scal_prod_global(  __restrict vec3_soa * const in_vect1, __restrict vec3_soa * const in_vect2 );
double l2norm2_global( __restrict  vec3_soa * const in_vect1);

void combine_in1xfactor_plus_in2(  __restrict vec3_soa * const in_vect1, const double factor, __restrict vec3_soa * const in_vect2, __restrict vec3_soa * const out);

void multiply_fermion_x_doublefactor(  __restrict vec3_soa * const in1, const double factor );

void combine_add_factor_x_in2_to_in1( __restrict vec3_soa * const in1,__restrict vec3_soa * const in2, double factor);

void combine_in1xferm_mass2_minus_in2_minus_in3(   __restrict vec3_soa * const in_vect1, double ferm_mass, __restrict vec3_soa * const in_vect2, __restrict vec3_soa * const in_vect3, __restrict vec3_soa * const out);

void combine_inside_loop( __restrict vec3_soa * const vect_out,  __restrict vec3_soa * const vect_r, const __restrict vec3_soa * const vect_s, const __restrict vec3_soa * const vect_p, const double omega);

void combine_in1xferm_mass_minus_in2(__restrict vec3_soa * const in_vect1,double ferm_mass2, __restrict vec3_soa * const in_vect2);

void combine_in1_minus_in2( __restrict vec3_soa * const in_vect1, __restrict vec3_soa * const in_vect2, __restrict vec3_soa * const out);

void assign_in_to_out(  __restrict vec3_soa * const in_vect1, __restrict vec3_soa * const out);

// Altre funzioni aggiunte per l'algebra lineare multishift "versatilizzato"
void set_vec3_soa_to_zero( __restrict vec3_soa* const fermion);
void multiple_combine_in1_minus_in2x_factor_back_into_in1( __restrict vec3_soa * const out, __restrict vec3_soa * const in, const int maxiter, __restrict const int * const flag, __restrict const double * const omegas);
void multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1( __restrict vec3_soa * const in1, int maxiter, __restrict const int * const flag, __restrict const double * const gammas, __restrict vec3_soa * const in2, __restrict const double * const zeta_iii);
void combine_in1_x_fact1_minus_in2_back_into_in2( __restrict vec3_soa * const in1, double fact1, __restrict vec3_soa * const in2 );

void combine_in1_minus_in2_allxfact( __restrict vec3_soa * const in1, __restrict vec3_soa * const in2, double fact, __restrict vec3_soa * const out ); 


#endif
