#ifndef FERMION_MATRIX_H
#define FERMION_MATRIX_H


#include "./struct_c_def.h"
#ifndef __GNUC__
 #include "openacc.h"
#endif
#include "./fermionic_utilities.h"
#include "../Include/fermion_parameters.h"
#include "../DbgTools/debug_macros_glvarcheck.h"


// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif



void acc_Deo( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in,__restrict ferm_param *pars,__restrict double_soa * backfield);
void acc_Doe( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in,__restrict ferm_param *pars,__restrict double_soa * backfield); 

void fermion_matrix_multiplication( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in, __restrict vec3_soa * const temp1, __restrict ferm_param *pars, __restrict double_soa * backfield);

void fermion_matrix_multiplication_shifted( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in, __restrict vec3_soa * const temp1, __restrict ferm_param *pars,__restrict double_soa * backfield, double shift);







#endif
