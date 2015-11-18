#ifndef FERMION_MATRIX_H
#define FERMION_MATRIX_H

#include "./struct_c_def.h"
#include "openacc.h"
#include "./fermionic_utilities.h"
#include "../DbgTools/debug_macros_glvarcheck.h"


void acc_Deo( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in,ferm_param *pars,double_soa * backfield);
void acc_Doe( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in,ferm_param *pars,double_soa * backfield); 
inline void fermion_matrix_multiplication( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in, __restrict vec3_soa * const temp1, ferm_param *pars,double_soa * backfield){
  SETREQUESTED(temp1);
  SETREQUESTED(out);
  acc_Doe(u,temp1,in,pars,backfield);
  acc_Deo(u,out,temp1,pars,backfield);
  combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass,out);// Nuova funzione in OpenAcc/fermionic_utilities.c
  SETFREE(temp1);
}
inline void fermion_matrix_multiplication_shifted( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in, __restrict vec3_soa * const temp1, ferm_param *pars,double_soa * backfield, double shift){
  SETREQUESTED(temp1);
  SETREQUESTED(out);
  acc_Doe(u,temp1,in,pars,backfield);
  acc_Deo(u,out,temp1,pars,backfield);
  combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass+shift,out);// Nuova funzione in OpenAcc/fermionic_utilities.c
  SETFREE(temp1);


}







#endif
