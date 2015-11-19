#ifndef INVERTER_SHIFT_MULTI_FULL_C_
#define INVERTER_SHIFT_MULTI_FULL_C_

#include "./struct_c_def.h"
#include "./fermion_matrix.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif


int multishift_invert(__restrict su3_soa * const u,
		      ferm_param * pars,
		      RationalApprox * approx,
		      double_soa * backfield,
		      __restrict vec3_soa * const out, // multi-fermion [nshifts]
		      __restrict vec3_soa * const in, // single ferm
		      double residuo,
		      __restrict vec3_soa * const loc_r,
		      __restrict vec3_soa * const loc_h,
		      __restrict vec3_soa * const loc_s,
		      __restrict vec3_soa * const loc_p,
		      __restrict vec3_soa * const shiftferm // multi-ferm [nshift]
		      );


void recombine_shifted_vec3_to_vec3(const __restrict vec3_soa* const in_shifted /*multi-fermion*/, 
				    const __restrict vec3_soa* const in, // [nshift]
				    __restrict vec3_soa * const out, // [1] 
				    const RationalApprox * const approx );


#endif

