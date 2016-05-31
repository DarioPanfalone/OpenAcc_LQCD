#ifndef INVERTER_SHIFT_MULTI_FULL_H_
#define INVERTER_SHIFT_MULTI_FULL_H_

#include "./struct_c_def.h"
#include "../RationalApprox/rationalapprox.h"
#include "../Include/fermion_parameters.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

extern int multishift_invert_iterations ; // global count of multishift CG iterations

int multishift_invert(__restrict su3_soa * const u,
		      __restrict ferm_param * pars,
		      RationalApprox * approx,
		      __restrict vec3_soa * out, // multi-fermion [nshifts]
		      __restrict vec3_soa * const  in, // single ferm
		      double residuo,
		      __restrict vec3_soa *  loc_r,
		      __restrict vec3_soa *  loc_h,
		      __restrict vec3_soa *  loc_s,
		      __restrict vec3_soa *  loc_p,
		      __restrict vec3_soa *  shiftferm, // multi-ferm [nshift]
              const int max_cg);


void recombine_shifted_vec3_to_vec3(const __restrict vec3_soa* const in_shifted /*multi-fermion*/, 
				    const __restrict vec3_soa* const in, // [nshift]
				    __restrict vec3_soa * const out, // [1] 
				    const RationalApprox * const approx );


#endif

