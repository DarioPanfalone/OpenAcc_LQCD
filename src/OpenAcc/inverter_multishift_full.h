#ifndef INVERTER_MULTISHIFT_FULL_H_
#define INVERTER_MULTISHIFT_FULL_H_

#include "./struct_c_def.h"
#include "../RationalApprox/rationalapprox.h"
#include "../Include/fermion_parameters.h"
#include "./inverter_wrappers.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif


int multishift_invert(__restrict const su3_soa * u,
											__restrict ferm_param * pars,
											RationalApprox * approx,
											__restrict vec3_soa * out, // multi-fermion [nshifts]
											__restrict const vec3_soa * in, // single ferm
											double residuo,
											__restrict vec3_soa *  loc_r,
											__restrict vec3_soa *  loc_h,
											__restrict vec3_soa *  loc_s,
											__restrict vec3_soa *  loc_p,
											__restrict vec3_soa *  shiftferm, // multi-ferm [nshift]
											const int max_cg,
											int * cg_return);


void recombine_shifted_vec3_to_vec3(const __restrict vec3_soa* in_shifted, // multi-fermion 
																		const __restrict vec3_soa* in, // [nshift]
																		__restrict vec3_soa *  out, // [1]
																		const RationalApprox * approx );


#endif
