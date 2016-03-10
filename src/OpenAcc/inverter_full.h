#ifndef INVERTER_FULL_H_
#define INVERTER_FULL_H_

#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif


#define max_cg 10000

#define DEBUG_INVERTER_FULL_OPENACC


int ker_invert_openacc(   __restrict su3_soa * const u,  // non viene aggiornata mai qui dentro
			  ferm_param *pars,
			  __restrict vec3_soa * const out,
			  __restrict vec3_soa * const in, // non viene aggiornato mai qui dentro
			  double res,
			  __restrict vec3_soa * const trialSolution, // non viene aggiornato mai qui dentro
			  __restrict vec3_soa * const loc_r,
			  __restrict vec3_soa * const loc_h,
			  __restrict vec3_soa * const loc_s,
			  __restrict vec3_soa * const loc_p
			  );

#endif

