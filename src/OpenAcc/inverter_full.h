#ifndef INVERTER_FULL_H_
#define INVERTER_FULL_H_

#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

#define INVERTER_SUCCESS 1
#define INVERTER_FAILURE 0



#define DEBUG_INVERTER_FULL_OPENACC


int ker_invert_openacc(   __restrict const su3_soa * u,  // non viene aggiornata mai qui dentro
			  ferm_param *pars,
			  __restrict vec3_soa * solution,
			  __restrict const vec3_soa * in, // non viene aggiornato mai qui dentro
			  double res,
			  __restrict vec3_soa * loc_r,
			  __restrict vec3_soa * loc_h,
			  __restrict vec3_soa * loc_s,
			  __restrict vec3_soa * loc_p,
              const int max_cg,
              double shift,
              int * cg_return );


#endif

