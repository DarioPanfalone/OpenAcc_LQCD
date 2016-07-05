#ifndef INVERTER_MIXEDP_H_
#define INVERTER_MIXEDP_H_

#include "./struct_c_def.h"
#include "./sp_struct_c_def.h"
#include "../Include/fermion_parameters.h"


// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

int inverter_mixed_precision( __restrict const su3_soa * ud, // non viene aggiornata mai qui dentro
              __restrict const su3_soa_f * u, // non viene aggiornata mai qui dentro
			  ferm_param *pars,
			  __restrict vec3_soa_f * solution,// single precision output
			  __restrict const vec3_soa * in, // non viene aggiornato mai qui dentro
			  double res,
			  __restrict vec3_soa   * loc_r,
			  __restrict vec3_soa_f * loc_r_f,   
			  __restrict vec3_soa_f * loc_h_f,   
			  __restrict vec3_soa_f * loc_s_f,   
			  __restrict vec3_soa_f * loc_p_f,   
			  __restrict vec3_soa   * t1,      //
			  __restrict vec3_soa   * t2,
              const int  max_cg,
              double shift  );





#endif

