#ifndef INVERTER_MIXEDP_H_
#define INVERTER_MIXEDP_H_

#include "./struct_c_def.h"
#include "./sp_struct_c_def.h"
#include "../Include/fermion_parameters.h"
#include "./inverter_package.h"


// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

int inverter_mixed_precision( inverter_package ip,
        ferm_param *pars,
        __restrict vec3_soa_f * solution,// single precision output
        __restrict const vec3_soa * in, // non viene aggiornato mai qui dentro
        double res,
        const int  max_cg,
        double shift  );





#endif

