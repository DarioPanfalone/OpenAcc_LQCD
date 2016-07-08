#ifndef INVERTER_MULTISHIFT_WRAPPER_H_
#define INVERTER_MULTISHIFT_WRAPPER_H_

#include "./inverter_package.h"
#include "../Include/fermion_parameters.h"
#include "../RationalApprox/rationalapprox.h"
#include "./struct_c_def.h"

void inverter_multishift_wrapper(inverter_package ip,
        ferm_param *pars,
        RationalApprox * approx,
        vec3_soa * out,
        vec3_soa * in,
        double res,
        int max_cg);


#endif
