#ifndef INVERTER_WRAPPERS_H_
#define INVERTER_WRAPPERS_H_

#include "./inverter_package.h"
#include "../Include/fermion_parameters.h"
#include "../RationalApprox/rationalapprox.h"
#include "./struct_c_def.h"



#define CONVERGENCE_CRITICAL 1
#define CONVERGENCE_NONCRITICAL 0

extern int multishift_invert_iterations ; // global count of multishift CG iterations

void convergence_messages(int conv_importance, int inverter_status);

int inverter_multishift_wrapper(inverter_package ip,
																ferm_param *pars,
																RationalApprox * approx,
																vec3_soa * out,
																const vec3_soa * in,
																double res,
																int max_cg,
																int convergence_importance);

int inverter_wrapper(inverter_package ip,
										 ferm_param *pars,
										 vec3_soa * out,
										 const vec3_soa * in,
										 double res,
										 int max_cg,
										 double shift,
										 int convergence_importance);


#endif
