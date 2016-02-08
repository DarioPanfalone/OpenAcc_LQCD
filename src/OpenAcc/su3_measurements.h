#ifndef SU3_MEASUREMENTS_H_
#define SU3_MEASUREMENTS_H_

#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"


void check_unitarity_device( __restrict su3_soa * const u, double * max_unitarity_deviation, double * avg_unitarity_deviation);
void check_unitarity_host( __restrict su3_soa * const u, double * max_unitarity_deviation, double * avg_unitarity_deviation);

double calc_momenta_action( const __restrict thmat_soa * const mom,
			    double_soa * tr_local, const int mu);

double  calc_plaquette_soloopenacc( __restrict  su3_soa * const tconf_acc, __restrict su3_soa * const local_plaqs, dcomplex_soa * const tr_local_plaqs);

double calc_force_norm(const __restrict tamat_soa * tipdot );




#endif
