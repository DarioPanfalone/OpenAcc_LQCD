#ifndef FIELD_CORR_H_
#define FIELD_CORR_H_

#include "./geometry.h"
#include "./su3_utilities.h"

void calc_field_corr(
										 __restrict const su3_soa * const u,
										 __restrict su3_soa * const field_corr,
										 __restrict su3_soa * const field_corr_aux,
										 __restrict su3_soa * const loc_plaq,
										 dcomplex_soa * const trace_local,
										 d_complex * const corr,
										 __restrict single_su3 * const closed_corr,
										 const int mu, const int nu, const int ro);

#endif
