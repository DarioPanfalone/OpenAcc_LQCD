#ifndef FIELD_CORR_H_
#define FIELD_CORR_H_

#include "./geometry.h"
#include "./su3_utilities.h"

void calc_field_corr(
        __restrict const su3_soa * const u,
		__restrict su3_soa * const field_corr,
		__restrict d_complex * const traccia,
        const int mu, const int nu, const int ro);


#endif
