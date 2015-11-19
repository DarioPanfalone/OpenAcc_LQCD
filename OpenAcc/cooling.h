#ifndef COOLING_C
#define COOLING_C

#include "./struct_c_def.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/su3_utilities.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

void compute_cooled_even_links(__restrict su3_soa   * const U, __restrict su3_soa   * const STAP);

void compute_cooled_odd_links(__restrict su3_soa   * const U,
			      __restrict su3_soa   * const STAP);

void cool_conf(__restrict su3_soa   * const U,
	       __restrict su3_soa   * const TMP);

#endif
