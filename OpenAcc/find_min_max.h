#ifndef FIND_MIN_MAX_H_
#define FIND_MIN_MAX_H_

#include "./struct_c_def.h"
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include "openacc.h"
#include <stdio.h>

// find the maximum eigenvalue of the fermion matrix
// use loc_h, loc_p, loc_r
double ker_find_max_eigenvalue_openacc(  __restrict su3_soa * const u,
					 double_soa * const backfield,
					 ferm_param *pars,
					 __restrict vec3_soa * const loc_r,
					 __restrict vec3_soa * const loc_h,
					 __restrict vec3_soa * const loc_p);

// find the minimum eigenvalue of the fermion matrix
// use loc_h, loc_p, loc_r
double ker_find_min_eigenvalue_openacc(  __restrict su3_soa * const u,
					 __restrict double_soa * const backfield,
					 __restrict ferm_param * const pars,
					 __restrict vec3_soa * const loc_r,
					 __restrict vec3_soa * const loc_h,
					 __restrict vec3_soa * const loc_p,
					 double max );


void find_min_max_eigenvalue_soloopenacc(  __restrict su3_soa * const u,
					   double_soa * const backfield,
					   ferm_param *pars,
					   __restrict vec3_soa * const loc_r,
					   __restrict vec3_soa * const loc_h,
					   __restrict vec3_soa * const loc_p1,
					   __restrict vec3_soa * const loc_p2,
					   double *minmax );



#endif
