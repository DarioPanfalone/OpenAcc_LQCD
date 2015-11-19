#ifndef CAYLEY_HAMILTON_H
#define CAYLEY_HAMILTON_H


// This is based on http://arXiv.org/abs/hep-lat/0311018v1,
// "Analytic Smearing of SU(3) Link Variables in Lattice QCD",
//  Morningstar & Peardon (2008)

#ifndef __GNUC__
 #include <accelmath.h>
#endif
#include <complex.h>
#include "./struct_c_def.h"
#include "./single_types.h"
#include "../DbgTools/debug_macros_glvarcheck.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif




/*
inline void Itamat_2ndDeg_poly_no3rdrow(d_complex f0, d_complex f1, d_complex f2, __restrict single_tamat * const QA, single_su3 * const out);

inline void Itamat_2ndDeg_poly(d_complex f0, d_complex f1, d_complex f2, single_tamat * QA, single_su3 * out);
*/




static inline double det_i_times_QA_soa( __restrict tamat_soa * const QA,const int idx);

static inline double Tr_i_times_QA_sq_soa( __restrict tamat_soa * const QA,const int idx);

static inline void CH_exponential_antihermitian_soa_nissalike(__restrict su3_soa * const exp_out,__restrict tamat_soa * const QA,const int idx);

static inline void CH_exponential_antihermitian_nissalike(single_su3 * const exp_out, __restrict single_tamat * const QA);
#endif
