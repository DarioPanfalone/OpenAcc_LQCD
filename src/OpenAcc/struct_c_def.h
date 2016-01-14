#ifndef STRUCT_C_DEF_H_
#define STRUCT_C_DEF_H_

#include "./geometry.h"
#include <complex.h>

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif


typedef double complex  d_complex;
typedef float  complex  f_complex;

typedef struct vec3_soa_t {
    int status;
    d_complex c0[sizeh];
    d_complex c1[sizeh];
    d_complex c2[sizeh];
} vec3_soa;
typedef struct dcomplex_soa_t {
    int status;
    d_complex c[sizeh];
} dcomplex_soa;
typedef struct double_soa_t {
    int status;
    double d[sizeh];
    int staggered_status; // 0 ok, 1 with staggered phases included.
} double_soa;
typedef struct vec3_t {
  d_complex c0;
  d_complex c1;
  d_complex c2;
} vec3;
typedef struct su3_soa_t {
  int status;
  int staggered_status; // 0 ok, 1 with staggered phases included.
  vec3_soa r0;
  vec3_soa r1;
  vec3_soa r2;
} su3_soa;
typedef struct tamat_soa_t {
  int status;
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double rc00[sizeh];   // Im(comp_00)
  double rc11[sizeh];   // Im(comp_11)
} tamat_soa;
typedef struct thmat_soa_t {
  int status;
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double rc00[sizeh];   // Re(comp_00)
  double rc11[sizeh];   // Re(comp_11)
} thmat_soa;



#endif


