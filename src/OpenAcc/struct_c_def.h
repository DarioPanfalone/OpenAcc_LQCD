#ifndef STRUCT_C_DEF_H_
#define STRUCT_C_DEF_H_

#include "./geometry.h"
#include <complex.h>

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

#include "./double_complex.h"

//TYPEDEF_FLOAT_COMPLEX  //leave this there for double_to_songle_transformer.py

typedef struct vec3_soa_t {
    d_complex c0[sizeh];
    d_complex c1[sizeh];
    d_complex c2[sizeh];
} vec3_soa;
typedef struct dcomplex_soa_t {
    d_complex c[sizeh];
} dcomplex_soa;

typedef struct double_soa_t {
    double d[sizeh];
} double_soa;

typedef struct vec3_t {
  d_complex c0;
  d_complex c1;
  d_complex c2;
} vec3;

typedef struct su3_soa_t {
  vec3_soa r0;
  vec3_soa r1;
  vec3_soa r2;
  double_soa K; //Adjoint vector K. This vector is part of the struct and its values directly modify the near defect links computation in the action. Obviously its lenght is sizeh.
  //int label; //the conf label -> moved to <rep_info> struct
} su3_soa;


typedef struct tamat_soa_t {
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double ic00[sizeh];   // Im(comp_00)
  double ic11[sizeh];   // Im(comp_11)
} tamat_soa;
typedef struct thmat_soa_t {
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double rc00[sizeh];   // Re(comp_00)
  double rc11[sizeh];   // Re(comp_11)
} thmat_soa;


typedef struct global_vec3_soa_t {
    d_complex c0[GL_SIZEH];
    d_complex c1[GL_SIZEH];
    d_complex c2[GL_SIZEH];
} global_vec3_soa;
typedef struct global_su3_soa_t {
  global_vec3_soa r0;
  global_vec3_soa r1;
  global_vec3_soa r2;
} global_su3_soa;



#endif


