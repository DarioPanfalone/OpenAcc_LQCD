#ifndef FLOAT_DOUBLE_CONV_H
#define FLOAT_DOUBLE_CONV_H

#include "../OpenAcc/struct_c_def.h"

typedef float  complex  f_complex;

typedef struct f_vec3_soa_t {
    int status;
    f_complex c0[sizeh];
    f_complex c1[sizeh];
    f_complex c2[sizeh];
} f_vec3_soa;
typedef struct f_vec3_t {
  f_complex c0;
  f_complex c1;
  f_complex c2;
} f_vec3;
typedef struct fcomplex_soa_t {
    int status;
    f_complex c[sizeh];
} fcomplex_soa;
typedef struct float_soa_t {
    int status;
    float d[sizeh];
} float_soa;
typedef struct f_su3_soa_t {
  int status;  
  f_vec3_soa r0;
  f_vec3_soa r1;
  f_vec3_soa r2;
} f_su3_soa;
typedef struct f_tamat_soa_t {
  int status;
  f_complex c01[sizeh]; // comp_01
  f_complex c02[sizeh]; // comp_02
  f_complex c12[sizeh]; // comp_12
  float rc00[sizeh];   // Im(comp_00)
  float rc11[sizeh];   // Im(comp_11)
} f_tamat_soa;
typedef struct f_thmat_soa_t {
  int status;
  f_complex c01[sizeh]; // comp_01
  f_complex c02[sizeh]; // comp_02
  f_complex c12[sizeh]; // comp_12
  float  rc00[sizeh];   // Re(comp_00)
  float  rc11[sizeh];   // Re(comp_11)
} thmat_soa;


////////////  VEC3_SOA    float <==> double conversions /////////////////////////
void convert_float_to_double_vec3_soa(__restrict f_vec3_soa * f_var,
				      __restrict   vec3_soa * d_var);
void convert_double_to_float_vec3_soa(__restrict   vec3_soa * d_var,
				      __restrict f_vec3_soa * f_var);
////////////  VEC3    float <==> double conversions /////////////////////////////
void convert_float_to_double_vec3(__restrict f_vec3 * f_var,
				  __restrict   vec3 * d_var);
void convert_double_to_float_vec3(__restrict   vec3 * d_var,
				  __restrict f_vec3 * f_var);
////////////  COMPLEX    float <==> double conversions //////////////////////////
void convert_float_to_double_complex_soa(__restrict fcomplex_soa * f_var,
					 __restrict dcomplex_soa * d_var);
void convert_double_to_float_complex_soa(__restrict dcomplex_soa * d_var,
					 __restrict fcomplex_soa * f_var);
////////////  REAL    float <==> double conversions /////////////////////////////
void convert_float_to_double_real_soa(__restrict  float_soa * f_var,
				      __restrict double_soa * d_var);
void convert_double_to_float_real_soa(__restrict double_soa * d_var,
				      __restrict  float_soa * f_var);
////////////  SU3_SOA    float <==> double conversions //////////////////////////
void convert_float_to_double_su3_soa(__restrict f_su3_soa * f_var,
				     __restrict   su3_soa * d_var);
void convert_double_to_float_su3_soa(__restrict   su3_soa * d_var,
				     __restrict f_su3_soa * f_var);
////////////  TAMAT_SOA    float <==> double conversions ////////////////////////
void convert_float_to_double_tamat_soa(__restrict f_tamat_soa * f_var,
				       __restrict   tamat_soa * d_var);
void convert_double_to_float_tamat_soa(__restrict   tamat_soa * d_var,
				       __restrict f_tamat_soa * f_var);
////////////  THMAT_SOA    float <==> double conversions ////////////////////////
void convert_float_to_double_thmat_soa(__restrict f_thmat_soa * f_var,
				       __restrict   thmat_soa * d_var);
void convert_double_to_float_thmat_soa(__restrict   thmat_soa * d_var,
				       __restrict f_thmat_soa * f_var);


#endif

