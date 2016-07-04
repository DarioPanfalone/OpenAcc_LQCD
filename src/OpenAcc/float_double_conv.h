#ifndef FLOAT_DOUBLE_CONV_H
#define FLOAT_DOUBLE_CONV_H

#include "./geometry.h"
#include "./struct_c_def.h"
#include "./sp_struct_c_def.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

////////////  VEC3_SOA    float <==> double conversions /////////////////////////
void convert_float_to_double_vec3_soa(__restrict const vec3_soa_f * f_var,
				      __restrict   vec3_soa * d_var);
void convert_double_to_float_vec3_soa(__restrict const   vec3_soa * d_var,
				      __restrict vec3_soa_f * f_var);
////////////  VEC3    float <==> double conversions /////////////////////////////
void convert_float_to_double_vec3(__restrict const vec3_f * f_var,
				  __restrict   vec3 * d_var);
void convert_double_to_float_vec3(__restrict  const  vec3 * d_var,
				  __restrict vec3_f * f_var);
////////////  COMPLEX    float <==> double conversions //////////////////////////
void convert_float_to_double_complex_soa(__restrict const  fcomplex_soa * f_var,
					 __restrict dcomplex_soa * d_var);
void convert_double_to_float_complex_soa(__restrict const  dcomplex_soa * d_var,
					 __restrict fcomplex_soa * f_var);
////////////  REAL    float <==> double conversions /////////////////////////////
void convert_float_to_double_real_soa(__restrict const   float_soa * f_var,
				      __restrict double_soa * d_var);
void convert_double_to_float_real_soa(__restrict const double_soa * d_var,
				      __restrict  float_soa * f_var);
////////////  SU3_SOA    float <==> double conversions //////////////////////////
void convert_float_to_double_su3_soa(__restrict const su3_soa_f * f_var,
				     __restrict   su3_soa * d_var);
void convert_double_to_float_su3_soa(__restrict const  su3_soa * d_var,
				     __restrict su3_soa_f * f_var);
////////////  TAMAT_SOA    float <==> double conversions ////////////////////////
void convert_float_to_double_tamat_soa(__restrict const tamat_soa_f * f_var,
				       __restrict   tamat_soa * d_var);
void convert_double_to_float_tamat_soa(__restrict const tamat_soa * d_var,
				       __restrict tamat_soa_f * f_var);
////////////  THMAT_SOA    float <==> double conversions ////////////////////////
void convert_float_to_double_thmat_soa(__restrict const thmat_soa_f * f_var,
				       __restrict   thmat_soa * d_var);
void convert_double_to_float_thmat_soa(__restrict const thmat_soa * d_var,
				       __restrict thmat_soa_f * f_var);


#endif

