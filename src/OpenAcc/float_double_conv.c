#ifndef FLOAT_DOUBLE_CONV_C
#define FLOAT_DOUBLE_CONV_C

#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/float_double_conv.h"

////////////  VEC3_SOA    float <==> double conversions /////////////////////////////
void convert_float_to_double_vec3_soa(__restrict f_vec3_soa * f_var,
				      __restrict   vec3_soa * d_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    d_var->c0[t] = ((double)crealf(f_var->c0[t])) + ((double)cimagf(f_var->c0[t]))*I;
    d_var->c1[t] = ((double)crealf(f_var->c1[t])) + ((double)cimagf(f_var->c1[t]))*I;
    d_var->c2[t] = ((double)crealf(f_var->c2[t])) + ((double)cimagf(f_var->c2[t]))*I;
  }
}
void convert_double_to_float_vec3_soa(__restrict   vec3_soa * d_var,
				      __restrict f_vec3_soa * f_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    f_var->c0[t] = ((float)creal(d_var->c0[t])) + ((float)cimag(d_var->c0[t]))*I;
    f_var->c1[t] = ((float)creal(d_var->c1[t])) + ((float)cimag(d_var->c1[t]))*I;
    f_var->c2[t] = ((float)creal(d_var->c2[t])) + ((float)cimag(d_var->c2[t]))*I;
  }
}
////////////  VEC3    float <==> double conversions /////////////////////////////
void convert_float_to_double_vec3(__restrict f_vec3 * f_var,
				  __restrict   vec3 * d_var){
  d_var->c0 = ((double)crealf(f_var->c0)) + ((double)cimagf(f_var->c0))*I;
  d_var->c1 = ((double)crealf(f_var->c1)) + ((double)cimagf(f_var->c1))*I;
  d_var->c2 = ((double)crealf(f_var->c2)) + ((double)cimagf(f_var->c2))*I;
}
void convert_double_to_float_vec3(__restrict   vec3 * d_var,
				  __restrict f_vec3 * f_var){
    f_var->c0 = ((float)creal(d_var->c0)) + ((float)cimag(d_var->c0))*I;
    f_var->c1 = ((float)creal(d_var->c1)) + ((float)cimag(d_var->c1))*I;
    f_var->c2 = ((float)creal(d_var->c2)) + ((float)cimag(d_var->c2))*I;
}

////////////  COMPLEX    float <==> double conversions /////////////////////////////
void convert_float_to_double_complex_soa(__restrict fcomplex_soa * f_var,
					 __restrict dcomplex_soa * d_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    d_var->c[t] = ((double)crealf(f_var->c[t])) + ((double)cimagf(f_var->c[t]))*I;
  }
}  
void convert_double_to_float_complex_soa(__restrict dcomplex_soa * d_var,
					 __restrict fcomplex_soa * f_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    f_var->c[t] = ((float)creal(d_var->c[t])) + ((float)cimag(d_var->c[t]))*I;
  }
}  

////////////  REAL    float <==> double conversions /////////////////////////////
void convert_float_to_double_real_soa(__restrict  float_soa * f_var,
				      __restrict double_soa * d_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    d_var->d[t] = ((double)(f_var->d[t]));
  }
}  
void convert_double_to_float_real_soa(__restrict double_soa * d_var,
				      __restrict  float_soa * f_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    f_var->d[t] = ((float)(d_var->d[t]));
  }
}  
////////////  SU3_SOA    float <==> double conversions /////////////////////////////
void convert_float_to_double_su3_soa(__restrict f_su3_soa * f_var,
				     __restrict   su3_soa * d_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    d_var->r0.c0[t] = ((double)crealf(f_var->r0.c0[t])) + ((double)cimagf(f_var->r0.c0[t]))*I;
    d_var->r0.c1[t] = ((double)crealf(f_var->r0.c1[t])) + ((double)cimagf(f_var->r0.c1[t]))*I;
    d_var->r0.c2[t] = ((double)crealf(f_var->r0.c2[t])) + ((double)cimagf(f_var->r0.c2[t]))*I;

    d_var->r1.c0[t] = ((double)crealf(f_var->r1.c0[t])) + ((double)cimagf(f_var->r1.c0[t]))*I;
    d_var->r1.c1[t] = ((double)crealf(f_var->r1.c1[t])) + ((double)cimagf(f_var->r1.c1[t]))*I;
    d_var->r1.c2[t] = ((double)crealf(f_var->r1.c2[t])) + ((double)cimagf(f_var->r1.c2[t]))*I;

    d_var->r2.c0[t] = ((double)crealf(f_var->r2.c0[t])) + ((double)cimagf(f_var->r2.c0[t]))*I;
    d_var->r2.c1[t] = ((double)crealf(f_var->r2.c1[t])) + ((double)cimagf(f_var->r2.c1[t]))*I;
    d_var->r2.c2[t] = ((double)crealf(f_var->r2.c2[t])) + ((double)cimagf(f_var->r2.c2[t]))*I;
  }
}
void convert_double_to_float_su3_soa(__restrict   su3_soa * d_var,
				     __restrict f_su3_soa * f_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    f_var->r0.c0[t] = ((float)creal(d_var->r0.c0[t])) + ((float)cimag(d_var->r0.c0[t]))*I;
    f_var->r0.c1[t] = ((float)creal(d_var->r0.c1[t])) + ((float)cimag(d_var->r0.c1[t]))*I;
    f_var->r0.c2[t] = ((float)creal(d_var->r0.c2[t])) + ((float)cimag(d_var->r0.c2[t]))*I;

    f_var->r1.c0[t] = ((float)creal(d_var->r1.c0[t])) + ((float)cimag(d_var->r1.c0[t]))*I;
    f_var->r1.c1[t] = ((float)creal(d_var->r1.c1[t])) + ((float)cimag(d_var->r1.c1[t]))*I;
    f_var->r1.c2[t] = ((float)creal(d_var->r1.c2[t])) + ((float)cimag(d_var->r1.c2[t]))*I;

    f_var->r2.c0[t] = ((float)creal(d_var->r2.c0[t])) + ((float)cimag(d_var->r2.c0[t]))*I;
    f_var->r2.c1[t] = ((float)creal(d_var->r2.c1[t])) + ((float)cimag(d_var->r2.c1[t]))*I;
    f_var->r2.c2[t] = ((float)creal(d_var->r2.c2[t])) + ((float)cimag(d_var->r2.c2[t]))*I;
  }
}
////////////  TAMAT_SOA    float <==> double conversions /////////////////////////////
void convert_float_to_double_tamat_soa(__restrict f_tamat_soa * f_var,
				       __restrict   tamat_soa * d_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    d_var->c01[t]  = ((double)crealf(f_var->c01[t]))  + ((double)cimagf(f_var->c01[t]))*I;
    d_var->c02[t]  = ((double)crealf(f_var->c02[t]))  + ((double)cimagf(f_var->c02[t]))*I;
    d_var->c12[t]  = ((double)crealf(f_var->c12[t]))  + ((double)cimagf(f_var->c12[t]))*I;
    d_var->rc00[t] = ((double)(f_var->rc00[t]));
    d_var->rc11[t] = ((double)(f_var->rc11[t]));
  }
}
void convert_double_to_float_tamat_soa(__restrict   tamat_soa * d_var,
				       __restrict f_tamat_soa * f_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    f_var->c01[t]  = ((float)creal(d_var->c01[t])) + ((float)cimag(d_var->c01[t]))*I;
    f_var->c02[t]  = ((float)creal(d_var->c02[t])) + ((float)cimag(d_var->c02[t]))*I;
    f_var->c12[t]  = ((float)creal(d_var->c12[t])) + ((float)cimag(d_var->c12[t]))*I;
    f_var->rc00[t] = ((float)(d_var->rc00[t]));
    f_var->rc11[t] = ((float)(d_var->rc11[t]));
  }
}
////////////  THMAT_SOA    float <==> double conversions /////////////////////////////
void convert_float_to_double_thmat_soa(__restrict f_thmat_soa * f_var,
				       __restrict   thmat_soa * d_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    d_var->c01[t]  = ((double)crealf(f_var->c01[t]))  + ((double)cimagf(f_var->c01[t]))*I;
    d_var->c02[t]  = ((double)crealf(f_var->c02[t]))  + ((double)cimagf(f_var->c02[t]))*I;
    d_var->c12[t]  = ((double)crealf(f_var->c12[t]))  + ((double)cimagf(f_var->c12[t]))*I;
    d_var->rc00[t] = ((double)(f_var->rc00[t]));
    d_var->rc11[t] = ((double)(f_var->rc11[t]));
  }
}
void convert_double_to_float_thmat_soa(__restrict   thmat_soa * d_var,
				       __restrict f_thmat_soa * f_var){
  int t;
#pragma acc kernels present(f_var) present(d_var)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    f_var->c01[t]  = ((float)creal(d_var->c01[t])) + ((float)cimag(d_var->c01[t]))*I;
    f_var->c02[t]  = ((float)creal(d_var->c02[t])) + ((float)cimag(d_var->c02[t]))*I;
    f_var->c12[t]  = ((float)creal(d_var->c12[t])) + ((float)cimag(d_var->c12[t]))*I;
    f_var->rc00[t] = ((float)(d_var->rc00[t]));
    f_var->rc11[t] = ((float)(d_var->rc11[t]));
  }
}


#endif

