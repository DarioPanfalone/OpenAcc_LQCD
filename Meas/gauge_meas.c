// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h
#ifndef GAUGE_MEAS_C
#define GAUGE_MEAS_C

#include "../OpenAcc/struct_c_def.c"
#include "../OpenAcc/alloc_vars.c"
#include "../OpenAcc/su3_utilities.c"



#pragma acc routine seq
static inline void comp_U_U_Udag_Udag(__restrict su3_soa * const mat1, int idx_mat1,
				      __restrict su3_soa * const mat2, int idx_mat2,
				      __restrict su3_soa * const mat3, int idx_mat3,
				      __restrict su3_soa * const mat4, int idx_mat4,
				      __restrict su3_soa * const out , int idx_out){

  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;
  d_complex C00,C01,C02,C10,C11,C12,C20,C21,C22;
  // LOAD A = MAT1 
  A00 = mat1->r0.c0[idx_mat1];
  A01 = mat1->r0.c1[idx_mat1];
  A02 = mat1->r0.c2[idx_mat1];
  A10 = mat1->r1.c0[idx_mat1];
  A11 = mat1->r1.c1[idx_mat1];
  A12 = mat1->r1.c2[idx_mat1];
  A20 = conj( ( A01 * A12 ) - ( A02 * A11) ) ;
  A21 = conj( ( A02 * A10 ) - ( A00 * A12) ) ;
  A22 = conj( ( A00 * A11 ) - ( A01 * A10) ) ;

  // LOAD B = MAT2
  B00 = mat2->r0.c0[idx_mat2];
  B01 = mat2->r0.c1[idx_mat2];
  B02 = mat2->r0.c2[idx_mat2];
  B10 = mat2->r1.c0[idx_mat2];
  B11 = mat2->r1.c1[idx_mat2];
  B12 = mat2->r1.c2[idx_mat2];
  B20 = conj( ( B01 * B12 ) - ( B02 * B11) ) ;
  B21 = conj( ( B02 * B10 ) - ( B00 * B12) ) ;
  B22 = conj( ( B00 * B11 ) - ( B01 * B10) ) ;

  // C = A * B = MAT1 * MAT2
  C00 = A00 * B00 + A01 * B10 + A02 * B20;
  C01 = A00 * B01 + A01 * B11 + A02 * B21;
  C02 = A00 * B02 + A01 * B12 + A02 * B22;
  C10 = A10 * B00 + A11 * B10 + A12 * B20;
  C11 = A10 * B01 + A11 * B11 + A12 * B21;
  C12 = A10 * B02 + A11 * B12 + A12 * B22;
  C20 = A20 * B00 + A21 * B10 + A22 * B20;
  C21 = A20 * B01 + A21 * B11 + A22 * B21;
  C22 = A20 * B02 + A21 * B12 + A22 * B22;

  // LOAD A = MAT3^DAG
  A00 = conj( mat3->r0.c0[idx_mat3] ) ;
  A10 = conj( mat3->r0.c1[idx_mat3] ) ;
  A20 = conj( mat3->r0.c2[idx_mat3] ) ;
  A01 = conj( mat3->r1.c0[idx_mat3] ) ;
  A11 = conj( mat3->r1.c1[idx_mat3] ) ;
  A21 = conj( mat3->r1.c2[idx_mat3] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // B = C * A = MAT1 * MAT2 * MAT3^DAG
  B00 = C00 * A00 + C01 * A10 + C02 * A20;
  B01 = C00 * A01 + C01 * A11 + C02 * A21;
  B02 = C00 * A02 + C01 * A12 + C02 * A22;
  B10 = C10 * A00 + C11 * A10 + C12 * A20;
  B11 = C10 * A01 + C11 * A11 + C12 * A21;
  B12 = C10 * A02 + C11 * A12 + C12 * A22;
  B20 = C20 * A00 + C21 * A10 + C22 * A20;
  B21 = C20 * A01 + C21 * A11 + C22 * A21;
  B22 = C20 * A02 + C21 * A12 + C22 * A22;

  // LOAD A = MAT4^DAG
  A00 = conj( mat4->r0.c0[idx_mat4] ) ;
  A10 = conj( mat4->r0.c1[idx_mat4] ) ;
  A20 = conj( mat4->r0.c2[idx_mat4] ) ;
  A01 = conj( mat4->r1.c0[idx_mat4] ) ;
  A11 = conj( mat4->r1.c1[idx_mat4] ) ;
  A21 = conj( mat4->r1.c2[idx_mat4] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // OUT = B * A = MAT1 * MAT2 * MAT3^DAG * MAT4^DAG
  out->r0.c0[idx_out] = B00 * A00 + B01 * A10 + B02 * A20;
  out->r0.c1[idx_out] = B00 * A01 + B01 * A11 + B02 * A21;
  out->r0.c2[idx_out] = B00 * A02 + B01 * A12 + B02 * A22;
  out->r1.c0[idx_out] = B10 * A00 + B11 * A10 + B12 * A20;
  out->r1.c1[idx_out] = B10 * A01 + B11 * A11 + B12 * A21;
  out->r1.c2[idx_out] = B10 * A02 + B11 * A12 + B12 * A22;
  out->r2.c0[idx_out] = B20 * A00 + B21 * A10 + B22 * A20;
  out->r2.c1[idx_out] = B20 * A01 + B21 * A11 + B22 * A21;
  out->r2.c2[idx_out] = B20 * A02 + B21 * A12 + B22 * A22;

}

#pragma acc routine seq
static inline void comp_and_add_U_Udag_Udag_U(__restrict su3_soa * const mat1, int idx_mat1,
					      __restrict su3_soa * const mat2, int idx_mat2,
					      __restrict su3_soa * const mat3, int idx_mat3,
					      __restrict su3_soa * const mat4, int idx_mat4,
					      __restrict su3_soa * const out , int idx_out){
  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;
  d_complex C00,C01,C02,C10,C11,C12,C20,C21,C22;
  // LOAD A = MAT1 
  A00 = mat1->r0.c0[idx_mat1];
  A01 = mat1->r0.c1[idx_mat1];
  A02 = mat1->r0.c2[idx_mat1];
  A10 = mat1->r1.c0[idx_mat1];
  A11 = mat1->r1.c1[idx_mat1];
  A12 = mat1->r1.c2[idx_mat1];
  A20 = conj( ( A01 * A12 ) - ( A02 * A11) ) ;
  A21 = conj( ( A02 * A10 ) - ( A00 * A12) ) ;
  A22 = conj( ( A00 * A11 ) - ( A01 * A10) ) ;

  // LOAD B = MAT2^DAG
  B00 = conj( mat2->r0.c0[idx_mat2] ) ;
  B10 = conj( mat2->r0.c1[idx_mat2] ) ;
  B20 = conj( mat2->r0.c2[idx_mat2] ) ;
  B01 = conj( mat2->r1.c0[idx_mat2] ) ;
  B11 = conj( mat2->r1.c1[idx_mat2] ) ;
  B21 = conj( mat2->r1.c2[idx_mat2] ) ;
  B02 = conj( ( B10 * B21 ) - ( B20 * B11) ) ;
  B12 = conj( ( B20 * B01 ) - ( B00 * B21) ) ;
  B22 = conj( ( B00 * B11 ) - ( B10 * B01) ) ;

  // C = A * B = MAT1 * MAT2^DAG
  C00 = A00 * B00 + A01 * B10 + A02 * B20;
  C01 = A00 * B01 + A01 * B11 + A02 * B21;
  C02 = A00 * B02 + A01 * B12 + A02 * B22;
  C10 = A10 * B00 + A11 * B10 + A12 * B20;
  C11 = A10 * B01 + A11 * B11 + A12 * B21;
  C12 = A10 * B02 + A11 * B12 + A12 * B22;
  C20 = A20 * B00 + A21 * B10 + A22 * B20;
  C21 = A20 * B01 + A21 * B11 + A22 * B21;
  C22 = A20 * B02 + A21 * B12 + A22 * B22;

  // LOAD A = MAT3^DAG
  A00 = conj( mat3->r0.c0[idx_mat3] ) ;
  A10 = conj( mat3->r0.c1[idx_mat3] ) ;
  A20 = conj( mat3->r0.c2[idx_mat3] ) ;
  A01 = conj( mat3->r1.c0[idx_mat3] ) ;
  A11 = conj( mat3->r1.c1[idx_mat3] ) ;
  A21 = conj( mat3->r1.c2[idx_mat3] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // B = C * A = MAT1 * MAT2^DAG * MAT3^DAG
  B00 = C00 * A00 + C01 * A10 + C02 * A20;
  B01 = C00 * A01 + C01 * A11 + C02 * A21;
  B02 = C00 * A02 + C01 * A12 + C02 * A22;
  B10 = C10 * A00 + C11 * A10 + C12 * A20;
  B11 = C10 * A01 + C11 * A11 + C12 * A21;
  B12 = C10 * A02 + C11 * A12 + C12 * A22;
  B20 = C20 * A00 + C21 * A10 + C22 * A20;
  B21 = C20 * A01 + C21 * A11 + C22 * A21;
  B22 = C20 * A02 + C21 * A12 + C22 * A22;

  // LOAD A = MAT4
  A00 = mat4->r0.c0[idx_mat4];
  A01 = mat4->r0.c1[idx_mat4];
  A02 = mat4->r0.c2[idx_mat4];
  A10 = mat4->r1.c0[idx_mat4];
  A11 = mat4->r1.c1[idx_mat4];
  A12 = mat4->r1.c2[idx_mat4];
  A20 = conj( ( A01 * A12 ) - ( A02 * A11) ) ;
  A21 = conj( ( A02 * A10 ) - ( A00 * A12) ) ;
  A22 = conj( ( A00 * A11 ) - ( A01 * A10) ) ;

  // OUT = B * A = MAT1 * MAT2^DAG * MAT3^DAG * MAT4
  out->r0.c0[idx_out] += B00 * A00 + B01 * A10 + B02 * A20;
  out->r0.c1[idx_out] += B00 * A01 + B01 * A11 + B02 * A21;
  out->r0.c2[idx_out] += B00 * A02 + B01 * A12 + B02 * A22;
  out->r1.c0[idx_out] += B10 * A00 + B11 * A10 + B12 * A20;
  out->r1.c1[idx_out] += B10 * A01 + B11 * A11 + B12 * A21;
  out->r1.c2[idx_out] += B10 * A02 + B11 * A12 + B12 * A22;
  out->r2.c0[idx_out] += B20 * A00 + B21 * A10 + B22 * A20;
  out->r2.c1[idx_out] += B20 * A01 + B21 * A11 + B22 * A21;
  out->r2.c2[idx_out] += B20 * A02 + B21 * A12 + B22 * A22;

}

#pragma acc routine seq
static inline void comp_and_add_Udag_Udag_U_U(__restrict su3_soa * const mat1, int idx_mat1,
					      __restrict su3_soa * const mat2, int idx_mat2,
					      __restrict su3_soa * const mat3, int idx_mat3,
					      __restrict su3_soa * const mat4, int idx_mat4,
					      __restrict su3_soa * const out , int idx_out){
  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;
  d_complex C00,C01,C02,C10,C11,C12,C20,C21,C22;
  // LOAD A = MAT1^DAG
  A00 = conj( mat1->r0.c0[idx_mat1] ) ;
  A10 = conj( mat1->r0.c1[idx_mat1] ) ;
  A20 = conj( mat1->r0.c2[idx_mat1] ) ;
  A01 = conj( mat1->r1.c0[idx_mat1] ) ;
  A11 = conj( mat1->r1.c1[idx_mat1] ) ;
  A21 = conj( mat1->r1.c2[idx_mat1] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // LOAD B = MAT2^DAG
  B00 = conj( mat2->r0.c0[idx_mat2] ) ;
  B10 = conj( mat2->r0.c1[idx_mat2] ) ;
  B20 = conj( mat2->r0.c2[idx_mat2] ) ;
  B01 = conj( mat2->r1.c0[idx_mat2] ) ;
  B11 = conj( mat2->r1.c1[idx_mat2] ) ;
  B21 = conj( mat2->r1.c2[idx_mat2] ) ;
  B02 = conj( ( B10 * B21 ) - ( B20 * B11) ) ;
  B12 = conj( ( B20 * B01 ) - ( B00 * B21) ) ;
  B22 = conj( ( B00 * B11 ) - ( B10 * B01) ) ;

  // C = A * B = MAT1^DAG * MAT2^DAG
  C00 = A00 * B00 + A01 * B10 + A02 * B20;
  C01 = A00 * B01 + A01 * B11 + A02 * B21;
  C02 = A00 * B02 + A01 * B12 + A02 * B22;
  C10 = A10 * B00 + A11 * B10 + A12 * B20;
  C11 = A10 * B01 + A11 * B11 + A12 * B21;
  C12 = A10 * B02 + A11 * B12 + A12 * B22;
  C20 = A20 * B00 + A21 * B10 + A22 * B20;
  C21 = A20 * B01 + A21 * B11 + A22 * B21;
  C22 = A20 * B02 + A21 * B12 + A22 * B22;

  // LOAD A = MAT3
  A00 = mat3->r0.c0[idx_mat3];
  A01 = mat3->r0.c1[idx_mat3];
  A02 = mat3->r0.c2[idx_mat3];
  A10 = mat3->r1.c0[idx_mat3];
  A11 = mat3->r1.c1[idx_mat3];
  A12 = mat3->r1.c2[idx_mat3];
  A20 = conj( ( A01 * A12 ) - ( A02 * A11) ) ;
  A21 = conj( ( A02 * A10 ) - ( A00 * A12) ) ;
  A22 = conj( ( A00 * A11 ) - ( A01 * A10) ) ;

  // B = C * A = MAT1^DAG * MAT2^DAG * MAT3
  B00 = C00 * A00 + C01 * A10 + C02 * A20;
  B01 = C00 * A01 + C01 * A11 + C02 * A21;
  B02 = C00 * A02 + C01 * A12 + C02 * A22;
  B10 = C10 * A00 + C11 * A10 + C12 * A20;
  B11 = C10 * A01 + C11 * A11 + C12 * A21;
  B12 = C10 * A02 + C11 * A12 + C12 * A22;
  B20 = C20 * A00 + C21 * A10 + C22 * A20;
  B21 = C20 * A01 + C21 * A11 + C22 * A21;
  B22 = C20 * A02 + C21 * A12 + C22 * A22;

  // LOAD A = MAT4
  A00 = mat4->r0.c0[idx_mat4];
  A01 = mat4->r0.c1[idx_mat4];
  A02 = mat4->r0.c2[idx_mat4];
  A10 = mat4->r1.c0[idx_mat4];
  A11 = mat4->r1.c1[idx_mat4];
  A12 = mat4->r1.c2[idx_mat4];
  A20 = conj( ( A01 * A12 ) - ( A02 * A11) ) ;
  A21 = conj( ( A02 * A10 ) - ( A00 * A12) ) ;
  A22 = conj( ( A00 * A11 ) - ( A01 * A10) ) ;

  // OUT = B * A = MAT1^DAG * MAT2^DAG * MAT3 * MAT4
  out->r0.c0[idx_out] += B00 * A00 + B01 * A10 + B02 * A20;
  out->r0.c1[idx_out] += B00 * A01 + B01 * A11 + B02 * A21;
  out->r0.c2[idx_out] += B00 * A02 + B01 * A12 + B02 * A22;
  out->r1.c0[idx_out] += B10 * A00 + B11 * A10 + B12 * A20;
  out->r1.c1[idx_out] += B10 * A01 + B11 * A11 + B12 * A21;
  out->r1.c2[idx_out] += B10 * A02 + B11 * A12 + B12 * A22;
  out->r2.c0[idx_out] += B20 * A00 + B21 * A10 + B22 * A20;
  out->r2.c1[idx_out] += B20 * A01 + B21 * A11 + B22 * A21;
  out->r2.c2[idx_out] += B20 * A02 + B21 * A12 + B22 * A22;

}

#pragma acc routine seq
static inline void comp_and_add_Udag_U_U_Udag(__restrict su3_soa * const mat1, int idx_mat1,
					      __restrict su3_soa * const mat2, int idx_mat2,
					      __restrict su3_soa * const mat3, int idx_mat3,
					      __restrict su3_soa * const mat4, int idx_mat4,
					      __restrict su3_soa * const out , int idx_out){
  
  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;
  d_complex C00,C01,C02,C10,C11,C12,C20,C21,C22;
  // LOAD A = MAT1^DAG
  A00 = conj( mat1->r0.c0[idx_mat1] ) ;
  A10 = conj( mat1->r0.c1[idx_mat1] ) ;
  A20 = conj( mat1->r0.c2[idx_mat1] ) ;
  A01 = conj( mat1->r1.c0[idx_mat1] ) ;
  A11 = conj( mat1->r1.c1[idx_mat1] ) ;
  A21 = conj( mat1->r1.c2[idx_mat1] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // LOAD B = MAT2
  B00 = mat2->r0.c0[idx_mat2];
  B01 = mat2->r0.c1[idx_mat2];
  B02 = mat2->r0.c2[idx_mat2];
  B10 = mat2->r1.c0[idx_mat2];
  B11 = mat2->r1.c1[idx_mat2];
  B12 = mat2->r1.c2[idx_mat2];
  B20 = conj( ( B01 * B12 ) - ( B02 * B11) ) ;
  B21 = conj( ( B02 * B10 ) - ( B00 * B12) ) ;
  B22 = conj( ( B00 * B11 ) - ( B01 * B10) ) ;

  // C = A * B = MAT1^DAG * MAT2
  C00 = A00 * B00 + A01 * B10 + A02 * B20;
  C01 = A00 * B01 + A01 * B11 + A02 * B21;
  C02 = A00 * B02 + A01 * B12 + A02 * B22;
  C10 = A10 * B00 + A11 * B10 + A12 * B20;
  C11 = A10 * B01 + A11 * B11 + A12 * B21;
  C12 = A10 * B02 + A11 * B12 + A12 * B22;
  C20 = A20 * B00 + A21 * B10 + A22 * B20;
  C21 = A20 * B01 + A21 * B11 + A22 * B21;
  C22 = A20 * B02 + A21 * B12 + A22 * B22;

  // LOAD A = MAT3
  A00 = mat3->r0.c0[idx_mat3];
  A01 = mat3->r0.c1[idx_mat3];
  A02 = mat3->r0.c2[idx_mat3];
  A10 = mat3->r1.c0[idx_mat3];
  A11 = mat3->r1.c1[idx_mat3];
  A12 = mat3->r1.c2[idx_mat3];
  A20 = conj( ( A01 * A12 ) - ( A02 * A11) ) ;
  A21 = conj( ( A02 * A10 ) - ( A00 * A12) ) ;
  A22 = conj( ( A00 * A11 ) - ( A01 * A10) ) ;

  // B = C * A = MAT1^DAG * MAT2 * MAT3
  B00 = C00 * A00 + C01 * A10 + C02 * A20;
  B01 = C00 * A01 + C01 * A11 + C02 * A21;
  B02 = C00 * A02 + C01 * A12 + C02 * A22;
  B10 = C10 * A00 + C11 * A10 + C12 * A20;
  B11 = C10 * A01 + C11 * A11 + C12 * A21;
  B12 = C10 * A02 + C11 * A12 + C12 * A22;
  B20 = C20 * A00 + C21 * A10 + C22 * A20;
  B21 = C20 * A01 + C21 * A11 + C22 * A21;
  B22 = C20 * A02 + C21 * A12 + C22 * A22;

  // LOAD A = MAT4^DAG
  A00 = conj( mat4->r0.c0[idx_mat4] ) ;
  A10 = conj( mat4->r0.c1[idx_mat4] ) ;
  A20 = conj( mat4->r0.c2[idx_mat4] ) ;
  A01 = conj( mat4->r1.c0[idx_mat4] ) ;
  A11 = conj( mat4->r1.c1[idx_mat4] ) ;
  A21 = conj( mat4->r1.c2[idx_mat4] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // OUT = B * A = MAT1^DAG * MAT2 * MAT3 * MAT4^DAG
  out->r0.c0[idx_out] += B00 * A00 + B01 * A10 + B02 * A20;
  out->r0.c1[idx_out] += B00 * A01 + B01 * A11 + B02 * A21;
  out->r0.c2[idx_out] += B00 * A02 + B01 * A12 + B02 * A22;
  out->r1.c0[idx_out] += B10 * A00 + B11 * A10 + B12 * A20;
  out->r1.c1[idx_out] += B10 * A01 + B11 * A11 + B12 * A21;
  out->r1.c2[idx_out] += B10 * A02 + B11 * A12 + B12 * A22;
  out->r2.c0[idx_out] += B20 * A00 + B21 * A10 + B22 * A20;
  out->r2.c1[idx_out] += B20 * A01 + B21 * A11 + B22 * A21;
  out->r2.c2[idx_out] += B20 * A02 + B21 * A12 + B22 * A22;


}

#pragma acc routine seq
static inline void combine_fourleaves_to_get_loc_q(__restrict su3_soa * const Qmn, int idx_Qmn,
						   __restrict su3_soa * const Qrs, int idx_Qrs,
						   int Eps,
						   double_soa * const locq      , int idx_locq){
  // Calcola  loc_q = Epsilon * Re[ Tr[ Qrs * (Qmn - Qmn^dag ) ] ] / (256 * pi * pi)
  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;
  d_complex C00,C01,C02,C10,C11,C12,C20,C21,C22;
  // LOAD A = Q_{MU,NU}
  A00 = Qmn->r0.c0[idx_Qmn];
  A01 = Qmn->r0.c1[idx_Qmn];
  A02 = Qmn->r0.c2[idx_Qmn];
  A10 = Qmn->r1.c0[idx_Qmn];
  A11 = Qmn->r1.c1[idx_Qmn];
  A12 = Qmn->r1.c2[idx_Qmn];
  A20 = Qmn->r2.c0[idx_Qmn];
  A21 = Qmn->r2.c1[idx_Qmn];
  A22 = Qmn->r2.c2[idx_Qmn];

  // LOAD B = Q_{RHO,SIGMA}
  B00 = Qrs->r0.c0[idx_Qrs];
  B01 = Qrs->r0.c1[idx_Qrs];
  B02 = Qrs->r0.c2[idx_Qrs];
  B10 = Qrs->r1.c0[idx_Qrs];
  B11 = Qrs->r1.c1[idx_Qrs];
  B12 = Qrs->r1.c2[idx_Qrs];
  B20 = Qrs->r2.c0[idx_Qrs];
  B21 = Qrs->r2.c1[idx_Qrs];
  B22 = Qrs->r2.c2[idx_Qrs];

  // COMPUTE C = A - A^DAG = Qmn - Qmn^DAG
  C00 = A00 - conj(A00);
  C01 = A01 - conj(A10);
  C02 = A02 - conj(A20);
  C10 = A10 - conj(A01);
  C11 = A11 - conj(A11);
  C12 = A12 - conj(A21);
  C20 = A20 - conj(A02);
  C21 = A21 - conj(A12);
  C22 = A22 - conj(A22);

  // COMPUTE THE DIAGONAL ELEMENTS OF THE PRODUCT A = B * C = Tr[ Qrs * (Qmn - Qmn^dag ) ]
  A00 = B00 * C00 + B01 * C10 + B02 * C20 +
        B10 * C01 + B11 * C11 + B12 * C21 +
        B20 * C02 + B21 * C12 + B22 * C22; 

  locq->d[idx_locq] += creal(A00) * Eps / (128.0 * M_PI * M_PI);

}





void compute_local_topological_charge(  __restrict su3_soa * const u,
					__restrict su3_soa * const quadri,
					double_soa * const loc_q,
					int mu, int nu){
  
  int x, y, z, t;
#pragma acc kernels present(u) present(quadri) present(loc_q) present(nnp_openacc) present(nnm_openacc)
#pragma acc loop independent gang 
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector 
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector 
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector 
        for(x=0; x < nx; x++) {
	  const int idxh = snum_acc(x,y,z,t);  // r
	  const int parity = (x+y+z+t) % 2;
	  loc_q[parity].d[idxh] = 0.0;
	  //	  for(nu=1; nu<4; nu++){
	    int b;
	    int a;
	    int rho;
	    int sigma;
	    int Epsilon;
	    
	    //ricorda che nu > mu sempre
	    //scelgo b (tra 0 e 3) che sia diverso da nu e da mu :
	    if(nu!=3)b=nu+1;
	    if(nu==3 && mu!=2)b=mu+1;
	    if(nu==3 && mu==2)b=1;
	    //scelgo a(tra 0 e 3) che sia diverso da mu, nu e b:
	    a=6-mu-nu-b; 
	    //ordino la direzione maggiore e quella minore per la plaquette:
	    if(a<b){rho=a;}else{rho=b;}
	    sigma=a+b-rho;
	    if((mu==0 && nu==2)||(mu==1 && nu==3)){
	      Epsilon=-1;
	    }else{
	      Epsilon=1;
	    }
	    
	    int idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
	    int idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
	    int idxmmu = nnm_openacc[idxh][mu][parity];// r-mu
	    int idxmnu = nnm_openacc[idxh][nu][parity];// r-nu
	    int idxmmupnu = nnp_openacc[idxmmu][nu][!parity]; // r-mu+nu
	    int idxmmumnu = nnm_openacc[idxmmu][nu][!parity]; // r-mu-nu
	    int idxpmumnu = nnm_openacc[idxpmu][nu][!parity]; // r+mu-nu
	    //piano mu-nu
	    comp_U_U_Udag_Udag(&u[2*mu+parity],   idxh,
			       &u[2*nu+!parity],  idxpmu,
			       &u[2*mu+!parity],  idxpnu,
			       &u[2*nu+parity],   idxh,
			       &quadri[parity],   idxh);
	    
	    comp_and_add_U_Udag_Udag_U(&u[2*nu+parity],   idxh,
				       &u[2*mu+parity],   idxmmupnu,
				       &u[2*nu+!parity],  idxmmu,
				       &u[2*mu+!parity],  idxmmu,
				       &quadri[parity],   idxh);
	    
	    comp_and_add_Udag_Udag_U_U(&u[2*mu+!parity],  idxmmu,
				       &u[2*nu+parity],   idxmmumnu,
				       &u[2*mu+parity],   idxmmumnu,
				       &u[2*nu+!parity],  idxmnu,
				       &quadri[parity],   idxh);
	    
	    comp_and_add_Udag_U_U_Udag(&u[2*nu+!parity],  idxmnu,
				       &u[2*mu+!parity],  idxmnu,
				       &u[2*nu+parity],   idxpmumnu,
				       &u[2*mu+parity],   idxh,
				       &quadri[parity],   idxh);
	    
	    int idxprho   = nnp_openacc[idxh][rho][parity];    // r+rho
	    int idxpsigma = nnp_openacc[idxh][sigma][parity];  // r+sigma
	    int idxmrho   = nnm_openacc[idxh][rho][parity];    // r-rho
	    int idxmsigma = nnm_openacc[idxh][sigma][parity];  // r-sigma
	    int idxmrhopsigma = nnp_openacc[idxmrho][sigma][!parity]; // r-rho+sigma
	    int idxmrhomsigma = nnm_openacc[idxmrho][sigma][!parity]; // r-rho-sigma
	    int idxprhomsigma = nnm_openacc[idxprho][sigma][!parity]; // r+rho-sigma
	    //piano rho-sigma
	    comp_U_U_Udag_Udag(&u[2*rho+parity],     idxh,
			       &u[2*sigma+!parity],  idxprho,
			       &u[2*rho+!parity],    idxpsigma,
			       &u[2*sigma+parity],   idxh,
			       &quadri[2+parity],    idxh);
	    
	    comp_and_add_U_Udag_Udag_U(&u[2*sigma+parity],   idxh,
				       &u[2*rho+parity],     idxmrhopsigma,
				       &u[2*sigma+!parity],  idxmrho,
				       &u[2*rho+!parity],    idxmrho,
				       &quadri[2+parity],    idxh);
	    
	    comp_and_add_Udag_Udag_U_U(&u[2*rho+!parity],    idxmrho,
				       &u[2*sigma+parity],   idxmrhomsigma,
				       &u[2*rho+parity],     idxmrhomsigma,
				       &u[2*sigma+!parity],  idxmsigma,
				       &quadri[2+parity],    idxh);
	    
	    comp_and_add_Udag_U_U_Udag(&u[2*sigma+!parity],  idxmsigma,
				       &u[2*rho+!parity],    idxmsigma,
				       &u[2*sigma+parity],   idxprhomsigma,
				       &u[2*rho+parity],     idxh,
				       &quadri[2+parity],    idxh);
	    // Calcola  loc_q = Epsilon * Re[ Tr[ Prs * (Pmn - Pmn^dag ) ] ] / (256 * pi * pi)
	    combine_fourleaves_to_get_loc_q(&quadri[  parity],   idxh,
					    &quadri[2+parity],   idxh,
					    Epsilon,
					    &loc_q[parity],      idxh);
	    //	      if((idxh==0)&&(parity==0)) printf("loc_q(mu=%d,nu=%d)= %.18lf\n",mu,nu,loc_q[0].d[0]);
	    
	    //	  } // closes nu
	} // closes x
      } // closes y
    } // closes z
  } // closes t
}


double reduce_loc_top_charge(double_soa * const loc_q){
  
  double result=0.0;
  int t;
#pragma acc kernels present(loc_q)
#pragma acc loop reduction(+:result)
  for(t=0; t<sizeh; t++) {
    result += loc_q[0].d[t];
    result += loc_q[1].d[t];
  }
  return result;
}


double compute_topological_charge(__restrict su3_soa * const u,
				  __restrict su3_soa * const quadri,
				  double_soa * const loc_q){  

  set_su3_soa_to_zero(quadri); // forse non serve a una mazza
  mult_conf_times_stag_phases(u);

  double temp_ch =0.0;
  compute_local_topological_charge(u,quadri,loc_q,0,1);// (x,y) - (z,t)
  temp_ch = reduce_loc_top_charge(loc_q);
  compute_local_topological_charge(u,quadri,loc_q,0,2);// (x,z) - (y,t)
  temp_ch += reduce_loc_top_charge(loc_q);
  compute_local_topological_charge(u,quadri,loc_q,0,3);// (x,t) - (y,z)
  temp_ch += reduce_loc_top_charge(loc_q);

  mult_conf_times_stag_phases(u);

  return  temp_ch;
}


#endif
