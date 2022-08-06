#ifndef SU3_UTILITIES_H_
#define SU3_UTILITIES_H_

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "./single_types.h"
#include "./cayley_hamilton.h"
#include "./action.h"
// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#include <math.h>
#else // assuming compilation with PGI for accelerators
#include <accelmath.h>
#include "../OpenAcc/deviceinit.h"
#endif

// multiply the whole configuration for the staggered phases field
// (only the first two lines)
// void mult_conf_times_stag_phases( __restrict su3_soa * const u);

// multiply the whole configuration for the staggered phases field
// (all three lines)
// void mult_gl3_soa_times_stag_phases( __restrict su3_soa * const u);

// multiply the whole configuration for the staggered phases field
// void mult_conf_times_stag_phases_nodev( __restrict su3_soa * const u);

void set_su3_soa_to_zero(__restrict su3_soa * const matrix);

//copy matrix in into matrix out, this has to happen on the host
void set_su3_soa_to_su3_soa(__restrict const su3_soa * const matrix_in,
														__restrict su3_soa * const matrix_out);

void set_su3_soa_to_su3_soa_trasl(__restrict const su3_soa * const matrix_in,
																	__restrict su3_soa * const matrix_out,int dir);

void set_su3_soa_to_su3_soa_device(__restrict const su3_soa * const matrix_in,
																	 __restrict su3_soa * const matrix_out);

void conf_times_staples_ta_part(__restrict const su3_soa * const u,
																__restrict const su3_soa * const loc_stap,
																__restrict tamat_soa * const tipdot);

void conf_times_staples_ta_part_addto_tamat(__restrict const su3_soa * const u,
																						__restrict const su3_soa * const loc_stap,
																						__restrict tamat_soa * const tipdot);

void RHO_times_conf_times_staples_ta_part(__restrict const su3_soa * const u,       
																					__restrict const su3_soa * const loc_stap,
																					__restrict tamat_soa * const tipdot, int istopo); // istopo = {0,1} -> rho = {fermrho,toporho}

void mom_sum_mult(__restrict thmat_soa * const mom,
									__restrict const tamat_soa * const ipdot, const double * factor,
									int id_factor);


void mom_exp_times_conf_soloopenacc(__restrict  su3_soa * const tconf_acc,
																		__restrict const thmat_soa * const tmomenta, 
																		const double * tdelta,  int id_delta);

// reunitarize the conf by brute force
void unitarize_conf(__restrict su3_soa * const u);



#pragma acc routine seq
static inline void loc_unitarize_conf(__restrict su3_soa * const cnf,
																			const int idx_cnf)
{
  d_complex A00,A01,A02,A10,A11,A12;
  //d_complex A20,A21,A22;
  A00 = cnf->r0.c0[idx_cnf];
  A01 = cnf->r0.c1[idx_cnf];
  A02 = cnf->r0.c2[idx_cnf];
  A10 = cnf->r1.c0[idx_cnf];
  A11 = cnf->r1.c1[idx_cnf];
  A12 = cnf->r1.c2[idx_cnf];

  // first row normalization
  double NORM = creal(A00)*creal(A00)+cimag(A00)*cimag(A00) 
		+ creal(A01)*creal(A01)+cimag(A01)*cimag(A01) 
		+ creal(A02)*creal(A02)+cimag(A02)*cimag(A02);
  NORM = 1.0/sqrt(NORM);
  A00 *= NORM;
  A01 *= NORM;
  A02 *= NORM;

	// scalar product with the second and subtraction (orthogonalize)
  d_complex SCAL_PROD = conj(A00) * A10 + conj(A01) * A11 + 
		conj(A02) * A12;

  A10 -= SCAL_PROD * A00;
  A11 -= SCAL_PROD * A01;
  A12 -= SCAL_PROD * A02;

  // second row normalization
  NORM = creal(A10)*creal(A10)+cimag(A10)*cimag(A10) + 
		creal(A11)*creal(A11)+cimag(A11)*cimag(A11) + 
		creal(A12)*creal(A12)+cimag(A12)*cimag(A12);
  NORM = 1.0/sqrt(NORM);
  A10 *= NORM;
  A11 *= NORM;
  A12 *= NORM;

  cnf->r0.c0[idx_cnf] = A00;
  cnf->r0.c1[idx_cnf] = A01;
  cnf->r0.c2[idx_cnf] = A02;
  cnf->r1.c0[idx_cnf] = A10;
  cnf->r1.c1[idx_cnf] = A11;
  cnf->r1.c2[idx_cnf] = A12;
  
}

// mat3 = mat1 * mat2
#pragma acc routine seq
static inline void mat1_times_mat2_into_mat3_absent_stag_phases(__restrict const su3_soa * const mat1,  const int idx_mat1,
																																__restrict const su3_soa * const mat2,  const int idx_mat2,
																																__restrict su3_soa * const mat3,   const int idx_mat3)
{
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];

  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];

  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_02 = mat2->r0.c2[idx_mat2];

  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];

  // compute 3rd mat2 row from the first two
  d_complex mat2_20 = conj(( mat2_01 * mat2_12 )-( mat2_02 * mat2_11));
  d_complex mat2_21 = conj(( mat2_02 * mat2_10 )-( mat2_00 * mat2_12));
  d_complex mat2_22 = conj(( mat2_00 * mat2_11 )-( mat2_01 * mat2_10));

  // compute the first two rows of the solution
  mat3->r0.c0[idx_mat3] = mat1_00 * mat2_00 + mat1_01 * mat2_10 
		+ mat1_02 * mat2_20 ;
  mat3->r0.c1[idx_mat3] = mat1_00 * mat2_01 + mat1_01 * mat2_11 
		+ mat1_02 * mat2_21 ;
  mat3->r0.c2[idx_mat3] = mat1_00 * mat2_02 + mat1_01 * mat2_12 
		+ mat1_02 * mat2_22 ;

  mat3->r1.c0[idx_mat3] = mat1_10 * mat2_00 + mat1_11 * mat2_10 
		+ mat1_12 * mat2_20 ;
  mat3->r1.c1[idx_mat3] = mat1_10 * mat2_01 + mat1_11 * mat2_11 
		+ mat1_12 * mat2_21 ;
  mat3->r1.c2[idx_mat3] = mat1_10 * mat2_02 + mat1_11 * mat2_12 
		+ mat1_12 * mat2_22 ;
}

// mat1 = mat1 * mat2 
#pragma acc routine seq
static inline void mat1_times_mat2_into_mat1_absent_stag_phases(__restrict su3_soa * const mat1,
																																const int idx_mat1,
																																__restrict  const su3_soa * const mat2,
																																const int idx_mat2)
{
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];

  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];

  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_02 = mat2->r0.c2[idx_mat2];

  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];

  // compute 3rd mat2 row from the first two                    
  d_complex mat2_20 = conj((mat2_01 * mat2_12 )-( mat2_02 * mat2_11));
  d_complex mat2_21 = conj((mat2_02 * mat2_10 )-( mat2_00 * mat2_12));
  d_complex mat2_22 = conj((mat2_00 * mat2_11 )-( mat2_01 * mat2_10));

  // compute the first two rows of the solution
  mat1->r0.c0[idx_mat1] = mat1_00 * mat2_00 + mat1_01 * mat2_10 
		+ mat1_02 * mat2_20 ;
  mat1->r0.c1[idx_mat1] = mat1_00 * mat2_01 + mat1_01 * mat2_11 
		+ mat1_02 * mat2_21 ;
  mat1->r0.c2[idx_mat1] = mat1_00 * mat2_02 + mat1_01 * mat2_12 
		+ mat1_02 * mat2_22 ;

  mat1->r1.c0[idx_mat1] = mat1_10 * mat2_00 + mat1_11 * mat2_10 
		+ mat1_12 * mat2_20 ;
  mat1->r1.c1[idx_mat1] = mat1_10 * mat2_01 + mat1_11 * mat2_11 
		+ mat1_12 * mat2_21 ;
  mat1->r1.c2[idx_mat1] = mat1_10 * mat2_02 + mat1_11 * mat2_12 
		+ mat1_12 * mat2_22 ;
}


#pragma acc routine seq
static inline void mat1_times_conj_mat2_into_mat3_absent_stag_phases(__restrict const su3_soa * const mat1, int idx_mat1,
																																		 __restrict const su3_soa * const mat2, int idx_mat2,
																																		 __restrict su3_soa * const mat3, int idx_mat3)
{

  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;

	// load A = mat1 
  A00 = mat1->r0.c0[idx_mat1];
  A01 = mat1->r0.c1[idx_mat1];
  A02 = mat1->r0.c2[idx_mat1];
  A10 = mat1->r1.c0[idx_mat1];
  A11 = mat1->r1.c1[idx_mat1];
  A12 = mat1->r1.c2[idx_mat1];
  A20 = conj( ( A01 * A12 ) - ( A02 * A11) ) ;
  A21 = conj( ( A02 * A10 ) - ( A00 * A12) ) ;
  A22 = conj( ( A00 * A11 ) - ( A01 * A10) ) ;

  // load B = mat2^DAG
  B00 = conj( mat2->r0.c0[idx_mat2] ) ;
  B10 = conj( mat2->r0.c1[idx_mat2] ) ;
  B20 = conj( mat2->r0.c2[idx_mat2] ) ;
  B01 = conj( mat2->r1.c0[idx_mat2] ) ;
  B11 = conj( mat2->r1.c1[idx_mat2] ) ;
  B21 = conj( mat2->r1.c2[idx_mat2] ) ;
  B02 = conj( ( B10 * B21 ) - ( B20 * B11) ) ;
  B12 = conj( ( B20 * B01 ) - ( B00 * B21) ) ;
  B22 = conj( ( B00 * B11 ) - ( B10 * B01) ) ;

  // mat3 = A * B = mat1 * mat3^DAG
  mat3->r0.c0[idx_mat3] = A00 * B00 + A01 * B10 + A02 * B20;
  mat3->r0.c1[idx_mat3] = A00 * B01 + A01 * B11 + A02 * B21;
  mat3->r0.c2[idx_mat3] = A00 * B02 + A01 * B12 + A02 * B22;
  mat3->r1.c0[idx_mat3] = A10 * B00 + A11 * B10 + A12 * B20;
  mat3->r1.c1[idx_mat3] = A10 * B01 + A11 * B11 + A12 * B21;
  mat3->r1.c2[idx_mat3] = A10 * B02 + A11 * B12 + A12 * B22;
  mat3->r2.c0[idx_mat3] = A20 * B00 + A21 * B10 + A22 * B20;
  mat3->r2.c1[idx_mat3] = A20 * B01 + A21 * B11 + A22 * B21;
  mat3->r2.c2[idx_mat3] = A20 * B02 + A21 * B12 + A22 * B22;

}


// mat1 = mat1 * hermitian_conjucate(mat2)
#pragma acc routine seq
static inline void mat1_times_conj_mat2_into_mat1_absent_stag_phases(__restrict su3_soa * const mat1,
																																		 const int idx_mat1,
																																		 __restrict const su3_soa * const mat2,
																																		 const int idx_mat2)
{
  
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];

  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];

  // construct (into the variables mat2_ij) the hermitian conjugate 
  // of the mat2 matrix
  d_complex mat2_00 = conj( mat2->r0.c0[idx_mat2] ) ;
  d_complex mat2_10 = conj( mat2->r0.c1[idx_mat2] ) ;
  d_complex mat2_20 = conj( mat2->r0.c2[idx_mat2] ) ;

  d_complex mat2_01 = conj( mat2->r1.c0[idx_mat2] ) ;
  d_complex mat2_11 = conj( mat2->r1.c1[idx_mat2] ) ;
  d_complex mat2_21 = conj( mat2->r1.c2[idx_mat2] ) ;

  // compute 3rd mat2 column from the first two
  d_complex mat2_02 = conj(( mat2_10 * mat2_21 )-( mat2_20 * mat2_11));
  d_complex mat2_12 = conj(( mat2_20 * mat2_01 )-( mat2_00 * mat2_21));
  d_complex mat2_22 = conj(( mat2_00 * mat2_11 )-( mat2_10 * mat2_01));

  // compute the first two rows of the solution
  mat1->r0.c0[idx_mat1] = mat1_00 * mat2_00 + mat1_01 * mat2_10 
		+ mat1_02 * mat2_20 ;
  mat1->r0.c1[idx_mat1] = mat1_00 * mat2_01 + mat1_01 * mat2_11 
		+ mat1_02 * mat2_21 ;
  mat1->r0.c2[idx_mat1] = mat1_00 * mat2_02 + mat1_01 * mat2_12 
		+ mat1_02 * mat2_22 ;

  mat1->r1.c0[idx_mat1] = mat1_10 * mat2_00 + mat1_11 * mat2_10 
		+ mat1_12 * mat2_20 ;
  mat1->r1.c1[idx_mat1] = mat1_10 * mat2_01 + mat1_11 * mat2_11 
		+ mat1_12 * mat2_21 ;
  mat1->r1.c2[idx_mat1] = mat1_10 * mat2_02 + mat1_11 * mat2_12 
		+ mat1_12 * mat2_22 ;
}


#pragma acc routine seq
static inline void conj_mat1_times_mat2_into_mat3_absent_stag_phases(
																																		 __restrict const su3_soa * const mat1, int idx_mat1,
																																		 __restrict const su3_soa * const mat2, int idx_mat2,
																																		 __restrict su3_soa * const mat3, int idx_mat3)
{

  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;

  // load A = mat1^DAG
  A00 = conj( mat1->r0.c0[idx_mat1] ) ;
  A10 = conj( mat1->r0.c1[idx_mat1] ) ;
  A20 = conj( mat1->r0.c2[idx_mat1] ) ;
  A01 = conj( mat1->r1.c0[idx_mat1] ) ;
  A11 = conj( mat1->r1.c1[idx_mat1] ) ;
  A21 = conj( mat1->r1.c2[idx_mat1] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // load B = mat2 
  B00 = mat2->r0.c0[idx_mat2];
  B01 = mat2->r0.c1[idx_mat2];
  B02 = mat2->r0.c2[idx_mat2];
  B10 = mat2->r1.c0[idx_mat2];
  B11 = mat2->r1.c1[idx_mat2];
  B12 = mat2->r1.c2[idx_mat2];
  B20 = conj( ( B01 * B12 ) - ( B02 * B11) ) ;
  B21 = conj( ( B02 * B10 ) - ( B00 * B12) ) ;
  B22 = conj( ( B00 * B11 ) - ( B01 * B10) ) ;

  // mat3 = A * B = mat1 * mat3^DAG
  mat3->r0.c0[idx_mat3] = A00 * B00 + A01 * B10 + A02 * B20;
  mat3->r0.c1[idx_mat3] = A00 * B01 + A01 * B11 + A02 * B21;
  mat3->r0.c2[idx_mat3] = A00 * B02 + A01 * B12 + A02 * B22;
  mat3->r1.c0[idx_mat3] = A10 * B00 + A11 * B10 + A12 * B20;
  mat3->r1.c1[idx_mat3] = A10 * B01 + A11 * B11 + A12 * B21;
  mat3->r1.c2[idx_mat3] = A10 * B02 + A11 * B12 + A12 * B22;
  mat3->r2.c0[idx_mat3] = A20 * B00 + A21 * B10 + A22 * B20;
  mat3->r2.c1[idx_mat3] = A20 * B01 + A21 * B11 + A22 * B21;
  mat3->r2.c2[idx_mat3] = A20 * B02 + A21 * B12 + A22 * B22;

}

#pragma acc routine seq
static inline void conj_mat1_times_conj_mat2_into_mat3_absent_stag_phases(
																																					__restrict su3_soa * const mat1, int idx_mat1,
																																					__restrict su3_soa * const mat2, int idx_mat2,
																																					__restrict su3_soa * const mat3, int idx_mat3)
{

  d_complex A00,A01,A02,A10,A11,A12,A20,A21,A22;
  d_complex B00,B01,B02,B10,B11,B12,B20,B21,B22;

  // load A = mat1^DAG
  A00 = conj( mat1->r0.c0[idx_mat1] ) ;
  A10 = conj( mat1->r0.c1[idx_mat1] ) ;
  A20 = conj( mat1->r0.c2[idx_mat1] ) ;
  A01 = conj( mat1->r1.c0[idx_mat1] ) ;
  A11 = conj( mat1->r1.c1[idx_mat1] ) ;
  A21 = conj( mat1->r1.c2[idx_mat1] ) ;
  A02 = conj( ( A10 * A21 ) - ( A20 * A11) ) ;
  A12 = conj( ( A20 * A01 ) - ( A00 * A21) ) ;
  A22 = conj( ( A00 * A11 ) - ( A10 * A01) ) ;

  // load B = mat2 
  B00 = conj(mat2->r0.c0[idx_mat2]);
  B01 = conj(mat2->r0.c1[idx_mat2]);
  B02 = conj(mat2->r0.c2[idx_mat2]);
  B10 = conj(mat2->r1.c0[idx_mat2]);
  B11 = conj(mat2->r1.c1[idx_mat2]);
  B12 = conj(mat2->r1.c2[idx_mat2]);
  B20 = conj( ( B01 * B12 ) - ( B02 * B11) ) ;
  B21 = conj( ( B02 * B10 ) - ( B00 * B12) ) ;
  B22 = conj( ( B00 * B11 ) - ( B01 * B10) ) ;

  // mat3 = A * B = mat1 * mat3^DAG
  mat3->r0.c0[idx_mat3] = A00 * B00 + A01 * B10 + A02 * B20;
  mat3->r0.c1[idx_mat3] = A00 * B01 + A01 * B11 + A02 * B21;
  mat3->r0.c2[idx_mat3] = A00 * B02 + A01 * B12 + A02 * B22;
  mat3->r1.c0[idx_mat3] = A10 * B00 + A11 * B10 + A12 * B20;
  mat3->r1.c1[idx_mat3] = A10 * B01 + A11 * B11 + A12 * B21;
  mat3->r1.c2[idx_mat3] = A10 * B02 + A11 * B12 + A12 * B22;
  mat3->r2.c0[idx_mat3] = A20 * B00 + A21 * B10 + A22 * B20;
  mat3->r2.c1[idx_mat3] = A20 * B01 + A21 * B11 + A22 * B21;
  mat3->r2.c2[idx_mat3] = A20 * B02 + A21 * B12 + A22 * B22;

}


#pragma acc routine seq
static inline void mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(__restrict const su3_soa * const mat1, const int idx_mat1,
																																						__restrict const su3_soa * const mat2, const int idx_mat2,
																																						double fact,
																																						__restrict su3_soa * const mat3,   const int idx_mat3)
{
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];
  
  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];
    
  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_20 = mat2->r2.c0[idx_mat2];

  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_21 = mat2->r2.c1[idx_mat2];

  d_complex mat2_02 = mat2->r0.c2[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];
  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  // compute the outcome
  mat3->r0.c0[idx_mat3] += fact*(mat1_00 * mat2_00 + mat1_01 * mat2_10 
																 + mat1_02 * mat2_20) ;
  mat3->r0.c1[idx_mat3] += fact*(mat1_00 * mat2_01 + mat1_01 * mat2_11 
																 + mat1_02 * mat2_21) ;
  mat3->r0.c2[idx_mat3] += fact*(mat1_00 * mat2_02 + mat1_01 * mat2_12 
																 + mat1_02 * mat2_22) ;
    
  mat3->r1.c0[idx_mat3] += fact*(mat1_10 * mat2_00 + mat1_11 * mat2_10 
																 + mat1_12 * mat2_20) ;
  mat3->r1.c1[idx_mat3] += fact*(mat1_10 * mat2_01 + mat1_11 * mat2_11 
																 + mat1_12 * mat2_21) ;
  mat3->r1.c2[idx_mat3] += fact*(mat1_10 * mat2_02 + mat1_11 * mat2_12 
																 + mat1_12 * mat2_22) ;
    
  mat3->r2.c0[idx_mat3] += fact*conj( ( mat1_01 * mat1_12 ) 
																			- ( mat1_02 * mat1_11) ) ;
  mat3->r2.c1[idx_mat3] += fact*conj( ( mat1_02 * mat1_10 ) 
																			- ( mat1_00 * mat1_12) ) ;
  mat3->r2.c2[idx_mat3] += fact*conj( ( mat1_00 * mat1_11 ) 
																			- ( mat1_01 * mat1_10) ) ;

}

#pragma acc routine seq
static inline void conj_mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(__restrict const su3_soa * const mat1, const int idx_mat1,
																																								 __restrict const su3_soa * const mat2, const int idx_mat2,
																																								 double fact,
																																								 __restrict su3_soa * const mat3,   const int idx_mat3)
{

	d_complex mat1_00 = conj(mat1->r0.c0[idx_mat1]);
	d_complex mat1_01 = conj(mat1->r1.c0[idx_mat1]);
	d_complex mat1_02 = conj(mat1->r2.c0[idx_mat1]);

	d_complex mat1_10 = conj(mat1->r0.c1[idx_mat1]);
	d_complex mat1_11 = conj(mat1->r1.c1[idx_mat1]);
	d_complex mat1_12 = conj(mat1->r2.c1[idx_mat1]);

	d_complex mat1_20 = conj(mat1->r0.c2[idx_mat1]);
	d_complex mat1_21 = conj(mat1->r1.c2[idx_mat1]);
	d_complex mat1_22 = conj(mat1->r2.c2[idx_mat1]);

	d_complex mat2_00 = mat2->r0.c0[idx_mat2];
	d_complex mat2_10 = mat2->r1.c0[idx_mat2];
	d_complex mat2_20 = mat2->r2.c0[idx_mat2];

	d_complex mat2_01 = mat2->r0.c1[idx_mat2];
	d_complex mat2_11 = mat2->r1.c1[idx_mat2];
	d_complex mat2_21 = mat2->r2.c1[idx_mat2];

	d_complex mat2_02 = mat2->r0.c2[idx_mat1];
	d_complex mat2_12 = mat2->r1.c2[idx_mat1];
	d_complex mat2_22 = mat2->r2.c2[idx_mat1];

	// compute the outcome
	mat3->r0.c0[idx_mat3] += fact*(mat1_00 * mat2_00 + mat1_01 * mat2_10 
																 + mat1_02 * mat2_20) ;
	mat3->r0.c1[idx_mat3] += fact*(mat1_00 * mat2_01 + mat1_01 * mat2_11 
																 + mat1_02 * mat2_21) ;
	mat3->r0.c2[idx_mat3] += fact*(mat1_00 * mat2_02 + mat1_01 * mat2_12 
																 + mat1_02 * mat2_22) ;
	
	mat3->r1.c0[idx_mat3] += fact*(mat1_10 * mat2_00 + mat1_11 * mat2_10 
																 + mat1_12 * mat2_20) ;
	mat3->r1.c1[idx_mat3] += fact*(mat1_10 * mat2_01 + mat1_11 * mat2_11 
																 + mat1_12 * mat2_21) ;
	mat3->r1.c2[idx_mat3] += fact*(mat1_10 * mat2_02 + mat1_11 * mat2_12 
																 + mat1_12 * mat2_22) ;
	
	mat3->r2.c0[idx_mat3] += fact*conj( ( mat1_01 * mat1_12 ) 
																			- ( mat1_02 * mat1_11) ) ;
	mat3->r2.c1[idx_mat3] += fact*conj( ( mat1_02 * mat1_10 ) 
																			- ( mat1_00 * mat1_12) ) ;
	mat3->r2.c2[idx_mat3] += fact*conj( ( mat1_00 * mat1_11 ) 
																			- ( mat1_01 * mat1_10) ) ;

}

#pragma acc routine seq
static inline void mat1_times_conj_mat2_times_fact_addto_mat3_absent_stag_phases(__restrict const su3_soa * const mat1, const int idx_mat1,
																																								 __restrict const su3_soa * const mat2, const int idx_mat2,
																																								 double fact,
																																								 __restrict su3_soa * const mat3,   const int idx_mat3)
{

	d_complex mat1_00 = mat1->r0.c0[idx_mat1];
	d_complex mat1_01 = mat1->r0.c1[idx_mat1];
	d_complex mat1_02 = mat1->r0.c2[idx_mat1];

	d_complex mat1_10 = mat1->r1.c0[idx_mat1];
	d_complex mat1_11 = mat1->r1.c1[idx_mat1];
	d_complex mat1_12 = mat1->r1.c2[idx_mat1];

	d_complex mat2_00 = conj(mat2->r0.c0[idx_mat2]);
	d_complex mat2_10 = conj(mat2->r0.c1[idx_mat2]);
	d_complex mat2_20 = conj(mat2->r0.c2[idx_mat2]);
						     					   
	d_complex mat2_01 = conj(mat2->r1.c0[idx_mat2]);
	d_complex mat2_11 = conj(mat2->r1.c1[idx_mat2]);
	d_complex mat2_21 = conj(mat2->r1.c2[idx_mat2]);

	d_complex mat2_02 = conj(mat2->r2.c0[idx_mat2]);
	d_complex mat2_12 = conj(mat2->r2.c1[idx_mat2]);
	d_complex mat2_22 = conj(mat2->r2.c2[idx_mat2]);

	// compute the outcome
	mat3->r0.c0[idx_mat3] += fact*(mat1_00 * mat2_00 + mat1_01 * mat2_10 
																 + mat1_02 * mat2_20) ;
	mat3->r0.c1[idx_mat3] += fact*(mat1_00 * mat2_01 + mat1_01 * mat2_11 
																 + mat1_02 * mat2_21) ;
	mat3->r0.c2[idx_mat3] += fact*(mat1_00 * mat2_02 + mat1_01 * mat2_12 
																 + mat1_02 * mat2_22) ;
	
	mat3->r1.c0[idx_mat3] += fact*(mat1_10 * mat2_00 + mat1_11 * mat2_10 
																 + mat1_12 * mat2_20) ;
	mat3->r1.c1[idx_mat3] += fact*(mat1_10 * mat2_01 + mat1_11 * mat2_11 
																 + mat1_12 * mat2_21) ;
	mat3->r1.c2[idx_mat3] += fact*(mat1_10 * mat2_02 + mat1_11 * mat2_12 
																 + mat1_12 * mat2_22) ;
	
	mat3->r2.c0[idx_mat3] += fact*conj( ( mat1_01 * mat1_12 ) 
																			- ( mat1_02 * mat1_11) ) ;
	mat3->r2.c1[idx_mat3] += fact*conj( ( mat1_02 * mat1_10 ) 
																			- ( mat1_00 * mat1_12) ) ;
	mat3->r2.c2[idx_mat3] += fact*conj( ( mat1_00 * mat1_11 ) 
																			- ( mat1_01 * mat1_10) ) ;

}


// routine for the computation of the 3 matrices which contributes 
// to the right part of the staple
// mat4 = mat1 * hermitian_conjucate(mat2)* hermitian_conjucate(mat3)

#pragma acc routine seq
static inline void mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(  
																																											__restrict const su3_soa * const matnu1, const int idx_mat_nu1,
																																											__restrict const su3_soa * const matmu2, const int idx_mat_mu2,
																																											__restrict const su3_soa * const matnu3, const int idx_mat_nu3,
																																											__restrict su3_soa * const mat4,   const int idx_mat4)
{

	d_complex mat1_00 = matnu1->r0.c0[idx_mat_nu1];
	d_complex mat1_01 = matnu1->r0.c1[idx_mat_nu1];
	d_complex mat1_02 = matnu1->r0.c2[idx_mat_nu1];

	d_complex mat1_10 = matnu1->r1.c0[idx_mat_nu1];
	d_complex mat1_11 = matnu1->r1.c1[idx_mat_nu1];
	d_complex mat1_12 = matnu1->r1.c2[idx_mat_nu1];

	// construct (into the variables mat2_ij) the hermitian conjugate of the mat2 matrix
	d_complex mat2_00 = conj( matmu2->r0.c0[idx_mat_mu2] ) ;
	d_complex mat2_10 = conj( matmu2->r0.c1[idx_mat_mu2] ) ;
	d_complex mat2_20 = conj( matmu2->r0.c2[idx_mat_mu2] ) ;

	d_complex mat2_01 = conj( matmu2->r1.c0[idx_mat_mu2] ) ;
	d_complex mat2_11 = conj( matmu2->r1.c1[idx_mat_mu2] ) ;
	d_complex mat2_21 = conj( matmu2->r1.c2[idx_mat_mu2] ) ;

	// compute 3rd mat2 column from the first two
	d_complex mat2_02 = conj((mat2_10 * mat2_21 )-( mat2_20 * mat2_11));
	d_complex mat2_12 = conj((mat2_20 * mat2_01 )-( mat2_00 * mat2_21));
	d_complex mat2_22 = conj((mat2_00 * mat2_11 )-( mat2_10 * mat2_01));

	// compute the first two rows of the result of m1 * ~m2 and assign to m1c2
	d_complex mat1c2_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 
		+ mat1_02 * mat2_20 ;
	d_complex mat1c2_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 
		+ mat1_02 * mat2_21 ;
	d_complex mat1c2_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 
		+ mat1_02 * mat2_22 ;

	d_complex mat1c2_10 = mat1_10 * mat2_00 + mat1_11 * mat2_10 
		+ mat1_12 * mat2_20 ;
	d_complex mat1c2_11 = mat1_10 * mat2_01 + mat1_11 * mat2_11 
		+ mat1_12 * mat2_21 ;
	d_complex mat1c2_12 = mat1_10 * mat2_02 + mat1_11 * mat2_12 
		+ mat1_12 * mat2_22 ;

	// construct (into the variables mat2_ij) the hermitian conjugate of the mat3 matrix
	mat2_00 = conj( matnu3->r0.c0[idx_mat_nu3] ) ;
	mat2_10 = conj( matnu3->r0.c1[idx_mat_nu3] ) ;
	mat2_20 = conj( matnu3->r0.c2[idx_mat_nu3] ) ;

	mat2_01 = conj( matnu3->r1.c0[idx_mat_nu3] ) ;
	mat2_11 = conj( matnu3->r1.c1[idx_mat_nu3] ) ;
	mat2_21 = conj( matnu3->r1.c2[idx_mat_nu3] ) ;

	// compute 3rd mat3 column from the first two
	mat2_02 = conj( ( mat2_10 * mat2_21 ) - ( mat2_20 * mat2_11) ) ;
	mat2_12 = conj( ( mat2_20 * mat2_01 ) - ( mat2_00 * mat2_21) ) ;
	mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_10 * mat2_01) ) ;

	// compute the first two rows of the result of m1c2 * ~m3 and assign to m1
	mat1_00 = mat1c2_00 * mat2_00 + mat1c2_01 * mat2_10 
		+ mat1c2_02 * mat2_20 ;
	mat1_01 = mat1c2_00 * mat2_01 + mat1c2_01 * mat2_11 
		+ mat1c2_02 * mat2_21 ;
	mat1_02 = mat1c2_00 * mat2_02 + mat1c2_01 * mat2_12 
		+ mat1c2_02 * mat2_22 ;

	mat1_10 = mat1c2_10 * mat2_00 + mat1c2_11 * mat2_10 
		+ mat1c2_12 * mat2_20 ;
	mat1_11 = mat1c2_10 * mat2_01 + mat1c2_11 * mat2_11 
		+ mat1c2_12 * mat2_21 ;
	mat1_12 = mat1c2_10 * mat2_02 + mat1c2_11 * mat2_12 
		+ mat1c2_12 * mat2_22 ;

#ifdef PAR_TEMP
	double K_mu_nu_right=(matnu1->K.d[idx_mat_nu1])*(matmu2->K.d[idx_mat_mu2])*(matnu3->K.d[idx_mat_nu3]);
		
	// write results inside mat4
	mat4->r0.c0[idx_mat4] += K_mu_nu_right*C_ZERO * mat1_00;
	mat4->r0.c1[idx_mat4] += K_mu_nu_right*C_ZERO * mat1_01;
	mat4->r0.c2[idx_mat4] += K_mu_nu_right*C_ZERO * mat1_02;

	mat4->r1.c0[idx_mat4] += K_mu_nu_right*C_ZERO * mat1_10;
	mat4->r1.c1[idx_mat4] += K_mu_nu_right*C_ZERO * mat1_11;
	mat4->r1.c2[idx_mat4] += K_mu_nu_right*C_ZERO * mat1_12;

	mat4->r2.c0[idx_mat4] += K_mu_nu_right*C_ZERO * conj( ( mat1_01 * mat1_12 )
																												- ( mat1_02 * mat1_11) ) ;
	mat4->r2.c1[idx_mat4] += K_mu_nu_right*C_ZERO * conj( ( mat1_02 * mat1_10 )
																												- ( mat1_00 * mat1_12) ) ;
	mat4->r2.c2[idx_mat4] += K_mu_nu_right*C_ZERO * conj( ( mat1_00 * mat1_11 )
																												- ( mat1_01 * mat1_10) ) ;
#else
	// write results inside mat4
	mat4->r0.c0[idx_mat4] += C_ZERO * mat1_00;
	mat4->r0.c1[idx_mat4] += C_ZERO * mat1_01;
	mat4->r0.c2[idx_mat4] += C_ZERO * mat1_02;

	mat4->r1.c0[idx_mat4] += C_ZERO * mat1_10;
	mat4->r1.c1[idx_mat4] += C_ZERO * mat1_11;
	mat4->r1.c2[idx_mat4] += C_ZERO * mat1_12;

	mat4->r2.c0[idx_mat4] += C_ZERO * conj( ( mat1_01 * mat1_12 )
																					- ( mat1_02 * mat1_11) ) ;
	mat4->r2.c1[idx_mat4] += C_ZERO * conj( ( mat1_02 * mat1_10 )
																					- ( mat1_00 * mat1_12) ) ;
	mat4->r2.c2[idx_mat4] += C_ZERO * conj( ( mat1_00 * mat1_11 )
																					- ( mat1_01 * mat1_10) ) ;
#endif
}

// routine for the computation of the 3 matrices which contributes to the left part of the staple
// mat4 = hermitian_conjucate(mat1)* hermitian_conjucate(mat2) * mat3
#pragma acc routine seq
static inline void conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(   
																																											__restrict const su3_soa * const matnu1, const int idx_mat_nu1,
																																											__restrict const su3_soa * const matmu2, const int idx_mat_mu2,
																																											__restrict const su3_soa * const matnu3, const int idx_mat_nu3,
																																											__restrict su3_soa * const mat4,   const int idx_mat4)
{
	// construct (into the variables mat1_ij) the hermitian conjugate 
	// of the mat1 matrix  
	d_complex mat1_00 = conj( matnu1->r0.c0[idx_mat_nu1] ) ;
	d_complex mat1_10 = conj( matnu1->r0.c1[idx_mat_nu1] ) ;
	d_complex mat1_20 = conj( matnu1->r0.c2[idx_mat_nu1] ) ;

	d_complex mat1_01 = conj( matnu1->r1.c0[idx_mat_nu1] ) ;
	d_complex mat1_11 = conj( matnu1->r1.c1[idx_mat_nu1] ) ;
	d_complex mat1_21 = conj( matnu1->r1.c2[idx_mat_nu1] ) ;

	// compute 3rd mat1 columN from the first two
	d_complex mat1_02 = conj(( mat1_10 * mat1_21 )-(mat1_20 * mat1_11));
	d_complex mat1_12 = conj(( mat1_20 * mat1_01 )-(mat1_00 * mat1_21));
	// d_complex mat1_22=conj(( mat1_00 * mat1_11 )-( mat1_10 * mat1_01));
	// not used

	// construct (into the variables mat2_ij) the hermitian conjugate 
	// of the mat2 matrix
	d_complex mat2_00 = conj( matmu2->r0.c0[idx_mat_mu2] ) ;
	d_complex mat2_10 = conj( matmu2->r0.c1[idx_mat_mu2] ) ;
	d_complex mat2_20 = conj( matmu2->r0.c2[idx_mat_mu2] ) ;

	d_complex mat2_01 = conj( matmu2->r1.c0[idx_mat_mu2] ) ;
	d_complex mat2_11 = conj( matmu2->r1.c1[idx_mat_mu2] ) ;
	d_complex mat2_21 = conj( matmu2->r1.c2[idx_mat_mu2] ) ;

	// compute 3rd mat2 column from the first two
	d_complex mat2_02 = conj(( mat2_10 * mat2_21 )-(mat2_20 * mat2_11));
	d_complex mat2_12 = conj(( mat2_20 * mat2_01 )-(mat2_00 * mat2_21));
	d_complex mat2_22 = conj(( mat2_00 * mat2_11 )-(mat2_10 * mat2_01));

	// compute the first two rows of the result of ~m1 * ~m2 
	// and assign to mc1c2
	d_complex matc1c2_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 
		+ mat1_02 * mat2_20 ;
	d_complex matc1c2_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 
		+ mat1_02 * mat2_21 ;
	d_complex matc1c2_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 
		+ mat1_02 * mat2_22 ;

	d_complex matc1c2_10 = mat1_10 * mat2_00 + mat1_11 * mat2_10 
		+ mat1_12 * mat2_20 ;
	d_complex matc1c2_11 = mat1_10 * mat2_01 + mat1_11 * mat2_11 
		+ mat1_12 * mat2_21 ;
	d_complex matc1c2_12 = mat1_10 * mat2_02 + mat1_11 * mat2_12 
		+ mat1_12 * mat2_22 ;

	// construct (into the variables mat2_ij)  the mat3 matrix
	mat2_00 = matnu3->r0.c0[idx_mat_nu3] ;
	mat2_01 = matnu3->r0.c1[idx_mat_nu3] ;
	mat2_02 = matnu3->r0.c2[idx_mat_nu3] ;

	mat2_10 = matnu3->r1.c0[idx_mat_nu3] ;
	mat2_11 = matnu3->r1.c1[idx_mat_nu3] ;
	mat2_12 = matnu3->r1.c2[idx_mat_nu3] ;

	// compute 3rd mat3 column from the first two
	mat2_20 = conj( ( mat2_01 * mat2_12 ) - ( mat2_02 * mat2_11) ) ;
	mat2_21 = conj( ( mat2_02 * mat2_10 ) - ( mat2_00 * mat2_12) ) ;
	mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_01 * mat2_10) ) ;

	// compute the first two rows of the result of mc1c2 * m3 
	// and assign to m1
	mat1_00 = matc1c2_00 * mat2_00 + matc1c2_01 * mat2_10 
		+ matc1c2_02 * mat2_20 ;
	mat1_01 = matc1c2_00 * mat2_01 + matc1c2_01 * mat2_11 
		+ matc1c2_02 * mat2_21 ;
	mat1_02 = matc1c2_00 * mat2_02 + matc1c2_01 * mat2_12 
		+ matc1c2_02 * mat2_22 ;

	mat1_10 = matc1c2_10 * mat2_00 + matc1c2_11 * mat2_10 
		+ matc1c2_12 * mat2_20 ;
	mat1_11 = matc1c2_10 * mat2_01 + matc1c2_11 * mat2_11 
		+ matc1c2_12 * mat2_21 ;
	mat1_12 = matc1c2_10 * mat2_02 + matc1c2_11 * mat2_12 
		+ matc1c2_12 * mat2_22 ;
    
#ifdef PAR_TEMP
	double K_mu_nu_left=1;
	K_mu_nu_left=(matnu1->K.d[idx_mat_nu1])*(matmu2->K.d[idx_mat_mu2])*(matnu3->K.d[idx_mat_nu3]);
    
	// write results inside mat4
	mat4->r0.c0[idx_mat4] += K_mu_nu_left*C_ZERO * mat1_00;
	mat4->r0.c1[idx_mat4] += K_mu_nu_left*C_ZERO * mat1_01;
	mat4->r0.c2[idx_mat4] += K_mu_nu_left*C_ZERO * mat1_02;

	mat4->r1.c0[idx_mat4] += K_mu_nu_left*C_ZERO * mat1_10;
	mat4->r1.c1[idx_mat4] += K_mu_nu_left*C_ZERO * mat1_11;
	mat4->r1.c2[idx_mat4] += K_mu_nu_left*C_ZERO * mat1_12;

	mat4->r2.c0[idx_mat4] += K_mu_nu_left*C_ZERO * conj( ( mat1_01 * mat1_12 )
																											 - ( mat1_02 * mat1_11) ) ;
	mat4->r2.c1[idx_mat4] += K_mu_nu_left*C_ZERO * conj( ( mat1_02 * mat1_10 )
																											 - ( mat1_00 * mat1_12) ) ;
	mat4->r2.c2[idx_mat4] += K_mu_nu_left*C_ZERO * conj( ( mat1_00 * mat1_11 )
																											 - ( mat1_01 * mat1_10) ) ;
#else
	// write results inside mat4
	mat4->r0.c0[idx_mat4] += C_ZERO * mat1_00;
	mat4->r0.c1[idx_mat4] += C_ZERO * mat1_01;
	mat4->r0.c2[idx_mat4] += C_ZERO * mat1_02;

	mat4->r1.c0[idx_mat4] += C_ZERO * mat1_10;
	mat4->r1.c1[idx_mat4] += C_ZERO * mat1_11;
	mat4->r1.c2[idx_mat4] += C_ZERO * mat1_12;

	mat4->r2.c0[idx_mat4] += C_ZERO * conj( ( mat1_01 * mat1_12 )
																					- ( mat1_02 * mat1_11) ) ;
	mat4->r2.c1[idx_mat4] += C_ZERO * conj( ( mat1_02 * mat1_10 )
																					- ( mat1_00 * mat1_12) ) ;
	mat4->r2.c2[idx_mat4] += C_ZERO * conj( ( mat1_00 * mat1_11 )
																					- ( mat1_01 * mat1_10) ) ;
#endif
}

#pragma acc routine seq
static inline void mat1_times_mat2_into_tamat3(
																							 __restrict const su3_soa * const mat1,
																							 const int idx_mat1,
																							 __restrict const su3_soa * const mat2,
																							 const int idx_mat2,
																							 __restrict tamat_soa * const mat3,
																							 const int idx_mat3)
{
  // load the first two rows of mat1 (that is a link variable)
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];
  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];
  // compute the 3rd row of mat1 (that is a link variable)
  d_complex mat1_20 = conj(( mat1_01 * mat1_12 )-( mat1_02 * mat1_11));
  d_complex mat1_21 = conj(( mat1_02 * mat1_10 )-( mat1_00 * mat1_12));
  d_complex mat1_22 = conj(( mat1_00 * mat1_11 )-( mat1_01 * mat1_10));

  // load all the rows of mat2 (that is a staple variable)
  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_02 = mat2->r0.c2[idx_mat2];
  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];
  d_complex mat2_20 = mat2->r2.c0[idx_mat2];
  d_complex mat2_21 = mat2->r2.c1[idx_mat2];
  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  // compute first row of the product mat1 * mat2
  d_complex mat3_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 
		+ mat1_02 * mat2_20;
  d_complex mat3_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 
		+ mat1_02 * mat2_21;
  d_complex mat3_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 
		+ mat1_02 * mat2_22;
  // compute second row of the product mat1 * mat2 and save it
  // into reusable variables
  mat1_00 = mat1_10*mat2_00+mat1_11*mat2_10+mat1_12*mat2_20; // mat3_10 
  mat1_01 = mat1_10*mat2_01+mat1_11*mat2_11+mat1_12*mat2_21; // mat3_11
  mat1_02 = mat1_10*mat2_02+mat1_11*mat2_12+mat1_12*mat2_22; // mat3_12
  // compute third row of the product mat1 * mat2 and save it 
  // into reusable variables
  mat1_10 = mat1_20*mat2_00+mat1_21*mat2_10+mat1_22*mat2_20; // mat3_20
  mat1_11 = mat1_20*mat2_01+mat1_21*mat2_11+mat1_22*mat2_21; // mat3_21
  mat1_12 = mat1_20*mat2_02+mat1_21*mat2_12+mat1_22*mat2_22; // mat3_22

  mat3->c01[idx_mat3]  = 0.5*(mat3_01-conj(mat1_00));
  mat3->c02[idx_mat3]  = 0.5*(mat3_02-conj(mat1_10));
  mat3->c12[idx_mat3]  = 0.5*(mat1_02-conj(mat1_11));
  mat3->ic00[idx_mat3] = cimag(mat3_00)-
		ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12));
  mat3->ic11[idx_mat3] = cimag(mat1_01)-
		ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12));
}

#pragma acc routine seq
static inline void mat1_times_mat2_addto_tamat3(
																								__restrict const su3_soa * const mat1,
																								const int idx_mat1,
																								__restrict const su3_soa * const mat2,
																								const int idx_mat2,
																								__restrict tamat_soa * const mat3,
																								const int idx_mat3)
{
  // load the first two rows of mat1 (that is a link variable)
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];
  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];
  // compute the 3rd row of mat1 (that is a link variable)
  d_complex mat1_20 = conj(( mat1_01 * mat1_12 )-( mat1_02 * mat1_11));
  d_complex mat1_21 = conj(( mat1_02 * mat1_10 )-( mat1_00 * mat1_12));
  d_complex mat1_22 = conj(( mat1_00 * mat1_11 )-( mat1_01 * mat1_10));

  // load all the rows of mat2 (that is a staple variable)
  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_02 = mat2->r0.c2[idx_mat2];
  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];
  d_complex mat2_20 = mat2->r2.c0[idx_mat2];
  d_complex mat2_21 = mat2->r2.c1[idx_mat2];
  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  // compute first row of the product mat1 * mat2
  d_complex mat3_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 
		+ mat1_02 * mat2_20;
  d_complex mat3_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 
		+ mat1_02 * mat2_21;
  d_complex mat3_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 
		+ mat1_02 * mat2_22;
  // compute second row of the product mat1 * mat2 and save it
  // into reusable variables
  mat1_00 = mat1_10*mat2_00+mat1_11*mat2_10+mat1_12*mat2_20; // mat3_10 
  mat1_01 = mat1_10*mat2_01+mat1_11*mat2_11+mat1_12*mat2_21; // mat3_11
  mat1_02 = mat1_10*mat2_02+mat1_11*mat2_12+mat1_12*mat2_22; // mat3_12
  // compute third row of the product mat1 * mat2 and save it 
  // into reusable variables
  mat1_10 = mat1_20*mat2_00+mat1_21*mat2_10+mat1_22*mat2_20; // mat3_20
  mat1_11 = mat1_20*mat2_01+mat1_21*mat2_11+mat1_22*mat2_21; // mat3_21
  mat1_12 = mat1_20*mat2_02+mat1_21*mat2_12+mat1_22*mat2_22; // mat3_22

  mat3->c01[idx_mat3] += 0.5*(mat3_01-conj(mat1_00));
  mat3->c02[idx_mat3] += 0.5*(mat3_02-conj(mat1_10));
  mat3->c12[idx_mat3] += 0.5*(mat1_02-conj(mat1_11));
  mat3->ic00[idx_mat3]+= cimag(mat3_00)-
		ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12));
  mat3->ic11[idx_mat3]+= cimag(mat1_01)-
		ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12));
}

#pragma acc routine seq
static inline void RHO_times_mat1_times_mat2_into_tamat3(
																												 __restrict const su3_soa * const mat1, const int idx_mat1,
																												 __restrict const su3_soa * const mat2, const int idx_mat2,
																												 __restrict tamat_soa * const mat3,     const int idx_mat3, int istopo)
{
	const double RHO=(istopo)?(double)gl_topo_rho:(double)gl_stout_rho;
  // load the first two rows of mat1 (that is a link variable)
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];
  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];
  // compute the 3rd row of mat1 (that is a link variable)
  d_complex mat1_20 = conj(( mat1_01 * mat1_12 )-( mat1_02 * mat1_11));
  d_complex mat1_21 = conj(( mat1_02 * mat1_10 )-( mat1_00 * mat1_12));
  d_complex mat1_22 = conj(( mat1_00 * mat1_11 )-( mat1_01 * mat1_10));

  // load all the rows of mat2 (that is a staple variable)
  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_02 = mat2->r0.c2[idx_mat2];
  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];
  d_complex mat2_20 = mat2->r2.c0[idx_mat2];
  d_complex mat2_21 = mat2->r2.c1[idx_mat2];
  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  // compute first row of the product mat1 * mat2
  d_complex mat3_00 = mat1_00*mat2_00+mat1_01*mat2_10+mat1_02*mat2_20;
  d_complex mat3_01 = mat1_00*mat2_01+mat1_01*mat2_11+mat1_02*mat2_21;
  d_complex mat3_02 = mat1_00*mat2_02+mat1_01*mat2_12+mat1_02*mat2_22;
  // compute second row of the product mat1 * mat2 and save it 
  // into reusable variables
  mat1_00 = mat1_10*mat2_00+mat1_11*mat2_10+mat1_12*mat2_20; // mat3_10 
  mat1_01 = mat1_10*mat2_01+mat1_11*mat2_11+mat1_12*mat2_21; // mat3_11
  mat1_02 = mat1_10*mat2_02+mat1_11*mat2_12+mat1_12*mat2_22; // mat3_12
  // compute third row of the product mat1 * mat2 and save it 
  // into reusable variables
  mat1_10 = mat1_20*mat2_00+mat1_21*mat2_10+mat1_22*mat2_20; // mat3_20
  mat1_11 = mat1_20*mat2_01+mat1_21*mat2_11+mat1_22*mat2_21; // mat3_21
  mat1_12 = mat1_20*mat2_02+mat1_21*mat2_12+mat1_22*mat2_22; // mat3_22

  // in addition to multiply by rho i have to divide by C_ZERO
  // because staples here are: staples * C_ZERO 
  // then, we have to remove it
  double tmp = RHO/C_ZERO;
  mat3->c01[idx_mat3]  = tmp*(0.5*(mat3_01-conj(mat1_00)));
  mat3->c02[idx_mat3]  = tmp*(0.5*(mat3_02-conj(mat1_10)));
  mat3->c12[idx_mat3]  = tmp*(0.5*(mat1_02-conj(mat1_11)));
  mat3->ic00[idx_mat3] = tmp*((cimag(mat3_00) -ONE_BY_THREE*
															 (cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12))));
  mat3->ic11[idx_mat3] = tmp*((cimag(mat1_01)-ONE_BY_THREE*
															 (cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12))));
}

// mat1 = mat1 * integer factor
#pragma acc routine seq
static inline void  mat1_times_int_factor( 
																					__restrict su3_soa * const mat1,
																					const int idx_mat1,
																					int factor)
{

  mat1->r0.c0[idx_mat1] *= factor;
  mat1->r0.c1[idx_mat1] *= factor;
  mat1->r0.c2[idx_mat1] *= factor;

  mat1->r1.c0[idx_mat1] *= factor;
  mat1->r1.c1[idx_mat1] *= factor;
  mat1->r1.c2[idx_mat1] *= factor;

  // third row needs not to be multiplied
  // mat1->r2.c0[idx_mat1] *= factor;
  // mat1->r2.c1[idx_mat1] *= factor;
  // mat1->r2.c2[idx_mat1] *= factor;

}

// mat1 = mat1 * integer factor
#pragma acc routine seq
static inline void  mat1_times_double_factor( 
																						 __restrict su3_soa * const mat1,
																						 const int idx_mat1,
																						 double factor)
{

  mat1->r0.c0[idx_mat1] *= factor;
  mat1->r0.c1[idx_mat1] *= factor;
  mat1->r0.c2[idx_mat1] *= factor;

  mat1->r1.c0[idx_mat1] *= factor;
  mat1->r1.c1[idx_mat1] *= factor;
  mat1->r1.c2[idx_mat1] *= factor;

  // third row needs not to be multiplied
  // mat1->r2.c0[idx_mat1] *= factor;
  // mat1->r2.c1[idx_mat1] *= factor;
  // mat1->r2.c2[idx_mat1] *= factor;

}


// mat1 = mat1 * integer factor
#pragma acc routine seq
static inline void  gl3_times_int_factor( 
																				 __restrict su3_soa * const mgl3,
																				 const int idx_mgl3,
																				 int factor){

  mgl3->r0.c0[idx_mgl3] *= factor;
  mgl3->r0.c1[idx_mgl3] *= factor;
  mgl3->r0.c2[idx_mgl3] *= factor;

  mgl3->r1.c0[idx_mgl3] *= factor;
  mgl3->r1.c1[idx_mgl3] *= factor;
  mgl3->r1.c2[idx_mgl3] *= factor;

  mgl3->r2.c0[idx_mgl3] *= factor;
  mgl3->r2.c1[idx_mgl3] *= factor;
  mgl3->r2.c2[idx_mgl3] *= factor;

}

// mat1 = mat1 * double factor
#pragma acc routine seq
static inline void  gl3_times_double_factor( 
																						__restrict su3_soa * const mgl3,
																						const int idx_mgl3,
																						double factor){

  mgl3->r0.c0[idx_mgl3] *= factor;
  mgl3->r0.c1[idx_mgl3] *= factor;
  mgl3->r0.c2[idx_mgl3] *= factor;

  mgl3->r1.c0[idx_mgl3] *= factor;
  mgl3->r1.c1[idx_mgl3] *= factor;
  mgl3->r1.c2[idx_mgl3] *= factor;

  mgl3->r2.c0[idx_mgl3] *= factor;
  mgl3->r2.c1[idx_mgl3] *= factor;
  mgl3->r2.c2[idx_mgl3] *= factor;

}

#pragma acc routine seq
static inline void  gl3_dag_times_double_factor( 
																								__restrict su3_soa * const mgl3,
																								const int idx_mgl3,
																								double factor){
  
  mgl3->r0.c0[idx_mgl3] = conj(mgl3->r0.c0[idx_mgl3])*factor;
  d_complex m01 = mgl3->r0.c1[idx_mgl3];
  d_complex m02 = mgl3->r0.c2[idx_mgl3];
  mgl3->r0.c1[idx_mgl3] = conj(mgl3->r1.c0[idx_mgl3])*factor;
  mgl3->r0.c2[idx_mgl3] = conj(mgl3->r2.c0[idx_mgl3])*factor;
  
  mgl3->r1.c0[idx_mgl3] = conj(m01)*factor;
  mgl3->r1.c1[idx_mgl3] = conj(mgl3->r1.c1[idx_mgl3])*factor;  
  d_complex m12 = mgl3->r1.c2[idx_mgl3];
  mgl3->r1.c2[idx_mgl3] = conj(mgl3->r2.c1[idx_mgl3])*factor;

  mgl3->r2.c0[idx_mgl3] = conj(m02)*factor;
  mgl3->r2.c1[idx_mgl3] = conj(m12)*factor;
  mgl3->r2.c2[idx_mgl3] = conj(mgl3->r2.c2[idx_mgl3])*factor;
  
}


// compute SU(3) matrix trace
#pragma acc routine seq
static inline d_complex matrix_trace_absent_stag_phase(__restrict su3_soa * const loc_plaq,
																											 const int idx)
{
  d_complex loc_plaq_00 = loc_plaq->r0.c0[idx];
  d_complex loc_plaq_01 = loc_plaq->r0.c1[idx];
  d_complex loc_plaq_10 = loc_plaq->r1.c0[idx];
  d_complex loc_plaq_11 = loc_plaq->r1.c1[idx];
  // 1) I have to reconstruct 22 entry: third row has never 
  // been computed, it is not stored in the r2.c2[idxh] member,
  // then I can not load it.
  // 2) I reconstruct it exploiting unitarity
  d_complex loc_plaq_22 =  conj( ( loc_plaq_00 * loc_plaq_11 ) -
																 ( loc_plaq_01 * loc_plaq_10) ) ;
  return (loc_plaq_00 + loc_plaq_11 + loc_plaq_22);
}

#pragma acc routine seq
static inline void set_traces_to_value(dcomplex_soa * const tr_local_plaqs,
																			 int idxh,
																			 double value_r,
																			 double value_i)
{
  tr_local_plaqs->c[idxh] = - value_r + I * value_i;
}

#pragma acc routine seq
static inline double half_tr_thmat_squared(const __restrict thmat_soa * const mom,
																					 int idx_mom)
{
  d_complex  C = mom->c01[idx_mom];
  d_complex  D = mom->c02[idx_mom];
  d_complex  E = mom->c12[idx_mom];
  double    A = mom->rc00[idx_mom];
  double    B = mom->rc11[idx_mom];
  return A*A + B*B + A*B + creal(C)*creal(C) + cimag(C)*cimag(C) 
		+ creal(D)*creal(D) + cimag(D)*cimag(D) 
		+ creal(E)*creal(E) + cimag(E)*cimag(E);
  
}

#pragma acc routine seq
static inline double half_tr_tamat_squared(const __restrict tamat_soa * const mom,
																					 int idx_mom)
{
  d_complex  C = mom->c01[idx_mom];
  d_complex  D = mom->c02[idx_mom];
  d_complex  E = mom->c12[idx_mom];
  double    A = mom->ic00[idx_mom];
  double    B = mom->ic11[idx_mom];
  return A*A + B*B + A*B 
		+ creal(C)*creal(C) + cimag(C)*cimag(C) 
		+ creal(D)*creal(D) + cimag(D)*cimag(D) 
		+ creal(E)*creal(E) + cimag(E)*cimag(E);
  
}



#pragma acc routine seq
static inline void assign_zero_to_su3_soa_component(__restrict su3_soa * const matrix_comp,
																										int idx)
{
  matrix_comp->r0.c0[idx]=0.0+I*0.0;
  matrix_comp->r0.c1[idx]=0.0+I*0.0;
  matrix_comp->r0.c2[idx]=0.0+I*0.0;

  matrix_comp->r1.c0[idx]=0.0+I*0.0;
  matrix_comp->r1.c1[idx]=0.0+I*0.0;
  matrix_comp->r1.c2[idx]=0.0+I*0.0;

  matrix_comp->r2.c0[idx]=0.0+I*0.0;
  matrix_comp->r2.c1[idx]=0.0+I*0.0;
  matrix_comp->r2.c2[idx]=0.0+I*0.0;
}

#pragma acc routine seq
static inline void assign_su3_soa_to_su3_soa_component(__restrict const su3_soa * const matrix_comp_in,
																											 __restrict su3_soa * const matrix_comp_out,
																											 int idx)
{
	
  matrix_comp_out->r0.c0[idx] =  matrix_comp_in->r0.c0[idx];
  matrix_comp_out->r0.c1[idx] =  matrix_comp_in->r0.c1[idx];
  matrix_comp_out->r0.c2[idx] =  matrix_comp_in->r0.c2[idx];

  matrix_comp_out->r1.c0[idx] =  matrix_comp_in->r1.c0[idx];
  matrix_comp_out->r1.c1[idx] =  matrix_comp_in->r1.c1[idx];
  matrix_comp_out->r1.c2[idx] =  matrix_comp_in->r1.c2[idx];

}


#pragma acc routine seq
static inline void assign_su3_soa_to_su3_soa_component_trasl(__restrict const su3_soa * const matrix_comp_in,
																														 __restrict su3_soa * const matrix_comp_out,
																														 int idx,int idx1)
{
    
	matrix_comp_out->r0.c0[idx1] =  matrix_comp_in->r0.c0[idx];
	matrix_comp_out->r0.c1[idx1] =  matrix_comp_in->r0.c1[idx];
	matrix_comp_out->r0.c2[idx1] =  matrix_comp_in->r0.c2[idx];
    
	matrix_comp_out->r1.c0[idx1] =  matrix_comp_in->r1.c0[idx];
	matrix_comp_out->r1.c1[idx1] =  matrix_comp_in->r1.c1[idx];
	matrix_comp_out->r1.c2[idx1] =  matrix_comp_in->r1.c2[idx];
    
	/*
    matrix_comp_out->r2.c0[idx1] =  matrix_comp_in->r2.c0[idx];
    matrix_comp_out->r2.c1[idx1] =  matrix_comp_in->r2.c1[idx];
    matrix_comp_out->r2.c2[idx1] =  matrix_comp_in->r2.c2[idx];
	*/
}

#pragma acc routine seq
static inline void thmat1_plus_tamat2_times_factor_into_thmat1(__restrict thmat_soa * const thm1,
																															 __restrict const tamat_soa * const tam2,
																															 int idx,
																															 const double fact)
{
  d_complex ifact = 0.0 + I*fact;
  thm1->c01[idx]  -= ifact * tam2->c01[idx]; // complex
  thm1->c02[idx]  -= ifact * tam2->c02[idx]; // complex
  thm1->c12[idx]  -= ifact * tam2->c12[idx]; // complex
  thm1->rc00[idx] += fact * tam2->ic00[idx]; // double
  thm1->rc11[idx] += fact * tam2->ic11[idx]; // double
  
}

#pragma acc routine seq
static inline void conf_left_exp_multiply(__restrict const su3_soa * const cnf, const int idx_cnf,
																					__restrict const single_su3 * const  EXP, 
																					__restrict single_su3 * AUX,
																					__restrict single_su3 * AUX_RIS)
{
  
  // multiply: U_new = exp(i*delta*H) * U_old =>   cnf = exp * cnf 

  // loading first two rows
  AUX->comp[0][0] = cnf->r0.c0[idx_cnf];
  AUX->comp[0][1] = cnf->r0.c1[idx_cnf];
  AUX->comp[0][2] = cnf->r0.c2[idx_cnf];
  AUX->comp[1][0] = cnf->r1.c0[idx_cnf];
  AUX->comp[1][1] = cnf->r1.c1[idx_cnf];
  AUX->comp[1][2] = cnf->r1.c2[idx_cnf];
  // reconstructing third row
  AUX->comp[2][0] = conj(AUX->comp[0][1] * AUX->comp[1][2] 
												 - AUX->comp[0][2] * AUX->comp[1][1]);
  AUX->comp[2][1] = conj(AUX->comp[0][2] * AUX->comp[1][0] 
												 - AUX->comp[0][0] * AUX->comp[1][2]);
  AUX->comp[2][2] = conj(AUX->comp[0][0] * AUX->comp[1][1] 
												 - AUX->comp[0][1] * AUX->comp[1][0]);


  for(int r=0;r<2;r++) // we do not need third tow in the product
    for(int c=0;c<3;c++)
      AUX_RIS->comp[r][c] = (EXP->comp[r][0] * AUX->comp[0][c] 
														 + EXP->comp[r][1] * AUX->comp[1][c] 
														 + EXP->comp[r][2] * AUX->comp[2][c]);      
	// AUX_RIS->comp[r][c] = (AUX->comp[r][0] * EXP->comp[0][c] 
	// + AUX->comp[r][1] * EXP->comp[1][c] 
	// + AUX->comp[r][2] * EXP->comp[2][c]);
	
}

#pragma acc routine seq
static inline void mom_exp_times_conf_soloopenacc_loc(const __restrict thmat_soa * const mom,
																											__restrict su3_soa *cnf,
																											int idx,
																											double delta)
{
  // constructing the matrix M = i*delta*moment
  // load first part
  // minus sign on M10,M20 e M21 components because
  // moltiplying by 1.0I makes the matrix become anti-hermitian
  single_su3 AUX, RES;

  single_tamat QA;
  QA.ic00 = -mom->rc00[idx]*delta;
  QA.ic11 = -mom->rc11[idx]*delta;
  QA.c01  = -mom->c01[idx] *delta*1.0I;
  QA.c02  = -mom->c02[idx] *delta*1.0I;
  QA.c12  = -mom->c12[idx] *delta*1.0I;

  CH_exponential_antihermitian_nissalike(&RES,&QA);

	/* 
		 single_su3 MOM;
		 MOM.comp[0][0] = mom->rc00[idx] * (delta * 1.0I);
		 MOM.comp[0][1] = mom->c01[idx] * (delta * 1.0I);
		 MOM.comp[0][2] = mom->c02[idx] * (delta * 1.0I);
		 MOM.comp[1][0] = - conj(MOM.comp[0][1]);
		 MOM.comp[1][1] = mom->rc11[idx] * (delta * 1.0I);
		 MOM.comp[1][2] = mom->c12[idx] * (delta * 1.0I);
		 MOM.comp[2][0] = - conj(MOM.comp[0][2]);
		 MOM.comp[2][1] = - conj(MOM.comp[1][2]);
		 MOM.comp[2][2] = - MOM.comp[0][0] - MOM.comp[1][1];


		 taylor_exponential_su3(&RES,&MOM,&AUX);
  */
  // multiply: U_new = exp(i*delta*H) * U_old =>   cnf = RES * cnf 
	single_su3  AUX_RIS;

  // load first two rows
  AUX.comp[0][0] = cnf->r0.c0[idx];
  AUX.comp[0][1] = cnf->r0.c1[idx];
  AUX.comp[0][2] = cnf->r0.c2[idx];
  AUX.comp[1][0] = cnf->r1.c0[idx];
  AUX.comp[1][1] = cnf->r1.c1[idx];
  AUX.comp[1][2] = cnf->r1.c2[idx];

	// reconstruct third row
  AUX.comp[2][0] = conj(AUX.comp[0][1] * AUX.comp[1][2] 
												- AUX.comp[0][2] * AUX.comp[1][1]);
  AUX.comp[2][1] = conj(AUX.comp[0][2] * AUX.comp[1][0] 
												- AUX.comp[0][0] * AUX.comp[1][2]);
  AUX.comp[2][2] = conj(AUX.comp[0][0] * AUX.comp[1][1] 
												- AUX.comp[0][1] * AUX.comp[1][0]);


  for(int r=0;r<2;r++) // we do not need third row in the product
    for(int c=0;c<3;c++)
      AUX_RIS.comp[r][c] = (RES.comp[r][0] * AUX.comp[0][c] 
														+ RES.comp[r][1] * AUX.comp[1][c] 
														+ RES.comp[r][2] * AUX.comp[2][c]);      
      
  // matrix reunitarization

  // normalizing first row
  double NORM = creal(AUX_RIS.comp[0][0])*creal(AUX_RIS.comp[0][0])
		+ cimag(AUX_RIS.comp[0][0])*cimag(AUX_RIS.comp[0][0]) 
		+ creal(AUX_RIS.comp[0][1])*creal(AUX_RIS.comp[0][1])
		+ cimag(AUX_RIS.comp[0][1])*cimag(AUX_RIS.comp[0][1]) 
		+ creal(AUX_RIS.comp[0][2])*creal(AUX_RIS.comp[0][2])
		+ cimag(AUX_RIS.comp[0][2])*cimag(AUX_RIS.comp[0][2]);
  NORM = 1.0/sqrt(NORM);

  for(int c=0;c<3;c++)
    AUX_RIS.comp[0][c] *= NORM;

  // scalar product by the second and subtraction (orthogonalize)
  d_complex SCAL_PROD = conj(AUX_RIS.comp[0][0]) * AUX_RIS.comp[1][0] 
		+ conj(AUX_RIS.comp[0][1]) * AUX_RIS.comp[1][1] 
		+ conj(AUX_RIS.comp[0][2]) * AUX_RIS.comp[1][2];

  for(int c=0;c<3;c++)
    AUX_RIS.comp[1][c] -= SCAL_PROD * AUX_RIS.comp[0][c];  

  // normalizing second row
  NORM = creal(AUX_RIS.comp[1][0])*creal(AUX_RIS.comp[1][0])
		+ cimag(AUX_RIS.comp[1][0])*cimag(AUX_RIS.comp[1][0]) 
		+ creal(AUX_RIS.comp[1][1])*creal(AUX_RIS.comp[1][1]) 
		+ cimag(AUX_RIS.comp[1][1])*cimag(AUX_RIS.comp[1][1]) 
		+ creal(AUX_RIS.comp[1][2])*creal(AUX_RIS.comp[1][2])
		+ cimag(AUX_RIS.comp[1][2])*cimag(AUX_RIS.comp[1][2]);
  NORM = 1.0/sqrt(NORM);

  AUX_RIS.comp[1][0] *= NORM;
  AUX_RIS.comp[1][1] *= NORM;
  AUX_RIS.comp[1][2] *= NORM;

  cnf->r0.c0[idx] = AUX_RIS.comp[0][0];
  cnf->r0.c1[idx] = AUX_RIS.comp[0][1];
  cnf->r0.c2[idx] = AUX_RIS.comp[0][2];
  cnf->r1.c0[idx] = AUX_RIS.comp[1][0];
  cnf->r1.c1[idx] = AUX_RIS.comp[1][1];
  cnf->r1.c2[idx] = AUX_RIS.comp[1][2];

}

#pragma acc routine seq
static inline void mom_exp_times_conf_soloopenacc_loc_split(const __restrict thmat_soa * const mom,
																														const __restrict su3_soa *cnf_in,
																														__restrict su3_soa *cnf_out,
																														int idx,
																														double delta)
{
  // constructing the matrix M = i*delta*moment
  // load first part
  // minus sign on M10,M20 e M21 components because
  // moltiplying by 1.0I makes the matrix become anti-hermitian
	single_su3 AUX, RES;
	single_tamat QA;
  QA.ic00 = -mom->rc00[idx]*delta;
  QA.ic11 = -mom->rc11[idx]*delta;
  QA.c01  = -mom->c01[idx] *delta*1.0I;
  QA.c02  = -mom->c02[idx] *delta*1.0I;
  QA.c12  = -mom->c12[idx] *delta*1.0I;

  CH_exponential_antihermitian_nissalike(&RES,&QA);
	/*
		single_su3 MOM;
		MOM.comp[0][0] = mom->rc00[idx] * (delta * 1.0I);
		MOM.comp[0][1] = mom->c01[idx] * (delta * 1.0I);
		MOM.comp[0][2] = mom->c02[idx] * (delta * 1.0I);
		MOM.comp[1][0] = - conj(MOM.comp[0][1]);
		MOM.comp[1][1] = mom->rc11[idx] * (delta * 1.0I);
		MOM.comp[1][2] = mom->c12[idx] * (delta * 1.0I);
		MOM.comp[2][0] = - conj(MOM.comp[0][2]);
		MOM.comp[2][1] = - conj(MOM.comp[1][2]);
		MOM.comp[2][2] = - MOM.comp[0][0] - MOM.comp[1][1];


		// EXPONENTIATION

		taylor_exponential_su3(&RES,&MOM,&AUX);
	*/  
  // multiply: U_new = exp(i*delta*H) * U_old =>   cnf = RES * cnf 
	single_su3  AUX_RIS;

  // loading first two rows
  AUX.comp[0][0] = cnf_in->r0.c0[idx];
  AUX.comp[0][1] = cnf_in->r0.c1[idx];
  AUX.comp[0][2] = cnf_in->r0.c2[idx];
  AUX.comp[1][0] = cnf_in->r1.c0[idx];
  AUX.comp[1][1] = cnf_in->r1.c1[idx];
  AUX.comp[1][2] = cnf_in->r1.c2[idx];
  // reconstructing third row
  AUX.comp[2][0] = conj(AUX.comp[0][1] * AUX.comp[1][2] 
												- AUX.comp[0][2] * AUX.comp[1][1]);
  AUX.comp[2][1] = conj(AUX.comp[0][2] * AUX.comp[1][0] 
												- AUX.comp[0][0] * AUX.comp[1][2]);
  AUX.comp[2][2] = conj(AUX.comp[0][0] * AUX.comp[1][1] 
												- AUX.comp[0][1] * AUX.comp[1][0]);


  for(int r=0;r<2;r++) // we do not need third tow in the product
    for(int c=0;c<3;c++)
      AUX_RIS.comp[r][c] = (RES.comp[r][0] * AUX.comp[0][c] 
														+ RES.comp[r][1] * AUX.comp[1][c] 
														+ RES.comp[r][2] * AUX.comp[2][c]);      
      
  // matrix reunitarization

  // normalizing first row
  double NORM = creal(AUX_RIS.comp[0][0])*creal(AUX_RIS.comp[0][0])
		+ cimag(AUX_RIS.comp[0][0])*cimag(AUX_RIS.comp[0][0]) 
		+ creal(AUX_RIS.comp[0][1])*creal(AUX_RIS.comp[0][1])
		+ cimag(AUX_RIS.comp[0][1])*cimag(AUX_RIS.comp[0][1]) 
		+ creal(AUX_RIS.comp[0][2])*creal(AUX_RIS.comp[0][2])
		+ cimag(AUX_RIS.comp[0][2])*cimag(AUX_RIS.comp[0][2]);
  NORM = 1.0/sqrt(NORM);

  for(int c=0;c<3;c++)
    AUX_RIS.comp[0][c] *= NORM;

  // scalar product by the second and subtraction (orthogonalize)
  d_complex SCAL_PROD = conj(AUX_RIS.comp[0][0]) * AUX_RIS.comp[1][0] 
		+ conj(AUX_RIS.comp[0][1]) * AUX_RIS.comp[1][1] 
		+ conj(AUX_RIS.comp[0][2]) * AUX_RIS.comp[1][2];

  for(int c=0;c<3;c++)
    AUX_RIS.comp[1][c] -= SCAL_PROD * AUX_RIS.comp[0][c];  

  // normalizing second row
  NORM = creal(AUX_RIS.comp[1][0])*creal(AUX_RIS.comp[1][0])
		+ cimag(AUX_RIS.comp[1][0])*cimag(AUX_RIS.comp[1][0]) 
		+ creal(AUX_RIS.comp[1][1])*creal(AUX_RIS.comp[1][1]) 
		+ cimag(AUX_RIS.comp[1][1])*cimag(AUX_RIS.comp[1][1]) 
		+ creal(AUX_RIS.comp[1][2])*creal(AUX_RIS.comp[1][2])
		+ cimag(AUX_RIS.comp[1][2])*cimag(AUX_RIS.comp[1][2]);
  NORM = 1.0/sqrt(NORM);

  AUX_RIS.comp[1][0] *= NORM;
  AUX_RIS.comp[1][1] *= NORM;
  AUX_RIS.comp[1][2] *= NORM;

  cnf_out->r0.c0[idx] = AUX_RIS.comp[0][0];
  cnf_out->r0.c1[idx] = AUX_RIS.comp[0][1];
  cnf_out->r0.c2[idx] = AUX_RIS.comp[0][2];
  cnf_out->r1.c0[idx] = AUX_RIS.comp[1][0];
  cnf_out->r1.c1[idx] = AUX_RIS.comp[1][1];
  cnf_out->r1.c2[idx] = AUX_RIS.comp[1][2];

}

#ifdef MULTIDEVICE

void set_su3_soa_to_zero_bulk(__restrict su3_soa * const matrix);

void conf_times_staples_ta_part_bulk(__restrict const su3_soa * const u,
																		 __restrict const su3_soa * const loc_stap,
																		 __restrict tamat_soa * const tipdot);

void mom_sum_mult_bulk(__restrict thmat_soa * const mom,
											 const __restrict tamat_soa * const ipdot,
											 const double * factor,
											 int id_factor);

void mom_exp_times_conf_soloopenacc_bulk(__restrict const su3_soa * const conf_old,
																				 __restrict su3_soa * const conf_new,
																				 __restrict const thmat_soa * const mom,
																				 const double * factor, 
																				 // this is the delta vector, where Omelyan dts are stored
																				 int id_factor);

void set_su3_soa_to_zero_d3c(__restrict su3_soa * const matrix,
														 int offset3, int thickness3);

void conf_times_staples_ta_part_d3c(__restrict const su3_soa * const u,
																		__restrict const su3_soa * const loc_stap,
																		__restrict tamat_soa * const tipdot,
																		int offset3, int thickness3);

void mom_sum_mult_d3c(__restrict thmat_soa * const mom,
											const __restrict tamat_soa * const ipdot,
											const double * factor,
											int id_factor,
											int offset3, int thickness3);

void mom_exp_times_conf_soloopenacc_d3c(__restrict const su3_soa * const conf_old,
																				__restrict su3_soa * const conf_new,
																				__restrict const thmat_soa * const mom,
																				const double * factor, 
																				// this is the delta vector, where Omelyan dts are stored
																				int id_factor,
																				int offset3, int thickness3);

#endif

#ifdef PAR_TEMP
void add_defect_coeffs_to_staple(__restrict const su3_soa * const u,
																 __restrict su3_soa * const loc_stap);

#ifdef MULTIDEVICE
void add_defect_coeffs_to_staple_bulk(__restrict const su3_soa * const u,
																			__restrict su3_soa * const loc_stap);

void add_defect_coeffs_to_staple_d3c(__restrict const su3_soa * const u,
																		 __restrict su3_soa * const loc_stap,
																		 int offset, int thickness);
#endif
#endif


#endif
