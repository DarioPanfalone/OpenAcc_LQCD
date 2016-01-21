#ifndef SU3_UTILITIES_H_
#define SU3_UTILITIES_H_

#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/single_types.h"

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
void mult_conf_times_stag_phases( __restrict su3_soa * const u);

// multiply the whole configuration for the staggered phases field
// (all three lines)
void mult_gl3_soa_times_stag_phases( __restrict su3_soa * const u);

// multiply the whole configuration for the staggered phases field
void mult_conf_times_stag_phases_nodev( __restrict su3_soa * const u);




void set_su3_soa_to_zero( __restrict su3_soa * const matrix);

//copy matrix in into matrix out, this has to happen on the host
void set_su3_soa_to_su3_soa( __restrict su3_soa * const matrix_in,
        __restrict su3_soa * const matrix_out);


void conf_times_staples_ta_part(__restrict su3_soa * const u, // constant --> is not updated
__restrict su3_soa * const loc_stap, // constant --> is not updated
__restrict tamat_soa * const tipdot);

void RHO_times_conf_times_staples_ta_part(__restrict su3_soa * const u,        // constant --> is not updated
  __restrict su3_soa * const loc_stap, // constant --> is not updated
  __restrict tamat_soa * const tipdot);

void mom_sum_mult( __restrict thmat_soa * const mom,
        const __restrict tamat_soa * const ipdot, double * factor,
        int id_factor);

void kernel_acc_mom_exp_times_conf( __restrict su3_soa * const conf,
thmat_soa * const mom, // e' costante e qui dentro non viene modificato
double * factor, // questo e' il vettore delta dove sono contenuti tutti i dt richiesti nell'omelyan
int id_factor);

// Ora che abbiamo tolto i prodotti con le fasi staggered questo e' diventato un wrapper inutile... lo togliamo poi ... 
void mom_exp_times_conf_soloopenacc(
        __restrict  su3_soa * const tconf_acc,
 thmat_soa * const tmomenta, // e' costante e qui dentro non viene modificata
 double * tdelta,  int id_delta);

// reunitarize the conf by brute force
void unitarize_conf( __restrict su3_soa * const u);



#pragma acc routine seq
static inline void loc_unitarize_conf(__restrict su3_soa * const cnf,
				      const int idx_cnf){
  d_complex A00,A01,A02,A10,A11,A12;
  //d_complex A20,A21,A22;
  A00 = cnf->r0.c0[idx_cnf];
  A01 = cnf->r0.c1[idx_cnf];
  A02 = cnf->r0.c2[idx_cnf];
  A10 = cnf->r1.c0[idx_cnf];
  A11 = cnf->r1.c1[idx_cnf];
  A12 = cnf->r1.c2[idx_cnf];

  //normalizzo la prima riga
  double NORM = creal(A00)*creal(A00)+cimag(A00)*cimag(A00) + creal(A01)*creal(A01)+cimag(A01)*cimag(A01) + creal(A02)*creal(A02)+cimag(A02)*cimag(A02);
  NORM = 1.0/sqrt(NORM);
  A00 *= NORM;
  A01 *= NORM;
  A02 *= NORM;
  //faccio il prodotto scalare con la seconda e sottraggo (ortogonalizzo)
  d_complex SCAL_PROD = conj(A00) * A10 + conj(A01) * A11 + conj(A02) * A12;

  A10 -= SCAL_PROD * A00;
  A11 -= SCAL_PROD * A01;
  A12 -= SCAL_PROD * A02;

  //normalizzo la seconda riga
  NORM = creal(A10)*creal(A10)+cimag(A10)*cimag(A10) + creal(A11)*creal(A11)+cimag(A11)*cimag(A11) + creal(A12)*creal(A12)+cimag(A12)*cimag(A12);
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
static inline void    mat1_times_mat2_into_mat3_absent_stag_phases( __restrict su3_soa * const mat1,
								    const int idx_mat1,
								    __restrict su3_soa * const mat2,
								    const int idx_mat2,
								    __restrict su3_soa * const mat3,
								    const int idx_mat3) {
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

  //Compute 3rd mat2 row from the first two                    
  d_complex mat2_20 = conj( ( mat2_01 * mat2_12 ) - ( mat2_02 * mat2_11) ) ;
  d_complex mat2_21 = conj( ( mat2_02 * mat2_10 ) - ( mat2_00 * mat2_12) ) ;
  d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_01 * mat2_10) ) ;

  //Compute the first two rows of the solution
  mat3->r0.c0[idx_mat3] = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  mat3->r0.c1[idx_mat3] = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  mat3->r0.c2[idx_mat3] = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  mat3->r1.c0[idx_mat3] = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  mat3->r1.c1[idx_mat3] = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  mat3->r1.c2[idx_mat3] = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;
}

// mat1 = mat1 * mat2 
#pragma acc routine seq
static inline void    mat1_times_mat2_into_mat1_absent_stag_phases( __restrict su3_soa * const mat1,
								    const int idx_mat1,
								    __restrict su3_soa * const mat2,
								    const int idx_mat2) {
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

  //Compute 3rd mat2 row from the first two                    
  d_complex mat2_20 = conj( ( mat2_01 * mat2_12 ) - ( mat2_02 * mat2_11) ) ;
  d_complex mat2_21 = conj( ( mat2_02 * mat2_10 ) - ( mat2_00 * mat2_12) ) ;
  d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_01 * mat2_10) ) ;

  //Compute the first two rows of the solution
  mat1->r0.c0[idx_mat1] = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  mat1->r0.c1[idx_mat1] = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  mat1->r0.c2[idx_mat1] = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  mat1->r1.c0[idx_mat1] = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  mat1->r1.c1[idx_mat1] = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  mat1->r1.c2[idx_mat1] = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;
}

// mat1 = mat1 * hermitian_conjucate(mat2)
#pragma acc routine seq
static inline void    mat1_times_conj_mat2_into_mat1_absent_stag_phases( __restrict su3_soa * const mat1,
									 const int idx_mat1,
									 __restrict su3_soa * const mat2,
									 const int idx_mat2){
  
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];

  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];

  // construct (into the variables mat2_ij) the hermitian conjugate of the mat2 matrix
  d_complex mat2_00 = conj( mat2->r0.c0[idx_mat2] ) ;
  d_complex mat2_10 = conj( mat2->r0.c1[idx_mat2] ) ;
  d_complex mat2_20 = conj( mat2->r0.c2[idx_mat2] ) ;

  d_complex mat2_01 = conj( mat2->r1.c0[idx_mat2] ) ;
  d_complex mat2_11 = conj( mat2->r1.c1[idx_mat2] ) ;
  d_complex mat2_21 = conj( mat2->r1.c2[idx_mat2] ) ;

  //Compute 3rd mat2 column from the first two
  d_complex mat2_02 = conj( ( mat2_10 * mat2_21 ) - ( mat2_20 * mat2_11) ) ;
  d_complex mat2_12 = conj( ( mat2_20 * mat2_01 ) - ( mat2_00 * mat2_21) ) ;
  d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_10 * mat2_01) ) ;

  //Compute the first two rows of the solution
  mat1->r0.c0[idx_mat1] = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  mat1->r0.c1[idx_mat1] = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  mat1->r0.c2[idx_mat1] = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  mat1->r1.c0[idx_mat1] = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  mat1->r1.c1[idx_mat1] = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  mat1->r1.c2[idx_mat1] = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;
}

// Routine for the computation of the 3 matrices which contributes to the right part of the staple
// mat4 = mat1 * hermitian_conjucate(mat2)* hermitian_conjucate(mat3)

#pragma acc routine seq
static inline void    mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(  __restrict su3_soa * const matnu1,
											   const int idx_mat_nu1,
											   __restrict su3_soa * const matmu2,
											   const int idx_mat_mu2,
											   __restrict su3_soa * const matnu3,
											   const int idx_mat_nu3,
											   __restrict su3_soa * const mat4,
											   const int idx_mat4){
  
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

  //Compute 3rd mat2 column from the first two
  d_complex mat2_02 = conj( ( mat2_10 * mat2_21 ) - ( mat2_20 * mat2_11) ) ;
  d_complex mat2_12 = conj( ( mat2_20 * mat2_01 ) - ( mat2_00 * mat2_21) ) ;
  d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_10 * mat2_01) ) ;

  //Compute the first two rows of the result of m1 * ~m2 and assign to m1c2
  d_complex mat1c2_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  d_complex mat1c2_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  d_complex mat1c2_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  d_complex mat1c2_10 = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  d_complex mat1c2_11 = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  d_complex mat1c2_12 = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;

  // construct (into the variables mat2_ij) the hermitian conjugate of the mat3 matrix
  mat2_00 = conj( matnu3->r0.c0[idx_mat_nu3] ) ;
  mat2_10 = conj( matnu3->r0.c1[idx_mat_nu3] ) ;
  mat2_20 = conj( matnu3->r0.c2[idx_mat_nu3] ) ;

  mat2_01 = conj( matnu3->r1.c0[idx_mat_nu3] ) ;
  mat2_11 = conj( matnu3->r1.c1[idx_mat_nu3] ) ;
  mat2_21 = conj( matnu3->r1.c2[idx_mat_nu3] ) ;

  //Compute 3rd mat3 column from the first two
  mat2_02 = conj( ( mat2_10 * mat2_21 ) - ( mat2_20 * mat2_11) ) ;
  mat2_12 = conj( ( mat2_20 * mat2_01 ) - ( mat2_00 * mat2_21) ) ;
  mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_10 * mat2_01) ) ;

  //Compute the first two rows of the result of m1c2 * ~m3 and assign to m1
  mat1_00 = mat1c2_00 * mat2_00 + mat1c2_01 * mat2_10 + mat1c2_02 * mat2_20 ;
  mat1_01 = mat1c2_00 * mat2_01 + mat1c2_01 * mat2_11 + mat1c2_02 * mat2_21 ;
  mat1_02 = mat1c2_00 * mat2_02 + mat1c2_01 * mat2_12 + mat1c2_02 * mat2_22 ;

  mat1_10 = mat1c2_10 * mat2_00 + mat1c2_11 * mat2_10 + mat1c2_12 * mat2_20 ;
  mat1_11 = mat1c2_10 * mat2_01 + mat1c2_11 * mat2_11 + mat1c2_12 * mat2_21 ;
  mat1_12 = mat1c2_10 * mat2_02 + mat1c2_11 * mat2_12 + mat1c2_12 * mat2_22 ;

  //Write results inside mat4
  mat4->r0.c0[idx_mat4] += C_ZERO * mat1_00;
  mat4->r0.c1[idx_mat4] += C_ZERO * mat1_01;
  mat4->r0.c2[idx_mat4] += C_ZERO * mat1_02;

  mat4->r1.c0[idx_mat4] += C_ZERO * mat1_10;
  mat4->r1.c1[idx_mat4] += C_ZERO * mat1_11;
  mat4->r1.c2[idx_mat4] += C_ZERO * mat1_12;

  mat4->r2.c0[idx_mat4] += C_ZERO * conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  mat4->r2.c1[idx_mat4] += C_ZERO * conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  mat4->r2.c2[idx_mat4] += C_ZERO * conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
}

// Routine for the computation of the 3 matrices which contributes to the left part of the staple
// mat4 = hermitian_conjucate(mat1)* hermitian_conjucate(mat2) * mat3
#pragma acc routine seq
static inline void    conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(   __restrict su3_soa * const matnu1,
											   const int idx_mat_nu1,
											   __restrict su3_soa * const matmu2,
											   const int idx_mat_mu2,
											   __restrict su3_soa * const matnu3,
											   const int idx_mat_nu3,
											   __restrict su3_soa * const mat4,
											   const int idx_mat4){
  // construct (into the variables mat1_ij) the hermitian conjugate of the mat1 matrix  
  d_complex mat1_00 = conj( matnu1->r0.c0[idx_mat_nu1] ) ;
  d_complex mat1_10 = conj( matnu1->r0.c1[idx_mat_nu1] ) ;
  d_complex mat1_20 = conj( matnu1->r0.c2[idx_mat_nu1] ) ;

  d_complex mat1_01 = conj( matnu1->r1.c0[idx_mat_nu1] ) ;
  d_complex mat1_11 = conj( matnu1->r1.c1[idx_mat_nu1] ) ;
  d_complex mat1_21 = conj( matnu1->r1.c2[idx_mat_nu1] ) ;

  //Compute 3rd mat1 column from the first two
  d_complex mat1_02 = conj( ( mat1_10 * mat1_21 ) - ( mat1_20 * mat1_11) ) ;
  d_complex mat1_12 = conj( ( mat1_20 * mat1_01 ) - ( mat1_00 * mat1_21) ) ;
  //d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_10 * mat1_01) ) ;//not used

  // construct (into the variables mat2_ij) the hermitian conjugate of the mat2 matrix
  d_complex mat2_00 = conj( matmu2->r0.c0[idx_mat_mu2] ) ;
  d_complex mat2_10 = conj( matmu2->r0.c1[idx_mat_mu2] ) ;
  d_complex mat2_20 = conj( matmu2->r0.c2[idx_mat_mu2] ) ;

  d_complex mat2_01 = conj( matmu2->r1.c0[idx_mat_mu2] ) ;
  d_complex mat2_11 = conj( matmu2->r1.c1[idx_mat_mu2] ) ;
  d_complex mat2_21 = conj( matmu2->r1.c2[idx_mat_mu2] ) ;

  //Compute 3rd mat2 column from the first two
  d_complex mat2_02 = conj( ( mat2_10 * mat2_21 ) - ( mat2_20 * mat2_11) ) ;
  d_complex mat2_12 = conj( ( mat2_20 * mat2_01 ) - ( mat2_00 * mat2_21) ) ;
  d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_10 * mat2_01) ) ;

  //Compute the first two rows of the result of ~m1 * ~m2 and assign to mc1c2
  d_complex matc1c2_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  d_complex matc1c2_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  d_complex matc1c2_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  d_complex matc1c2_10 = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  d_complex matc1c2_11 = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  d_complex matc1c2_12 = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;

  // construct (into the variables mat2_ij)  the mat3 matrix
  mat2_00 = matnu3->r0.c0[idx_mat_nu3] ;
  mat2_01 = matnu3->r0.c1[idx_mat_nu3] ;
  mat2_02 = matnu3->r0.c2[idx_mat_nu3] ;

  mat2_10 = matnu3->r1.c0[idx_mat_nu3] ;
  mat2_11 = matnu3->r1.c1[idx_mat_nu3] ;
  mat2_12 = matnu3->r1.c2[idx_mat_nu3] ;

  //Compute 3rd mat3 column from the first two
  mat2_20 = conj( ( mat2_01 * mat2_12 ) - ( mat2_02 * mat2_11) ) ;
  mat2_21 = conj( ( mat2_02 * mat2_10 ) - ( mat2_00 * mat2_12) ) ;
  mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_01 * mat2_10) ) ;

  //Compute the first two rows of the result of mc1c2 * m3 and assign to m1
  mat1_00 = matc1c2_00 * mat2_00 + matc1c2_01 * mat2_10 + matc1c2_02 * mat2_20 ;
  mat1_01 = matc1c2_00 * mat2_01 + matc1c2_01 * mat2_11 + matc1c2_02 * mat2_21 ;
  mat1_02 = matc1c2_00 * mat2_02 + matc1c2_01 * mat2_12 + matc1c2_02 * mat2_22 ;

  mat1_10 = matc1c2_10 * mat2_00 + matc1c2_11 * mat2_10 + matc1c2_12 * mat2_20 ;
  mat1_11 = matc1c2_10 * mat2_01 + matc1c2_11 * mat2_11 + matc1c2_12 * mat2_21 ;
  mat1_12 = matc1c2_10 * mat2_02 + matc1c2_11 * mat2_12 + matc1c2_12 * mat2_22 ;

  //Write results inside mat4
  mat4->r0.c0[idx_mat4] += C_ZERO * mat1_00;
  mat4->r0.c1[idx_mat4] += C_ZERO * mat1_01;
  mat4->r0.c2[idx_mat4] += C_ZERO * mat1_02;

  mat4->r1.c0[idx_mat4] += C_ZERO * mat1_10;
  mat4->r1.c1[idx_mat4] += C_ZERO * mat1_11;
  mat4->r1.c2[idx_mat4] += C_ZERO * mat1_12;

  mat4->r2.c0[idx_mat4] += C_ZERO * conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  mat4->r2.c1[idx_mat4] += C_ZERO * conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  mat4->r2.c2[idx_mat4] += C_ZERO * conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
}

#pragma acc routine seq
static inline void mat1_times_mat2_into_tamat3(__restrict su3_soa * const mat1,
					       const int idx_mat1,
					       __restrict su3_soa * const mat2,
					       const int idx_mat2,
					       __restrict tamat_soa * const mat3,
					       const int idx_mat3){
  //Load the first two rows of mat1 (that is a link variable)
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];
  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];
  //Compute the 3rd row of mat1 (that is a link variable)
  d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;

  //Load all the rows of mat2 (that is a staple variable)
  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_02 = mat2->r0.c2[idx_mat2];
  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];
  d_complex mat2_20 = mat2->r2.c0[idx_mat2];
  d_complex mat2_21 = mat2->r2.c1[idx_mat2];
  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  // Compute first row of the product mat1 * mat2
  d_complex mat3_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20;
  d_complex mat3_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21;
  d_complex mat3_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22;
  // Compute second row of the product mat1 * mat2 and save it into reusable variables
  mat1_00 = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20; // mat3_10 
  mat1_01 = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21; // mat3_11
  mat1_02 = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22; // mat3_12
  // Compute third row of the product mat1 * mat2 and save it into reusable variables
  mat1_10 = mat1_20 * mat2_00 + mat1_21 * mat2_10 + mat1_22 * mat2_20; // mat3_20
  mat1_11 = mat1_20 * mat2_01 + mat1_21 * mat2_11 + mat1_22 * mat2_21; // mat3_21
  mat1_12 = mat1_20 * mat2_02 + mat1_21 * mat2_12 + mat1_22 * mat2_22; // mat3_22

  mat3->c01[idx_mat3]  = 0.5*(mat3_01-conj(mat1_00));
  mat3->c02[idx_mat3]  = 0.5*(mat3_02-conj(mat1_10));
  mat3->c12[idx_mat3]  = 0.5*(mat1_02-conj(mat1_11));
  mat3->rc00[idx_mat3] = cimag(mat3_00)-ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12));
  mat3->rc11[idx_mat3] = cimag(mat1_01)-ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12));
}

#pragma acc routine seq
static inline void RHO_times_mat1_times_mat2_into_tamat3(__restrict su3_soa * const mat1,
							 const int idx_mat1,
							 __restrict su3_soa * const mat2,
							 const int idx_mat2,
							 __restrict tamat_soa * const mat3,
							 const int idx_mat3){
  //Load the first two rows of mat1 (that is a link variable)
  d_complex mat1_00 = mat1->r0.c0[idx_mat1];
  d_complex mat1_01 = mat1->r0.c1[idx_mat1];
  d_complex mat1_02 = mat1->r0.c2[idx_mat1];
  d_complex mat1_10 = mat1->r1.c0[idx_mat1];
  d_complex mat1_11 = mat1->r1.c1[idx_mat1];
  d_complex mat1_12 = mat1->r1.c2[idx_mat1];
  //Compute the 3rd row of mat1 (that is a link variable)
  d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;

  //Load all the rows of mat2 (that is a staple variable)
  d_complex mat2_00 = mat2->r0.c0[idx_mat2];
  d_complex mat2_01 = mat2->r0.c1[idx_mat2];
  d_complex mat2_02 = mat2->r0.c2[idx_mat2];
  d_complex mat2_10 = mat2->r1.c0[idx_mat2];
  d_complex mat2_11 = mat2->r1.c1[idx_mat2];
  d_complex mat2_12 = mat2->r1.c2[idx_mat2];
  d_complex mat2_20 = mat2->r2.c0[idx_mat2];
  d_complex mat2_21 = mat2->r2.c1[idx_mat2];
  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  // Compute first row of the product mat1 * mat2
  d_complex mat3_00 = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20;
  d_complex mat3_01 = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21;
  d_complex mat3_02 = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22;
  // Compute second row of the product mat1 * mat2 and save it into reusable variables
  mat1_00 = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20; // mat3_10 
  mat1_01 = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21; // mat3_11
  mat1_02 = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22; // mat3_12
  // Compute third row of the product mat1 * mat2 and save it into reusable variables
  mat1_10 = mat1_20 * mat2_00 + mat1_21 * mat2_10 + mat1_22 * mat2_20; // mat3_20
  mat1_11 = mat1_20 * mat2_01 + mat1_21 * mat2_11 + mat1_22 * mat2_21; // mat3_21
  mat1_12 = mat1_20 * mat2_02 + mat1_21 * mat2_12 + mat1_22 * mat2_22; // mat3_22

  // oltre a moltiplicare per RHO devo anche dividere per C_ZERO
  // perche' le staples che entrano qui dentro sono le staples * C_ZERO --> lo devo togliere!!
  double tmp = RHO/C_ZERO;
  mat3->c01[idx_mat3]  = tmp*(0.5*(mat3_01-conj(mat1_00)));
  mat3->c02[idx_mat3]  = tmp*(0.5*(mat3_02-conj(mat1_10)));
  mat3->c12[idx_mat3]  = tmp*(0.5*(mat1_02-conj(mat1_11)));
  mat3->rc00[idx_mat3] = tmp*((cimag(mat3_00)-ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12))));
  mat3->rc11[idx_mat3] = tmp*((cimag(mat1_01)-ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12))));
}

// mat1 = mat1 * integer factor
#pragma acc routine seq
static inline void   mat1_times_int_factor( __restrict su3_soa * const mat1,
					   const int idx_mat1,
					   int factor){

  mat1->r0.c0[idx_mat1] *= factor;
  mat1->r0.c1[idx_mat1] *= factor;
  mat1->r0.c2[idx_mat1] *= factor;

  mat1->r1.c0[idx_mat1] *= factor;
  mat1->r1.c1[idx_mat1] *= factor;
  mat1->r1.c2[idx_mat1] *= factor;

  //Third row needs not to be multiplied
  //  mat1->r2.c0[idx_mat1] *= factor;
  //  mat1->r2.c1[idx_mat1] *= factor;
  //  mat1->r2.c2[idx_mat1] *= factor;

}

// mat1 = mat1 * integer factor
#pragma acc routine seq
static inline void   gl3_times_int_factor( __restrict su3_soa * const mgl3,
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

// calcola la traccia della matrice di su3
#pragma acc routine seq
static inline d_complex matrix_trace_absent_stag_phase(__restrict su3_soa * const loc_plaq,
						       const int idx){
  d_complex loc_plaq_00 = loc_plaq->r0.c0[idx];
  d_complex loc_plaq_01 = loc_plaq->r0.c1[idx];
  d_complex loc_plaq_10 = loc_plaq->r1.c0[idx];
  d_complex loc_plaq_11 = loc_plaq->r1.c1[idx];
  // 1) devo ricostruire per forza l'elemento 22: la terza riga non la ho calcolata mai, quindi dove sembrerebbe doverci essere in realta' non c'e niente, quindi non la posso caricare
  // 2) carico quello che mi serve e ricostruisco
  d_complex loc_plaq_22 =  conj( ( loc_plaq_00 * loc_plaq_11 ) - ( loc_plaq_01 * loc_plaq_10) ) ;
  return (loc_plaq_00 + loc_plaq_11 + loc_plaq_22);
}

#pragma acc routine seq
static inline void set_traces_to_value( dcomplex_soa * const tr_local_plaqs,
				        int idxh,
					double value_r,
					double value_i){

  tr_local_plaqs->c[idxh] = - value_r + I * value_i;

}

#pragma acc routine seq
static inline double half_tr_thmat_squared( const __restrict thmat_soa * const mom,
					    int idx_mom){
  d_complex  C = mom->c01[idx_mom];
  d_complex  D = mom->c02[idx_mom];
  d_complex  E = mom->c12[idx_mom];
  double    A = mom->rc00[idx_mom];
  double    B = mom->rc11[idx_mom];
  return A*A + B*B + A*B + creal(C)*creal(C) + cimag(C)*cimag(C) + creal(D)*creal(D) + cimag(D)*cimag(D) + creal(E)*creal(E) + cimag(E)*cimag(E);
  
}

#pragma acc routine seq
static inline void assign_zero_to_su3_soa_component(__restrict su3_soa * const matrix_comp,
						    int idx){
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
static inline void assign_su3_soa_to_su3_soa_component(__restrict su3_soa * const matrix_comp_in,
						       __restrict su3_soa * const matrix_comp_out,
						       int idx){

  matrix_comp_out->r0.c0[idx] =  matrix_comp_in->r0.c0[idx];
  matrix_comp_out->r0.c1[idx] =  matrix_comp_in->r0.c1[idx];
  matrix_comp_out->r0.c2[idx] =  matrix_comp_in->r0.c2[idx];

  matrix_comp_out->r1.c0[idx] =  matrix_comp_in->r1.c0[idx];
  matrix_comp_out->r1.c1[idx] =  matrix_comp_in->r1.c1[idx];
  matrix_comp_out->r1.c2[idx] =  matrix_comp_in->r1.c2[idx];

}

#pragma acc routine seq
static inline void thmat1_plus_tamat2_times_factor_into_thmat1(__restrict thmat_soa * const thm1,
							       const __restrict tamat_soa * const tam2,
							       int idx,
							       const double fact){
  d_complex ifact = 0.0 + I*fact;
  thm1->c01[idx]  -= ifact * tam2->c01[idx]; // complex
  thm1->c02[idx]  -= ifact * tam2->c02[idx]; // complex
  thm1->c12[idx]  -= ifact * tam2->c12[idx]; // complex
  thm1->rc00[idx] += fact * tam2->rc00[idx];  // double
  thm1->rc11[idx] += fact * tam2->rc11[idx];  // double
  
}

#pragma acc routine seq
static inline void extract_mom(const __restrict thmat_soa * const mom,
			       int idx_mom,
			       double delta,
			       single_su3 * M){
  //COSTRUISCO LA MATRICE M = i*delta*momento
  //leggo la prima parte
  // il segno meno sulle componenti M10,M20 e M21 c'e' perche' dopo aver
  // moltiplicato per 1.0I la matrice diventa anti-hermitiana
  M->comp[0][0] = mom->rc00[idx_mom] * (delta * 1.0I);
  M->comp[0][1] = mom->c01[idx_mom] * (delta * 1.0I);
  M->comp[0][2] = mom->c02[idx_mom] * (delta * 1.0I);
  M->comp[1][0] = - conj(M->comp[0][1]);
  M->comp[1][1] = mom->rc11[idx_mom] * (delta * 1.0I);
  M->comp[1][2] = mom->c12[idx_mom] * (delta * 1.0I);
  M->comp[2][0] = - conj(M->comp[0][2]);
  M->comp[2][1] = - conj(M->comp[1][2]);
  M->comp[2][2] = - M->comp[0][0] - M->comp[1][1];
}

#pragma acc routine seq
static inline void matrix_exp_openacc(const __restrict single_su3 * const MOM,
				      __restrict single_su3 * AUX,
				      __restrict single_su3 * RES){
  // exp x = 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5))))
  // first iteration
  // ris=1+x/5
  for(int r=0;r<3;r++){
    for(int c=0;c<3;c++)
      RES->comp[r][c] = MOM->comp[r][c] * 0.2;
    RES->comp[r][r] = RES->comp[r][r] + 1.0;
  }
  // second iteration
  // ris=1.0+x/4*(1+x/5)
  for(int r=0;r<3;r++){
    for(int c=0;c<3;c++)
      AUX->comp[r][c] = (MOM->comp[r][0] * RES->comp[0][c] + MOM->comp[r][1] * RES->comp[1][c] + MOM->comp[r][2] * RES->comp[2][c]) * 0.25;
    AUX->comp[r][r] = AUX->comp[r][r] + 1.0;
  }
  // third iteration
  // ris=1.0+x/3.0*(1.0+x/4*(1+x/5))
  for(int r=0;r<3;r++){
    for(int c=0;c<3;c++)
      RES->comp[r][c] = (MOM->comp[r][0] * AUX->comp[0][c] + MOM->comp[r][1] * AUX->comp[1][c] + MOM->comp[r][2] * AUX->comp[2][c]) * ONE_BY_THREE;
    RES->comp[r][r] = RES->comp[r][r] + 1.0;
  }
  // fourth iteration
  // ris=1.0+x/2.0*(1.0+x/3.0*(1.0+x/4*(1+x/5)))
  for(int r=0;r<3;r++){
    for(int c=0;c<3;c++)
      AUX->comp[r][c] = (MOM->comp[r][0] * RES->comp[0][c] + MOM->comp[r][1] * RES->comp[1][c] + MOM->comp[r][2] * RES->comp[2][c]) * 0.5;
    AUX->comp[r][r] = AUX->comp[r][r] + 1.0;
  }
  // fifth iteration
  // ris=1.0+x*(1.0+x/2.0*(1.0+x/3.0*(1.0+x/4*(1+x/5))))
  for(int r=0;r<3;r++){
    for(int c=0;c<3;c++)
      RES->comp[r][c] = (MOM->comp[r][0] * AUX->comp[0][c] + MOM->comp[r][1] * AUX->comp[1][c] + MOM->comp[r][2] * AUX->comp[2][c]);
    RES->comp[r][r] = RES->comp[r][r] + 1.0;
  }
}

#pragma acc routine seq
static inline void conf_left_exp_multiply(__restrict su3_soa * const cnf,  // e' costante e qui dentro non viene modificata
					  const int idx_cnf,
					  const __restrict single_su3 * const  EXP, 
					  __restrict single_su3 * AUX,
					  __restrict single_su3 * AUX_RIS){
  
  //Multiply: U_new = exp(i*delta*H) * U_old =>   cnf = EXP * cnf 

  // leggo le prime due righe
  AUX->comp[0][0] = cnf->r0.c0[idx_cnf];
  AUX->comp[0][1] = cnf->r0.c1[idx_cnf];
  AUX->comp[0][2] = cnf->r0.c2[idx_cnf];
  AUX->comp[1][0] = cnf->r1.c0[idx_cnf];
  AUX->comp[1][1] = cnf->r1.c1[idx_cnf];
  AUX->comp[1][2] = cnf->r1.c2[idx_cnf];
  // ricostruisco la terza
  AUX->comp[2][0] = conj(AUX->comp[0][1] * AUX->comp[1][2] - AUX->comp[0][2] * AUX->comp[1][1]);
  AUX->comp[2][1] = conj(AUX->comp[0][2] * AUX->comp[1][0] - AUX->comp[0][0] * AUX->comp[1][2]);
  AUX->comp[2][2] = conj(AUX->comp[0][0] * AUX->comp[1][1] - AUX->comp[0][1] * AUX->comp[1][0]);


  for(int r=0;r<2;r++) // qui il loop va fino solo a 2 perche non mi serve calcolare anche la terza riga del prodotto
    for(int c=0;c<3;c++)
      AUX_RIS->comp[r][c] = (EXP->comp[r][0] * AUX->comp[0][c] + EXP->comp[r][1] * AUX->comp[1][c] + EXP->comp[r][2] * AUX->comp[2][c]);      
//      AUX_RIS->comp[r][c] = (AUX->comp[r][0] * EXP->comp[0][c] + AUX->comp[r][1] * EXP->comp[1][c] + AUX->comp[r][2] * EXP->comp[2][c]);
      
}

#pragma acc routine seq
static inline void project_on_su3(__restrict su3_soa * const cnf,
				  const int idx_cnf,
				  __restrict single_su3 *  AUX){
  
  //normalizzo la prima riga
  double NORM = creal(AUX->comp[0][0])*creal(AUX->comp[0][0])+cimag(AUX->comp[0][0])*cimag(AUX->comp[0][0]) + creal(AUX->comp[0][1])*creal(AUX->comp[0][1])+cimag(AUX->comp[0][1])*cimag(AUX->comp[0][1]) + creal(AUX->comp[0][2])*creal(AUX->comp[0][2])+cimag(AUX->comp[0][2])*cimag(AUX->comp[0][2]);
  NORM = 1.0/sqrt(NORM);

  for(int c=0;c<3;c++)
    AUX->comp[0][c] *= NORM;

  //faccio il prodotto scalare con la seconda e sottraggo (ortogonalizzo)
  d_complex SCAL_PROD = conj(AUX->comp[0][0]) * AUX->comp[1][0] + conj(AUX->comp[0][1]) * AUX->comp[1][1] + conj(AUX->comp[0][2]) * AUX->comp[1][2];

  for(int c=0;c<3;c++)
    AUX->comp[1][c] -= SCAL_PROD * AUX->comp[0][c];  

  //normalizzo la seconda riga
  NORM = creal(AUX->comp[1][0])*creal(AUX->comp[1][0])+cimag(AUX->comp[1][0])*cimag(AUX->comp[1][0]) + creal(AUX->comp[1][1])*creal(AUX->comp[1][1])+cimag(AUX->comp[1][1])*cimag(AUX->comp[1][1]) + creal(AUX->comp[1][2])*creal(AUX->comp[1][2])+cimag(AUX->comp[1][2])*cimag(AUX->comp[1][2]);
  NORM = 1.0/sqrt(NORM);

  AUX->comp[1][0] *= NORM;
  AUX->comp[1][1] *= NORM;
  AUX->comp[1][2] *= NORM;

  cnf->r0.c0[idx_cnf] = AUX->comp[0][0];
  cnf->r0.c1[idx_cnf] = AUX->comp[0][1];
  cnf->r0.c2[idx_cnf] = AUX->comp[0][2];
  cnf->r1.c0[idx_cnf] = AUX->comp[1][0];
  cnf->r1.c1[idx_cnf] = AUX->comp[1][1];
  cnf->r1.c2[idx_cnf] = AUX->comp[1][2];

  // temporaneo --> poi togliere!!
  //   cnf->r2.c0[idx_cnf]= conj(AUX->comp[0][1] * AUX->comp[1][2] - AUX->comp[0][2] * AUX->comp[1][1]);
  //   cnf->r2.c1[idx_cnf]= conj(AUX->comp[0][2] * AUX->comp[1][0] - AUX->comp[0][0] * AUX->comp[1][2]);
  //   cnf->r2.c2[idx_cnf]= conj(AUX->comp[0][0] * AUX->comp[1][1] - AUX->comp[0][1] * AUX->comp[1][0]);
  
}


#endif

