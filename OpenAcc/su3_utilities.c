#ifndef SU3_UTILITIES_C_
#define SU3_UTILITIES_C_

// ROUTINE TO CHOOSE AND INITIALIZE THE OPENACC DEVICE
void SELECT_INIT_ACC_DEVICE(acc_device_t my_device_type, int dev_index) {

  // Initialize context for this device type
  acc_init(my_device_type);

  // Get available devices of this type
  int num_devices = acc_get_num_devices(my_device_type);
  printf("Number of OpenAcc exploitable devices found: %d \n", num_devices);

  // Pick the device number dev_index 
  acc_set_device_num(dev_index, my_device_type);
  printf("Selected device number: %d \n", dev_index);

}

void SHUTDOWN_ACC_DEVICE(acc_device_t my_device_type) {

  // Close context for this device type
  acc_shutdown(my_device_type);

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
  d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_10 * mat1_01) ) ;

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

  mat3->c01[idx_mat3]  = RHO*(0.5*(mat3_01-conj(mat1_00)));
  mat3->c02[idx_mat3]  = RHO*(0.5*(mat3_02-conj(mat1_10)));
  mat3->c12[idx_mat3]  = RHO*(0.5*(mat1_02-conj(mat1_11)));
  mat3->rc00[idx_mat3] = RHO*((cimag(mat3_00)-ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12))));
  mat3->rc11[idx_mat3] = RHO*((cimag(mat1_01)-ONE_BY_THREE*(cimag(mat3_00)+cimag(mat1_01)+cimag(mat1_12))));
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

// multiply the whole configuration for the staggered phases field
void mult_conf_times_stag_phases( __restrict su3_soa * const u){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u)
#pragma acc loop independent gang // gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector // gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector // //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector // vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

	  //even sites
	  x = 2*hx + ((y+z+t) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  // dir  2  =  y even
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[2], idxh, eta);
	  // dir  4  =  z even
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[4], idxh, eta);
	  // dir  6  =  t even
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[6], idxh, eta);

	  //odd sites
	  x = 2*hx + ((y+z+t+1) & 0x1);
	  idxh = snum_acc(x,y,z,t);	  
	  // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  // dir  3  =  y odd
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[3], idxh, eta);
	  // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[5], idxh, eta);
	  // dir  7  =  t odd
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[7], idxh, eta);	  

	  
	}
      }
    }
  }

}




// multiply the whole configuration for the staggered phases field
void mult_conf_times_stag_phases_nodev( __restrict su3_soa * const u){
  int hx,y,z,t,idxh;
  for(t=0; t<nt; t++) {
    for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

	  //even sites
	  x = 2*hx + ((y+z+t) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  // dir  2  =  y even
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[2], idxh, eta);
	  // dir  4  =  z even
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[4], idxh, eta);
	  // dir  6  =  t even
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[6], idxh, eta);

	  //odd sites
	  x = 2*hx + ((y+z+t+1) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  // dir  3  =  y odd
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[3], idxh, eta);
	  // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[5], idxh, eta);
	  // dir  7  =  t odd
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[7], idxh, eta);
	  
	}
      }
    }
  }

}









// routine for the computation of the average of the plaquettes computed on the plane mu-nu
// 1) all the plaquettes on the plane mu-nu are computed and saved locally
// 2) finally the reduction of the traces is performed
double calc_loc_plaquettes_removing_stag_phases_nnptrick(   __restrict su3_soa * const u,
							    __restrict su3_soa * const loc_plaq,
							    dcomplex_soa * const tr_local_plaqs,
							    const int mu,
							    const int nu){

  int x, y, z, t;
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh,idxpmu,idxpnu;
	  int parity;
	  int dir_muA,dir_nuB;
	  int dir_muC,dir_nuD;

	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;

	  dir_muA = 2*mu +  parity;
	  dir_muC = 2*mu + !parity;
	  idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
	    
	  dir_nuB = 2*nu + !parity;
	  dir_nuD = 2*nu +  parity;
	  idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
	  //       r+nu (C)  r+mu+nu
	  //          +<---+
	  // nu       |    ^
	  // ^    (D) V    | (B)
	  // |        +--->+
	  // |       r  (A)  r+mu
	  // +---> mu

	  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
	  
	  tr_local_plaqs[parity].c[idxh] = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
	  

	}  // x
      }  // y
    }  // z
  }  // t

  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=0; t<sizeh; t++) {
    res_R_p += creal(tr_local_plaqs[0].c[t]);
    res_R_p += creal(tr_local_plaqs[1].c[t]);
    res_I_p += cimag(tr_local_plaqs[0].c[t]);
    res_I_p += cimag(tr_local_plaqs[1].c[t]);
  }

  //  plaqs[0] = res_R_p;
  //  plaqs[1] = res_I_p;

  return res_R_p;
}// closes routine


static inline double half_tr_thmat_squared( const __restrict thmat_soa * const mom,
					    int idx_mom){
  d_complex  C = mom->c01[idx_mom];
  d_complex  D = mom->c02[idx_mom];
  d_complex  E = mom->c12[idx_mom];
  double    A = mom->rc00[idx_mom];
  double    B = mom->rc11[idx_mom];
  return A*A + B*B + A*B + creal(C)*creal(C) + cimag(C)*cimag(C) + creal(D)*creal(D) + cimag(D)*cimag(D) + creal(E)*creal(E) + cimag(E)*cimag(E);
  
}

double calc_momenta_action( const __restrict thmat_soa * const mom,
			    double_soa * tr_local,
			    const int mu){
  int t;

#pragma acc kernels present(mom) present(tr_local)
#pragma acc loop independent //gang(nt)
  for(t=0; t<sizeh; t++) {
	  tr_local[0].d[t] = half_tr_thmat_squared(&mom[mu],t);
  }  // t


  double result=0.0;
#pragma acc kernels present(tr_local)
#pragma acc loop reduction(+:result)
  for(t=0; t<sizeh; t++) {
    result += tr_local[0].d[t];
  }
  
  return result;


}// closes routine






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


void set_su3_soa_to_zero( __restrict su3_soa * const matrix){
  int hx, y, z, t;
#pragma acc kernels present(matrix)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(hx=0; hx < nxh; hx++) {
	  int x,idxh;
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);
	  assign_zero_to_su3_soa_component(&matrix[0],idxh);
	  assign_zero_to_su3_soa_component(&matrix[1],idxh);
	  assign_zero_to_su3_soa_component(&matrix[2],idxh);
	  assign_zero_to_su3_soa_component(&matrix[3],idxh);
	  assign_zero_to_su3_soa_component(&matrix[4],idxh);
	  assign_zero_to_su3_soa_component(&matrix[5],idxh);
	  assign_zero_to_su3_soa_component(&matrix[6],idxh);
	  assign_zero_to_su3_soa_component(&matrix[7],idxh);
	}
      }
    }
  }
}

//copy matrix in into matrix out, this has to happen on the host
void set_su3_soa_to_su3_soa( __restrict su3_soa * const matrix_in,
			     __restrict su3_soa * const matrix_out){
  int hx, y, z, t;
  for(t=0; t<nt; t++) {
    for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
	for(hx=0; hx < nxh; hx++) {
	  int x,idxh;
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[0],&matrix_out[0],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[1],&matrix_out[1],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[2],&matrix_out[2],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[3],&matrix_out[3],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[4],&matrix_out[4],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[5],&matrix_out[5],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[6],&matrix_out[6],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[7],&matrix_out[7],idxh);
	}
      }
    }
  }
}

// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
void calc_loc_staples_removing_stag_phases_nnptrick(  __restrict su3_soa * const u,
						      __restrict su3_soa * const loc_stap,
						      const int mu,
						      const int nu){
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      

  int x, y, z, t;
#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh,idx_pmu,idx_pnu,idx_pmu_mnu,idx_mnu;
	  int parity;
	  int dir_link;
	  int dir_nu_1R,dir_nu_1L;
	  int dir_mu_2R,dir_mu_2L;
	  int dir_nu_3R,dir_nu_3L;

	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;



	  dir_link = 2*mu + parity;
	  dir_nu_1R = 2*nu + !parity;
	  dir_mu_2R = 2*mu + !parity;
	  dir_nu_3R = 2*nu +  parity;

	  dir_nu_1L = 2*nu +  parity;
	  dir_mu_2L = 2*mu + !parity;
	  dir_nu_3L = 2*nu + !parity;

	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine



void calc_loc_staples_removing_stag_phases_nnptrick_all(  __restrict su3_soa * const u,
							  __restrict su3_soa * const loc_stap ){
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int x, y, z, t, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
 #pragma acc loop independent gang 
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4) 
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32) 
	for(x=0; x < nx; x++) {

     #pragma acc loop seq 
	  for(mu=0; mu<4; mu++){
      #pragma acc loop seq
	    for(iter=0; iter<3; iter++){

	      int nu;
	      if (mu==0) { nu = iter + 1; }
	      else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
	      else if (mu==2) { nu = iter + (iter >> 1); }
	      else if (mu==3) { nu = iter; }
	      else { //error 
	      }

	      const int idxh = snum_acc(x,y,z,t);  // r 
	      const int parity = (x+y+z+t) % 2;
#pragma acc cache (nnp_openacc[idxh:8])

	      const int dir_link = 2*mu + parity;
	      const int dir_mu_2R = 2*mu + !parity;
	      const int dir_mu_2L = 2*mu + !parity;
	      const int idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])

	      const int dir_nu_1R = 2*nu + !parity;
	      const int dir_nu_3R = 2*nu +  parity;
	      const int dir_nu_1L = 2*nu +  parity;
	      const int dir_nu_3L = 2*nu + !parity;

	      const int idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu

	      //computation of the Right part of the staple
	      mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1R],       idx_pmu,
										 &u[dir_mu_2R],       idx_pnu,
										 &u[dir_nu_3R],       idxh,
										 &loc_stap[dir_link], idxh);

	      const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
	      const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu

	      //computation of the Left  part of the staple
	      conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases(&u[dir_nu_1L],       idx_pmu_mnu,
										 &u[dir_mu_2L],       idx_mnu,
										 &u[dir_nu_3L],       idx_mnu,
										 &loc_stap[dir_link], idxh);

	    }  // mu
	  }  // iter

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine


// tamattamat
void conf_times_staples_ta_part(__restrict su3_soa * const u,        // constant --> is not updated
			        __restrict su3_soa * const loc_stap, // constant --> is not updated
				__restrict tamat_soa * const tipdot){

  int x, y, z, t;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine

// tamattamat
void RHO_times_conf_times_staples_ta_part(__restrict su3_soa * const u,        // constant --> is not updated
					  __restrict su3_soa * const loc_stap, // constant --> is not updated
					  __restrict tamat_soa * const tipdot){

  int x, y, z, t;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    RHO_times_mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine

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

void mom_sum_mult( __restrict thmat_soa * const mom,
		   const __restrict tamat_soa * const ipdot,
		   double * factor,
		   int id_factor){
  // !!!!!!!!!!!!!!!  factor  is  equal to   -beta/3.0*timestep !!!!!!!!!!!!!!!!!!!!
  int x, y, z, t;
#pragma acc kernels present(mom) present(ipdot) present(factor)
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
            thmat1_plus_tamat2_times_factor_into_thmat1(&mom[dir_link],&ipdot[dir_link],idxh,factor[id_factor]);
	  }
	}  // x
      }  // y
    }  // z
  }  // t
}// closes routine


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



void kernel_acc_mom_exp_times_conf( __restrict su3_soa * const conf,
				    thmat_soa * const mom, // e' costante e qui dentro non viene modificato
				    double * factor, // questo e' il vettore delta dove sono contenuti tutti i dt richiesti nell'omelyan
				    int id_factor){

  int x, y, z, t;
#pragma acc kernels present(mom) present(conf) present(factor)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
        for(x=0; x < nx; x++) {
          int idxh;
          int parity;
          int dir_link;
          int mu;
	  single_su3 mom_aux[1];
	  single_su3 expo[1];
	  single_su3 aux[1];
          idxh = snum_acc(x,y,z,t);  // r
          parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){
	    dir_link = 2*mu + parity;

	    extract_mom(&mom[dir_link],idxh,factor[id_factor],&mom_aux[0]);
	    matrix_exp_openacc(&mom_aux[0],&aux[0],&expo[0]);
	    conf_left_exp_multiply(&conf[dir_link],idxh,&expo[0],&aux[0],&mom_aux[0]);
	    project_on_su3(&conf[dir_link],idxh,&mom_aux[0]);
	  }
	  
        }  // x
      }  // y
    }  // z
  }  // t

}




// Ora che abbiamo tolto i prodotti con le fasi staggered questo e' diventato un wrapper inutile... lo togliamo poi ... 
void mom_exp_times_conf_soloopenacc( __restrict  su3_soa * const tconf_acc,
				     thmat_soa * const tmomenta, // e' costante e qui dentro non viene modificata
				     double * tdelta,
				     int id_delta){
  //    mult_conf_times_stag_phases(tconf_acc);  // toglie le fasi staggered dalla conf
  kernel_acc_mom_exp_times_conf(tconf_acc,tmomenta,tdelta, id_delta);
  //    mult_conf_times_stag_phases(tconf_acc); // rimette le fasi staggered nella conf
}


double  calc_plaquette_soloopenacc( __restrict  su3_soa * const tconf_acc, __restrict su3_soa * const local_plaqs, dcomplex_soa * const tr_local_plaqs){

  double tempo=0.0;
  // tolgo le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);
  // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)
  for(int mu=0;mu<3;mu++){
    for(int nu=mu+1;nu<4;nu++){
      // sommo i 6 risultati in tempo
      tempo  += calc_loc_plaquettes_removing_stag_phases_nnptrick(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
    }
  }
  // rimetto le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);

  return tempo;

  }



#endif





