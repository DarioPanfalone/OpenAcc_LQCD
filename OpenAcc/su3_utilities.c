#ifndef SU3_UTILITIES_C_
#define SU3_UTILITIES_C_

// ROUTINE TO CHOOSE THE ACTIVE GPU
//   should be moved elsewhere!
void select_working_gpu_homemade(int dev_index){
  acc_device_t my_device_type;
  my_device_type = acc_device_nvidia;
  int num_devices = acc_get_num_devices(my_device_type);
  printf("Number of OpenAcc exploitable GPUs found: %d \n",num_devices);

  // Pick the device number dev_index 
  acc_set_device_num(dev_index,my_device_type);
  printf("Selected GPU number: %d \n",dev_index);
}

// mat3 = mat1 * mat2 
static inline void    mat1_times_mat2_into_mat3_absent_stag_phases( const __restrict su3_soa * const mat1,
								  const int idx_mat1,
								  const __restrict su3_soa * const mat2,
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

// mat1 = mat1 * hermitian_conjucate(mat2)
static inline void    mat1_times_conj_mat2_into_mat1_absent_stag_phases( __restrict su3_soa * const mat1,
								       const int idx_mat1,
								       const __restrict su3_soa * const mat2,
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
static inline void    mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases( const __restrict su3_soa * const matnu1,
											  const int idx_mat_nu1,
											  const __restrict su3_soa * const matmu2,
											  const int idx_mat_mu2,
											  const __restrict su3_soa * const matnu3,
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
  mat4->r0.c0[idx_mat4] += mat1_00;
  mat4->r0.c1[idx_mat4] += mat1_01;
  mat4->r0.c2[idx_mat4] += mat1_02;

  mat4->r1.c0[idx_mat4] += mat1_10;
  mat4->r1.c1[idx_mat4] += mat1_11;
  mat4->r1.c2[idx_mat4] += mat1_12;

  mat4->r2.c0[idx_mat4] += conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  mat4->r2.c1[idx_mat4] += conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  mat4->r2.c2[idx_mat4] += conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
}

// Routine for the computation of the 3 matrices which contributes to the left part of the staple
// mat4 = hermitian_conjucate(mat1)* hermitian_conjucate(mat2) * mat3
static inline void    conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases( const __restrict su3_soa * const matnu1,
											  const int idx_mat_nu1,
											  const __restrict su3_soa * const matmu2,
											  const int idx_mat_mu2,
											  const __restrict su3_soa * const matnu3,
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
  mat4->r0.c0[idx_mat4] += mat1_00;
  mat4->r0.c1[idx_mat4] += mat1_01;
  mat4->r0.c2[idx_mat4] += mat1_02;

  mat4->r1.c0[idx_mat4] += mat1_10;
  mat4->r1.c1[idx_mat4] += mat1_11;
  mat4->r1.c2[idx_mat4] += mat1_12;

  mat4->r2.c0[idx_mat4] += conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  mat4->r2.c1[idx_mat4] += conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  mat4->r2.c2[idx_mat4] += conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
}


static inline void mat1_times_mat2_into_tamat3(const __restrict su3_soa * const mat1,
					       const int idx_mat1,
					       const __restrict su3_soa * const mat2,
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



// mat1 = mat1 * integer factor
static inline void   mat1_times_int_factor( __restrict su3_soa * const mat1,
					   const int idx_mat1,
					   int factor){

  mat1->r0.c0[idx_mat1] *= factor;
  mat1->r0.c1[idx_mat1] *= factor;
  mat1->r0.c2[idx_mat1] *= factor;

  mat1->r1.c0[idx_mat1] *= factor;
  mat1->r1.c1[idx_mat1] *= factor;
  mat1->r1.c2[idx_mat1] *= factor;

  //Third row is not multiplied
  /*
  mat1->r2.c0[idx_mat1] *= factor;
  mat1->r2.c1[idx_mat1] *= factor;
  mat1->r2.c2[idx_mat1] *= factor;
  */
}

// calcola la traccia della matrice di su3
static inline d_complex matrix_trace_absent_stag_phase(const __restrict su3_soa * const loc_plaq,
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
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
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
							    double * const plaqs,
							    const int mu,
							    const int nu){

  int x, y, z, t;
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
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

  plaqs[0] = res_R_p;
  plaqs[1] = res_I_p;

  return res_R_p;
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

void set_su3_soa_to_zero( __restrict su3_soa * const matrix){
  int hx, y, z, t;
#pragma acc kernels present(matrix)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
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

// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
void calc_loc_staples_removing_stag_phases_nnptrick(   const __restrict su3_soa * const u,
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
#pragma acc kernels present(u) present(loc_stap)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
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
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine



// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
void calc_loc_staples_removing_stag_phases_nnptrick_all(   const __restrict su3_soa * const u,
							   __restrict su3_soa * const loc_stap
							   ){
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
  int mu,nu,iter;
#pragma acc kernels present(u) present(loc_stap)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int looping_directions[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
	  int idxh,idx_pmu,idx_pnu,idx_pmu_mnu,idx_mnu;
	  int parity;
	  int dir_link;
	  int dir_nu_1R,dir_nu_1L;
	  int dir_mu_2R,dir_mu_2L;
	  int dir_nu_3R,dir_nu_3L;

	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  /*

	  ////////////////////////////////// DIR X  //////////////////////////////////////////////////
	  mu=0;
	  iter=0;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);


	  iter=1;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);


	  iter=2;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  ////////////////////////////////// DIR Y  //////////////////////////////////////////////////
	  mu=1;
	  iter=0;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  iter=1;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  iter=2;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  ////////////////////////////////// DIR Z  //////////////////////////////////////////////////
	  mu=2;
	  iter=0;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  iter=1;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  iter=2;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  ////////////////////////////////// DIR T  //////////////////////////////////////////////////
	  mu=3;
	  iter=0;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  iter=1;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	  iter=2;
	  nu = looping_directions[mu][iter];
	  dir_link = 2*mu + parity;	  dir_nu_1R = 2*nu + !parity;	  dir_mu_2R = 2*mu + !parity;	  dir_nu_3R = 2*nu +  parity;
	  dir_nu_1L = 2*nu +  parity;	  dir_mu_2L = 2*mu + !parity;	  dir_nu_3L = 2*nu + !parity;
	  idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
	  idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
	  idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu 
	  idx_mnu= nnm_openacc[idxh][nu][parity] ;          // r-nu
	  //computation of the Right part of the staple
	  mat1_times_conj_mat2__times_conj_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_into_mat4_absent_stag_phases(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);


	  */
	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine





// tamattamat
void conf_times_staples_ta_part(const __restrict su3_soa * const u,
				const __restrict su3_soa * const loc_stap,
				__restrict tamat_soa * const ipdot){

  int x, y, z, t;
#pragma acc kernels present(u) present(loc_stap) present(ipdot)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&ipdot[dir_link],idxh);

	  }

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine

static inline void tamat1_plus_tamat2_times_factor_into_tamat1(__restrict tamat_soa * const tam1,
							       const __restrict tamat_soa * const tam2,
							       int idx,
							       const double fact){

  tam1->c01[idx]  += fact * tam2->c01[idx]; // complex
  tam1->c02[idx]  += fact * tam2->c02[idx]; // complex
  tam1->c12[idx]  += fact * tam2->c12[idx]; // complex
  tam1->rc00[idx] += fact * tam2->rc00[idx];  // double
  tam1->rc11[idx] += fact * tam2->rc11[idx];  // double
  
}

void mom_sum_mult( __restrict tamat_soa * const mom,
		   const __restrict tamat_soa * const ipdot,
		   const double factor){
  // !!!!!!!!!!!!!!!  factor  must be equal to   -beta/3.0*timestep !!!!!!!!!!!!!!!!!!!!
  int x, y, z, t;
#pragma acc kernels present(mom) present(ipdot)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
            tamat1_plus_tamat2_times_factor_into_tamat1(&mom[dir_link],&ipdot[dir_link],idxh,factor);
	  }
	}  // x
      }  // y
    }  // z
  }  // t
}// closes routine



static inline void exponetiate_acc(const __restrict thmat_soa * const mom,
				   const int idx_mom,
				   __restrict su3_soa * const aux,
				   const int idx_aux){

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


}

void kernel_acc_mom_exp_times_conf(su3_soa *conf,su3_soa *aux_conf,const thmat_soa * mom){

  int x, y, z, t;
#pragma acc kernels present(mom) present(conf) present(aux_conf)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(x=0; x < nx; x++) {
          int idxh;
          int parity;
          int dir_link;
          int mu;
          idxh = snum_acc(x,y,z,t);  // r
          parity = (x+y+z+t) % 2;
          for(mu=0;mu<4;mu++){
            dir_link = 2*mu + parity;
	    /////////////////////////////////////////////mettere l-istruzione giusta
          }

        }  // x
      }  // y
    }  // z
  }  // t

}


void mom_exp_times_conf_openacc(su3COM_soa *conf,const thmatCOM_soa * com_mom){
  su3_soa  * conf_acc,*aux_conf_acc;
  thmat_soa * momenta;
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));    // --> 4*size
  posix_memalign((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa));    // --> 4*size
  posix_memalign((void **)&ipdot, ALIGN, 8*sizeof(thmat_soa)); // --> 4*size

  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);
  for(dir=0;dir<8;dir++)  convert_thmatCOM_soa_to_thmat_soa(&com_mom[dir],&momenta[dir]);

  select_working_gpu_homemade(0);

  struct timeval t0, t1,t2,t3;

  gettimeofday ( &t0, NULL );
#pragma acc data copy(conf_acc[0:8]) copyin(momenta[0:8]) create(aux_conf_acc[0:8])
  {
    gettimeofday ( &t1, NULL );
    mult_conf_times_stag_phases(conf_acc);

    kernel_acc_mom_exp_times_conf(conf_acc,aux_conf_acc,momenta)

    mult_conf_times_stag_phases(conf_acc);
    gettimeofday ( &t2, NULL );
  }
    gettimeofday ( &t3, NULL );


  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("FULL STAPLES CALC OPENACC                       Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);

  free(conf_acc);
  free(ipdot);


}





void  calc_plaquette_openacc(const su3COM_soa *conf){
  su3_soa  * conf_acc, * local_plaqs;
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));    // --> 4*size
  posix_memalign((void **)&local_plaqs, ALIGN, 2*sizeof(su3_soa)); // --> size --> 1 plaquetta per sito (a fissato piano mu-nu)

  compute_nnp_and_nnm_openacc();

  dcomplex_soa * tr_local_plaqs;
  posix_memalign((void **)&tr_local_plaqs, ALIGN, 2*sizeof(dcomplex_soa)); // --> size complessi --> vettore per sommare tracce di plaquette locali

  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);


  double *plaq;
  posix_memalign((void **)&plaq, ALIGN, 2*sizeof(double));
  plaq[0] = 0.00;
  plaq[1] = 0.00;

  double tempo=0.0;
  select_working_gpu_homemade(0);
  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

#pragma acc data copy(conf_acc[0:8]) copyout(plaq[0:2]) create(local_plaqs[0:2])   create(tr_local_plaqs[0:2]) copyin(nnp_openacc)
  {
    gettimeofday ( &t1, NULL );
    // rimuovo le fasi staggered
    mult_conf_times_stag_phases(conf_acc);  
    // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)
    for(int mu=0;mu<3;mu++){
      for(int nu=mu+1;nu<4;nu++){
	// sommo i 6 risultati in tempo
	tempo  += calc_loc_plaquettes_removing_stag_phases_nnptrick(conf_acc,local_plaqs,tr_local_plaqs,plaq,mu,nu);
      }
    }
    // rimetto le fasi staggered
    mult_conf_times_stag_phases(conf_acc);
    gettimeofday ( &t2, NULL );
    
  }
  gettimeofday ( &t3, NULL );

  // stampo il valore della plaquette
  printf("Plaquette OPENACC = %.18lf \n",tempo/((double)(nx*ny*nz*nt*6.0*3.0)));

  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("FULL PLAQUETTE CALC OPENACC                     Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);


  free(conf_acc);
  free(local_plaqs);
  free(plaq);
  free(tr_local_plaqs);
  }


//QUESTA ROUTINE NON FUNZIONA ANCORA CORRETTAMENTE ==> non da i risultati corretti
// TUTTAVIA COMPILANDO CON O SENZA I PRGMA OPENACC I RISULTATI SONO UGUALI 
void  calc_staples_openacc(const su3COM_soa *conf,su3COM_soa *COM_staples){
  su3_soa  * conf_acc, * local_staples;
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));      // --> 4*size
  posix_memalign((void **)&local_staples, ALIGN, 8*sizeof(su3_soa)); // --> 4*size

  compute_nnp_and_nnm_openacc();

  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);

  double tempo=0.0;
  select_working_gpu_homemade(0);
  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

  int looping_directions[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};

#pragma acc data copy(conf_acc[0:8]) copyout(local_staples[0:8])
  {
    gettimeofday ( &t1, NULL );

    mult_conf_times_stag_phases(conf_acc);
    set_su3_soa_to_zero(local_staples);
    for(int mu=0;mu<4;mu++){
      for(int iter=0;iter<3;iter++){
	calc_loc_staples_removing_stag_phases_nnptrick(conf_acc,local_staples,mu,looping_directions[mu][iter]);

      }
    }

    //-> in this routine the mu and nu loops are inside the kernel
    //         ANYHOW STILL NOT WORKING
    //    calc_loc_staples_removing_stag_phases_nnptrick_all(conf_acc,local_staples); 

    mult_conf_times_stag_phases(conf_acc);
    gettimeofday ( &t2, NULL );
  }
  gettimeofday ( &t3, NULL );

  printf("Staple00 = ( %.18lf , %.18lf )\n", creal(local_staples[0].r0.c0[0]) , cimag(local_staples[0].r0.c0[0]));
  printf("Staple33 = ( %.18lf , %.18lf )\n", creal(local_staples[0].r2.c2[0]) , cimag(local_staples[0].r2.c2[0]));

  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("FULL STAPLES CALC OPENACC                       Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);



  for(dir=0;dir<8;dir++)  convert_su3_soa_to_su3COM_soa(&local_staples[dir],&COM_staples[dir]);

  free(conf_acc);
  free(local_staples);

  }




void  calc_ipdot_gauge_openacc(const su3COM_soa *conf,tamatCOM_soa * com_ipdot){
//void  calc_ipdot_gauge_openacc(const su3COM_soa *conf){
  su3_soa  * conf_acc, * local_staples;
  tamat_soa * ipdot;
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));      // --> 4*size
  posix_memalign((void **)&local_staples, ALIGN, 8*sizeof(su3_soa)); // --> 4*size
  posix_memalign((void **)&ipdot, ALIGN, 8*sizeof(tamat_soa)); // --> 4*size

  compute_nnp_and_nnm_openacc();

  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);

  double tempo=0.0;
  select_working_gpu_homemade(0);
  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

  int looping_directions[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};

#pragma acc data copyin(conf_acc[0:8]) create(local_staples[0:8]) copyout(ipdot[0:8])
  {
    gettimeofday ( &t1, NULL );

    mult_conf_times_stag_phases(conf_acc);
    set_su3_soa_to_zero(local_staples);
    for(int mu=0;mu<4;mu++){
      for(int iter=0;iter<3;iter++){
	calc_loc_staples_removing_stag_phases_nnptrick(conf_acc,local_staples,mu,looping_directions[mu][iter]);
	
      }
    }
    conf_times_staples_ta_part(conf_acc,local_staples,ipdot);
    
    mult_conf_times_stag_phases(conf_acc);
    gettimeofday ( &t2, NULL );
  }
  gettimeofday ( &t3, NULL );

  int index = 19;

  printf("Ipdot 00 = ( %.18lf )\n", ipdot[0].rc00[index]);
  printf("Ipdot 11 = ( %.18lf )\n", ipdot[0].rc11[index]);
  printf("Ipdot 01 = ( %.18lf , %.18lf )\n", creal(ipdot[0].c01[index]) , cimag(ipdot[0].c01[index]));
  printf("Ipdot 02 = ( %.18lf , %.18lf )\n", creal(ipdot[0].c02[index]) , cimag(ipdot[0].c02[index]));
  printf("Ipdot 12 = ( %.18lf , %.18lf )\n", creal(ipdot[0].c12[index]) , cimag(ipdot[0].c12[index]));

  // scrivere le routines di conversione convert_tamat_soa_to_tamatCOM_soa
    for(dir=0;dir<8;dir++)  convert_tamat_soa_to_tamatCOM_soa(&ipdot[dir],&com_ipdot[dir]);

    /*
  printf("Ipdot 00 = ( %.18lf )\n", com_ipdot[0].rc00[index]);
  printf("Ipdot 11 = ( %.18lf )\n", com_ipdot[0].rc11[index]);
  printf("Ipdot 01 = ( %.18lf , %.18lf )\n", com_ipdot[0].c01[index].Re , com_ipdot[0].c01[index].Im);
  printf("Ipdot 02 = ( %.18lf , %.18lf )\n", com_ipdot[0].c02[index].Re , com_ipdot[0].c02[index].Im);
  printf("Ipdot 12 = ( %.18lf , %.18lf )\n", com_ipdot[0].c12[index].Re , com_ipdot[0].c12[index].Im);
    */


  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("FULL STAPLES CALC OPENACC                       Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);




  free(conf_acc);
  free(local_staples);
  free(ipdot);
  }


#endif





