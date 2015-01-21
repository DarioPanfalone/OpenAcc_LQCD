
#ifndef SU3_UTILITIES_C_
#define SU3_UTILITIES_C_

/*
#ifndef INCLUDE_STRUCT_C_DEF
#define INCLUDE_STRUCT_C_DEF
#include "./struct_c_def.c"
#include "openacc.h"
#endif
*/

static inline void   mat1_times_mat2_into_mat3_present_stag_phases( const __restrict su3_soa * const mat1,
								   const int idx_mat1,
								   const __restrict su3_soa * const mat2,
								   const int idx_mat2,
								   const int eta2,
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

  // Load 3rd mat2 row from global memory                      
  //  d_complex mat2_20 = mat2->r2.c0[idx_mat2];                  
  //  d_complex mat2_21 = mat2->r2.c1[idx_mat2];                  
  //  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  //Compute 3rd mat2 row from the first two                    
  d_complex mat2_20 = conj( ( mat2_01 * mat2_12 ) - ( mat2_02 * mat2_11) ) ;
  d_complex mat2_21 = conj( ( mat2_02 * mat2_10 ) - ( mat2_00 * mat2_12) ) ;
  d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_01 * mat2_10) ) ;
  //Multiply 3rd row by eta               
  mat2_20 = (mat2_20)*eta2;
  mat2_21 = (mat2_21)*eta2;
  mat2_22 = (mat2_22)*eta2;

  //Compute the first two rows of the solution
  mat3->r0.c0[idx_mat3] = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  mat3->r0.c1[idx_mat3] = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  mat3->r0.c2[idx_mat3] = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  mat3->r1.c0[idx_mat3] = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  mat3->r1.c1[idx_mat3] = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  mat3->r1.c2[idx_mat3] = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;
}

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

  // Load 3rd mat2 row from global memory                      
  //  d_complex mat2_20 = mat2->r2.c0[idx_mat2];                  
  //  d_complex mat2_21 = mat2->r2.c1[idx_mat2];                  
  //  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

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

static inline void    mat1_times_mat2_into_mat1_present_stag_phases( __restrict su3_soa * const mat1,
								   const int idx_mat1,
								   const __restrict su3_soa * const mat2,
								   const int idx_mat2,
								   const int eta2){
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

  // Load 3rd mat2 row from global memory                      
  //  d_complex mat2_20 = mat2->r2.c0[idx_mat2];                  
  //  d_complex mat2_21 = mat2->r2.c1[idx_mat2];                  
  //  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

  //Compute 3rd mat2 row from the first two                    
  d_complex mat2_20 = conj( ( mat2_01 * mat2_12 ) - ( mat2_02 * mat2_11) ) ;
  d_complex mat2_21 = conj( ( mat2_02 * mat2_10 ) - ( mat2_00 * mat2_12) ) ;
  d_complex mat2_22 = conj( ( mat2_00 * mat2_11 ) - ( mat2_01 * mat2_10) ) ;
  //Multiply 3rd row by eta               
  mat2_20 = (mat2_20)*eta2;
  mat2_21 = (mat2_21)*eta2;
  mat2_22 = (mat2_22)*eta2;

  //Compute the first two rows of the solution
  mat1->r0.c0[idx_mat1] = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  mat1->r0.c1[idx_mat1] = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  mat1->r0.c2[idx_mat1] = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  mat1->r1.c0[idx_mat1] = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  mat1->r1.c1[idx_mat1] = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  mat1->r1.c2[idx_mat1] = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;
}

static inline void    mat1_times_mat2_into_mat1_absent_stag_phases( __restrict su3_soa * const mat1,
								  const int idx_mat1,
								  const __restrict su3_soa * const mat2,
								  const int idx_mat2){
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

  // Load 3rd mat2 row from global memory                      
  //  d_complex mat2_20 = mat2->r2.c0[idx_mat2];                  
  //  d_complex mat2_21 = mat2->r2.c1[idx_mat2];                  
  //  d_complex mat2_22 = mat2->r2.c2[idx_mat2];

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


static inline void    mat1_times_conj_mat2_into_mat1_present_stag_phases( __restrict su3_soa * const mat1,
									const int idx_mat1,
									const __restrict su3_soa * const mat2,
									const int idx_mat2,
									const int eta2){
  
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
  //Multiply 3rd column by eta               
  mat2_02 = (mat2_02)*eta2;
  mat2_12 = (mat2_12)*eta2;
  mat2_22 = (mat2_22)*eta2;

  //Compute the first two rows of the solution
  mat1->r0.c0[idx_mat1] = mat1_00 * mat2_00 + mat1_01 * mat2_10 + mat1_02 * mat2_20 ;
  mat1->r0.c1[idx_mat1] = mat1_00 * mat2_01 + mat1_01 * mat2_11 + mat1_02 * mat2_21 ;
  mat1->r0.c2[idx_mat1] = mat1_00 * mat2_02 + mat1_01 * mat2_12 + mat1_02 * mat2_22 ;

  mat1->r1.c0[idx_mat1] = mat1_10 * mat2_00 + mat1_11 * mat2_10 + mat1_12 * mat2_20 ;
  mat1->r1.c1[idx_mat1] = mat1_10 * mat2_01 + mat1_11 * mat2_11 + mat1_12 * mat2_21 ;
  mat1->r1.c2[idx_mat1] = mat1_10 * mat2_02 + mat1_11 * mat2_12 + mat1_12 * mat2_22 ;
}

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


void multiply_conf_times_stag_phases(__restrict su3_soa * const u){
  int x, y, z, t;
  int idxh,eta;
#pragma acc kernels present(u) 
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(x=0; x < nx; x++) {

	  idxh = snum_acc(x,y,z,t);
	  // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  // dir  1  =  x odd    --> eta = 1 , no multiplication needed

	  // dir  2  =  y even
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[2], idxh, eta);
	  // dir  3  =  y odd
	  mat1_times_int_factor(&u[3], idxh, eta);

	  // dir  4  =  z even
	  eta = 1 - ( 2*((x+y) & 0x1) );
          mat1_times_int_factor(&u[4], idxh, eta);
	  // dir  5  =  z odd
          mat1_times_int_factor(&u[5], idxh, eta);

	  // dir  6  =  t even
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
          mat1_times_int_factor(&u[6], idxh, eta);
	  // dir  7  =  t odd
          mat1_times_int_factor(&u[7], idxh, eta);	  
	}
      }
    }
  }
}


static inline d_complex matrix_trace(const __restrict su3_soa * const loc_plaq,
		       const int idx){
  d_complex loc_plaq_00 = loc_plaq->r0.c0[idx];
  d_complex loc_plaq_01 = loc_plaq->r0.c1[idx];
  d_complex loc_plaq_10 = loc_plaq->r1.c0[idx];
  d_complex loc_plaq_11 = loc_plaq->r1.c1[idx];
  // 1) devo ricostruire per forza l'elemento 22: la terza riga non la ho calcolata mai, quindi dove sembrerebbe doverci essere in realta' non c'e niente
  // 2) carico quello che mi serve e ricostruisco
  // 3) il segno meno che c'e' e' per le fasi staggered 
  d_complex loc_plaq_22 = - conj( ( loc_plaq_00 * loc_plaq_11 ) - ( loc_plaq_01 * loc_plaq_10) ) ;
  return -(loc_plaq_00 + loc_plaq_11 + loc_plaq_22); // anche il meno che c'e' qui e' per le fasi staggered
}
static inline void compute_matrix_trace_and_add_to_reducible_vector(const __restrict su3_soa * const loc_plaq,
								    const int idx,
								    __restrict dcomplex_soa * const tr_local_plaqs
								    ){
  d_complex loc_plaq_00 = loc_plaq->r0.c0[idx];
  d_complex loc_plaq_01 = loc_plaq->r0.c1[idx];
  d_complex loc_plaq_10 = loc_plaq->r1.c0[idx];
  d_complex loc_plaq_11 = loc_plaq->r1.c1[idx];
  // 1) devo ricostruire per forza l'elemento 22: la terza riga non la ho calcolata mai, quindi dove sembrerebbe doverci essere in realta' non c'e niente
  // 2) carico quello che mi serve e ricostruisco
  // 3) il segno meno che c'e' e' per le fasi staggered 
  d_complex loc_plaq_22 = - conj( ( loc_plaq_00 * loc_plaq_11 ) - ( loc_plaq_01 * loc_plaq_10) ) ;
  tr_local_plaqs->c[idx] +=   -(loc_plaq_00 + loc_plaq_11 + loc_plaq_22); // anche il meno che c'e' qui e' per le fasi staggered
}


void calc_loc_plaquettes_present_stag_phases(   const __restrict su3_soa * const u,
						__restrict su3_soa * const loc_plaq,
						dcomplex_soa * const tr_local_plaqs,
						double * const plaqs){
  int x, y, z, t;

  //  bool V[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent  //#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent    //#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent      //#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent	//#pragma acc loop independent vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh,idxpmu,idxpnu;
	  int parity;
	  int mu,nu,dir_muA,dir_nuB,dir_muC,dir_nuD;
	  int coordinates[4];
	  int coordinates_neigh_nnp_mu[4];
	  int coordinates_neigh_nnp_nu[4];

	  int dimensions[4]={nx,ny,nz,nt};
	  int V[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
	  coordinates[0] = x;
	  coordinates[1] = y;
	  coordinates[2] = z;
	  coordinates[3] = t;
	  coordinates_neigh_nnp_mu[0] = x;  
	  coordinates_neigh_nnp_mu[1] = y;  
	  coordinates_neigh_nnp_mu[2] = z;  
	  coordinates_neigh_nnp_mu[3] = t;  

	  coordinates_neigh_nnp_nu[0] = x;  
	  coordinates_neigh_nnp_nu[1] = y;  
	  coordinates_neigh_nnp_nu[2] = z;  
	  coordinates_neigh_nnp_nu[3] = t;  

	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;

      	  tr_local_plaqs[parity].c[idxh] = 0.0 + I * 0.0;
	  
	  for(mu=0; mu<3; mu++){
	    dir_muA = 2*mu +  parity;
	    dir_muC = 2*mu + !parity;
	    
	    coordinates_neigh_nnp_mu[mu]  = coordinates[mu] + 1;
	    coordinates_neigh_nnp_mu[mu] *= (((coordinates_neigh_nnp_mu[mu]-dimensions[mu]) >> 31) & 0x1);
	    idxpmu = snum_acc(coordinates_neigh_nnp_mu[0],coordinates_neigh_nnp_mu[1],coordinates_neigh_nnp_mu[2],coordinates_neigh_nnp_mu[3]); // r+mu
	    // alla fine del loop, come ultima cosa mettere di nuovo coordinates_neigh_nnp_mu[mu] = coordinates[mu]; !!!
	    
	    for(nu=mu+1; nu<4; nu++){
	      dir_nuB = 2*nu + !parity;
	      dir_nuD = 2*nu +  parity;

	      coordinates_neigh_nnp_nu[nu]  = coordinates[nu] + 1;
	      coordinates_neigh_nnp_nu[nu] *= (((coordinates_neigh_nnp_nu[nu]-dimensions[nu]) >> 31) & 0x1);
	      idxpnu = snum_acc(coordinates_neigh_nnp_nu[0],coordinates_neigh_nnp_nu[1],coordinates_neigh_nnp_nu[2],coordinates_neigh_nnp_nu[3]); // r+nu
	      // alla fine del loop, come ultima cosa mettere di nuovo coordinates_neigh_nnp_nu[nu] = coordinates[nu]; !!! 
	      int  idx, eta_B, eta_C, eta_D;
	      //	      printf("idxh   = %i          x= %i   y= %i   z= %i   t= %i\n",idxh,coordinates[0],coordinates[1],coordinates[2],coordinates[3]);
	      //	      printf("idxh   = %i  mu= %i  x= %i   y= %i   z= %i   t= %i\n",idxh,mu,coordinates_neigh_nnp_mu[0],coordinates_neigh_nnp_mu[1],coordinates_neigh_nnp_mu[2],coordinates_neigh_nnp_mu[3]);
	      //	      printf("idxh   = %i  nu= %i  x= %i   y= %i   z= %i   t= %i\n",idxh,nu,coordinates_neigh_nnp_nu[0],coordinates_neigh_nnp_nu[1],coordinates_neigh_nnp_nu[2],coordinates_neigh_nnp_nu[3]);
	      

	      //       r+nu (C)  r+mu+nu
	      //          +<---+
	      // nu       |    ^
	      // ^    (D) V    | (B)
	      // |        +--->+
	      // |       r  (A)  r+mu
	      // +---> mu
	      
	      //	      //eta_A non e' necessaria per il calcolo perche' la terza riga di A non viene mai ricostruita
	      //	      //eta_A = V[mu][0] + (1-(2*(x & 0x1) ))*V[mu][1] + (1-(2*((x+y) & 0x1) ))*V[mu][2]+ (1-(2*((x+y+z) & 0x1) ))*V[mu][3]*(1- 2*(int)(t/(nt-1)));

	      eta_B = V[nu][0] + (1-(2*(coordinates_neigh_nnp_mu[0] & 0x1) ))*V[nu][1] + (1-(2*((coordinates_neigh_nnp_mu[0]+coordinates_neigh_nnp_mu[1]) & 0x1) ))*V[nu][2]+ (1-(2*((coordinates_neigh_nnp_mu[0]+coordinates_neigh_nnp_mu[1]+coordinates_neigh_nnp_mu[2]) & 0x1) ))*V[nu][3]*(1- 2*(int)(coordinates_neigh_nnp_mu[3]/(nt-1)));
	      
	      eta_C = V[mu][0] + (1-(2*(coordinates_neigh_nnp_nu[0] & 0x1) ))*V[mu][1] + (1-(2*((coordinates_neigh_nnp_nu[0]+coordinates_neigh_nnp_nu[1]) & 0x1) ))*V[mu][2]+ (1-(2*((coordinates_neigh_nnp_nu[0]+coordinates_neigh_nnp_nu[1]+coordinates_neigh_nnp_nu[2]) & 0x1) ))*V[mu][3]*(1- 2*(int)(coordinates_neigh_nnp_nu[3]/(nt-1)));
	      
	      eta_D = V[nu][0] + (1-(2*(x & 0x1) ))*V[nu][1] + (1-(2*((x+y) & 0x1) ))*V[nu][2]+ (1-(2*((x+y+z) & 0x1) ))*V[nu][3]*(1- 2*(int)(t/(nt-1)));
	      
	      mat1_times_mat2_into_mat3_present_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,eta_B,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
	      mat1_times_conj_mat2_into_mat1_present_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu,eta_C);              // LOC_PLAQ = LOC_PLAQ * C
	      mat1_times_conj_mat2_into_mat1_present_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh,eta_D);                // LOC_PLAQ = LOC_PLAQ * D
	      
	      //	      tr_local_plaqs[parity].c[idxh] += matrix_trace(&loc_plaq[parity],idxh);

	      compute_matrix_trace_and_add_to_reducible_vector(&loc_plaq[parity],idxh,&tr_local_plaqs[parity]);
	      
	      coordinates_neigh_nnp_nu[nu] = coordinates[nu];	      
	    }  // nu
	    
	    coordinates_neigh_nnp_mu[mu] = coordinates[mu];
	  }  // mu
	  tr_local_plaqs[parity].c[idxh] = 3.0 - I * 3.0;


	}  // x
      }  // y
    }  // z
  }  // t

  printf("LOCAL00 = %f \n",creal(tr_local_plaqs[0].c[0]));


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

}// closes routine




void calc_loc_plaquettes_removing_stag_phases(   __restrict su3_soa * const u,
						 __restrict su3_soa * const loc_plaq,
						 dcomplex_soa * const tr_local_plaqs,
						 double * const plaqs){
  int x, y, z, t;
  int mu,nu;
  int idxh,idxpmu,idxpnu;
  int parity;
  int dir_muA,dir_nuB;
  int dir_muC,dir_nuD;
  int coordinates[4];
  int coordinates_neigh_nnp_mu[4];
  int coordinates_neigh_nnp_nu[4];
  int dimensions[4];
  dimensions[0]=nx;
  dimensions[1]=ny;
  dimensions[2]=nz;
  dimensions[3]=nt;

  // toglie le fasi staggered (ed anche le condizioni al contorno antiperiodiche)
  multiply_conf_times_stag_phases(u);


#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent
  for(t=0; t<sizeh; t++) {
    tr_local_plaqs[0].c[idxh] = 0.0 + I * 0.0;
    tr_local_plaqs[1].c[idxh] = 0.0 + I * 0.0;
  }



  for(mu=0; mu<3; mu++){
    for(nu=mu+1; nu<4; nu++){
      
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  coordinates[0] = x;
	  coordinates[1] = y;
	  coordinates[2] = z;
	  coordinates[3] = t;
	  coordinates_neigh_nnp_mu[0] = x;
          coordinates_neigh_nnp_mu[1] = y;
          coordinates_neigh_nnp_mu[2] = z;
          coordinates_neigh_nnp_mu[3] = t;

          coordinates_neigh_nnp_nu[0] = x;
          coordinates_neigh_nnp_nu[1] = y;
          coordinates_neigh_nnp_nu[2] = z;
          coordinates_neigh_nnp_nu[3] = t;

	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;

	  //	  tr_local_plaqs[parity].c[idxh] = 0.0 + I * 0.0;

	  //	  for(mu=0; mu<3; mu++){
	    dir_muA = 2*mu +  parity;
	    dir_muC = 2*mu + !parity;
	    
	    coordinates_neigh_nnp_mu[mu]  = coordinates[mu] + 1;
	    coordinates_neigh_nnp_mu[mu] *= (((coordinates_neigh_nnp_mu[mu]-dimensions[mu]) >> 31) & 0x1);
	    idxpmu = snum_acc(coordinates_neigh_nnp_mu[0],coordinates_neigh_nnp_mu[1],coordinates_neigh_nnp_mu[2],coordinates_neigh_nnp_mu[3]); // r+mu
	    // alla fine del loop, come ultima cosa mettere di nuovo coordinates_neigh_nnp_mu[mu] = coordinates[mu]; !!!
	    
	    //	    for(nu=mu+1; nu<4; nu++){
	      dir_nuB = 2*nu + !parity;
	      dir_nuD = 2*nu +  parity;

	      coordinates_neigh_nnp_nu[nu]  = coordinates[nu] + 1;
	      coordinates_neigh_nnp_nu[nu] *= (((coordinates_neigh_nnp_nu[nu]-dimensions[nu]) >> 31) & 0x1);
	      idxpnu = snum_acc(coordinates_neigh_nnp_nu[0],coordinates_neigh_nnp_nu[1],coordinates_neigh_nnp_nu[2],coordinates_neigh_nnp_nu[3]); // r+nu
	      // alla fine del loop, come ultima cosa mettere di nuovo coordinates_neigh_nnp_nu[nu] = coordinates[nu]; !!! 
	      int  idx;
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
	      
	      tr_local_plaqs[parity].c[idxh] += matrix_trace(&loc_plaq[parity],idxh);

	      coordinates_neigh_nnp_nu[nu] = coordinates[nu];	      
	      //	    }  // nu
	    
	    coordinates_neigh_nnp_mu[mu] = coordinates[mu];
	    //	  }  // mu


	}  // x
      }  // y
    }  // z
  }  // t

	  }  // nu
	  }  // mu

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

  // rimette le fasi staggered
  multiply_conf_times_stag_phases(u);
  plaqs[1] = 23.4;
}// closes routine




void  calc_plaquette_openacc(const su3COM_soa *conf){
  su3_soa  * conf_acc, * local_plaqs;
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));    // --> 4*size
  posix_memalign((void **)&local_plaqs, ALIGN, 2*sizeof(su3_soa)); // --> size --> 1 plaquetta per sito (a fissato piano mu-nu)

  dcomplex_soa * tr_local_plaqs;
  posix_memalign((void **)&tr_local_plaqs, ALIGN, 2*sizeof(dcomplex_soa)); // --> size complessi --> vettore per sommare tracce di plaquette locali

  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);


  double *plaq;
  posix_memalign((void **)&plaq, ALIGN, 2*sizeof(double));
  plaq[0] = 0.00;
  plaq[1] = 0.00;

  //#pragma acc data copyin(conf_acc[0:8]) copy(plaq[0:2]) create(local_plaqs[0:2]) create(tr_local_plaqs[0:2])
#pragma acc data copyin(conf_acc[0:8]) copy(plaq[0:2]) copy(local_plaqs[0:2]) copy(tr_local_plaqs[0:2])
  {
    calc_loc_plaquettes_present_stag_phases(conf_acc,local_plaqs,tr_local_plaqs,plaq);
    //    calc_loc_plaquettes_removing_stag_phases(conf_acc,local_plaqs,tr_local_plaqs,plaq);

  }

  printf("Plaq RE = %.18lf \n",plaq[0]/((double)(nx*ny*nz*nt*6.0*3.0)));
  printf("Plaq IM = %.18lf \n",plaq[1]/((double)(nx*ny*nz*nt*6.0*3.0)));


  free(conf_acc);
  free(local_plaqs);
  free(plaq);
  free(tr_local_plaqs);
}


#endif


