#ifndef SU3_UTILITIES_V2_C_
#define SU3_UTILITIES_V2_C_

//#pragma acc routine seq
static inline void read_mat_from_su3_soa_tworows( __restrict su3_soa * const mat,
						  const int idx_mat,
						  d_complex Omat[3][3]){
  Omat[0][0] = mat->r0.c0[idx_mat];
  Omat[0][1] = mat->r0.c1[idx_mat];
  Omat[0][2] = mat->r0.c2[idx_mat];
  Omat[1][0] = mat->r1.c0[idx_mat];
  Omat[1][1] = mat->r1.c1[idx_mat];
  Omat[1][2] = mat->r1.c2[idx_mat];
}
//#pragma acc routine seq
static inline void write_su3_soa_tworows_from_mat( __restrict su3_soa * const mat,
						   const int idx_mat,
						   d_complex Imat[3][3]){
  mat->r0.c0[idx_mat] = Imat[0][0];
  mat->r0.c1[idx_mat] = Imat[0][1];
  mat->r0.c2[idx_mat] = Imat[0][2];
  mat->r1.c0[idx_mat] = Imat[1][0];
  mat->r1.c1[idx_mat] = Imat[1][1];
  mat->r1.c2[idx_mat] = Imat[1][2];
}

//#pragma acc routine seq
static inline void write_tamat_soa_taking_ta_from_mat( __restrict tamat_soa * const mat,
						       const int idx_mat,
						       d_complex Imat[3][3]){


  mat->c01[idx_mat]  = 0.5*(Imat[0][1]-conj(Imat[1][0]));
  mat->c02[idx_mat]  = 0.5*(Imat[0][2]-conj(Imat[2][0]));
  mat->c12[idx_mat]  = 0.5*(Imat[1][2]-conj(Imat[2][1]));
  mat->rc00[idx_mat] = cimag(Imat[0][0])-ONE_BY_THREE*(cimag(Imat[0][0])+cimag(Imat[1][1])+cimag(Imat[2][2]));
  mat->rc11[idx_mat] = cimag(Imat[1][1])-ONE_BY_THREE*(cimag(Imat[0][0])+cimag(Imat[1][1])+cimag(Imat[2][2]));
}

//#pragma acc routine seq
static inline void write_tamat_soa_taking_ta_from_mat_times_dfact( __restrict tamat_soa * const mat,
								   const int idx_mat,
								   d_complex Imat[3][3],
								   double factor){


  mat->c01[idx_mat]  = factor*0.5*(Imat[0][1]-conj(Imat[1][0]));
  mat->c02[idx_mat]  = factor*0.5*(Imat[0][2]-conj(Imat[2][0]));
  mat->c12[idx_mat]  = factor*0.5*(Imat[1][2]-conj(Imat[2][1]));
  mat->rc00[idx_mat] = factor*(cimag(Imat[0][0])-ONE_BY_THREE*(cimag(Imat[0][0])+cimag(Imat[1][1])+cimag(Imat[2][2])));
  mat->rc11[idx_mat] = factor*(cimag(Imat[1][1])-ONE_BY_THREE*(cimag(Imat[0][0])+cimag(Imat[1][1])+cimag(Imat[2][2])));
}

//#pragma acc routine seq
static inline void write_su3_soa_reconstr_from_mat( __restrict su3_soa * const mat,
						    const int idx_mat,
						    d_complex Imat[3][3]){
  mat->r0.c0[idx_mat] = Imat[0][0];
  mat->r0.c1[idx_mat] = Imat[0][1];
  mat->r0.c2[idx_mat] = Imat[0][2];
  mat->r1.c0[idx_mat] = Imat[1][0];
  mat->r1.c1[idx_mat] = Imat[1][1];
  mat->r1.c2[idx_mat] = Imat[1][2];
  mat->r2.c0[idx_mat] = conj(Imat[0][1]*Imat[1][2]-Imat[0][2]*Imat[1][1]);
  mat->r2.c1[idx_mat] = conj(Imat[0][2]*Imat[1][0]-Imat[0][0]*Imat[1][2]);
  mat->r2.c2[idx_mat] = conj(Imat[0][0]*Imat[1][1]-Imat[0][1]*Imat[1][0]);
}

//#pragma acc routine seq
static inline void write_su3_soa_reconstr_from_mat_times_dfact( __restrict su3_soa * const mat,
							       const int idx_mat,
							       d_complex Imat[3][3],
							       double factor){
  mat->r0.c0[idx_mat] = factor * Imat[0][0];
  mat->r0.c1[idx_mat] = factor * Imat[0][1];
  mat->r0.c2[idx_mat] = factor * Imat[0][2];
  mat->r1.c0[idx_mat] = factor * Imat[1][0];
  mat->r1.c1[idx_mat] = factor * Imat[1][1];
  mat->r1.c2[idx_mat] = factor * Imat[1][2];
  mat->r2.c0[idx_mat] = factor * conj(Imat[0][1]*Imat[1][2]-Imat[0][2]*Imat[1][1]);
  mat->r2.c1[idx_mat] = factor * conj(Imat[0][2]*Imat[1][0]-Imat[0][0]*Imat[1][2]);
  mat->r2.c2[idx_mat] = factor * conj(Imat[0][0]*Imat[1][1]-Imat[0][1]*Imat[1][0]);
}

//#pragma acc routine seq
static inline void read_mat_from_su3_soa_reconstr( __restrict su3_soa * const mat,
						   const int idx_mat,
						   d_complex Omat[3][3]){
  Omat[0][0] = mat->r0.c0[idx_mat];
  Omat[0][1] = mat->r0.c1[idx_mat];
  Omat[0][2] = mat->r0.c2[idx_mat];
  Omat[1][0] = mat->r1.c0[idx_mat];
  Omat[1][1] = mat->r1.c1[idx_mat];
  Omat[1][2] = mat->r1.c2[idx_mat];
  Omat[2][0] = conj( ( Omat[0][1] * Omat[1][2] ) - ( Omat[0][2] * Omat[1][1]) ) ;
  Omat[2][1] = conj( ( Omat[0][2] * Omat[1][0] ) - ( Omat[0][0] * Omat[1][2]) ) ;
  Omat[2][2] = conj( ( Omat[0][0] * Omat[1][1] ) - ( Omat[0][1] * Omat[1][0]) ) ;
}

//#pragma acc routine seq
static inline void read_mat_from_su3_soa_totally( __restrict su3_soa * const mat,
						  const int idx_mat,
						  d_complex Omat[3][3]){
  Omat[0][0] = mat->r0.c0[idx_mat];
  Omat[0][1] = mat->r0.c1[idx_mat];
  Omat[0][2] = mat->r0.c2[idx_mat];
  Omat[1][0] = mat->r1.c0[idx_mat];
  Omat[1][1] = mat->r1.c1[idx_mat];
  Omat[1][2] = mat->r1.c2[idx_mat];
  Omat[2][0] = mat->r2.c0[idx_mat];
  Omat[2][1] = mat->r2.c1[idx_mat];
  Omat[2][2] = mat->r2.c2[idx_mat];
}
//#pragma acc routine seq
static inline void read_conjmat_from_su3_soa_reconstr( __restrict su3_soa * const mat,
						       const int idx_mat,
						       d_complex Omat[3][3]){
  Omat[0][0] = conj(mat->r0.c0[idx_mat]);
  Omat[1][0] = conj(mat->r0.c1[idx_mat]);
  Omat[2][0] = conj(mat->r0.c2[idx_mat]);
  Omat[0][1] = conj(mat->r1.c0[idx_mat]);
  Omat[1][1] = conj(mat->r1.c1[idx_mat]);
  Omat[2][1] = conj(mat->r1.c2[idx_mat]);
  Omat[0][2] = conj( ( Omat[1][0] * Omat[2][1] ) - ( Omat[2][0] * Omat[1][1]) ) ;
  Omat[1][2] = conj( ( Omat[2][0] * Omat[0][1] ) - ( Omat[0][0] * Omat[2][1]) ) ;
  Omat[2][2] = conj( ( Omat[0][0] * Omat[1][1] ) - ( Omat[1][0] * Omat[0][1]) ) ;
}

//#pragma acc routine seq
static inline void M3_eq_M1_times_M2_2rows( d_complex M1[3][3],
					    d_complex M2[3][3],
					    d_complex M3[3][3]){
  for(int i=0;i<2;i++)
    for(int j=0;j<3;j++)
      M3[i][j] = M1[i][0] * M2[0][j] +  M1[i][1] * M2[1][j] +  M1[i][2] * M2[2][j];
}
//#pragma acc routine seq
static inline void M3_eq_M1_times_M2_allrows( d_complex M1[3][3],
					      d_complex M2[3][3],
					      d_complex M3[3][3]){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      M3[i][j] = M1[i][0] * M2[0][j] +  M1[i][1] * M2[1][j] +  M1[i][2] * M2[2][j];
}

// mat3 = mat1 * mat2 
#pragma acc routine seq
static inline void    mat1_times_mat2_into_mat3_absent_stag_phases_V2( __restrict su3_soa * const mat1,
								       const int idx_mat1,
								       __restrict su3_soa * const mat2,
								       const int idx_mat2,
								       __restrict su3_soa * const mat3,
								       const int idx_mat3) {
  d_complex MAT1[3][3];
  d_complex MAT2[3][3];
  d_complex MAT3[3][3];
  read_mat_from_su3_soa_tworows(mat1,idx_mat1,&MAT1);
  read_mat_from_su3_soa_reconstr(mat2,idx_mat2,MAT2);
  M3_eq_M1_times_M2_2rows(MAT1,MAT2,MAT3);
  write_su3_soa_tworows_from_mat(mat3,idx_mat3,MAT3);
}

// mat1 = mat1 * mat2 
#pragma acc routine seq
static inline void    mat1_times_mat2_into_mat1_absent_stag_phases_V2( __restrict su3_soa * const mat1,
								       const int idx_mat1,
								       __restrict su3_soa * const mat2,
								       const int idx_mat2) {
  d_complex MAT1[3][3];
  d_complex MAT2[3][3];
  d_complex MAT3[3][3];
  read_mat_from_su3_soa_tworows(mat1,idx_mat1,MAT1);
  read_mat_from_su3_soa_reconstr(mat2,idx_mat2,MAT2);
  M3_eq_M1_times_M2_2rows(MAT1,MAT2,MAT3);
  write_su3_soa_tworows_from_mat(mat1,idx_mat1,MAT3);
}

// mat1 = mat1 * hermitian_conjucate(mat2)
#pragma acc routine seq
static inline void    mat1_times_conj_mat2_into_mat1_absent_stag_phases_V2( __restrict su3_soa * const mat1,
									    const int idx_mat1,
									    __restrict su3_soa * const mat2,
									    const int idx_mat2){
  d_complex MAT1[3][3];
  d_complex MAT2[3][3];
  d_complex MAT3[3][3];
  read_mat_from_su3_soa_tworows(mat1,idx_mat1,MAT1);
  read_conjmat_from_su3_soa_reconstr(mat2,idx_mat2,MAT2);
  M3_eq_M1_times_M2_2rows(MAT1,MAT2,MAT3);
  write_su3_soa_tworows_from_mat(mat1,idx_mat1,MAT3);
}


// Routine for the computation of the 3 matrices which contributes to the right part of the staple
// mat4 = mat1 * hermitian_conjucate(mat2)* hermitian_conjucate(mat3)
#pragma acc routine seq
static inline void    mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases_V2(  __restrict su3_soa * const mat1,
											      const int idx_mat1,
											      __restrict su3_soa * const mat2,
											      const int idx_mat2,
											      __restrict su3_soa * const mat3,
											      const int idx_mat3,
											      __restrict su3_soa * const mat4,
											      const int idx_mat4){
  d_complex MAT1[3][3];
  d_complex MAT2[3][3];
  d_complex MAT3[3][3];
  read_mat_from_su3_soa_tworows(mat1,idx_mat1,MAT1);
  read_conjmat_from_su3_soa_reconstr(mat2,idx_mat2,MAT2);
  M3_eq_M1_times_M2_2rows(MAT1,MAT2,MAT3);
  read_conjmat_from_su3_soa_reconstr(mat3,idx_mat3,MAT1);
  M3_eq_M1_times_M2_2rows(MAT3,MAT1,MAT2);
  write_su3_soa_reconstr_from_mat_times_dfact(mat4,idx_mat4,MAT2,C_ZERO);
}

// Routine for the computation of the 3 matrices which contributes to the left part of the staple
// mat4 = hermitian_conjucate(mat1)* hermitian_conjucate(mat2) * mat3
#pragma acc routine seq
static inline void    conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases_V2(   __restrict su3_soa * const mat1,
											       const int idx_mat1,
											       __restrict su3_soa * const mat2,
											       const int idx_mat2,
											       __restrict su3_soa * const mat3,
											       const int idx_mat3,
											       __restrict su3_soa * const mat4,
											       const int idx_mat4){
  d_complex MAT1[3][3];
  d_complex MAT2[3][3];
  d_complex MAT3[3][3];
  read_conjmat_from_su3_soa_reconstr(mat1,idx_mat1,MAT1);
  read_conjmat_from_su3_soa_reconstr(mat2,idx_mat2,MAT2);
  M3_eq_M1_times_M2_2rows(MAT1,MAT2,MAT3);
  read_mat_from_su3_soa_reconstr(mat3,idx_mat3,MAT1);
  M3_eq_M1_times_M2_2rows(MAT3,MAT1,MAT2);
  write_su3_soa_reconstr_from_mat_times_dfact(mat4,idx_mat4,MAT2,C_ZERO);
}

#pragma acc routine seq
static inline void mat1_times_mat2_into_tamat3_V2(__restrict su3_soa * const mat1,
						  const int idx_mat1,
						  __restrict su3_soa * const mat2,
						  const int idx_mat2,
						  __restrict tamat_soa * const mat3,
						  const int idx_mat3){
  d_complex MAT1[3][3];
  d_complex MAT2[3][3];
  d_complex MAT3[3][3];
  read_mat_from_su3_soa_reconstr(mat1,idx_mat1,MAT1);
  //Load all the rows of mat2 (that is a staple variable)
  read_mat_from_su3_soa_totally(mat2,idx_mat2,MAT2);
  M3_eq_M1_times_M2_allrows(MAT1,MAT2,MAT3);
  write_tamat_soa_taking_ta_from_mat(mat3,idx_mat3,MAT3);
}

#pragma acc routine seq
static inline void RHO_times_mat1_times_mat2_into_tamat3_V2(__restrict su3_soa * const mat1,
							    const int idx_mat1,
							    __restrict su3_soa * const mat2,
							    const int idx_mat2,
							    __restrict tamat_soa * const mat3,
							    const int idx_mat3){
  d_complex MAT1[3][3];
  d_complex MAT2[3][3];
  d_complex MAT3[3][3];
  read_mat_from_su3_soa_reconstr(mat1,idx_mat1,MAT1);
  //Load all the rows of mat2 (that is a staple variable)
  read_mat_from_su3_soa_totally(mat2,idx_mat2,MAT2);
  M3_eq_M1_times_M2_allrows(MAT1,MAT2,MAT3);
  // oltre a moltiplicare per RHO devo anche dividere per C_ZERO
  // perche' le staples che entrano qui dentro sono le staples * C_ZERO --> lo devo togliere!!
  double tmp = RHO/C_ZERO;
  write_tamat_soa_taking_ta_from_mat_times_dfact(mat3,idx_mat3,MAT3,tmp);
}




// routine for the computation of the average of the plaquettes computed on the plane mu-nu
// 1) all the plaquettes on the plane mu-nu are computed and saved locally
// 2) finally the reduction of the traces is performed
double calc_loc_plaquettes_removing_stag_phases_nnptrick_V2(   __restrict su3_soa * const u,
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

	  mat1_times_mat2_into_mat3_absent_stag_phases_V2(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases_V2(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases_V2(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
	  
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


// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
void calc_loc_staples_removing_stag_phases_nnptrick_V2(  __restrict su3_soa * const u,
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
	  mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases_V2(&u[dir_nu_1R],idx_pmu,&u[dir_mu_2R],idx_pnu,&u[dir_nu_3R],idxh,&loc_stap[dir_link],idxh);
	  //computation of the Left  part of the staple
	  conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases_V2(&u[dir_nu_1L],idx_pmu_mnu,&u[dir_mu_2L],idx_mnu,&u[dir_nu_3L],idx_mnu,&loc_stap[dir_link],idxh);

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine



void calc_loc_staples_removing_stag_phases_nnptrick_all_V2(  __restrict su3_soa * const u,
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
	      mat1_times_conj_mat2_times_conj_mat3_addto_mat4_absent_stag_phases_V2(&u[dir_nu_1R],       idx_pmu,
										    &u[dir_mu_2R],       idx_pnu,
										    &u[dir_nu_3R],       idxh,
										    &loc_stap[dir_link], idxh);

	      const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
	      const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu

	      //computation of the Left  part of the staple
	      conj_mat1_times_conj_mat2_times_mat3_addto_mat4_absent_stag_phases_V2(&u[dir_nu_1L],       idx_pmu_mnu,
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
void conf_times_staples_ta_part_V2(__restrict su3_soa * const u,        // constant --> is not updated
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
	    mat1_times_mat2_into_tamat3_V2(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);
	    
	  }

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine

// tamattamat
void RHO_times_conf_times_staples_ta_part_V2(__restrict su3_soa * const u,        // constant --> is not updated
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
	    RHO_times_mat1_times_mat2_into_tamat3_V2(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine




double  calc_plaquette_soloopenacc_V2( __restrict  su3_soa * const tconf_acc, __restrict su3_soa * const local_plaqs, dcomplex_soa * const tr_local_plaqs){
  
  double tempo=0.0;
  // tolgo le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);
  // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)
  for(int mu=0;mu<3;mu++){
    for(int nu=mu+1;nu<4;nu++){
      // sommo i 6 risultati in tempo
      tempo  += calc_loc_plaquettes_removing_stag_phases_nnptrick_V2(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
    }
  }
  // rimetto le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);

  return tempo;
  
}



#endif





