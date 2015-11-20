#ifndef SU3_UTILITIES_C_
#define SU3_UTILITIES_C_

#include <math.h>
#include "su3_utilities.h"


// multiply the whole configuration for the staggered phases field
// (only the first two lines)
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
// (all three lines)
void mult_gl3_soa_times_stag_phases( __restrict su3_soa * const u){
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
	  gl3_times_int_factor(&u[2], idxh, eta);
	  // dir  4  =  z even
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  gl3_times_int_factor(&u[4], idxh, eta);
	  // dir  6  =  t even
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  gl3_times_int_factor(&u[6], idxh, eta);

	  //odd sites
	  x = 2*hx + ((y+z+t+1) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  // dir  3  =  y odd
	  eta = 1 - ( 2*(x & 0x1) );
	  gl3_times_int_factor(&u[3], idxh, eta);
	  // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  gl3_times_int_factor(&u[5], idxh, eta);
	  // dir  7  =  t odd
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  gl3_times_int_factor(&u[7], idxh, eta);	  

	  
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


// reunitarize the conf by brute force
void unitarize_conf( __restrict su3_soa * const u){
  int idxh, dirindex;
#pragma acc kernels present(u)
#pragma acc loop independent
  for(dirindex = 0 ; dirindex < 8 ; dirindex++){
#pragma acc loop independent
    for( idxh = 0 ; idxh < sizeh; idxh++){
      loc_unitarize_conf(&u[dirindex],idxh);
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
	  
	  d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
	  tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
	  

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




void set_su3_soa_to_zero( __restrict su3_soa * const matrix){
  SETINUSE(matrix);
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
    SETINUSE(loc_stap);
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







void calc_loc_staples_removing_stag_phases_nnptrick_all_only_even(  __restrict su3_soa * const u,
								    __restrict su3_soa * const loc_stap ){
  SETINUSE(loc_stap);
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int hx, y, z, t, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
 #pragma acc loop independent gang 
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4) 
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32) 
	for(hx=0; hx < nxh; hx++) {

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
	      const int x = 2*hx + ((y+z+t) & 0x1); //           x = 2*hx + ((y+z+t+1) & 0x1); (for the odd case)
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





void calc_loc_staples_removing_stag_phases_nnptrick_all_only_odd(
        __restrict su3_soa * const u,__restrict su3_soa * const loc_stap ){
  SETINUSE(loc_stap);
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int hx, y, z, t, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
 #pragma acc loop independent gang 
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4) 
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32) 
	for(hx=0; hx < nxh; hx++) {

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
	      const int x = 2*hx + ((y+z+t+1) & 0x1);
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

    SETINUSE(tipdot);
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



void mom_exp_times_conf_soloopenacc(
        __restrict  su3_soa * const tconf_acc,
 thmat_soa * const tmomenta, // e' costante e qui dentro non viene modificata
 double * tdelta,  int id_delta){
  kernel_acc_mom_exp_times_conf(tconf_acc,tmomenta,tdelta, id_delta);
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





