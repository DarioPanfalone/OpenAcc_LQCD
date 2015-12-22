#ifndef PLAQUETTES_C_
#define PLAQUETTES_C_


#include "./plaquettes.h"
#include "./su3_utilities.h"
#include "../DbgTools/debug_macros_glvarcheck.h"

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





#endif
