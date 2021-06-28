#ifndef PLAQUETTES_C_
#define PLAQUETTES_C_

#include "./geometry.h"
#include "./plaquettes.h"
#include "./su3_utilities.h"

// routine for the computation of the average of the plaquettes computed on the plane mu-nu
// 1) all the plaquettes on the plane mu-nu are computed and saved locally
// 2) finally the reduction of the traces is performed
double calc_loc_plaquettes_nnptrick(
        __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
        __restrict su3_soa * const loc_plaq, //la placchetta locale.
        dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
        const int mu, const int nu)
{
    double K_mu_nu; //MOD.
    
  int d0, d1, d2, d3;

    
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
          for(d0=0; d0 < nd0; d0++) {
              
              int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
              int parity; //parity
              int dir_muA,dir_nuB; //mu and nu directions.
              int dir_muC,dir_nuD;

	  idxh = snum_acc(d0,d1,d2,d3);  // the site on the  half-lattice.
	  parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
	  
	  dir_muA = 2*mu +  parity;
	  dir_muC = 2*mu + !parity;
	  idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
	    
	  dir_nuB = 2*nu + !parity;
	  dir_nuD = 2*nu +  parity;
	  idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
	  //       r+nu (C)  r+mu+nu
	  //          +<---+
	  // nu       |    ^
	  // ^    (D) V    | (B)
	  // |        +--->+
	  // |       r  (A)  r+mu
	  // +---> mu

      //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
              
	  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
	  
	  d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
	  tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
              
             /* printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);*/
           //MOD****************************************//
              
              //K_mu_nu computation;
              K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);

              
              tr_local_plaqs[parity].c[idxh]=K_mu_nu*tr_local_plaqs[parity].c[idxh];
          //*****************************************//

              /*printf("%f +i%f : (%d,%d,%d,%d) \n ",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I,d0,d1,d2,d3);*/
              /*printf("(%d,%d,%d,%d)\n",d0,d1,d2,d3);*/
              
	}  // d0
      }  // d1
    }  // d2
  }  // d3

  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
  int t;  // ONLY GOOD FOR 1D CUT
    

    
#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
    res_R_p += creal(tr_local_plaqs[0].c[t]); //even sites plaquettes
      
    res_R_p += creal(tr_local_plaqs[1].c[t]); //odd sites plaquettes
  }


  #pragma acc kernels present(tr_local_plaqs)
   
    printf("ecco1 %f(%d)  %d %d  \n",creal(tr_local_plaqs[1].c[snum_acc(31,6,6,6)]),snum_acc(31,6,6,6),mu,nu);
    #pragma acc kernels present(tr_local_plaqs)

    printf("ecco2 %f(%d)  %d %d \n",creal(tr_local_plaqs[0].c[snum_acc(31,6,6,6)]),snum_acc(31,6,6,6),mu,nu);
  
  
    
  return res_R_p;
}// closes routine
// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples

//This function doesn't have a mu, nu parameters. It just compute all staples.
void calc_loc_staples_nnptrick_all(  
        __restrict const su3_soa * const u,
        __restrict su3_soa * const loc_stap )
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int d0, d1, d2, d3, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) { //for for all lattice sites.
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {

     #pragma acc loop seq 
	  for(mu=0; mu<4; mu++){ //for directions.
      #pragma acc loop seq
	    for(iter=0; iter<3; iter++){ //for dell'iter. The iteration along all possible directions.
            

	      int nu;
	      if (mu==0) { nu = iter + 1; }
	      else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
	      else if (mu==2) { nu = iter + (iter >> 1); }
	      else if (mu==3) { nu = iter; }
          else {  printf("NU ERROR!\n");
	      }
            
            
	      const int idxh = snum_acc(d0,d1,d2,d3);  // r  //the site.
	      const int parity = (d0+d1+d2+d3) % 2;
#pragma acc cache (nnp_openacc[idxh:8])

	      const int dir_link = 2*mu + parity;
	      const int dir_mu_2R = 2*mu + !parity; //obvious these links are on the neighbour, whose parity is different.
	      const int dir_mu_2L = 2*mu + !parity;
	      const int idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])

	      const int dir_nu_1R = 2*nu + !parity;
	      const int dir_nu_3R = 2*nu +  parity;
	      const int dir_nu_1L = 2*nu +  parity;
	      const int dir_nu_3L = 2*nu + !parity;

	      const int idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu

	      //computation of the Right part of the staple

            
            //The computation is
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
            
            
            //Adding K_mu(x) to the staple.
            int K_mu;
            
  
            K_mu=u[dir_link].K.d[idxh];
         
            
        
            if(idxh==snum_acc(31,6,6,6) & parity==1 ){
                /*printf("%d %d\n",iter,nu);*/
                printf("%d %d\n",mu,nu);
                 printf("CONTROLLO STAPLE[idxh] K_mu:%f[%d] \n",K_mu,idxh);
            }
            
            
            loc_stap[dir_link].r0.c0[idxh] *= K_mu;
            loc_stap[dir_link].r0.c1[idxh] *= K_mu;
            loc_stap[dir_link].r0.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r1.c0[idxh] *= K_mu;
            loc_stap[dir_link].r1.c1[idxh] *= K_mu;
            loc_stap[dir_link].r1.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r2.c0[idxh] *= K_mu;
            loc_stap[dir_link].r2.c1[idxh] *= K_mu;
            loc_stap[dir_link].r2.c2[idxh] *= K_mu;
            
            
            

	    }  // mu
	  }  // iter

	}  // d0
      }  // d1
    }  // d2
  }  // d3
  #pragma acc kernels present(loc_stap)
    printf("ecco (31,6,6,6): %f\n",creal(loc_stap[1].r0.c0[snum_acc(31,6,6,6)]));
    
    
}// closes routine

#ifdef MULTIDEVICE
void calc_loc_staples_nnptrick_all_bulk(  
        __restrict const su3_soa * const u,
        __restrict su3_soa * const loc_stap )
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int d0, d1, d2, d3, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
#pragma acc loop independent gang(STAPGANG3) 
  for(d3=D3_HALO+GAUGE_HALO; d3<nd3-D3_HALO-GAUGE_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {

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

	      const int idxh = snum_acc(d0,d1,d2,d3);  // r 
	      const int parity = (d0+d1+d2+d3) % 2;
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
            
            //Adding K_mu(x) to the staple.
            int K_mu=u[dir_link].K.d[idxh];
            
            loc_stap[dir_link].r0.c0[idxh] *= K_mu;
            loc_stap[dir_link].r0.c1[idxh] *= K_mu;
            loc_stap[dir_link].r0.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r1.c0[idxh] *= K_mu;
            loc_stap[dir_link].r1.c1[idxh] *= K_mu;
            loc_stap[dir_link].r1.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r2.c0[idxh] *= K_mu;
            loc_stap[dir_link].r2.c1[idxh] *= K_mu;
            loc_stap[dir_link].r2.c2[idxh] *= K_mu;
            
            
            
          

	    }  // mu
	  }  // iter

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine

void calc_loc_staples_nnptrick_all_d3c(  
        __restrict const su3_soa * const u,
        __restrict su3_soa * const loc_stap,
       int offset, int thickness)
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int d0, d1, d2, d3, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
#pragma acc loop independent gang
  for(d3=offset; d3<offset+thickness; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {

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

	      const int idxh = snum_acc(d0,d1,d2,d3);  // r 
	      const int parity = (d0+d1+d2+d3) % 2;
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
            
            //Adding K_mu(x) to the staple.
            int K_mu=u[dir_link].K.d[idxh];
            
            loc_stap[dir_link].r0.c0[idxh] *= K_mu;
            loc_stap[dir_link].r0.c1[idxh] *= K_mu;
            loc_stap[dir_link].r0.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r1.c0[idxh] *= K_mu;
            loc_stap[dir_link].r1.c1[idxh] *= K_mu;
            loc_stap[dir_link].r1.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r2.c0[idxh] *= K_mu;
            loc_stap[dir_link].r2.c1[idxh] *= K_mu;
            loc_stap[dir_link].r2.c2[idxh] *= K_mu;

	    }  // mu
	  }  // iter

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine


#endif


void calc_loc_staples_nnptrick_all_only_even( 
        __restrict const su3_soa * const u,
        __restrict su3_soa * const loc_stap )
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int hd0, d1, d2, d3, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
#pragma acc loop independent gang(STAPGANG3) 
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(hd0=0; hd0 < nd0h; hd0++) {

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
	      const int d0 = 2*hd0 + ((d1+d2+d3) & 0x1); //           d0 = 2*hd0 + ((d1+d2+d3+1) & 0x1); (for the odd case)
	      const int idxh = snum_acc(d0,d1,d2,d3);  // r 
	      const int parity = (d0+d1+d2+d3) % 2;
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

            //Adding K_mu(x) to the staple.
            int K_mu=u[dir_link].K.d[idxh];
            
            loc_stap[dir_link].r0.c0[idxh] *= K_mu;
            loc_stap[dir_link].r0.c1[idxh] *= K_mu;
            loc_stap[dir_link].r0.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r1.c0[idxh] *= K_mu;
            loc_stap[dir_link].r1.c1[idxh] *= K_mu;
            loc_stap[dir_link].r1.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r2.c0[idxh] *= K_mu;
            loc_stap[dir_link].r2.c1[idxh] *= K_mu;
            loc_stap[dir_link].r2.c2[idxh] *= K_mu;
            
            
	    }  // mu
	  }  // iter

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine
void calc_loc_staples_nnptrick_all_only_odd(
        __restrict const su3_soa * const u,
        __restrict su3_soa * const loc_stap )
{
  //       r+mu-nu  r+mu   r+mu+nu
  //          +<-----+----->+
  //          |  1L  ^  1R  |
  // mu    2L |      |      | 2R
  // ^        V  3L  |  3R  V
  // |        +----->+<-----+
  // |       r-nu    r     r+nu
  // +---> nu       
  //            r is idxh in the following      


  int hd0, d1, d2, d3, mu, iter;

#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
#pragma acc loop independent gang(STAPGANG3) 
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
     for(d1=0; d1<nd1; d1++) {
	for(hd0=0; hd0 < nd0h; hd0++) {

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
	      const int d0 = 2*hd0 + ((d1+d2+d3+1) & 0x1);
	      const int idxh = snum_acc(d0,d1,d2,d3);  // r 
	      const int parity = (d0+d1+d2+d3) % 2;
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
            
            //Adding K_mu(x) to the staple.
            int K_mu=u[dir_link].K.d[idxh];
            
            loc_stap[dir_link].r0.c0[idxh] *= K_mu;
            loc_stap[dir_link].r0.c1[idxh] *= K_mu;
            loc_stap[dir_link].r0.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r1.c0[idxh] *= K_mu;
            loc_stap[dir_link].r1.c1[idxh] *= K_mu;
            loc_stap[dir_link].r1.c2[idxh] *= K_mu;
            
            loc_stap[dir_link].r2.c0[idxh] *= K_mu;
            loc_stap[dir_link].r2.c1[idxh] *= K_mu;
            loc_stap[dir_link].r2.c2[idxh] *= K_mu;

	    }  // mu
	  }  // iter

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine





#endif
