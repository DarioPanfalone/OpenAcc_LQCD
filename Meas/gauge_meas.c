// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h
#ifndef GAUGE_MEAS_C
#define GAUGE_MEAS_C

#include "../OpenAcc/cooling.c"


#pragma acc routine seq
static inline void comp_U_U_Udag_Udag(){

}
#pragma acc routine seq
static inline void comp_and_add_U_Udag_Udag_U(){
}
#pragma acc routine seq
static inline void comp_and_add_Udag_Udag_U_U(){

}
#pragma acc routine seq
static inline void comp_and_add_Udag_U_U_Udag(){

}

#pragma acc routine seq
static inline void combine_fourleaves_to_get_loc_q(){

}





void compute_local_topological_charge(  __restrict su3_soa * const u,
					__restrict su3_soa * const quadri,
					double_soa * const loc_q){
  
  int x, y, z, t;
  int mu,iter;
#pragma acc kernels present(u) present(quadri) present(loc_q)
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

#pragma acc loop seq
          for(mu=0; mu<4; mu++){
#pragma acc loop seq
            for(iter=0; iter<3; iter++){
	      int nu;
              if (mu==0) { nu = iter + 1; }
              else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
              else if (mu==2) { nu = iter + (iter >> 1); }
              else if (mu==3) { nu = iter; }
              else { 
              }
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
					      &loc_q[parity],idxh);

	    }// closes iter (i.e. nu)
	  } // closes mu
	} // closes x
      } // closes y
    } // closes z
  } // closes t
}

#endif
