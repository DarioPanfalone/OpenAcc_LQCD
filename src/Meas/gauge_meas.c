// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h
#ifndef GAUGE_MEAS_C
#define GAUGE_MEAS_C

#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/su3_utilities.h"
#include "./gauge_meas.h"


#ifdef __GNUC__
 #include "math.h"
 #ifndef M_PI
  #define M_PI 3.14159265358979323846
 #endif
#endif



char gauge_outfilename[50];
char gauge_outfile_header[100];

void compute_local_topological_charge(  __restrict su3_soa * const u,
					__restrict su3_soa * const quadri,
					double_soa * const loc_q,
					int mu, int nu){
  
  int x, y, z, t;
#pragma acc kernels present(u) present(quadri) present(loc_q) present(nnp_openacc) present(nnm_openacc)
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
	  //	  for(nu=1; nu<4; nu++){
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
					    &loc_q[parity],      idxh);
	    //	      if((idxh==0)&&(parity==0)) printf("loc_q(mu=%d,nu=%d)= %.18lf\n",mu,nu,loc_q[0].d[0]);
	    
	    //	  } // closes nu
	} // closes x
      } // closes y
    } // closes z
  } // closes t
}


double reduce_loc_top_charge(double_soa * const loc_q){
  
  double result=0.0;
  int t;
#pragma acc kernels present(loc_q)
#pragma acc loop reduction(+:result)
  for(t=0; t<sizeh; t++) {
    result += loc_q[0].d[t];
    result += loc_q[1].d[t];
  }
  return result;
}


double compute_topological_charge(__restrict su3_soa * const u,
				  __restrict su3_soa * const quadri,
				  double_soa * const loc_q){  

  set_su3_soa_to_zero(quadri); // forse non serve a una mazza

  double temp_ch =0.0;
  compute_local_topological_charge(u,quadri,loc_q,0,1);// (x,y) - (z,t)
  temp_ch = reduce_loc_top_charge(loc_q);
  compute_local_topological_charge(u,quadri,loc_q,0,2);// (x,z) - (y,t)
  temp_ch += reduce_loc_top_charge(loc_q);
  compute_local_topological_charge(u,quadri,loc_q,0,3);// (x,t) - (y,z)
  temp_ch += reduce_loc_top_charge(loc_q);


  return  temp_ch;
}


#endif
