#ifndef FERMION_FORCE_C
#define FERMION_FORCE_C

#include "../Include/common_defines.h"
#include "./stouting.h"
#include "./struct_c_def.h"
#include "./inverter_multishift_full.h"
#include "./alloc_vars.h"
#include "./backfield.h"
#include "./fermion_force.h"
#include "./fermion_force_utilities.h"
#include "../Include/fermion_parameters.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "./su3_utilities.h"
#include "./plaquettes.h"
#include "./action.h"

#ifndef __GNUC__
 #define TIMING_FERMION_FORCE
#endif

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif





void compute_sigma_from_sigma_prime_backinto_sigma_prime(  __restrict su3_soa    * Sigma, // la var globale e' auxbis_conf_acc [sia input che ouptput]
							   __restrict thmat_soa  * Lambda, // la var globale e' aux_th
							   __restrict tamat_soa  * QA, // la var globale e' aux_ta
							   __restrict su3_soa    * const U, // la var globale e' .... per adesso conf_acc
							   __restrict su3_soa    * const TMP// la var globale e' aux_conf_acc //PARCHEGGIO??
							   ){




  if(verbosity_lv > 2) printf("\t\tSIGMA_PRIME --> SIGMA\n");
  if(verbosity_lv > 5){// printing stuff
#pragma acc update host(Sigma[0:8])
  printf("-------------Sigma[old]------------------\n");                                                                                             
  printf("Sigma[old]00 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c0[0]),cimag(Sigma[0].r0.c0[0]));                                               
  printf("Sigma[old]01 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c1[0]),cimag(Sigma[0].r0.c1[0]));                                               
  printf("Sigma[old]02 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c2[0]),cimag(Sigma[0].r0.c2[0]));                                               
  printf("Sigma[old]10 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c0[0]),cimag(Sigma[0].r1.c0[0]));                                               
  printf("Sigma[old]11 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c1[0]),cimag(Sigma[0].r1.c1[0]));                                               
  printf("Sigma[old]12 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c2[0]),cimag(Sigma[0].r1.c2[0]));                                               
  printf("Sigma[old]20 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c0[0]),cimag(Sigma[0].r2.c0[0]));                                               
  printf("Sigma[old]21 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c1[0]),cimag(Sigma[0].r2.c1[0]));                                               
  printf("Sigma[old]22 = %.18lf + (%.18lf)*I\n\n",creal(Sigma[0].r2.c2[0]),cimag(Sigma[0].r2.c2[0]));                
 }

  SETREQUESTED(TMP);
  set_su3_soa_to_zero(TMP);

   if(verbosity_lv > 5)printf("         Removed stag phases  \n");
  calc_loc_staples_removing_stag_phases_nnptrick_all(U,TMP);
   if(verbosity_lv > 5)printf("         computed staples  \n");


  RHO_times_conf_times_staples_ta_part(U,TMP,QA);
  // check: TMP = local staples.
  SETFREE(TMP);
   if(verbosity_lv > 5) printf("         computed Q  \n");
 if(verbosity_lv > 5) {// printing stuff
#pragma acc update host(QA[0:8])
  printf("-------------Q------------------\n");
  printf("Q00 = %.18lf\n",QA[0].rc00[0]);
  printf("Q00 = %.18lf\n",QA[0].rc11[0]);
  printf("Q01 = %.18lf + (%.18lf)*I\n",creal(QA[0].c01[0]),cimag(QA[0].c01[0]));
  printf("Q02 = %.18lf + (%.18lf)*I\n",creal(QA[0].c02[0]),cimag(QA[0].c02[0]));
  printf("Q12 = %.18lf + (%.18lf)*I\n\n",creal(QA[0].c12[0]),cimag(QA[0].c12[0]));
  }
  SETREQUESTED(Lambda);
  compute_lambda(Lambda,Sigma,U,QA,TMP);
 if(verbosity_lv > 5)   printf("         computed Lambda  \n");

 if(verbosity_lv > 5) {// printing stuff
#pragma acc update host(Lambda[0:8])
  printf("-------------LAMBDA------------------\n");
  printf("Lambda00 = %.18lf\n",Lambda[0].rc00[0]);
  printf("Lambda00 = %.18lf\n",Lambda[0].rc11[0]);
  printf("Lambda01 = %.18lf + (%.18lf)*I\n",creal(Lambda[0].c01[0]),cimag(Lambda[0].c01[0]));
  printf("Lambda02 = %.18lf + (%.18lf)*I\n",creal(Lambda[0].c02[0]),cimag(Lambda[0].c02[0]));
  printf("Lambda12 = %.18lf + (%.18lf)*I\n\n",creal(Lambda[0].c12[0]),cimag(Lambda[0].c12[0]));
  }
  SETREQUESTED(Sigma);
  compute_sigma(Lambda,U,Sigma,QA,TMP);
 if(verbosity_lv > 5)   printf("         computed Sigma  \n");

 if(verbosity_lv > 5) {// printing stuff
#pragma acc update host(Sigma[0:8])
  printf("-------------Sigma[new]------------------\n");                                                                                             
  printf("Sigma[new]00 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c0[0]),cimag(Sigma[0].r0.c0[0]));                                               
  printf("Sigma[new]01 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c1[0]),cimag(Sigma[0].r0.c1[0]));                                               
  printf("Sigma[new]02 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r0.c2[0]),cimag(Sigma[0].r0.c2[0]));                                               
  printf("Sigma[new]10 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c0[0]),cimag(Sigma[0].r1.c0[0]));                                               
  printf("Sigma[new]11 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c1[0]),cimag(Sigma[0].r1.c1[0]));                                               
  printf("Sigma[new]12 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r1.c2[0]),cimag(Sigma[0].r1.c2[0]));                                               
  printf("Sigma[new]20 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c0[0]),cimag(Sigma[0].r2.c0[0]));                                               
  printf("Sigma[new]21 = %.18lf + (%.18lf)*I\n",creal(Sigma[0].r2.c1[0]),cimag(Sigma[0].r2.c1[0]));                                               
  printf("Sigma[new]22 = %.18lf + (%.18lf)*I\n\n",creal(Sigma[0].r2.c2[0]),cimag(Sigma[0].r2.c2[0]));                
  }

  SETFREE(Lambda);
 if(verbosity_lv > 5)   printf("         Restored stag phases  \n");
  
  
}







//STANDARD VERSION OF THE FERMIONIC FORCE
void fermion_force_soloopenacc(__restrict su3_soa    * tconf_acc, // la configurazione qui dentro e' costante e non viene modificata           
#ifdef STOUT_FERMIONS        
			       __restrict su3_soa * tstout_conf_acc_arr,// parking
			       __restrict su3_soa * gl3_aux, // gl(3) parking
#endif
			       __restrict tamat_soa  * tipdot_acc,
			       __restrict ferm_param * tfermion_parameters,// [nflavs] 
			       int tNDiffFlavs,
			       __restrict vec3_soa * ferm_in_acc, // [NPS_tot]         
			       double res,
			       __restrict su3_soa  * taux_conf_acc,
			       __restrict vec3_soa * tferm_shiftmulti_acc,//parking variable [max_ps*max_approx_order]           
			       __restrict vec3_soa * tkloc_r, // parking 
			       __restrict vec3_soa * tkloc_h, // parking 
			       __restrict vec3_soa * tkloc_s, // parking 
			       __restrict vec3_soa * tkloc_p, // parking 
			       __restrict vec3_soa * tk_p_shiftferm//parking variable [max_approx_order]           
			       ){

  if(verbosity_lv > 2){
    printf("\tInside fermion force soloopenacc\n");
  }

#ifdef TIMING_FERMION_FORCE
  struct timeval t1,t2;
  gettimeofday ( &t1, NULL );
#endif
  
  __restrict su3_soa * conf_to_use; // CONF TO USE IN CALCULATION OF 
  // FERMION FORCE
#ifdef STOUT_FERMIONS
  stout_wrapper(tconf_acc,tstout_conf_acc_arr);// calcolo 
  conf_to_use =  &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);
  set_su3_soa_to_zero(gl3_aux); // pseudo ipdot
#else
  conf_to_use = tconf_acc;
#endif
  
  set_tamat_soa_to_zero(tipdot_acc);

  for(int iflav = 0; iflav < tNDiffFlavs; iflav++) {
    set_su3_soa_to_zero(taux_conf_acc);
    int ifps = tfermion_parameters[iflav].index_of_the_first_ps;
    for(int ips = 0 ; ips < tfermion_parameters[iflav].number_of_ps ; ips++){
      multishift_invert(conf_to_use, &tfermion_parameters[iflav], 
			&(tfermion_parameters[iflav].approx_md), 
			tferm_shiftmulti_acc, &(ferm_in_acc[ifps+ips]), res, 
			tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);

      ker_openacc_compute_fermion_force(conf_to_use, taux_conf_acc, tferm_shiftmulti_acc, tkloc_s, tkloc_h, &(tfermion_parameters[iflav]));

    }
    
    // JUST MULTIPLY BY STAGGERED PHASES,
   // BACK FIELD AND/OR CHEMICAL POTENTIAL 
    multiply_backfield_times_force(&(tfermion_parameters[iflav]),taux_conf_acc,gl3_aux);

   
  }
#ifdef STOUT_FERMIONS

  for(int stout_level = act_params.stout_steps ; stout_level > 1 ; stout_level--){
    if(verbosity_lv > 2) printf("\t\tSigma' to Sigma [lvl %d to lvl %d]\n",stout_level,stout_level-1);
    conf_to_use = &(tstout_conf_acc_arr[8*(stout_level-2)]);
    compute_sigma_from_sigma_prime_backinto_sigma_prime(gl3_aux, aux_th,aux_ta,conf_to_use, taux_conf_acc );
  }
   if(verbosity_lv > 2)  printf("\t\tSigma' to Sigma [lvl 1 to lvl 0]\n");
  compute_sigma_from_sigma_prime_backinto_sigma_prime(gl3_aux, aux_th,aux_ta,tconf_acc, taux_conf_acc );

#endif

  multiply_conf_times_force_and_take_ta_nophase(tconf_acc, gl3_aux,tipdot_acc);



  /*
#pragma acc update host(tipdot_acc[0:8])
  printf("-------------FFORCE------------------\n");
       printf("F00 = %.18lf\n",tipdot_acc[0].rc00[0]);
       printf("F11 = %.18lf\n",tipdot_acc[0].rc11[0]);
       printf("F01 = %.18lf + (%.18lf)*I\n",creal(tipdot_acc[0].c01[0]),cimag(tipdot_acc[0].c01[0]));
       printf("F02 = %.18lf + (%.18lf)*I\n",creal(tipdot_acc[0].c02[0]),cimag(tipdot_acc[0].c02[0]));
       printf("F12 = %.18lf + (%.18lf)*I\n\n",creal(tipdot_acc[0].c12[0]),cimag(tipdot_acc[0].c12[0]));

  */

#ifdef TIMING_FERMION_FORCE
    gettimeofday ( &t2, NULL );
    double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("\t\tFULL FERMION FORCE COMPUTATION  PreKer->PostKer :%f sec  \n",dt_preker_to_postker);
#endif
 if(verbosity_lv > 2){
     printf("\t\tCompleted fermion force openacc\n");
 }
}


#endif
