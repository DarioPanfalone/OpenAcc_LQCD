#ifndef FERMION_FORCE_C
#define FERMION_FORCE_C

#include "./stouting.c"
#include "./struct_c_def.c"
#include "./inverter_multishift_full.c"

#define TIMING_FERMION_FORCE


//STANDARD VERSION OF THE FERMIONIC FORCE
void fermion_force_soloopenacc(__restrict su3_soa    * tconf_acc, // la configurazione qui dentro e' costante e non viene modificata           
#ifdef STOUT_FERMIONS        
                   __restrict su3_soa * tstout_conf_acc_arr,// parking
                   __restrict su3_soa * gl3_aux, // gl(3) parking
#endif
			       __restrict double_soa * backfield,
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

    //  printf("############################################ \n");
    //  printf("#### Inside fermion force soloopenacc ###### \n");
    //  printf("############################################ \n");

#ifdef TIMING_FERMION_FORCE
    struct timeval t1,t2;
    gettimeofday ( &t1, NULL );
#endif

    su3_soa * conf_to_use; // CONF TO USE IN CALCULATION OF 
    // FERMION FORCE
#ifdef STOUT_FERMIONS
    stout_wrapper(tconf_acc,tstout_conf_acc_arr);// calcolo 
    conf_to_use =  &(tstout_conf_acc_arr[8*(STOUT_STEPS-1)]);
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
                    &(tfermion_parameters[iflav].approx_md), backfield,
                    tferm_shiftmulti_acc, &(ferm_in_acc[ifps+ips]), res, 
                    tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
            ker_openacc_compute_fermion_force(conf_to_use, backfield, taux_conf_acc, tferm_shiftmulti_acc, tkloc_s, tkloc_h, &(tfermion_parameters[iflav]));
        }

#ifdef STOUT_FERMIONS
 #if defined(IMCHEMPOT) || defined(BACKFIELD)
        // JUST MULTIPLY BY BACK FIELD AND/OR CHEMICAL POTENTIAL
        multiply_backfield_times_force(&(tfermion_parameters[iflav],backfield,taux_conf_acc,gl3_aux);
 #else               
       accumulate_gl3soa_into_gl3soa(taux_conf_acc,gl3_aux); 
 #endif 
#else
       multiply_conf_times_force_and_take_ta_even(tconf_acc,&(tfermion_parameters[iflav]),backfield, taux_conf_acc,tipdot_acc);
       multiply_conf_times_force_and_take_ta_odd(tconf_acc,&(tfermion_parameters[iflav]),backfield, taux_conf_acc,tipdot_acc);
#endif
                }

#ifdef STOUT_FERMIONS
  for(int stout_level = STOUT_STEPS ; stout_level > 1 ; stout_level--){
       conf_to_use = &(tstout_conf_acc_arr[8*(stout_level-2)]);
       compute_sigma_from_sigma_prime_backinto_sigma_prime(gl3_aux, aux_th,aux_ta,conf_to_use, aux_conf_acc );
       }
       compute_sigma_from_sigma_prime_backinto_sigma_prime(gl3_aux, aux_th,aux_ta,tconf_acc, aux_conf_acc );
       multiply_conf_times_force_and_take_ta_even_nophase(tconf_acc,&(tfermion_parameters[iflav]), taux_conf_acc,tipdot_acc);
       multiply_conf_times_force_and_take_ta_odd_nophase(tconf_acc,&(tfermion_parameters[iflav]), taux_conf_acc,tipdot_acc);
#endif

#ifdef TIMING_FERMION_FORCE
    gettimeofday ( &t2, NULL );
    double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("FULL FERMION FORCE COMPUTATION                  PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

    //  printf("########################################### \n");
    //  printf("#### Completed fermion force openacc ###### \n");
    //  printf("########################################### \n");

}






void compute_sigma_from_sigma_prime_backinto_sigma_prime(  __restrict su3_soa    * Sigma, // la var globale e' auxbis_conf_acc [sia input che ouptput]
							   __restrict thmat_soa  * Lambda, // la var globale e' aux_th
							   __restrict tamat_soa  * QA, // la var globale e' aux_ta
							   __restrict su3_soa    * const U, // la var globale e' .... per adesso conf_acc
							   __restrict su3_soa    * const TMP// la var globale e' aux_conf_acc //PARCHEGGIO??
							   ){

  printf("INSIDE fermion_force_soloopenacc_stout \n");
  set_su3_soa_to_zero(TMP);
  mult_conf_times_stag_phases(U);
  printf("         Removed stag phases  \n");
  calc_loc_staples_removing_stag_phases_nnptrick_all(U,TMP);
  printf("         computed staples  \n");
  RHO_times_conf_times_staples_ta_part(U,TMP,QA);
  printf("         computed Q  \n");
  compute_lambda(Lambda,Sigma,U,QA,TMP);
  printf("         computed Lambda  \n");
  compute_sigma(Lambda,U,Sigma,QA,TMP);
  
  mult_conf_times_stag_phases(U);
  printf("         Restored stag phases  \n");
  
  
}






#endif
