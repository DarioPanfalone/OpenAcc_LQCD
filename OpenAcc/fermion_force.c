#ifndef FERMION_FORCE_C
#define FERMION_FORCE_C

#include "./stouting.c"
#include "./struct_c_def.c"

#define TIMING_FERMION_FORCE


//STANDARD VERSION OF THE FERMIONIC FORCE
void fermion_force_soloopenacc(__restrict su3_soa    * tconf_acc, // la configurazione qui dentro e' costante e non viene modificata           
                   __restrict su3_soa * tstout_conf_acc_arr,
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

  stout_wrapper(tconf_acc,tstout_conf_acc_arr);

  set_tamat_soa_to_zero(tipdot_acc);
  for(int iflav = 0; iflav < tNDiffFlavs; iflav++) {
    set_su3_soa_to_zero(taux_conf_acc);
    int ifps = tfermion_parameters[iflav].index_of_the_first_ps;
    for(int ips = 0 ; ips < tfermion_parameters[iflav].number_of_ps ; ips++){
        // USING STOUTED GAUGE MATRIX
        multishift_invert(&(tstout_conf_acc_arr[8*(STOUT_STEPS-1)]), &tfermion_parameters[iflav], &(tfermion_parameters[iflav].approx_md),
                backfield, tferm_shiftmulti_acc, &(ferm_in_acc[ifps+ips]), res, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
        // USING STOUTED GAUGE MATRIX
        ker_openacc_compute_fermion_force(&(tstout_conf_acc_arr[8*(STOUT_STEPS-1)]), backfield, taux_conf_acc, tferm_shiftmulti_acc, tkloc_s, tkloc_h, &(tfermion_parameters[iflav]));
    }
    // JUST MULTIPLY BY BACK FIELD...
    multiply_conf_times_force_and_take_ta_even(tconf_acc,&(tfermion_parameters[iflav]),backfield, taux_conf_acc,tipdot_acc);// WRONG
    multiply_conf_times_force_and_take_ta_odd(tconf_acc,&(tfermion_parameters[iflav]),backfield, taux_conf_acc,tipdot_acc);// WRONG
  }

  // BACK CHAIN TO OBTAIN THE FORCE AT STOUTING LEVEL 0
  // ....
  //
  //






  /*
  //////TEMPORANEISSIMO
  print_tamat_soa(tipdot_acc,"FERM_FORCE_BF"); 
  mult_conf_times_stag_phases(tconf_acc);
  calc_loc_staples_removing_stag_phases_nnptrick_all(tconf_acc,taux_conf_acc);
  RHO_times_conf_times_staples_ta_part(tconf_acc,taux_conf_acc,aux_ta);

  compute_lambda(tipdot_acc,aux_th,aux_ta,tconf_acc,auxbis_conf_acc);


  mult_conf_times_stag_phases(tconf_acc);
  print_tamat_soa(tipdot_acc,"FERM_FORCE_AF");
  print_thmat_soa(aux_th,"LAMBDA"); 
  ///////////
  */




#ifdef TIMING_FERMION_FORCE
    gettimeofday ( &t2, NULL );
    double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("FULL FERMION FORCE COMPUTATION                  PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

    //  printf("########################################### \n");
    //  printf("#### Completed fermion force openacc ###### \n");
    //  printf("########################################### \n");

}












/*

//STANDARD VERSION OF THE FERMIONIC FORCE
void fermion_force_soloopenacc_stout(__restrict su3_soa    * tconf_acc, // la configurazione qui dentro e' costante e non viene modificata           
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
  
#ifdef TIMING_FERMION_FORCE
  struct timeval t1,t2;
  gettimeofday ( &t1, NULL );
#endif


  set_tamat_soa_to_zero(tipdot_acc);
  for(int iflav = 0; iflav < tNDiffFlavs; iflav++) {
    set_su3_soa_to_zero(taux_conf_acc);
    int ifps = tfermion_parameters[iflav].index_of_the_first_ps;
    for(int ips = 0 ; ips < tfermion_parameters[iflav].number_of_ps ; ips++){
      multishift_invert(tconf_acc, &tfermion_parameters[iflav], &(tfermion_parameters[iflav].approx_md),
			backfield, tferm_shiftmulti_acc, &(ferm_in_acc[ifps+ips]), res, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
      ker_openacc_compute_fermion_force(tconf_acc, backfield, taux_conf_acc, tferm_shiftmulti_acc, tkloc_s, tkloc_h, &(tfermion_parameters[iflav]));
    }
    multiply_conf_times_force_and_take_ta_even(tconf_acc,&(tfermion_parameters[iflav]),backfield, taux_conf_acc,tipdot_acc);
    multiply_conf_times_force_and_take_ta_odd(tconf_acc,&(tfermion_parameters[iflav]),backfield, taux_conf_acc,tipdot_acc);
  }
  




#ifdef TIMING_FERMION_FORCE
    gettimeofday ( &t2, NULL );
    double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("FULL FERMION FORCE COMPUTATION                  PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

    //  printf("########################################### \n");
    //  printf("#### Completed fermion force openacc ###### \n");
    //  printf("########################################### \n");

}

*/



void fermion_force_soloopenacc_stout(  __restrict su3_soa  * ferm_force, // la var globale e' auxbis_conf_acc
				       __restrict thmat_soa  * Lambda, // la var globale e' aux_th
				       __restrict tamat_soa  * QA, // la var globale e' aux_ta
				       __restrict su3_soa   * const U, // la var globale e' .... per adesso conf_acc
				       __restrict su3_soa   * const TMP// la var globale e' aux_conf_acc
				       ){

  printf("INSIDE fermion_force_soloopenacc_stout \n");
  set_su3_soa_to_zero(TMP);
  mult_conf_times_stag_phases(U);
  printf("         Removed stag phases  \n");
  calc_loc_staples_removing_stag_phases_nnptrick_all(U,TMP);
  printf("         computed staples  \n");
  RHO_times_conf_times_staples_ta_part(U,TMP,QA);
  printf("         computed Q  \n");
  compute_lambda(Lambda,ferm_force,U,QA,TMP);
  printf("         computed Lambda  \n");
  compute_sigma(Lambda,U,ferm_force,QA,TMP);

  mult_conf_times_stag_phases(U);
  printf("         Restored stag phases  \n");


}






#endif
