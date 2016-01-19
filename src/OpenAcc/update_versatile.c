// se metro==1 allora fa il test di metropolis
// se metro==0 allora non fa il test di metropolis --> termalizzazione
#ifndef UPDATE_VERSATILE_C_
#define UPDATE_VERSATILE_C_

#include "../Include/common_defines.h"
#include "./update_versatile.h"
#include "./random_assignement.h"
#include "./find_min_max.h"
#include "./rettangoli.h"
#include "./alloc_vars.h"
#include "./stouting.h"
#include "./md_integrator.h"
#include "./su3_utilities.h"
#include "./su3_measurements.h"
#include "./inverter_multishift_full.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "./fermionic_utilities.h"
#include "./action.h"
#include "../Rand/random.h"


//#define NORANDOM  // FOR debug, check also main.c 
#ifdef NORANDOM
#include "../DbgTools/dbgtools.h"
#endif

#ifdef __GNUC__
#include "sys/time.h"
#endif

action_param act_params;

int UPDATE_SOLOACC_UNOSTEP_VERSATILE(su3_soa *tconf_acc,
        double res_metro, double res_md, int id_iter,int acc,int metro){
  
#ifdef STOUT_FERMIONS        
    su3_soa *tstout_conf_acc_arr = gstout_conf_acc_arr;
#endif

  
  printf("UPDATE_SOLOACC_UNOSTEP_VERSATILE_TLSM_STDFERM: OK \n");
  // DEFINIZIONE DI TUTTI I dt NECESSARI PER L'INTEGRATORE OMELYAN
  int iterazioni = id_iter+1;
  double dt_tot;
  double dt_pretrans_to_preker;
  double dt_preker_to_postker;
  double dt_postker_to_posttrans;

  int mu;
  double **minmaxeig; 
  minmaxeig = (double**)malloc(NDiffFlavs*sizeof(double*));
  for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++)
      minmaxeig[iflav] = (double*) malloc(2*sizeof(double));

  double p1;
  double p2;
  int accettata;
  double delta_S;
  double action_in,action_fin,action_mom_in,action_mom_fin,action_ferm_in,action_ferm_fin;

  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

  if(metro==1){
    // store old conf   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
    set_su3_soa_to_su3_soa(tconf_acc,conf_acc_bkp);
    if(verbosity_lv > 2) printf("Backup copy of the initial gauge conf : OK \n");

  }

//#pragma acc data copyin(delta[0:7]) // should not be needed?
  {
  
    gettimeofday ( &t1, NULL );

    // ESTRAZIONI RANDOM
#ifdef NORANDOM
    printf("NORANDOM mode, loading momenta from memory.\n");
    if(read_thmat_soa(momenta,"momenta_norndtest")){
        printf("GENERATING MOMENTA FILE FOR YOUR CONVENIENCE, RE-RUN THIS TEST\n");
        generate_Momenta_gauss(momenta);
        print_thmat_soa(momenta,"momenta_norndtest");
    }
#else
    generate_Momenta_gauss(momenta);
#endif

    printf("Momenta generated/read : OK \n");

    //    read_thmat_soa(momenta,"momenta");
#pragma acc update device(momenta[0:8])
    
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
      int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	 
#ifdef NORANDOM
      char psferm_filename[20];

      char ps_index_str[5];
      sprintf(ps_index_str,"%d",ps_index);
      strcpy(psferm_filename,"fermion_norndtest");
      strcat(psferm_filename,ps_index_str);
      if(read_vec3_soa(&ferm_phi_acc[ps_index],psferm_filename)){
          printf("GENERATING FERMION FILE FOR YOUR CONVENIENCE, RE-RUN THIS TEST\n");
	      generate_vec3_soa_gauss(&ferm_phi_acc[ps_index]);
          print_vec3_soa(&ferm_phi_acc[ps_index],psferm_filename);

      }
#else
      if(verbosity_lv > 3 )printf("Ferm generation (flav=%d,ps=%d,ps_index=%d) : OK \n",iflav,ips,ps_index);
	  generate_vec3_soa_gauss(&ferm_phi_acc[ps_index]);
#endif
      }
    }// end for iflav
#pragma acc update device(ferm_phi_acc[0:NPS_tot])

#ifdef STOUT_FERMIONS 
    // USO DELLA VERSIONE STOUTATA GIA' PER LO STIRACCHIAMENTO
    // STOUTING...(ALREADY ON DEVICE)
    stout_wrapper(tconf_acc,tstout_conf_acc_arr);
    gconf_as_fermionmatrix = &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);
#else
    gconf_as_fermionmatrix = tconf_acc;
#endif




    // STIRACCHIAMENTO DELL'APPROX RAZIONALE FIRST_INV
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
       if(verbosity_lv > 2 ) printf("Rat approx rescale (flav=%d)\n",iflav);

      SETREQUESTED(kloc_p);
      SETREQUESTED(kloc_s);
#ifdef NORANDOM
      if(read_vec3_soa(kloc_p,"kloc_p_norndtest")){
          generate_vec3_soa_gauss(kloc_p);
          print_vec3_soa(kloc_p,"kloc_p_norndtest");
          printf("GENERATED kloc_p_norndtest FOR NORANDOM TEST, RE-RUN THIS TEST\n");
      }
      if(read_vec3_soa(kloc_s,"kloc_s_norndtest")){
          generate_vec3_soa_gauss(kloc_s);
          print_vec3_soa(kloc_s,"kloc_s_norndtest");
          printf("GENERATED kloc_s_norndtest FOR NORANDOM TEST, RE-RUN THIS TEST\n");
      }
#else
      // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
#endif
      // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      // USING STOUTED GAUGE MATRIX
      //printf("    before min and max eig comp : OK \n");
      SETREQUESTED(kloc_r);
      SETREQUESTED(kloc_h);
      find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig[iflav]);
      if(verbosity_lv > 3 ) printf("    find min and max eig : OK \n");
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig[iflav]);
      if(verbosity_lv > 4 ) printf("    rat approx rescaled : OK \n\n");

#pragma acc update device(approx_fi[0:1])
    }//end for iflav
    
    if(metro==1){
      /////////////// INITIAL ACTION COMPUTATION ////////////////////////////////////////////
      if(GAUGE_ACTION == 0)// Standard gauge action 
	action_in = BETA_BY_THREE*calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      if(GAUGE_ACTION == 1){ //Tlsym gauge action
	action_in = C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
	action_in += C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      }
      action_mom_in = 0.0;
      for(mu =0;mu<8;mu++)  action_mom_in += calc_momenta_action(momenta,d_local_sums,mu);
      action_ferm_in=0;
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
	for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	  
	  int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	  action_ferm_in += real_scal_prod_global(&ferm_phi_acc[ps_index],&ferm_phi_acc[ps_index]);
	}
      }// end for iflav
      ///////////////////////////////////////////////////////////////////////////////////////
    }
    printf(" Initial Action Computed : OK \n");

    // FIRST INV APPROX CALC --> calcolo del fermione CHI

    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	
        int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
        // USING STOUTED GAUGE MATRIX
        multishift_invert(gconf_as_fermionmatrix, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
        recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), &(ferm_chi_acc[ps_index]),&(fermions_parameters[iflav].approx_fi));
	
      }
    }// end for iflav
    if(verbosity_lv > 3) printf(" Computed the fermion CHI : OK \n");
    
    // STIRACCHIAMENTO DELL'APPROX RAZIONALE MD
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      // recupero gli autovalori già calcolati che sono salvati in approx_fi
      RationalApprox *approx_md = &(fermions_parameters[iflav].approx_md);
      RationalApprox *approx_md_mother = &(fermions_parameters[iflav].approx_md_mother);
      rescale_rational_approximation(approx_md_mother,approx_md,minmaxeig[iflav]);
#pragma acc update device(approx_md[0:1])
    }//end for iflav
    
    // DINAMICA MOLECOLARE (stouting implicitamente usato in calcolo forza fermionica)
    multistep_2MN_SOLOOPENACC(ipdot_acc,tconf_acc,
#ifdef STOUT_FERMIONS
			      tstout_conf_acc_arr,
			      auxbis_conf_acc, // globale
#endif
			      u1_back_field_phases,aux_conf_acc,fermions_parameters,NDiffFlavs,
			      ferm_chi_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,
			      k_p_shiftferm,momenta,local_sums,res_md);
    
    if(verbosity_lv > 1) printf(" Molecular Dynamics Completed \n");


#ifdef STOUT_FERMIONS
    // STOUTING...(ALREADY ON DEVICE)
    stout_wrapper(tconf_acc,tstout_conf_acc_arr);
    gconf_as_fermionmatrix = &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);
#else
    gconf_as_fermionmatrix = tconf_acc;
#endif

    if(metro==1){
      // STIRACCHIAMENTO DELL'APPROX RAZIONALE LAST_INV
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
	// generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
	generate_vec3_soa_gauss(kloc_p);
	generate_vec3_soa_gauss(kloc_s);
	// update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
    // USING STOUTED CONF
	find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig[iflav]);
	//#pragma acc update device(minmaxeig[0:2])
	RationalApprox *approx_li = &(fermions_parameters[iflav].approx_li);
	RationalApprox *approx_li_mother = &(fermions_parameters[iflav].approx_li_mother);
	rescale_rational_approximation(approx_li_mother,approx_li,minmaxeig[iflav]);
#pragma acc update device(approx_li[0:1])
      }

      // LAST INV APPROX CALC 
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
	for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
        int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
        // USING STOUTED CONF
        multishift_invert(gconf_as_fermionmatrix, &fermions_parameters[iflav], 
                &(fermions_parameters[iflav].approx_li), u1_back_field_phases,
                ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), res_metro, 
                kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
        recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), &(ferm_phi_acc[ps_index]),&(fermions_parameters[iflav].approx_li));
	}
      }
      printf(" Final Action Computed : OK \n");
      
      ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
      if(GAUGE_ACTION == 0) // Standard gauge action
      action_fin = BETA_BY_THREE * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      if(GAUGE_ACTION == 1){ // Tlsym gauge action
      action_fin = C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      action_fin += C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      }

      action_mom_fin = 0.0;
      for(mu =0;mu<8;mu++)    action_mom_fin += calc_momenta_action(momenta,d_local_sums,mu);
      
      action_ferm_fin=0;
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
          for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
              int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
              action_ferm_fin += real_scal_prod_global(&ferm_chi_acc[ps_index],&ferm_phi_acc[ps_index]);
          }
      } // end for iflav


      ////////////////////////////////////////////////////////////////////////////////////////

      // delta_S = action_new - action_old
      delta_S  = - (-action_in+action_mom_in+action_ferm_in) + (-action_fin+action_mom_fin+action_ferm_fin);
      if(verbosity_lv > 2){
      printf("iterazione %i:  Gauge_ACTION  (in and out) = %.18lf , %.18lf\n",iterazioni,-action_in,-action_fin);
      printf("iterazione %i:  Momen_ACTION  (in and out) = %.18lf , %.18lf\n",iterazioni,action_mom_in,action_mom_fin);
      printf("iterazione %i:  Fermi_ACTION  (in and out) = %.18lf , %.18lf\n",iterazioni,action_ferm_in,action_ferm_fin);
      }
      printf("iterazione %i:  DELTA_ACTION = %.18lf, ",iterazioni,delta_S);
     


#ifdef NORANDOM      
      printf("Always accept in NORANDOM MODE!!!\n");
#endif

      if(delta_S<0){
	accettata=1;
      }
      else
	{
	  p1=exp(-delta_S);
#ifdef NORANDOM      
      p2=0;
#else
	  p2=casuale();
#endif
	  if(p2<p1)
	    {
	      accettata=1;
	    }
	  else
	    {
	      accettata=0;
	      // configuration reject
	    }
	}
    }
    gettimeofday ( &t2, NULL );
  }// end pragma acc data copyin(delta[0:7])
  
gettimeofday ( &t3, NULL );
  
  if(metro==1){
    if(accettata==1){
      acc++;
      printf("ACCEPTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
      // configuration accepted   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
      set_su3_soa_to_su3_soa(tconf_acc,conf_acc_bkp);
    }else{
      printf("REJECTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
      // configuration rejected   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
      set_su3_soa_to_su3_soa(conf_acc_bkp,tconf_acc);
#pragma acc update device(tconf_acc[0:8])
      // sul device aggiorniamo la conf rimettendo quella del passo precedente
    }
  }

  if(metro==0){// accetta sempre in fase di termalizzazione
    acc++;
  }

  dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  if(metro==0){
    printf("   FULL UPDATE COMPUTATION TIME NOMETRO - Tot time : %f sec \n",dt_tot);
    printf("\t\tPreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
    printf("\t\tPreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
    printf("\t\tPostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);
  }
  if(metro==1){
    printf("   FULL UPDATE COMPUTATION TIME SIMETRO - Tot time : %f sec  \n",dt_tot);
    printf("\t\tPreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
    printf("\t\tPreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
    printf("\t\tPostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);
  }

  return acc;

}

#endif
