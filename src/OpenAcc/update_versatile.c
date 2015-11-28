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
#include "./inverter_multishift_full.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "./fermionic_utilities.h"
#include "sys/time.h"

extern double casuale();



#define PRINT_DETAILS_INSIDE_UPDATE

int UPDATE_SOLOACC_UNOSTEP_VERSATILE(su3_soa *tconf_acc,
        double res_metro, double res_md, int id_iter,int acc,int metro){
  
#ifdef STOUT_FERMIONS        
    su3_soa *tstout_conf_acc_arr = gstout_conf_acc_arr;
#endif



  
  printf("UPDATE_SOLOACC_UNOSTEP_VERSATILE_TLSM_STDFERM: OK \n");
  // DEFINIZIONE DI TUTTI I dt NECESSARI PER L'INTEGRATORE OMELYAN
  const double lambda=0.1931833275037836; // Omelyan Et Al.
  const double gs=0.5/(double) gauge_scale_acc;
  double delta[7];
  delta[0]= -cimag(ieps_acc) * lambda;
  delta[1]= -cimag(ieps_acc) * (1.0-2.0*lambda);
  delta[2]= -cimag(ieps_acc) * 2.0*lambda;
  delta[3]= -cimag(ieps_acc) * gs*lambda * beta_by_three;
  delta[4]=  cimag(iepsh_acc)* gs;
  delta[5]= -cimag(ieps_acc) * gs*(1.0-2.0*lambda)*beta_by_three;
  delta[6]= -cimag(ieps_acc) * gs*2.0*lambda*beta_by_three;

  int iterazioni = id_iter+1;
  double dt_tot;
  double dt_pretrans_to_preker;
  double dt_preker_to_postker;
  double dt_postker_to_posttrans;

  int mu;
  double minmaxeig[2]; 

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
#ifdef PRINT_DETAILS_INSIDE_UPDATE
    printf("Backup copy of the initial gauge conf : OK \n");
#endif

  }

#pragma acc data copyin(delta[0:7])
  {
  
    gettimeofday ( &t1, NULL );

    // ESTRAZIONI RANDOM
    generate_Momenta_gauss(momenta);
#ifdef PRINT_DETAILS_INSIDE_UPDATE
    printf("Momenta generated : OK \n");
#endif

    //    read_thmat_soa(momenta,"momenta");
#pragma acc update device(momenta[0:8])
    
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
      int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
#ifdef PRINT_DETAILS_INSIDE_UPDATE
	  printf("Ferm generation (flav=%d,ps=%d,ps_index=%d) : OK \n",iflav,ips,ps_index);
#endif
	  generate_vec3_soa_gauss(&ferm_phi_acc[ps_index]);
      }
    }// end for iflav
#pragma acc update device(ferm_phi_acc[0:NPS_tot])

#ifdef STOUT_FERMIONS 
    // USO DELLA VERSIONE STOUTATA GIA' PER LO STIRACCHIAMENTO
    // STOUTING...(ALREADY ON DEVICE)
    stout_wrapper(tconf_acc,tstout_conf_acc_arr);
    gconf_as_fermionmatrix = &(tstout_conf_acc_arr[8*(STOUT_STEPS-1)]);
#else
    gconf_as_fermionmatrix = tconf_acc;
#endif




    // STIRACCHIAMENTO DELL'APPROX RAZIONALE FIRST_INV
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
#ifdef PRINT_DETAILS_INSIDE_UPDATE
      printf("Rat approx rescale (flav=%d) : OK \n",iflav);
#endif

      // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
      SETREQUESTED(kloc_p);
      generate_vec3_soa_gauss(kloc_p);
      SETREQUESTED(kloc_s);
      generate_vec3_soa_gauss(kloc_s);
      // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      // USING STOUTED GAUGE MATRIX
      //printf("    before min and max eig comp : OK \n");
      SETREQUESTED(kloc_r);
      SETREQUESTED(kloc_h);
      find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
#ifdef PRINT_DETAILS_INSIDE_UPDATE
      printf("    find min and max eig : OK \n");
#endif
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig);
#ifdef PRINT_DETAILS_INSIDE_UPDATE
      printf("    rat approx rescaled : OK \n\n");
#endif

#pragma acc update device(approx_fi[0:1])
    }//end for iflav
    
    if(metro==1){
      /////////////// INITIAL ACTION COMPUTATION ////////////////////////////////////////////
      if(GAUGE_ACTION == 0)// Standard gauge action 
	action_in = beta_by_three*calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      if(GAUGE_ACTION == 1){ //Tlsym gauge action
	action_in = C_ZERO * beta_by_three * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
	action_in += C_ONE * beta_by_three * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
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
#ifdef PRINT_DETAILS_INSIDE_UPDATE
    printf(" Initial Action Computed : OK \n");
#endif

    // FIRST INV APPROX CALC --> calcolo del fermione CHI

    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	
        int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
        // USING STOUTED GAUGE MATRIX
        multishift_invert(gconf_as_fermionmatrix, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
        recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), &(ferm_chi_acc[ps_index]),&(fermions_parameters[iflav].approx_fi));
	
      }
    }// end for iflav
#ifdef PRINT_DETAILS_INSIDE_UPDATE
    printf(" Computed the fermion CHI : OK \n");
#endif
    
    // STIRACCHIAMENTO DELL'APPROX RAZIONALE MD
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      // recupero gli autovalori già calcolati che sono salvati in approx_fi
      minmaxeig[0] = (fermions_parameters[iflav].approx_fi.lambda_min);
      minmaxeig[1] = (fermions_parameters[iflav].approx_fi.lambda_max);
      minmaxeig[0] = minmaxeig[0] / minmaxeig[1];
      minmaxeig[0] = minmaxeig[0] / 0.95;
      minmaxeig[1] = minmaxeig[1] / 1.05;
      RationalApprox *approx_md = &(fermions_parameters[iflav].approx_md);
      RationalApprox *approx_md_mother = &(fermions_parameters[iflav].approx_md_mother);
      rescale_rational_approximation(approx_md_mother,approx_md,minmaxeig);
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
			      k_p_shiftferm,momenta,local_sums,delta,res_md);
    
#ifdef PRINT_DETAILS_INSIDE_UPDATE
    printf(" Molecular Dynamics Completed : OK \n");
#endif


#ifdef STOUT_FERMIONS
    // STOUTING...(ALREADY ON DEVICE)
    stout_wrapper(tconf_acc,tstout_conf_acc_arr);
    gconf_as_fermionmatrix = &(tstout_conf_acc_arr[8*(STOUT_STEPS-1)]);
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
	find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
	//#pragma acc update device(minmaxeig[0:2])
	RationalApprox *approx_li = &(fermions_parameters[iflav].approx_li);
	RationalApprox *approx_li_mother = &(fermions_parameters[iflav].approx_li_mother);
	rescale_rational_approximation(approx_li_mother,approx_li,minmaxeig);
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
#ifdef PRINT_DETAILS_INSIDE_UPDATE
      printf(" Final Action Computed : OK \n");
#endif
      
      ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
      if(GAUGE_ACTION == 0) // Standard gauge action
      action_fin = beta_by_three * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      if(GAUGE_ACTION == 1){ // Tlsym gauge action
      action_fin = C_ZERO * beta_by_three * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      action_fin += C_ONE * beta_by_three * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
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
#ifdef PRINT_DETAILS_INSIDE_UPDATE
      printf("iterazione %i:  Gauge_ACTION  (in and out) = %.18lf , %.18lf\n",iterazioni,-action_in,-action_fin);
      printf("iterazione %i:  Momen_ACTION  (in and out) = %.18lf , %.18lf\n",iterazioni,action_mom_in,action_mom_fin);
      printf("iterazione %i:  Fermi_ACTION  (in and out) = %.18lf , %.18lf\n",iterazioni,action_ferm_in,action_ferm_fin);
#endif
      printf("iterazione %i:  DELTA_ACTION               = %.18lf\n",iterazioni,delta_S);
      
      if(delta_S<0){
	accettata=1;
      }
      else
	{
	  p1=exp(-delta_S);
	  p2=casuale();
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
    printf("   FULL UPDATE COMPUTATION TIME NOMETRO            Tot time          : %f sec  \n",dt_tot);
    printf("                                                   PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
    printf("                                                   PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
    printf("                                                   PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);
  }
  if(metro==1){
    printf("   FULL UPDATE COMPUTATION TIME SIMETRO            Tot time          : %f sec  \n",dt_tot);
    printf("                                                   PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
    printf("                                                   PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
    printf("                                                   PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);
  }

  return acc;

}

#endif
