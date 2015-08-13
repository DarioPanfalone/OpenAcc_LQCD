// se metro==1 allora fa il test di metropolis
// se metro==0 allora non fa il test di metropolis --> termalizzazione
int UPDATE_SOLOACC_UNOSTEP_VERSATILE(su3_soa *tconf_acc,double res_metro, double res_md, int id_iter,int acc,int metro){

  
  printf("############################################### \n");
  printf("## Inside UPDATE_SOLOACC_UNOSTEP_VERSATILE() ## \n");
  printf("############################################### \n");



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

  int iterazioni = id_iter;
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
    printf("Backup copy of the initial gauge conf : OK \n");
  }

#pragma acc data copyin(delta[0:7])
  {
    gettimeofday ( &t1, NULL );

    // ESTRAZIONI RANDOM
    generate_Momenta_gauss(momenta);
    printf("Momenta generated : OK \n");
#pragma acc update device(momenta[0:8])
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	  printf("Ferm generation (flav=%d,ps=%d) : OK \n",iflav,ips);
          int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	  printf("   find the index: OK \n"); 
	  vec3_soa *temp = &ferm_phi_acc[ps_index];
	  printf("   copy the pointer: OK \n"); 
	  generate_vec3_soa_gauss(temp);
	  printf("   generate the ferm: OK \n\n"); 
#pragma acc update device(temp[0:1])
      }
    }// end for iflav
      
    // STIRACCHIAMENTO DELL'APPROX RAZIONALE FIRST_INV
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      printf("Rat approx rescale (flav=%d) : OK \n",iflav);
      // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
      generate_vec3_soa_gauss(kloc_p);
      printf("    generate kloc_p: OK \n");
      generate_vec3_soa_gauss(kloc_s);
      printf("    generate kloc_s: OK \n");
      // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      find_min_max_eigenvalue_soloopenacc(tconf_acc,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
      printf("    find min and max eig : OK \n");
#pragma acc update device(minmaxeig[0:2])
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      printf("    get approx_fi pointer : OK \n");
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      printf("    get approx_fi_mother pointer : OK \n");
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig);
      printf("    rat approx rescaled : OK \n\n");
#pragma acc update device(approx_fi[0:1])
    }//end for iflav
    
    if(metro==1){
      /////////////// INITIAL ACTION COMPUTATION ////////////////////////////////////////////
      action_in = beta_by_three*calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
      printf("Gauge action computed : OK \n");
      action_mom_in = 0.0;
      for(mu =0;mu<8;mu++)  action_mom_in += calc_momenta_action(momenta,d_local_sums,mu);
      printf("Momenta action computed : OK \n");
      action_ferm_in=0;
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
          for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){

              int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
              action_ferm_in += real_scal_prod_global(&ferm_phi_acc[ps_index],&ferm_phi_acc[ps_index]);
          }
      }// end for iflav
      printf("Ferm action computed : OK \n\n");
      ///////////////////////////////////////////////////////////////////////////////////////
    }

    // FIRST INV APPROX CALC --> calcolo del fermione CHI
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	  printf("First inv approx calc (flav=%d,ps=%d) : OK \n",iflav,ips);
          int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	  printf("    determined the index : OK \n");
	  multishift_invert(tconf_acc, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
	  printf("    computed the inverse : OK \n");
	  recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), &(ferm_chi_acc[ps_index]),&(fermions_parameters[iflav].approx_fi));
	  printf("    recombined the fermions : OK \n");
      }
    }// end for iflav
    
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

    // DINAMICA MOLECOLARE
    multistep_2MN_SOLOOPENACC(ipdot_acc,tconf_acc,u1_back_field_phases,aux_conf_acc,fermions_parameters,NDiffFlavs,ferm_chi_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm,momenta,local_sums,delta,res_md);


    if(metro==1){
      // STIRACCHIAMENTO DELL'APPROX RAZIONALE LAST_INV
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
	// generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
	generate_vec3_soa_gauss(kloc_p);
	generate_vec3_soa_gauss(kloc_s);
	// update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
	find_min_max_eigenvalue_soloopenacc(tconf_acc,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
	RationalApprox *approx_li = &(fermions_parameters[iflav].approx_li);
	RationalApprox *approx_li_mother = &(fermions_parameters[iflav].approx_li_mother);
	rescale_rational_approximation(approx_li_mother,approx_li,minmaxeig);
#pragma acc update device(approx_li[0:1])
      }

      // LAST INV APPROX CALC 
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
	for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
        int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
        multishift_invert(tconf_acc, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_li), u1_back_field_phases, ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
        recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), &(ferm_phi_acc[ps_index]),&(fermions_parameters[iflav].approx_li));
	}
      }
      
      ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
      action_fin = beta_by_three * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
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
      printf("iterazione %i:  DELTA_ACTION                       = %.18lf\n",iterazioni,delta_S);
      
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
    printf("Id_iter %i   FULL UPDATE COMPUTATION TIME NOMETRO            Tot time          : %f sec  \n",id_iter,dt_tot);
    printf("Id_iter %i                                                   PreTrans->Preker  : %f sec  \n",id_iter,dt_pretrans_to_preker);
    printf("Id_iter %i                                                   PreKer->PostKer   : %f sec  \n",id_iter,dt_preker_to_postker);
    printf("Id_iter %i                                                   PostKer->PostTrans: %f sec  \n",id_iter,dt_postker_to_posttrans);
  }
  if(metro==1){
    printf("Id_iter %i   FULL UPDATE COMPUTATION TIME SIMETRO            Tot time          : %f sec  \n",id_iter,dt_tot);
    printf("Id_iter %i                                                   PreTrans->Preker  : %f sec  \n",id_iter,dt_pretrans_to_preker);
    printf("Id_iter %i                                                   PreKer->PostKer   : %f sec  \n",id_iter,dt_preker_to_postker);
    printf("Id_iter %i                                                   PostKer->PostTrans: %f sec  \n",id_iter,dt_postker_to_posttrans);
  }

  return acc;

}


/*
void THERM_UPDATE_SOLOACC_NOMETRO(su3_soa *tconf_acc,double res_metro, double res_md, int id_iter){

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

  int iterazioni = id_iter;
  double dt_tot;
  double dt_pretrans_to_preker;
  double dt_preker_to_postker;
  double dt_postker_to_posttrans;

  int mu;
  double minmaxeig[2]; 

  
  
  
  
  
  
  
  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );




#pragma acc data copyin(delta[0:7])
  {
    gettimeofday ( &t1, NULL );

    // ESTRAZIONI RANDOM
    generate_Momenta_gauss(momenta);
#pragma acc update device(momenta[0:8])
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	vec3_soa *temp = &ferm_phi_acc[iflav][ips];
	generate_vec3_soa_gauss(temp);
#pragma acc update device(temp[0:1])
      }
    }//end for iflav

    // STIRACCHIAMENTO DELL'APPROX RAZIONALE FIRST_INV
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
      // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      find_min_max_eigenvalue_soloopenacc(tconf_acc,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig);
#pragma acc update device(approx_fi[0:1])
    }//end for iflav

   
   
   
   
   
   
   
   
   
   
   
   
    // FIRST INV APPROX CALC --> calcolo del fermione CHI
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	multishift_invert(tconf_acc, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
	recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), &(ferm_chi_acc[iflav][ips]),&(fermions_parameters[iflav].approx_fi));
      }
    }// end for iflav

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
    
    // DINAMICA MOLECOLARE
    multistep_2MN_SOLOOPENACC(ipdot_acc,tconf_acc,u1_back_field_phases,aux_conf_acc,fermions_parameters,NDiffFlavs,ferm_chi_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm,momenta,local_sums,delta,res_md);
    


















    // LAST INV APPROX CALC
    












    
    
    
    
    
    
    
    
    
    
    
    
    gettimeofday ( &t2, NULL );

























  }//end pragma acc data copyin(delta[0:7])

  gettimeofday ( &t3, NULL );
 
















  dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("Id_iter %i   FULL UPDATE COMPUTATION TIME                    Tot time          : %f sec  \n",id_iter,dt_tot);
  printf("Id_iter %i                                                   PreTrans->Preker  : %f sec  \n",id_iter,dt_pretrans_to_preker);
  printf("Id_iter %i                                                   PreKer->PostKer   : %f sec  \n",id_iter,dt_preker_to_postker);



}



int UPDATE_SOLOACC_UNOSTEP_METRO(su3_soa *tconf_acc,double res_metro, double res_md, int id_iter,int acc){

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

  int iterazioni = id_iter;
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

  // store old conf   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
  set_su3_soa_to_su3_soa(tconf_acc,conf_acc_bkp);

#pragma acc data copyin(delta[0:7])
  {
    gettimeofday ( &t1, NULL );

    // ESTRAZIONI RANDOM
    generate_Momenta_gauss(momenta);
#pragma acc update device(momenta[0:8])
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	vec3_soa *temp = &ferm_phi_acc[iflav][ips];
	generate_vec3_soa_gauss(temp);
#pragma acc update device(temp[0:1])
      }
    }// end for iflav
      
    // STIRACCHIAMENTO DELL'APPROX RAZIONALE FIRST_INV
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
      // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      find_min_max_eigenvalue_soloopenacc(tconf_acc,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig);
#pragma acc update device(approx_fi[0:1])
    }//end for iflav
    
    /////////////// INITIAL ACTION COMPUTATION ////////////////////////////////////////////
    action_in = beta_by_three*calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
    action_mom_in = 0.0;
    for(mu =0;mu<8;mu++)  action_mom_in += calc_momenta_action(momenta,d_local_sums,mu);
    action_ferm_in=scal_prod_between_multiferm(ferm_phi_acc,ferm_phi_acc);
    ///////////////////////////////////////////////////////////////////////////////////////


    // FIRST INV APPROX CALC --> calcolo del fermione CHI
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	multishift_invert(tconf_acc, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
	recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), &(ferm_chi_acc[iflav][ips]),&(fermions_parameters[iflav].approx_fi));
      }
    }// end for iflav
    
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

    // DINAMICA MOLECOLARE
    multistep_2MN_SOLOOPENACC(ipdot_acc,tconf_acc,u1_back_field_phases,aux_conf_acc,fermions_parameters,NDiffFlavs,ferm_chi_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm,momenta,local_sums,delta,res_md);
    // STIRACCHIAMENTO DELL'APPROX RAZIONALE LAST_INV
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
      // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      find_min_max_eigenvalue_soloopenacc(tconf_acc,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
      RationalApprox *approx_li = &(fermions_parameters[iflav].approx_li);
      RationalApprox *approx_li_mother = &(fermions_parameters[iflav].approx_li_mother);
      rescale_rational_approximation(approx_li_mother,approx_li,minmaxeig);
#pragma acc update device(approx_li[0:1])
    }
    // LAST INV APPROX CALC 
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	multishift_invert(tconf_acc, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_li), u1_back_field_phases, ferm_shiftmulti_acc[ips], &(ferm_chi_acc[iflav][ips]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
	recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc[ips], &(ferm_chi_acc[iflav][ips]), &(ferm_phi_acc[iflav][ips]),&(fermions_parameters[iflav].approx_li));
      }
    }
    
    ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
    action_fin = beta_by_three * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
    action_mom_fin = 0.0;
    for(mu =0;mu<8;mu++)    action_mom_fin += calc_momenta_action(momenta,d_local_sums,mu);
    action_ferm_fin=scal_prod_between_multiferm(ferm_chi_acc,ferm_phi_acc);
    ////////////////////////////////////////////////////////////////////////////////////////

    gettimeofday ( &t2, NULL );

    // delta_S = action_new - action_old
    delta_S  = - (-action_in+action_mom_in+action_ferm_in) + (-action_fin+action_mom_fin+action_ferm_fin);
    printf("iterazione %i:  DELTA_ACTION                       = %.18lf\n",iterazioni,delta_S);
    
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
    
  }// end pragma acc data copyin(delta[0:7])

  gettimeofday ( &t3, NULL );

  if(accettata==1){
    acc = 1;
    printf("ACCEPTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
    // configuration accepted   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
    set_su3_soa_to_su3_soa(tconf_acc,conf_acc_bkp);
  }else{
    acc = 0;
    printf("REJECTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
    // configuration rejected   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
    set_su3_soa_to_su3_soa(conf_acc_bkp,tconf_acc);
#pragma acc update device(tconf_acc[0:8])
    // sul device aggiorniamo la conf rimettendo quella del passo precedente
  }
  
  dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("Id_iter %i   FULL UPDATE COMPUTATION TIME                    Tot time          : %f sec  \n",id_iter,dt_tot);
  printf("Id_iter %i                                                   PreTrans->Preker  : %f sec  \n",id_iter,dt_pretrans_to_preker);
  printf("Id_iter %i                                                   PreKer->PostKer   : %f sec  \n",id_iter,dt_preker_to_postker);
  printf("Id_iter %i                                                   PostKer->PostTrans: %f sec  \n",id_iter,dt_postker_to_posttrans);

  return acc;

}
*/

