void THERM_UPDATE_SOLOACC_NOMETRO(su3_soa *conf_acc,double res_metro, double res_md, int id_iter){

  const double lambda=0.1931833275037836; // Omelyan Et Al.
  const double gs=0.5/(double) gauge_scale_acc;
  double delta[7];

  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

  double fattore_arbitrario = 1.0;
  delta[0]= -cimag(ieps_acc) * lambda;
  delta[1]= -cimag(ieps_acc) * (1.0-2.0*lambda);
  delta[2]= -cimag(ieps_acc) * 2.0*lambda;
  delta[3]= -cimag(ieps_acc) * gs*lambda * beta_by_three;
  delta[4]=  cimag(iepsh_acc)* gs;
  delta[5]= -cimag(ieps_acc) * gs*(1.0-2.0*lambda)*beta_by_three;
  delta[6]= -cimag(ieps_acc) * gs*2.0*lambda*beta_by_three;

  int index = 0;

  double placchetta;
  int usestoredeigen;
  int iterazioni = id_iter;
  double dt_tot;
  double dt_pretrans_to_preker;
  double dt_preker_to_postker;
  double dt_postker_to_posttrans;


  double action_in;
  double action_fin;
  double action_mom_in;
  double action_mom_fin;
  double action_ferm_in;
  double action_ferm_fin;
  int mu;
  double minmaxeig[2]; 











#pragma acc data copyin(delta[0:7])
  {
    gettimeofday ( &t1, NULL );

    placchetta =  calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
    printf("iterazione %i   -iniziale-    PLACCHETTA CALCOLATA SU OPENACC  %.18lf  \n",iterazioni,placchetta/(2.0*sizeh*6.0*3.0));

    generate_Momenta_gauss(momenta);
#pragma acc update device(momenta[0:8])
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	vec3_soa *temp = &ferm_phi_acc[iflav][ips];
	generate_vec3_soa_gauss(temp);
    //    generate_MultiFermion_gauss(ferm_phi_acc);
#pragma acc update device(temp[0:1])
      }
    }

    // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
    generate_vec3_soa_gauss(kloc_p);
    generate_vec3_soa_gauss(kloc_s);
    // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])

    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      usestoredeigen = 1; // quindi li calcola
      find_min_max_eigenvalue_soloopenacc(conf_acc,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,usestoredeigen,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
      usestoredeigen = 0;
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig);
#pragma acc update device(approx_fi[0:1])
    }

















      // first_inv_approx_calc
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++)
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++)
	{
	  multishift_invert(conf_acc, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
	  recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), &(ferm_chi_acc[iflav][ips]),&(fermions_parameters[iflav].approx_fi));
	}
    // rescale the approx for the md
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      RationalApprox *approx_md = &(fermions_parameters[iflav].approx_md);
      RationalApprox *approx_md_mother = &(fermions_parameters[iflav].approx_md_mother);
     rescale_rational_approximation(approx_md_mother,approx_md,minmaxeig);
#pragma acc update device(approx_md[0:1])
    }
    
    // dinamica molecolare
    multistep_2MN_SOLOOPENACC(ipdot_acc,conf_acc,u1_back_field_phases,aux_conf_acc,fermions_parameters,NDiffFlavs,ferm_chi_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm,momenta,local_sums,delta,res_md);
    
   
    placchetta =  calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
    printf("iterazione %i  -finale-     PLACCHETTA CALCOLATA SU OPENACC  %.18lf  \n",iterazioni,placchetta/(2.0*sizeh*6.0*3.0));
    
    gettimeofday ( &t2, NULL );
    dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("iterazione %i       AccUpdateTime   : %f sec  \n",iterazioni,dt_preker_to_postker);
    
  }
  
  gettimeofday ( &t3, NULL );
  
  dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("Id_iter %i   FULL UPDATE COMPUTATION TIME                    Tot time          : %f sec  \n",id_iter,dt_tot);
  printf("Id_iter %i                                                   PreTrans->Preker  : %f sec  \n",id_iter,dt_pretrans_to_preker);
  printf("Id_iter %i                                                   PreKer->PostKer   : %f sec  \n",id_iter,dt_preker_to_postker);
  printf("Id_iter %i                                                   PostKer->PostTrans: %f sec  \n",id_iter,dt_postker_to_posttrans);
}

// modificato per versatilita' fermionica fino qui.




int UPDATE_SOLOACC_UNOSTEP_versatile(su3_soa *conf_acc,double res_metro, double res_md, int id_iter,int acc, int metro){

  const double lambda=0.1931833275037836; // Omelyan Et Al.
  const double gs=0.5/(double) gauge_scale_acc;
  double delta[7];

  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

  double fattore_arbitrario = 1.0;
  delta[0]= -cimag(ieps_acc) * lambda;
  delta[1]= -cimag(ieps_acc) * (1.0-2.0*lambda);
  delta[2]= -cimag(ieps_acc) * 2.0*lambda;
  delta[3]= -cimag(ieps_acc) * gs*lambda * beta_by_three;
  delta[4]=  cimag(iepsh_acc)* gs;
  delta[5]= -cimag(ieps_acc) * gs*(1.0-2.0*lambda)*beta_by_three;
  delta[6]= -cimag(ieps_acc) * gs*2.0*lambda*beta_by_three;

  int index = 0;

  double placchetta;
  int usestoredeigen;
  int iterazioni = id_iter + 1;
  double dt_tot;
  double dt_pretrans_to_preker;
  double dt_preker_to_postker;
  double dt_postker_to_posttrans;


  double action_in;
  double action_fin;
  double action_mom_in;
  double action_mom_fin;
  double action_ferm_in;
  double action_ferm_fin;
  int mu;
  double minmaxeig[2]; 

  double p1;
  double p2;
  int accettata;
  double delta_S;

  if(metro==1){
  // store old conf   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
  set_su3_soa_to_su3_soa(conf_acc,conf_acc_bkp);
  }

#pragma acc data copyin(delta[0:7])
  {
    gettimeofday ( &t1, NULL );
    if(metro==0){
    placchetta =  calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
    printf("iterazione %i   -iniziale-    PLACCHETTA CALCOLATA SU OPENACC  %.18lf  \n",iterazioni,placchetta/(2.0*sizeh*6.0*3.0));
    }
    generate_Momenta_gauss(momenta);
#pragma acc update device(momenta[0:8])
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	vec3_soa *temp = &ferm_phi_acc[iflav][ips];
	generate_vec3_soa_gauss(temp);
	//    generate_MultiFermion_gauss(ferm_phi_acc);
#pragma acc update device(temp[0:1])
      }
    }
      
    // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
    generate_vec3_soa_gauss(kloc_p);
    generate_vec3_soa_gauss(kloc_s);
    // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])

    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      usestoredeigen = 1; // quindi li calcola
      find_min_max_eigenvalue_soloopenacc(conf_acc,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,usestoredeigen,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
      usestoredeigen = 0;
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig);
#pragma acc update device(approx_fi[0:1])
      }
    if(metro==1){
      /////////////// INITIAL ACTION COMPUTATION ////////////////////////////////////////////
      action_in = beta_by_three*calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
      printf("iterazione %i  -iniziale-   update  PLACCHETTA CALCOLATA SU OPENACC  %.18lf  \n",iterazioni,action_in/(2.0*sizeh*6.0*3.0*beta_by_three));
      action_mom_in = 0.0;
      for(mu =0;mu<8;mu++){
	action_mom_in += calc_momenta_action(momenta,d_local_sums,mu);
      }
      action_ferm_in=scal_prod_between_multiferm(ferm_phi_acc,ferm_phi_acc);
      ///////////////////////////////////////////////////////////////////////////////////////
    }
      // first_inv_approx_calc 
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++)
          for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++)
          {
              multishift_invert(conf_acc, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), res_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
              recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc[ips], &(ferm_phi_acc[iflav][ips]), &(ferm_chi_acc[iflav][ips]),&(fermions_parameters[iflav].approx_fi));
          }
      // rescale the approx for the md
      for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
	RationalApprox *approx_md = &(fermions_parameters[iflav].approx_md);
	RationalApprox *approx_md_mother = &(fermions_parameters[iflav].approx_md_mother);
	rescale_rational_approximation(approx_md_mother,approx_md,minmaxeig);
#pragma acc update device(approx_md[0:1])
      }
      
      // dinamica molecolare
    multistep_2MN_SOLOOPENACC(ipdot_acc,conf_acc,u1_back_field_phases,aux_conf_acc,fermions_parameters,NDiffFlavs,ferm_chi_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm,momenta,local_sums,delta,res_md);
     
    if(metro==0){   
    placchetta =  calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
    printf("iterazione %i  -finale-     PLACCHETTA CALCOLATA SU OPENACC  %.18lf  \n",iterazioni,placchetta/(2.0*sizeh*6.0*3.0));
    
    gettimeofday ( &t2, NULL );
    dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("iterazione %i       AccUpdateTime   : %f sec  \n",iterazioni,dt_preker_to_postker);
    }
    else{
     
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      // compute the highest and lowest eigenvalues of (M^dag M)
      usestoredeigen = 1; // quindi li calcola
      find_min_max_eigenvalue_soloopenacc(conf_acc,kloc_r,kloc_h,kloc_p,kloc_s,usestoredeigen,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
      usestoredeigen = 0;
      rescale_rational_approximation(approx_mother3,approx3,minmaxeig,-(1.0/4.0));
#pragma acc update device(approx3[0:1])
      // last_inv_approx_calc
      ker_invert_openacc_shiftmulti(conf_acc,ferm_shiftmulti_acc,ferm_chi_acc,res_metro,approx3,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
      ker_openacc_recombine_shiftmulti_to_multi(ferm_shiftmulti_acc,ferm_chi_acc,ferm_phi_acc,approx3);


      ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
      action_fin = beta_by_three * calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
      printf("iterazione %i  -finale-   update  PLACCHETTA CALCOLATA SU OPENACC  %.18lf  \n",iterazioni,action_fin/(2.0*sizeh*6.0*3.0*beta_by_three));
      action_mom_fin = 0.0;
      for(mu =0;mu<8;mu++){
	action_mom_fin += calc_momenta_action(momenta,d_local_sums,mu);
      }
      action_ferm_fin=scal_prod_between_multiferm(ferm_chi_acc,ferm_phi_acc);
      ////////////////////////////////////////////////////////////////////////////////////////

      gettimeofday ( &t2, NULL );

      // delta_S = action_new - action_old
      delta_S  = - (-action_in+action_mom_in+action_ferm_in) + (-action_fin+action_mom_fin+action_ferm_fin);

      printf("iterazione %i   INSIDE PRAGMA DATA - DELTA_ACTION  = %.18lf\n",iterazioni,delta_S);

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
  }
    
  gettimeofday ( &t3, NULL );
  
  if(metro==0){
  
  dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("Id_iter %i   FULL UPDATE COMPUTATION TIME                    Tot time          : %f sec  \n",id_iter,dt_tot);
  printf("Id_iter %i                                                   PreTrans->Preker  : %f sec  \n",id_iter,dt_pretrans_to_preker);
  printf("Id_iter %i                                                   PreKer->PostKer   : %f sec  \n",id_iter,dt_preker_to_postker);
  printf("Id_iter %i                                                   PostKer->PostTrans: %f sec  \n",id_iter,dt_postker_to_posttrans);
  }else{

  if(accettata==1){
    acc +=1;
    printf("ACCEPTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
    // configuration accepted   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
    set_su3_soa_to_su3_soa(conf_acc,conf_acc_bkp);
  }else{
    printf("REJECTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
    // configuration rejected   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
    set_su3_soa_to_su3_soa(conf_acc_bkp,conf_acc);
#pragma acc update device(conf_acc[0:8])
    // sul device aggiorniamo la conf rimettendo quella del passo precedente
  }

  //  for(dir=0;dir<8;dir++)  convert_su3_soa_to_su3COM_soa(&conf_acc[dir],&conf[dir]);


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

}    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
      // compute the highest and lowest eigenvalues of (M^dag M)
      usestoredeigen = 1; // quindi li calcola
      find_min_max_eigenvalue_soloopenacc(conf_acc,kloc_r,kloc_h,kloc_p,kloc_s,usestoredeigen,minmaxeig);
#pragma acc update device(minmaxeig[0:2])
      usestoredeigen = 0;
      rescale_rational_approximation(approx_mother3,approx3,minmaxeig,-(1.0/4.0));
#pragma acc update device(approx3[0:1])
      // last_inv_approx_calc
      ker_invert_openacc_shiftmulti(conf_acc,ferm_shiftmulti_acc,ferm_chi_acc,res_metro,approx3,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
      ker_openacc_recombine_shiftmulti_to_multi(ferm_shiftmulti_acc,ferm_chi_acc,ferm_phi_acc,approx3);


      ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
      action_fin = beta_by_three * calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
      printf("iterazione %i  -finale-   update  PLACCHETTA CALCOLATA SU OPENACC  %.18lf  \n",iterazioni,action_fin/(2.0*sizeh*6.0*3.0*beta_by_three));
      action_mom_fin = 0.0;
      for(mu =0;mu<8;mu++){
	action_mom_fin += calc_momenta_action(momenta,d_local_sums,mu);
      }
      action_ferm_fin=scal_prod_between_multiferm(ferm_chi_acc,ferm_phi_acc);
      ////////////////////////////////////////////////////////////////////////////////////////

      gettimeofday ( &t2, NULL );

      // delta_S = action_new - action_old
      delta_S  = - (-action_in+action_mom_in+action_ferm_in) + (-action_fin+action_mom_fin+action_ferm_fin);

      printf("iterazione %i   INSIDE PRAGMA DATA - DELTA_ACTION  = %.18lf\n",iterazioni,delta_S);

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
  gettimeofday ( &t3, NULL );

  if(accettata==1){
    acc +=1;
    printf("ACCEPTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
    // configuration accepted   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
    set_su3_soa_to_su3_soa(conf_acc,conf_acc_bkp);
  }else{
    printf("REJECTED   ---> [acc/iter] = [%i/%i] \n",acc,iterazioni);
    // configuration rejected   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
    set_su3_soa_to_su3_soa(conf_acc_bkp,conf_acc);
#pragma acc update device(conf_acc[0:8])
    // sul device aggiorniamo la conf rimettendo quella del passo precedente
  }




  //  for(dir=0;dir<8;dir++)  convert_su3_soa_to_su3COM_soa(&conf_acc[dir],&conf[dir]);


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
