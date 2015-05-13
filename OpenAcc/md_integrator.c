// MULTISTEP VERSION of minimum_norm2_B

// 2nd order Minimum Norm integrator (2MN)
// See reference hep-lat/0505020 Takaishi, De Forcrand
//
// Scheme (needs three force calculations per step):
// 2MN(dt) = exp (-l * dt * dS/dq) exp (dt/2 * p) * ...
//           exp (-(1- 2l) * dt * dS/dq) exp (dt/2 * p)  exp (-l * dt * dS/dq)
// p       : momenta
// - dS/sq : force contibution
// dt      : integration step
// l       : lambda parameter (1/6 Sexton Weingarten, 0.193... Omelyan et al)
//    
// Total scheme:
// [2MN]^N
// N       : number of Molecular Dynamics steps
//
// See reference hep-lat/0506011 Urbach et al. for the multiple time scale 
//
// p->p+a dS/dq=p+ia(-i dS/dq)=p+ia*ipdot

// the various dt involved are stored into the delta array which is allocated into UPDATE_ACC(args)
void multistep_2MN_gauge(su3_soa *conf_acc,su3_soa *local_staples,tamat_soa *ipdot,thmat_soa *momenta,double * delta)
 {
 int md;
 // Step for the P
 // P' = P - l*dt*dS/dq
 // delta[3]=-cimag(ieps_acc)*scale*lambda;
 calc_ipdot_gauge_soloopenacc(conf_acc,local_staples,ipdot);
 mom_sum_mult(momenta,ipdot,delta,3);
 for(md=1; md<gauge_scale; md++){
   // Step for the Q
   // Q' = exp[dt/2 *i P] Q
   // delta[4]=cimag(iepsh_acc)*scale;
   mom_exp_times_conf_soloopenacc(conf_acc,momenta,delta,4);
   // Step for the P
   // P' = P - (1-2l)*dt*dS/dq
   // delta[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
   calc_ipdot_gauge_soloopenacc(conf_acc,local_staples,ipdot);
   mom_sum_mult(momenta,ipdot,delta,5);
   // Step for the Q
   // Q' = exp[dt/2 *i P] Q
   // delta[4]=cimag(iepsh_acc)*scale;
   mom_exp_times_conf_soloopenacc(conf_acc,momenta,delta,4);
   // Step for the P
   // P' = P - 2l*dt*dS/dq
   // delta[6]=-cimag(ieps_acc)*2.0*lambda*scale;
   calc_ipdot_gauge_soloopenacc(conf_acc,local_staples,ipdot);
   mom_sum_mult(momenta,ipdot,delta,6);
 }
 
 // Step for the Q
 // Q' = exp[dt/2 *i P] Q
 // delta[4]=cimag(iepsh_acc)*scale;
 mom_exp_times_conf_soloopenacc(conf_acc,momenta,delta,4);
 // Step for the P
 // P' = P - (1-2l)*dt*dS/dq
 calc_ipdot_gauge_soloopenacc(conf_acc,local_staples,ipdot);
 // calc_ipdot_gauge();
 // delta[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
 mom_sum_mult(momenta,ipdot,delta,5);
 // Step for the Q
 // Q' = exp[dt/2 *i P] Q
 // delta[4]=cimag(iepsh_acc)*scale;
 mom_exp_times_conf_soloopenacc(conf_acc,momenta,delta,4);
 // Step for the P
 // P' = P - l*dt*dS/dq
 // delta[3]=-cimag(ieps_acc)*lambda*scale;
 calc_ipdot_gauge_soloopenacc(conf_acc,local_staples,ipdot);
 mom_sum_mult(momenta,ipdot,delta,3);
 
 }
void multistep_2MN_SOLOOPENACC( tamat_soa * ipdot_acc,
				su3_soa  * conf_acc,
				su3_soa  * aux_conf_acc,
				ACC_MultiFermion * ferm_in_acc,
				ACC_MultiFermion * ferm_out_acc,
				ACC_ShiftMultiFermion * ferm_shiftmulti_acc,
				vec3_soa * kloc_r,
				vec3_soa * kloc_h,
				vec3_soa * kloc_s,
				vec3_soa * kloc_p,
				ACC_ShiftFermion *k_p_shiftferm,
				thmat_soa * momenta,
				dcomplex_soa * local_sums,
				double * delta,
				double res,
				COM_RationalApprox *approx)
{
  
  int md;
  
  // Step for the P
  // P' = P - l*dt*dS/dq
  //    delta[0]=-cimag(ieps_acc)*lambda;
  fermion_force_soloopenacc(conf_acc,ipdot_acc,ferm_in_acc,res,approx,ferm_out_acc,aux_conf_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
  mom_sum_mult(momenta,ipdot_acc,delta,0);
  
  for(md=1; md<no_md; md++){
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    multistep_2MN_gauge(conf_acc,aux_conf_acc,ipdot_acc,momenta,delta);
    // Step for the P
    // P' = P - (1-2l)*dt*dS/dq
    // delta[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
    fermion_force_soloopenacc(conf_acc,ipdot_acc,ferm_in_acc,res,approx,ferm_out_acc,aux_conf_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
    mom_sum_mult(momenta,ipdot_acc,delta,1);
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    multistep_2MN_gauge(conf_acc,aux_conf_acc,ipdot_acc,momenta,delta);
    // Step for the P
    // P' = P - 2l*dt*dS/dq
    // delta[2]=-cimag(ieps_acc)*(2.0*lambda);
    fermion_force_soloopenacc(conf_acc,ipdot_acc,ferm_in_acc,res,approx,ferm_out_acc,aux_conf_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
    mom_sum_mult(momenta,ipdot_acc,delta,2);
  }  
  // Step for the Q
  // Q' = exp[dt/2 *i P] Q
  multistep_2MN_gauge(conf_acc,aux_conf_acc,ipdot_acc,momenta,delta);
  // Step for the P
  // P' = P - (1-2l)*dt*dS/dq
  // delta[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
  fermion_force_soloopenacc(conf_acc,ipdot_acc,ferm_in_acc,res,approx,ferm_out_acc,aux_conf_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
  mom_sum_mult(momenta,ipdot_acc,delta,1);
  // Step for the Q
  // Q' = exp[dt/2 *i P] Q
  multistep_2MN_gauge(conf_acc,aux_conf_acc,ipdot_acc,momenta,delta);
  // Step for the P
  // P' = P - l*dt*dS/dq
  // delta[0]=-cimag(ieps_acc)*lambda;
  fermion_force_soloopenacc(conf_acc,ipdot_acc,ferm_in_acc,res,approx,ferm_out_acc,aux_conf_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
  mom_sum_mult(momenta,ipdot_acc,delta,0);
    
}





void UPDATE_ACC(su3COM_soa *conf,double residue_metro,double residue_md,const COM_RationalApprox *approx1,const COM_RationalApprox *approx2,const COM_RationalApprox *approx3,const COM_MultiFermion *in,thmatCOM_soa * com_mom){

  // approx1   --> approssimazione per la prima inversione
  // approx2   --> approssimazione per la dinamica molecolare
  // approx3   --> approssimazione per l'ultima inversione

  // defined in struct_c_def.c
  // contains the assignement of no_md_acc, ieps, ...
  initialize_global_variables();
  compute_nnp_and_nnm_openacc();
  
  tamat_soa * ipdot_acc;
  su3_soa  * conf_acc;
  su3_soa  * aux_conf_acc;
  ACC_MultiFermion * ferm_phi_acc;
  ACC_MultiFermion * ferm_chi_acc;
  ACC_MultiFermion * ferm_aux_acc;
  ACC_ShiftMultiFermion * ferm_shiftmulti_acc;
  vec3_soa * kloc_r;
  vec3_soa * kloc_h;
  vec3_soa * kloc_s;
  vec3_soa * kloc_p;
  ACC_ShiftFermion *k_p_shiftferm;
  thmat_soa * momenta;
  dcomplex_soa * local_sums;
  double_soa * d_local_sums;
  double delta[7];


  posix_memalign((void **)&momenta, ALIGN, 8*sizeof(thmat_soa));   //  -->  4*size
  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&k_p_shiftferm, ALIGN, sizeof(ACC_ShiftFermion));
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&ipdot_acc, ALIGN, 8*sizeof(tamat_soa));
  posix_memalign((void **)&ferm_phi_acc  , ALIGN, sizeof(ACC_MultiFermion));
  posix_memalign((void **)&ferm_chi_acc  , ALIGN, sizeof(ACC_MultiFermion));
  posix_memalign((void **)&ferm_aux_acc  , ALIGN, sizeof(ACC_MultiFermion));
  posix_memalign((void **)&ferm_shiftmulti_acc, ALIGN, sizeof(ACC_ShiftMultiFermion));
  posix_memalign((void **)&local_sums, ALIGN, 2*sizeof(dcomplex_soa));  // --> size complessi --> vettore per sommare cose locali
  posix_memalign((void **)&d_local_sums, ALIGN, sizeof(double_soa));  // --> sizeh double
  
  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);
  convert_COM_MultiFermion_to_ACC_MultiFermion(in,ferm_phi_acc);
  for(dir=0;dir<8;dir++)  convert_thmatCOM_soa_to_thmat_soa(&com_mom[dir],&momenta[dir]);
  
  const double lambda=0.1931833275037836; // Omelyan Et Al.
  const double gs=0.5/(double) gauge_scale_acc;

  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

  double fattore_arbitrario = 1.0;
  delta[0]= -cimag(ieps_acc) * lambda;
  delta[1]= -cimag(ieps_acc) * (1.0-2.0*lambda);
  delta[2]= -cimag(ieps_acc) * 2.0*lambda;
  delta[3]= -cimag(ieps_acc) * gs*lambda;
  delta[4]=  cimag(iepsh_acc)* gs;
  delta[5]= -cimag(ieps_acc) * gs*(1.0-2.0*lambda);
  delta[6]= -cimag(ieps_acc) * gs*2.0*lambda;

  int index = 0;

  printf("MOMENTA AND CONF - INSIDE OPENACC MD INTEGRATOR - BEFORE \n");
  printf("Momenta 00 = ( %.18lf )\n", momenta[0].rc00[index]);
  printf("Momenta 11 = ( %.18lf )\n", momenta[0].rc11[index]);
  printf("Momenta 01 = ( %.18lf , %.18lf )\n", creal(momenta[0].c01[index]) , cimag(momenta[0].c01[index]));
  printf("Momenta 02 = ( %.18lf , %.18lf )\n", creal(momenta[0].c02[index]) , cimag(momenta[0].c02[index]));
  printf("Momenta 12 = ( %.18lf , %.18lf )\n", creal(momenta[0].c12[index]) , cimag(momenta[0].c12[index]));

  printf("Conf 00 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r0.c0[index]) , cimag(conf_acc[0].r0.c0[index]));
  printf("Conf 01 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r0.c1[index]) , cimag(conf_acc[0].r0.c1[index]));
  printf("Conf 02 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r0.c2[index]) , cimag(conf_acc[0].r0.c2[index]));
  printf("Conf 10 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r1.c0[index]) , cimag(conf_acc[0].r1.c0[index]));
  printf("Conf 11 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r1.c1[index]) , cimag(conf_acc[0].r1.c1[index]));
  printf("Conf 12 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r1.c2[index]) , cimag(conf_acc[0].r1.c2[index]));
  printf("Conf 20 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r2.c0[index]) , cimag(conf_acc[0].r2.c0[index]));
  printf("Conf 21 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r2.c1[index]) , cimag(conf_acc[0].r2.c1[index]));
  printf("Conf 22 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r2.c2[index]) , cimag(conf_acc[0].r2.c2[index]));


  double action_in;
  double action_fin;
  double action_mom_in;
  double action_mom_fin;
  double action_ferm_in;
  double action_ferm_fin;
  int mu;
  double temppp[8];

#pragma acc data copy(conf_acc[0:8]) copy(momenta[0:8]) create(aux_conf_acc[0:8]) copy(ferm_phi_acc[0:1]) copy(ferm_aux_acc[0:1]) copy(approx1[0:1]) copy(approx2[0:1]) copy(approx3[0:1])  copy(ferm_chi_acc[0:1])  create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1])  create(k_p_shiftferm[0:1]) create(ferm_shiftmulti_acc[0:1]) create(ipdot_acc[0:8]) copyin(delta[0:7])  copyin(nnp_openacc) copyin(nnm_openacc) create(local_sums[0:2]) create(d_local_sums[0:1])
  {

    gettimeofday ( &t1, NULL );

    // ora come ora la risuzione delle local sums viene fatta ogni volta.
    // eventualmente si puo' rendere piu efficiente l'algoritmo facendogli fare la riduzione solo alla fine
    //    --> la cosa puo' essere rilevante soprattutto per la versione salamica
    action_in = beta_by_three * calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
    action_mom_in = 0.0;
    for(mu =0;mu<8;mu++){
      //      temppp[mu] = calc_momenta_action(momenta,d_local_sums,mu);
      action_mom_in += calc_momenta_action(momenta,d_local_sums,mu);
    }
    action_ferm_in=scal_prod_between_multiferm(ferm_phi_acc,ferm_phi_acc);
    
    // first inv approx calc
    ker_invert_openacc_shiftmulti(conf_acc,ferm_shiftmulti_acc,ferm_phi_acc,residue_metro,approx1,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
    ker_openacc_recombine_shiftmulti_to_multi(ferm_shiftmulti_acc,ferm_phi_acc,ferm_chi_acc,approx1);

    multistep_2MN_SOLOOPENACC(ipdot_acc,conf_acc,aux_conf_acc,ferm_chi_acc,ferm_aux_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm,momenta,local_sums,delta,residue_md,approx2);



    action_fin = beta_by_three * calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
    action_mom_fin = 0.0;
    for(mu =0;mu<8;mu++){
      action_mom_fin += calc_momenta_action(momenta,d_local_sums,mu);
    }

    gettimeofday ( &t2, NULL );    
  } 

  gettimeofday ( &t3, NULL );

  printf("PLAQ IN  = %.18lf \n",action_in);
  printf("PLAQ OUT = %.18lf \n",action_fin);
  printf("MOM ACT IN  = %.18lf \n",action_mom_in);
  printf("MOM ACT OUT = %.18lf \n",action_mom_fin);
  printf("FERM ACT IN  = %.18lf \n",action_ferm_in);

  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);


  for(dir=0;dir<8;dir++)  convert_su3_soa_to_su3COM_soa(&conf_acc[dir],&conf[dir]);
  for(dir=0;dir<8;dir++)  convert_thmat_soa_to_thmatCOM_soa(&momenta[dir],&com_mom[dir]);

  printf("MOMENTA AND CONF - INSIDE OPENACC MD INTEGRATOR - AFTER \n");
  printf("Momenta 00 = ( %.18lf )\n", momenta[0].rc00[index]);
  printf("Momenta 11 = ( %.18lf )\n", momenta[0].rc11[index]);
  printf("Momenta 01 = ( %.18lf , %.18lf )\n", creal(momenta[0].c01[index]) , cimag(momenta[0].c01[index]));
  printf("Momenta 02 = ( %.18lf , %.18lf )\n", creal(momenta[0].c02[index]) , cimag(momenta[0].c02[index]));
  printf("Momenta 12 = ( %.18lf , %.18lf )\n", creal(momenta[0].c12[index]) , cimag(momenta[0].c12[index]));

  printf("Conf 00 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r0.c0[index]) , cimag(conf_acc[0].r0.c0[index]));
  printf("Conf 01 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r0.c1[index]) , cimag(conf_acc[0].r0.c1[index]));
  printf("Conf 02 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r0.c2[index]) , cimag(conf_acc[0].r0.c2[index]));
  printf("Conf 10 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r1.c0[index]) , cimag(conf_acc[0].r1.c0[index]));
  printf("Conf 11 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r1.c1[index]) , cimag(conf_acc[0].r1.c1[index]));
  printf("Conf 12 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r1.c2[index]) , cimag(conf_acc[0].r1.c2[index]));
  printf("Conf 20 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r2.c0[index]) , cimag(conf_acc[0].r2.c0[index]));
  printf("Conf 21 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r2.c1[index]) , cimag(conf_acc[0].r2.c1[index]));
  printf("Conf 22 = ( %.18lf , %.18lf )\n", creal(conf_acc[0].r2.c2[index]) , cimag(conf_acc[0].r2.c2[index]));

  printf("Ipdot 00 = ( %.18lf )\n", ipdot_acc[0].rc00[index]);
  printf("Ipdot 11 = ( %.18lf )\n", ipdot_acc[0].rc11[index]);
  printf("Ipdot 01 = ( %.18lf , %.18lf )\n", creal(ipdot_acc[0].c01[index]) , cimag(ipdot_acc[0].c01[index]));
  printf("Ipdot 02 = ( %.18lf , %.18lf )\n", creal(ipdot_acc[0].c02[index]) , cimag(ipdot_acc[0].c02[index]));
  printf("Ipdot 12 = ( %.18lf , %.18lf )\n", creal(ipdot_acc[0].c12[index]) , cimag(ipdot_acc[0].c12[index]));



  printf("FULL UPDATE COMPUTATION TIME                    Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);

  
  
  free(conf_acc);
  free(aux_conf_acc);
  free(ipdot_acc);
  free(ferm_phi_acc);
  free(ferm_chi_acc);
  free(ferm_aux_acc);
  free(ferm_shiftmulti_acc);
  free(kloc_r);
  free(kloc_s);
  free(kloc_h);
  free(kloc_p);
  free(k_p_shiftferm);
  
  
  
  
  
}
