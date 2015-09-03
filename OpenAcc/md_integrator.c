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
void multistep_2MN_gauge(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta,double * delta)
 {
 int md;
 // Step for the P
 // P' = P - l*dt*dS/dq
 // delta[3]=-cimag(ieps_acc)*scale*lambda;
 calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
 mom_sum_mult(tmomenta,tipdot,delta,3);
 for(md=1; md<gauge_scale; md++){
   // Step for the Q
   // Q' = exp[dt/2 *i P] Q
   // delta[4]=cimag(iepsh_acc)*scale;
   mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,delta,4);
   // Step for the P
   // P' = P - (1-2l)*dt*dS/dq
   // delta[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
   calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
   mom_sum_mult(tmomenta,tipdot,delta,5);
   // Step for the Q
   // Q' = exp[dt/2 *i P] Q
   // delta[4]=cimag(iepsh_acc)*scale;
   mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,delta,4);
   // Step for the P
   // P' = P - 2l*dt*dS/dq
   // delta[6]=-cimag(ieps_acc)*2.0*lambda*scale;
   calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
   mom_sum_mult(tmomenta,tipdot,delta,6);
 }
 
 // Step for the Q
 // Q' = exp[dt/2 *i P] Q
 // delta[4]=cimag(iepsh_acc)*scale;
 mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,delta,4);
 // Step for the P
 // P' = P - (1-2l)*dt*dS/dq
 calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
 // calc_ipdot_gauge();
 // delta[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
 mom_sum_mult(tmomenta,tipdot,delta,5);
 // Step for the Q
 // Q' = exp[dt/2 *i P] Q
 // delta[4]=cimag(iepsh_acc)*scale;
 mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,delta,4);
 // Step for the P
 // P' = P - l*dt*dS/dq
 // delta[3]=-cimag(ieps_acc)*lambda*scale;
 calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
 mom_sum_mult(tmomenta,tipdot,delta,3);
 
 }
void multistep_2MN_SOLOOPENACC( tamat_soa * tipdot_acc,
				su3_soa  * tconf_acc,
				double_soa * backfield,
				su3_soa  * taux_conf_acc,
				ferm_param * tfermions_parameters,// [nflavs]
				int tNDiffFlavs,
				vec3_soa * ferm_in_acc, //[NPS_tot], will be ferm_chi_acc
				//ACC_MultiFermion * ferm_in_acc,
				vec3_soa * tferm_shiftmulti_acc,// parking variable [max_ps*max_approx_order]
				//ACC_ShiftMultiFermion * ferm_shiftmulti_acc,
				vec3_soa * tkloc_r, // parking
				vec3_soa * tkloc_h, // parking
				vec3_soa * tkloc_s, // parking
				vec3_soa * tkloc_p, // parking
				vec3_soa * tk_p_shiftferm, // parking, [max_nshift]
				//ACC_ShiftFermion *k_p_shiftferm,
				thmat_soa * tmomenta,
				dcomplex_soa * local_sums,
				double * delta,
				//COM_RationalApprox *approx), // included in ferm_param
				double res)
{

 
  printf("######################################## \n");
  printf("## Inside multistep_2MN_SOLOOPENACC() ## \n");
  printf("######################################## \n");






  int md;
  
  // Step for the P
  // P' = P - l*dt*dS/dq
  //    delta[0]=-cimag(ieps_acc)*lambda;
  fermion_force_soloopenacc(tconf_acc, backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
  //  fermion_force_soloopenacc(conf_acc,ipdot_acc,ferm_in_acc,res,approx,aux_conf_acc,ferm_shiftmulti_acc,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm); // OLD

#pragma acc update host(tmomenta[0:8])
  print_thmat_soa(tmomenta,"momenta_before");

  mom_sum_mult(tmomenta,tipdot_acc,delta,0);

#pragma acc update host(tmomenta[0:8])
  print_thmat_soa(tmomenta,"momenta_after");

#pragma acc update host(tipdot_acc[0:8])
  print_tamat_soa(tipdot_acc,"tipdot");

  
  for(md=1; md<no_md; md++){
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta,delta);
    // Step for the P
    // P' = P - (1-2l)*dt*dS/dq
    // delta[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
  fermion_force_soloopenacc(tconf_acc, backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
    mom_sum_mult(tmomenta,tipdot_acc,delta,1);
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta,delta);
    // Step for the P
    // P' = P - 2l*dt*dS/dq
    // delta[2]=-cimag(ieps_acc)*(2.0*lambda);
  fermion_force_soloopenacc(tconf_acc, backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
    mom_sum_mult(tmomenta,tipdot_acc,delta,2);
  }  
  // Step for the Q
  // Q' = exp[dt/2 *i P] Q
  multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta,delta);
  // Step for the P
  // P' = P - (1-2l)*dt*dS/dq
  // delta[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
  fermion_force_soloopenacc(tconf_acc, backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
  mom_sum_mult(tmomenta,ipdot_acc,delta,1);
  // Step for the Q
  // Q' = exp[dt/2 *i P] Q
  multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta,delta);
  // Step for the P
  // P' = P - l*dt*dS/dq
  // delta[0]=-cimag(ieps_acc)*lambda;
  fermion_force_soloopenacc(tconf_acc, backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
  mom_sum_mult(tmomenta,tipdot_acc,delta,0);
    
}// end multistep_2MN_SOLOOPENACC()

