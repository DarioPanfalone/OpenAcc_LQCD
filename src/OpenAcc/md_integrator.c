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

#ifndef MD_INTEGRATOR_C
#define MD_INTEGRATOR_C

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "./fermion_force.h"
#include "./md_integrator.h"
#include "./alloc_vars.h"
#include "./ipdot_gauge.h"
#include "./su3_utilities.h"
#include "../Include/common_defines.h"
#include "../Include/fermion_parameters.h"

int no_md;// number of MD steps
int gauge_scale;   // Update fermions every gauge_scale gauge updates

int no_md_acc,gauge_scale_acc;
double epsilon_acc;
d_complex ieps_acc,iepsh_acc;
double deltas_Omelyan[7];


void initialize_md_global_variables(void )
{
  gauge_scale = 5;
  no_md = 11;

  epsilon_acc = 1.0/((double)(no_md));
  ieps_acc  = 0.0 + (epsilon_acc) * 1.0I;
  iepsh_acc = 0.0 + (epsilon_acc) * 0.5 * 1.0I;

  const double lambda=0.1931833275037836; // Omelyan Et Al.
  const double gs=0.5/(double) gauge_scale;

  deltas_Omelyan[0]= -cimag(ieps_acc) * lambda;
  deltas_Omelyan[1]= -cimag(ieps_acc) * (1.0-2.0*lambda);
  deltas_Omelyan[2]= -cimag(ieps_acc) * 2.0*lambda;
  deltas_Omelyan[3]= -cimag(ieps_acc) * gs*lambda * beta_by_three;
  deltas_Omelyan[4]=  cimag(iepsh_acc)* gs;
  deltas_Omelyan[5]= -cimag(ieps_acc) * gs*(1.0-2.0*lambda)*beta_by_three;
  deltas_Omelyan[6]= -cimag(ieps_acc) * gs*2.0*lambda*beta_by_three;

}

void multistep_2MN_gauge(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta)
 {
 int md;
 // Step for the P
 // P' = P - l*dt*dS/dq
 // deltas_Omelyan[3]=-cimag(ieps_acc)*scale*lambda;

 mult_conf_times_stag_phases(tconf_acc);

 calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
 mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);
 for(md=1; md<gauge_scale; md++){
   if(verbosity_lv > 2) printf("Gauge step %d of %d...\n",md,gauge_scale);
   // Step for the Q
   // Q' = exp[dt/2 *i P] Q
   // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
   mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
   // Step for the P
   // P' = P - (1-2l)*dt*dS/dq
   // deltas_Omelyan[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
   calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
   mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,5);
   // Step for the Q
   // Q' = exp[dt/2 *i P] Q
   // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
   mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
   // Step for the P
   // P' = P - 2l*dt*dS/dq
   // deltas_Omelyan[6]=-cimag(ieps_acc)*2.0*lambda*scale;
   calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
   mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,6);
 }
 
 // Step for the Q
 // Q' = exp[dt/2 *i P] Q
 // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
 mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
 // Step for the P
 // P' = P - (1-2l)*dt*dS/dq
 calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
 // calc_ipdot_gauge();
 // deltas_Omelyan[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
 mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,5);
 // Step for the Q
 // Q' = exp[dt/2 *i P] Q
 // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
 mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
 // Step for the P
 // P' = P - l*dt*dS/dq
 // deltas_Omelyan[3]=-cimag(ieps_acc)*lambda*scale;
 calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
 mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);


 mult_conf_times_stag_phases(tconf_acc);
 
 }
void multistep_2MN_SOLOOPENACC( tamat_soa * tipdot_acc,
				su3_soa  * tconf_acc,
#ifdef STOUT_FERMIONS
                su3_soa  * tstout_conf_acc_arr, // huge parking for stouting
                su3_soa  * tauxbis_conf_acc, 
#endif
				double_soa * backfield,
				su3_soa  * taux_conf_acc,
				ferm_param * tfermions_parameters,// [nflavs]
				int tNDiffFlavs,
				vec3_soa * ferm_in_acc, //[NPS_tot], will be ferm_chi_acc
				vec3_soa * tferm_shiftmulti_acc,// parking variable [max_ps*max_approx_order]
				vec3_soa * tkloc_r, // parking
				vec3_soa * tkloc_h, // parking
				vec3_soa * tkloc_s, // parking
				vec3_soa * tkloc_p, // parking
				vec3_soa * tk_p_shiftferm, // parking, [max_nshift]
				thmat_soa * tmomenta,
				dcomplex_soa * local_sums,
				double res)
{

 
  int md;
  
  // Step for the P
  // P' = P - l*dt*dS/dq
  //    deltas_Omelyan[0]=-cimag(ieps_acc)*lambda;
  //  DEOTT_fermion_force_soloopenacc(tconf_acc, 
  fermion_force_soloopenacc(tconf_acc, 
#ifdef STOUT_FERMIONS
          tstout_conf_acc_arr, tauxbis_conf_acc, // parkeggio
#endif
          backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, 
          ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r,
	  tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);

  mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,0);
  
  for(md=1; md<no_md; md++){
      printf("\n\n\t\tRUNNING MD STEP %d OF %d...\n", md, no_md);
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
    // Step for the P
    // P' = P - (1-2l)*dt*dS/dq
    // deltas_Omelyan[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
    //    DEOTT_fermion_force_soloopenacc(tconf_acc, 
    fermion_force_soloopenacc(tconf_acc, 
#ifdef STOUT_FERMIONS
          tstout_conf_acc_arr, tauxbis_conf_acc, // parkeggio
#endif
          backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs,
          ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc,
          tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
    mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,1);
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
    // Step for the P
    // P' = P - 2l*dt*dS/dq
    // deltas_Omelyan[2]=-cimag(ieps_acc)*(2.0*lambda);
    //    DEOTT_fermion_force_soloopenacc(tconf_acc,
    fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
          tstout_conf_acc_arr, tauxbis_conf_acc, // parkeggio
#endif
          backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
    mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,2);
  }  
  // Step for the Q
  // Q' = exp[dt/2 *i P] Q
  multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
  // Step for the P
  // P' = P - (1-2l)*dt*dS/dq
  // deltas_Omelyan[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
  //  DEOTT_fermion_force_soloopenacc(tconf_acc,
  fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
          tstout_conf_acc_arr, tauxbis_conf_acc, // parkeggio
#endif
          backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
  mom_sum_mult(tmomenta,ipdot_acc,deltas_Omelyan,1);
  // Step for the Q
  // Q' = exp[dt/2 *i P] Q
  multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
  // Step for the P
  // P' = P - l*dt*dS/dq
  // deltas_Omelyan[0]=-cimag(ieps_acc)*lambda;
  //DEOTT_fermion_force_soloopenacc(tconf_acc,
  fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
          tstout_conf_acc_arr, tauxbis_conf_acc, // parkeggio
#endif
          backfield, tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
  mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,0);
    
}// end multistep_2MN_SOLOOPENACC()

#endif

