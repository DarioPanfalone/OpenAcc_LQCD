#ifndef FERMION_FORCE_H
#define FERMION_FORCE_H


#include "./stouting.h"
#include "./struct_c_def.h"
#include "./inverter_multishift_full.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif



void compute_sigma_from_sigma_prime_backinto_sigma_prime(  __restrict su3_soa    * Sigma, // la var globale e' auxbis_conf_acc [sia input che ouptput]
							   __restrict thmat_soa  * Lambda, // la var globale e' aux_th
							   __restrict tamat_soa  * QA, // la var globale e' aux_ta
							   __restrict su3_soa    * const U, // la var globale e' .... per adesso conf_acc
							   __restrict su3_soa    * const TMP// la var globale e' aux_conf_acc //PARCHEGGIO??
							   );

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
			       )



#endif
