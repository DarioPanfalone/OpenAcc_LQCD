#ifndef STOUTING_H
#define STOUTING_H

#include "./cayley_hamilton.h"
#include "./struct_c_def.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/su3_utilities.h"


void exp_minus_QA_times_conf(__restrict su3_soa * const tu,
			     __restrict tamat_soa * const QA,
                 __restrict su3_soa * const tu_out,
			     __restrict su3_soa * const exp_aux);

inline void stout_isotropic(
 __restrict su3_soa * const u,               // --> input conf
 __restrict su3_soa * const uprime,          // --> output conf [stouted]
 __restrict su3_soa * const local_staples,   // --> parking variable
 __restrict su3_soa * const auxiliary,       // --> parking variable
 __restrict tamat_soa * const tipdot){       // --> parking variable

    SETREQUESTED(uprime);
    SETREQUESTED(local_staples);
    SETREQUESTED(auxiliary);
    SETREQUESTED(tipdot);


  set_su3_soa_to_zero(local_staples);

  mult_conf_times_stag_phases(u);

  calc_loc_staples_removing_stag_phases_nnptrick_all(u,local_staples);

  RHO_times_conf_times_staples_ta_part(u,local_staples,tipdot);
  SETFREE(local_staples);


  exp_minus_QA_times_conf(u,tipdot,uprime,auxiliary);
  SETFREE(tipdot);


  mult_conf_times_stag_phases(u);
  mult_conf_times_stag_phases(uprime);

}

void compute_lambda(__restrict thmat_soa * const L, // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
__restrict su3_soa   * const SP, // Sigma primo --> input (forza fermionica del passo precedente)
__restrict su3_soa   * const U,    // la configurazione di gauge --> input
__restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
__restrict su3_soa   * const TMP  // variabile di parcheggio
		    );

#ifdef STOUT_FERMIONS
inline void stout_wrapper(su3_soa * tconf_acc, su3_soa * tstout_conf_acc_arr){


  for(int mu = 0; mu < 8*STOUT_STEPS; mu ++) SETREQUESTED((&tstout_conf_acc_arr[mu]));
    stout_isotropic(tconf_acc, tstout_conf_acc_arr, auxbis_conf_acc, glocal_staples, gipdot );
    for(int stoutlevel=1;stoutlevel < STOUT_STEPS; stoutlevel++)
        stout_isotropic(&(tstout_conf_acc_arr[8*(stoutlevel-1)]),&(tstout_conf_acc_arr[8*stoutlevel]),auxbis_conf_acc, glocal_staples,  gipdot );

}
#endif

void compute_sigma(__restrict thmat_soa * const L,  // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
 __restrict su3_soa   * const U,  // la configurazione di gauge --> input
 __restrict su3_soa   * const S,  // entra Sigma primo (input: fermforce del passo precedente) ED esce Sigma --> sia input che ouput
 __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
 __restrict su3_soa   * const TMP // variabile di parcheggio
		   );


#endif
