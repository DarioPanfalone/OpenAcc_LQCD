#ifndef IPDOT_GAUGE_H
#define IPDOT_GAUGE_H

#include "./su3_utilities.h"
#include "./rettangoli.h"

void calc_ipdot_gauge_soloopenacc_std( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);

void calc_ipdot_gauge_soloopenacc_tlsm( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);


// VERSATILE WRAPPER WHICH CHOOSES BETWEEN STD GAUGE ACTION OR TLSM GAUGE ACTION
inline void calc_ipdot_gauge_soloopenacc( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot){
  if(GAUGE_ACTION==0){
    calc_ipdot_gauge_soloopenacc_std(tconf_acc,local_staples,tipdot);
  }
  if(GAUGE_ACTION==1){
    calc_ipdot_gauge_soloopenacc_tlsm(tconf_acc,local_staples,tipdot);
  }


}


#endif
