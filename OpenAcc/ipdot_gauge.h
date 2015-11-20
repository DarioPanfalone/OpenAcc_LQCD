#ifndef IPDOT_GAUGE_H
#define IPDOT_GAUGE_H

#include "./su3_utilities.h"
#include "./rettangoli.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif


void calc_ipdot_gauge_soloopenacc_std( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);

void calc_ipdot_gauge_soloopenacc_tlsm( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);


// VERSATILE WRAPPER WHICH CHOOSES BETWEEN STD GAUGE ACTION OR TLSM GAUGE ACTION
void calc_ipdot_gauge_soloopenacc( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);


#endif
