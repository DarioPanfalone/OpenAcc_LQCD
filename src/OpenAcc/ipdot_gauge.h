#ifndef IPDOT_GAUGE_H
#define IPDOT_GAUGE_H

#include "struct_c_def.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif


void calc_ipdot_gauge_soloopenacc_std( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);

void calc_ipdot_gauge_soloopenacc_tlsm( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);


// VERSATILE WRAPPER WHICH CHOOSES BETWEEN STD GAUGE ACTION OR TLSM GAUGE ACTION
void calc_ipdot_gauge_soloopenacc( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot);

#ifdef MULTIDEVICE
void calc_ipdot_gauge_soloopenacc_std_bulk( 
        __restrict  su3_soa * const tconf_acc, 
        __restrict su3_soa * const local_staples,
        __restrict tamat_soa * const tipdot);

void calc_ipdot_gauge_soloopenacc_tlsm_bulk( 
        __restrict  su3_soa * const tconf_acc,  
        __restrict su3_soa * const local_staples,
        __restrict tamat_soa * const tipdot);

void calc_ipdot_gauge_soloopenacc_bulk( 
        __restrict  su3_soa * const tconf_acc,  
        __restrict su3_soa * const local_staples,
        __restrict tamat_soa * const tipdot);

void calc_ipdot_gauge_soloopenacc_std_d3c( 
        __restrict  su3_soa * const tconf_acc, 
        __restrict su3_soa * const local_staples,
        __restrict tamat_soa * const tipdot,
        int offset3, int thickness3);

void calc_ipdot_gauge_soloopenacc_tlsm_d3c( 
        __restrict  su3_soa * const tconf_acc,  
        __restrict su3_soa * const local_staples,
        __restrict tamat_soa * const tipdot,
        int offset3, int thickness3);

void calc_ipdot_gauge_soloopenacc_d3c( 
        __restrict  su3_soa * const tconf_acc,  
        __restrict su3_soa * const local_staples,
        __restrict tamat_soa * const tipdot,
        int offset3, int thickness3);

#endif



#endif
