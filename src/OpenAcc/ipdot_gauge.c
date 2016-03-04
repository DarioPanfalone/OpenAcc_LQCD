#ifndef IPDOT_GAUGE_C
#define IPDOT_GAUGE_C

#include "struct_c_def.h"
#include "./su3_utilities.h"
#include "./su3_measurements.h"
#include "./plaquettes.h"
#include "./rettangoli.h"
#include "./ipdot_gauge.h"
#include "../Include/markowchain.h"

#include "sys/time.h"
//#define TIMING_STAPLES

extern int verbosity_lv;

extern tamat_soa * ipdot_g_old; // see alloc_vars.c

// FOR ASYNC TRANSFERS-MULTIDEVICE: SPLIT BORDERS-BULK
void calc_ipdot_gauge_soloopenacc_std( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot){

#ifdef TIMING_STAPLES
    struct timeval t1,t2;
    gettimeofday ( &t1, NULL );
#endif

    set_su3_soa_to_zero(local_staples);
    calc_loc_staples_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);
    conf_times_staples_ta_part(tconf_acc,local_staples,tipdot);

#ifdef TIMING_STAPLES
    gettimeofday ( &t2, NULL );
    double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("FULL STAPLES CALC OPENACC                       PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

}

// FOR ASYNC TRANSFERS-MULTIDEVICE: SPLIT BORDERS-BULK
void calc_ipdot_gauge_soloopenacc_tlsm( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot){

#ifdef TIMING_STAPLES
    struct timeval t1,t2;
    gettimeofday ( &t1, NULL );
#endif

    set_su3_soa_to_zero(local_staples);
    calc_loc_staples_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);

    // QUESTA CHE FA TUTTO IN UNA BOTTA SEMBRA ANDARE PIU' PIANO
    //    calc_loc_improved_staples_typeABC_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);

    calc_loc_improved_staples_typeA_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);
    calc_loc_improved_staples_typeB_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);
    calc_loc_improved_staples_typeC_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);

    conf_times_staples_ta_part(tconf_acc,local_staples,tipdot);

#ifdef TIMING_STAPLES
    gettimeofday ( &t2, NULL );
    double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
    printf("FULL STAPLES CALC OPENACC                       PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

}

// FOR ASYNC TRANSFERS-MULTIDEVICE: SPLIT BORDERS-BULK
void calc_ipdot_gauge_soloopenacc( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot){
    if(GAUGE_ACTION==0){
        calc_ipdot_gauge_soloopenacc_std(tconf_acc,local_staples,tipdot);
    }
    if(GAUGE_ACTION==1){
        calc_ipdot_gauge_soloopenacc_tlsm(tconf_acc,local_staples,tipdot);
    }

    if(mkwch_pars.save_diagnostics == 1){
        double  force_norm, diff_force_norm;
        force_norm = calc_force_norm(tipdot);
        diff_force_norm = calc_diff_force_norm(tipdot,ipdot_g_old);
        copy_ipdot_into_old(tipdot,ipdot_g_old);

        FILE *foutfile = 
            fopen(mkwch_pars.diagnostics_filename,"at");
        fprintf(foutfile,"GFHN %e \tDGFHN %e \t",force_norm,diff_force_norm);
        fclose(foutfile);



        if(verbosity_lv > 1)
            printf("\t\t\tGauge Force Half Norm: %e, Diff with previous: %e \n", 
                    force_norm, diff_force_norm);
    } 




}


#endif
