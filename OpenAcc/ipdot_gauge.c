#ifndef IPDOT_GAUGE_C
#define IPDOT_GAUGE_C

#include "./su3_utilities.c"
#include "./rettangoli.c"



void calc_ipdot_gauge_soloopenacc_std( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot){

#ifdef TIMING_ALL
  struct timeval t1,t2;
  gettimeofday ( &t1, NULL );
#endif

  set_su3_soa_to_zero(local_staples);
  calc_loc_staples_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);
  conf_times_staples_ta_part(tconf_acc,local_staples,tipdot);

#ifdef TIMING_ALL
  gettimeofday ( &t2, NULL );
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  printf("FULL STAPLES CALC OPENACC                       PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

}

void calc_ipdot_gauge_soloopenacc_tlsm( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot){

#ifdef TIMING_ALL
  struct timeval t1,t2;
  gettimeofday ( &t1, NULL );
#endif

  set_su3_soa_to_zero(local_staples);
  calc_loc_staples_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);

  // QUESTA CHE FA TUTTO IN UNA BOTTA SEMBRA ANDARE PIU' PIANO
  //  calc_loc_improved_staples_typeABC_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);

  calc_loc_improved_staples_typeA_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);
  calc_loc_improved_staples_typeB_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);
  calc_loc_improved_staples_typeC_removing_stag_phases_nnptrick_all(tconf_acc,local_staples);
  
  conf_times_staples_ta_part(tconf_acc,local_staples,tipdot);

#ifdef TIMING_ALL
  gettimeofday ( &t2, NULL );
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  printf("FULL STAPLES CALC OPENACC                       PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
#endif

}



// VERSATILE WRAPPER WHICH CHOOSES BETWEEN STD GAUGE ACTION OR TLSM GAUGE ACTION
void calc_ipdot_gauge_soloopenacc( __restrict  su3_soa * const tconf_acc,  __restrict su3_soa * const local_staples,__restrict tamat_soa * const tipdot){
  if(GAUGE_ACTION==0){
    calc_ipdot_gauge_soloopenacc_std(tconf_acc,local_staples,tipdot);
  }
  if(GAUGE_ACTION==1){
    calc_ipdot_gauge_soloopenacc_tlsm(tconf_acc,local_staples,tipdot);
  }


}


#endif
