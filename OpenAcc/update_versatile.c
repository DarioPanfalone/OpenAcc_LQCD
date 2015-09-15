// WRAPPER TO VARIOUS UPDATES 
#ifndef UPDATE_VERSATILE_C
#define UPDATE_VERSATILE_C

#include "./update_standard_action.c"
#include "./update_tlsm_stdferm.c"


int UPDATE_SOLOACC_UNOSTEP_VERSATILE(su3_soa *tconf_acc,double res_metro, double res_md, int id_iter,int acc,int metro){
  int output;
  if(GAUGE_ACTION == 0){
    output = UPDATE_SOLOACC_UNOSTEP_VERSATILE_STDGAUGE_STDFERM(tconf_acc,res_metro,res_md,id_iter,acc,metro);
  }
  if(GAUGE_ACTION == 1){
    output = UPDATE_SOLOACC_UNOSTEP_VERSATILE_TLSM_STDFERM(tconf_acc,res_metro,res_md,id_iter,acc,metro);
  }
  return output;
}


#endif
