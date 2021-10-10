//
//  Mod_functions.h
//  
//
//  Created by Luca Parente on 02/06/21.
//

#ifndef Mod_functions_h
#define Mod_functions_h

#include <stdio.h>
#include "Mod_functions.h"
#include "./struct_c_def.h"
#include "./alloc_settings.h"
#include "./alloc_vars.h"
#include "../Include/common_defines.h"
#include "./action.h"


typedef struct defect_info_t{

    int  defect_swap_min[4][4];
    int  defect_swap_max[4][4];
    #ifdef GAUGE_ACT_TLSM
    int  defect_swap_min_TLSM[4][4][2]; //cause mind that in symanzik action you have 2 different plaquette : 1x2 & 2x1
    int  defect_swap_max_TLSM[4][4][2];
    #endif
    int def_axis_mapped;
    int def_mapped_perp_dir[3];



}defect_info;
extern defect_info *def;
extern defect_info *pef;

int init_k_test(su3_soa *conf_acc,double c_r);
int n_replicas_reader(const char* input_filename);
void counter_size_function(int d0,int d1,int d2,int d3);
int init_k(su3_soa * conf,double c_r,int def_axis,int * def_vet,defect_info * def,int defect_info_config);
void printing_k_mu(su3_soa * conf);
void trasl_conf( __restrict const su3_soa *  const tconf_acc,
                 __restrict const su3_soa *  const taux_conf);

int replicas_swap_1(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet );
int replicas_swap(su3_soa * conf1,su3_soa * conf2);
int label_print(su3_soa ** conf_hasen, int replicas_number,FILE *file,int step_number);



double calc_Delta_S_Symanzik_SWAP(
                                           __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                           __restrict const su3_soa * const w,
                                           __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                           dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                           const int mu, const int nu, defect_info * def);

double calc_Delta_S_Wilson_SWAP(  __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                             __restrict const su3_soa * const w,
                                         __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                         dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                         const int mu, const int nu, defect_info * def);




double  calc_Delta_S_soloopenacc_SWAP( __restrict  su3_soa * const tconf_acc,
                                          __restrict  su3_soa * const tconf_acc2,
                                        __restrict su3_soa * const local_plaqs,
                                        dcomplex_soa * const tr_local_plaqs,defect_info * def, int improved );

int metro_SWAP(su3_soa ** conf_hasenbusch,
               __restrict su3_soa * const loc_plaq, //la placchetta locale.
               dcomplex_soa * const tr_local_plaqs,
               int rep_indx1, int rep_indx2,defect_info * def);

void All_Conf_SWAP(su3_soa ** conf_hasenbusch,
                   __restrict su3_soa * const loc_plaq, //la placchetta locale.
                   dcomplex_soa * const tr_local_plaqs,
                   
                   int replicas_number, defect_info * def, int *swap_num,int * all_swap_vet,int * acceptance_vet);


#endif /* Mod_functions_h */


