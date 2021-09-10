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



int init_k_test(su3_soa *conf_acc,double c_r);
int n_replicas_reader(const char* input_filename);
void counter_size_function(int d0,int d1,int d2,int d3);
int init_k(su3_soa * conf,double c_r,int def_axis,int * def_vet);
void printing_k_mu(su3_soa * conf);
void trasl_conf(su3_soa * tconf_acc,su3_soa * taux_conf);

int replicas_swap_1(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet );
int replicas_swap(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet );
int label_print(su3_soa ** conf_hasen, int replicas_number,FILE *file,int step_number);



double calc_Delta_S_Symanzik_SWAP(
                                           __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                           __restrict const su3_soa * const w,
                                           __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                           dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                           const int mu, const int nu, int def_axis, int *def_vet);

double calc_Delta_S_Wilson_SWAP(  __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                             __restrict const su3_soa * const w,
                                         __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                         dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                         const int mu, const int nu, int def_axis, int *def_vet);




double  calc_Delta_soloopenacc_SWAP( __restrict  su3_soa * const tconf_acc,
                                          __restrict  su3_soa * const tconf_acc2,
                                        __restrict su3_soa * const local_plaqs,
                                        dcomplex_soa * const tr_local_plaqs,int def_axis, int *def_vet, int improved );

int metro_SWAP(su3_soa ** conf_hasenbusch,
               __restrict su3_soa * const loc_plaq, //la placchetta locale.
               dcomplex_soa * const tr_local_plaqs,
               int rep_indx1, int rep_indx2,int defect_axis,int * defect_coordinates);

void All_Conf_SWAP(su3_soa ** conf_hasenbusch,
                   __restrict su3_soa * const loc_plaq, //la placchetta locale.
                   dcomplex_soa * const tr_local_plaqs,
                   
                   int replicas_number, int defect_axis, int * defect_coordinates, int *swap_num,int * all_swap_vet,int * acceptance_vet);


#endif /* Mod_functions_h */


