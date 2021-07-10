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



int init_k_test(su3_soa *conf_acc,int c_r);
int n_replicas_reader(const char* input_filename);
void counter_size_function(int d0,int d1,int d2,int d3);
int init_k(su3_soa * conf,double c_r,int def_axis,int * def_vet);
void printing_k_mu(su3_soa * conf);

int replicas_swap_1(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet );
int replicas_swap(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet );
int label_print(su3_soa ** conf_hasen, int replicas_number,FILE *file,int step_number);


double calc_loc_plaquettes_rectangles_SWAP(
                                           __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                           __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                           dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                           const int mu, const int nu, int def_axis, int *def_vet);

double calc_loc_plaquettes_nnptrick_SWAP(  __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                         __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                         dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                         const int mu, const int nu, int def_axis, int *def_vet);
double  calc_plaquette_soloopenacc_SWAP( __restrict  su3_soa * const tconf_acc,
                                        __restrict su3_soa * const local_plaqs,
                                        dcomplex_soa * const tr_local_plaqs,int def_axis, int *def_vet, int improved );

#endif /* Mod_functions_h */

