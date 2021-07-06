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

#endif /* Mod_functions_h */

