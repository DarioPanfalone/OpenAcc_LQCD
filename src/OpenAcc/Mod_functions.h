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


#pragma acc routine seq
void mat_times_value(su3_soa * mat,int idx,double value){
    
    mat->r0.c0[idx] =value*mat->r0.c0[idx] ;
    mat->r0.c1[idx] =value*mat->r0.c1[idx] ;
    mat->r0.c2[idx] =value*mat->r0.c2[idx] ;
    mat->r1.c0[idx] =value*mat->r1.c0[idx] ;
    mat->r1.c1[idx] =value*mat->r1.c1[idx] ;
    mat->r1.c2[idx] =value*mat->r1.c2[idx] ;
    mat->r2.c0[idx] =value*mat->r2.c0[idx] ;
    mat->r2.c1[idx] =value*mat->r2.c1[idx] ;
    mat->r2.c2[idx] =value*mat->r2.c2[idx] ;
    
    return;
    
    
}

#endif /* Mod_functions_h */

