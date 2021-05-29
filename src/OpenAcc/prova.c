
#include "prova.h"
#include <complex.h>
#include <stdlib.h>

#include "./double_complex.h"
#ifdef __GNUC__
#define _POSIX_C_SOURCE 200809L // not to have warning on posix memalign
#endif

#include <stdio.h>
#include <stdlib.h>




#define sizeh 30
#define ALIGN 128
typedef struct double_soa_t {
    double d[sizeh];
} double_soa;


typedef struct vec3_soa_t {
    d_complex c0[sizeh];
    d_complex c1[sizeh];
    d_complex c2[sizeh];
} vec3_soa;


typedef struct su3_soa_t {
    vec3_soa r0;
    vec3_soa r1;
    vec3_soa r2;
    double_soa K; //Adjoint vector K. This vector is part of the struct and its values directly modify the near defect links computation in the action. Obviously its lenght is sizeh.
} su3_soa;



void init_k(su3_soa * conf){
    int mu,i;
    for(mu=0;mu<8;mu++){
        for(i=0; i<sizeh; i++ )
            conf[mu].K.d[i]=1;
    }
    
}




void init_vec3(vec3_soa vec3,int n){
    int i;
    for(i=0;i<sizeh;i++){
        vec3.c0[i]=n;
        vec3.c1[i]=n;
        vec3.c2[i]=n;
    }
    
    
}


/*
void print_vec3(){
    for(i=0)
    
}
*/



int main(){
    

    
   
    
    
    
    
    return 0;
}
