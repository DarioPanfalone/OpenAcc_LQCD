//
//  Mod_functions.c
//  
//
//  Created by Luca Parente on 02/06/21.
//

#include "Mod_functions.h"
#include "./struct_c_def.h"
#include "./alloc_settings.h"
#include "./alloc_vars.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/////funzione aggiunta//////////////////////
//per ora la metto qui poi la sposto in struct_def.h

void init_k(su3_soa * conf,int c_r){
    int mu,i;
    for(mu=0;mu<8;mu++){
        for(i=0; i<sizeh; i++ )
            conf[mu].K.d[i]=c_r;
    }
}


int init_k_test(su3_soa *conf_acc){
    int kk2=0;
    int mu1=0;
    for(mu1=0;mu1<8;mu1++){
        for(kk2=0;kk2<sizeh;kk2++){
            if(conf_acc[mu1].K.d[kk2]!=(1 && c_r)){
                return 1;
            }
        printf("%d:ku[%d]:%f\n",mu1,kk2, conf_acc[mu1].K.d[kk2]);
        }
    }
    return 0;
}


