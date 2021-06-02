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

/////funzione aggiunta//////////////////////
//per ora la metto qui poi la sposto in struct_def.h

void init_k(su3_soa * conf,int c_r){
    int mu,i;
    for(mu=0;mu<8;mu++){
        for(i=0; i<sizeh; i++ )
            conf[mu].K.d[i]=c_r;
    }
}


