#ifndef IO_H_
#define IO_H_

#include "./struct_c_def.h"

void print_su3_soa_ASCII(su3_soa * conf, const char* nomefile,int conf_id_iter);
void print_su3_soa_ildg_binary(su3_soa * conf, const char* nomefile,int conf_id_iter);



//WRAPPER
inline void save_conf(su3_soa * const conf, const char* nomefile,int conf_id_iter,int use_ildg){

    if(use_ildg){
        printf("Using ILDG format.\n");
        print_su3_soa_ildg_binary(conf,nomefile,conf_id_iter);
    }
    else 
        print_su3_soa_ASCII(conf,nomefile,conf_id_iter);

}

int read_su3_soa_ASCII(su3_soa * conf, const char* nomefile,int * conf_id_iter );
int read_su3_soa_ildg_binary(su3_soa * conf, const char* nomefile,int * conf_id_iter );

//WRAPPER
inline int read_conf(su3_soa * conf, const char* nomefile,int * conf_id_iter, int use_ildg ){


    printf("Reading configuration. ");
    if(use_ildg){ 
        printf("Using ILDG format.\n");
        return read_su3_soa_ildg_binary(conf,nomefile,conf_id_iter);
    }
    else {
        printf("Using ascii format.\n");
        return read_su3_soa_ASCII(conf,nomefile,conf_id_iter );
    }


}

#endif 
