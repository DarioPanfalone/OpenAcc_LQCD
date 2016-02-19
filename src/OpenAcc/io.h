#ifndef IO_H_
#define IO_H_

#include "./struct_c_def.h"


void print_su3_soa_ASCII(su3_soa * const conf, const char* nomefile,int conf_id_iter);
int read_su3_soa_ASCII(su3_soa * conf, const char* nomefile,int * conf_id_iter );

void print_su3_soa_ildg_binary(su3_soa * const conf, const char* nomefile,int conf_id_iter);
int read_su3_soa_ildg_binary(su3_soa * conf, const char* nomefile,int * conf_id_iter );



#endif 
