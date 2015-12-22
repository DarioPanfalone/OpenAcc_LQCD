#ifndef IO_C_
#define IO_C_

#include "./struct_c_def.h"

void print_su3_soa_ASCII(su3_soa * const conf, const char* nomefile,int conf_id_iter);
int read_su3_soa_ASCII(su3_soa * conf, const char* nomefile,int * conf_id_iter );


#endif 
