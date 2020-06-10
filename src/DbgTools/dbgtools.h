#ifndef DBGTOOLS_H_
#define DBGTOOLS_H_

#include "../OpenAcc/struct_c_def.h"

extern int ipdot_f_reset; // flag to be set to 1 when md starts, 
                   // and to be set to zero when first force calculation happens
                   // if it is 1 force difference is not calculated.
extern int ipdot_g_reset; // same

void print_vec3_soa_wrapper(vec3_soa * const fermion,
        const char* nomefile);
int read_vec3_soa_wrapper(vec3_soa * fermion, 
        const char* nomefile);

void print_vec3_soa(vec3_soa * const fermion, const char* nomefile);
int read_vec3_soa(vec3_soa * fermion, const char* nomefile);
void print_su3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter);
void print_1su3_soa(su3_soa * const conf, const char* nomefile);
void read_su3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter );
void print_tamat_soa(tamat_soa * const ipdot, const char* nomefile);
int read_tamat_soa(tamat_soa * ipdot, const char* nomefile);
void print_thmat_soa(thmat_soa * const ipdot, const char* nomefile);
void print_1thmat_soa(thmat_soa * const ipdot, const char* nomefile);
int read_thmat_soa(thmat_soa * ipdot, const char* nomefile);
void print_double_soa(double_soa * const backfield, const char* nomefile);
void print_1double_soa(double_soa * const vettore, const char* nomefile);
void dbg_print_su3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter);
int dbgread_su3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter );
void calc_loc_abelian_plaquettes(const double_soa* phases,
        double_soa * loc_abelian_plaquettes, const int mu, const int nu );

void print_all_abelian_plaquettes(const double_soa* phases, const char * filename);


void dbgprint_gl3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter);
int dbgread_gl3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter );



#endif
