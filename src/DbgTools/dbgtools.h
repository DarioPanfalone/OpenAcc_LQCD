#ifndef DBGTOOLS_H_
#define DBGTOOLS_H_

#include "../OpenAcc/struct_c_def.h"

extern int ipdot_f_reset; // flag to be set to 1 when md starts, 
                   // and to be set to zero when first force calculation happens
                   // if it is 1 force difference is not calculated.
extern int ipdot_g_reset; // same


// multi rank data structure read/write functions - for debug/testing
// multi rank wrappers

// gl3 / su3
void print_gl3_soa_wrapper(const su3_soa * gl3, const char* nomefile);
int read_gl3_soa_wrapper(su3_soa * gl3, const char* nomefile);

// fermions
void print_vec3_soa_wrapper(const vec3_soa * fermion, const char* nomefile);
int read_vec3_soa_wrapper(vec3_soa * fermion, const char* nomefile);

// tamat
void print_tamat_soa_wrapper(const tamat_soa * tamat, const char* nomefile);
int read_tamat_soa_wrapper(tamat_soa * tamat, const char* nomefile);

// thmat
void print_thmat_soa_wrapper(const thmat_soa * thmat, const char* nomefile);
int read_thmat_soa_wrapper(thmat_soa * thmat, const char* nomefile);

// dcomplex
void print_dcomplex_soa_wrapper(const dcomplex_soa *dcarr,
        const char* nomefile);
int read_dcomplex_soa_wrapper(dcomplex_soa * dcarr, const char* nomefile);


// double
void print_double_soa_wrapper(const double_soa * darr, const char* nomefile);
int read_double_soa_wrapper(double_soa * darr, const char* nomefile);


// single rank functions
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
