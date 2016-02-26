#ifndef _MARKOWCHAIN_H_
#define _MARKOWCHAIN_H_

typedef struct mc_param_t{
    int ntraj                  ;
    int therm_ntraj            ;
    int storeconfinterval       ;
    int saveconfinterval;
    int use_ildg;
    double residue_metro;
    double expected_max_eigenvalue; // verbosity level
    char store_conf_name[200];
    char save_conf_name[200];
    int seed;
    double eps_gen;
    int input_vbl; // verbosity level
    int save_diagnostics;
    char diagnostics_filename[200];
}mc_param;

extern mc_param mkwch_pars;

#endif
