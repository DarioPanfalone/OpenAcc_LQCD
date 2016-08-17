#ifndef _MONTECARLO_PARAMS_H_
#define _MONTECARLO_PARAMS_H_

typedef struct MC_PARAM_T{
    int ntraj                  ;
    int MaxConfIdIter; 
    double MaxRunTimeS ;
    int therm_ntraj            ;
    int storeconfinterval       ;
    int saveconfinterval;
    char RandGenStatusFilename[200];
    char store_conf_name[200];
    char save_conf_name[200];
    int seed;
    double eps_gen;
    int JarzynskiMode;
} mc_params_t;

extern mc_params_t mc_params;

#endif
