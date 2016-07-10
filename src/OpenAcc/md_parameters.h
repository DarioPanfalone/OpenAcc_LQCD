#ifndef MD_PARAMETERS_H
#define MD_PARAMETERS_H


typedef struct md_param_t{

    int no_md;// number of MD steps
    int gauge_scale;   // Update fermions every gauge_scale gauge updates
    double t ;
    double residue_metro;
    double expected_max_eigenvalue;
    int singlePrecMD ; 
    double residue_md;
    int max_cg_iterations;
    int recycleInvsForce;

} md_param; 

extern double deltas_Omelyan[7];// must be declared here to copy it in the device in the main
extern float deltas_Omelyan_f[7];// must be declared here to copy it in the device in the main

void initialize_md_global_variables(md_param);


extern md_param md_parameters;

#endif
