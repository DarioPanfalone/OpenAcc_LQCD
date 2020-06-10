#ifndef MD_PARAMETERS_C
#define MD_PARAMETERS_C

#include "./md_parameters.h"
#include "./struct_c_def.h"
#include "./action.h"

double deltas_Omelyan[7];
float deltas_Omelyan_f[7];

void initialize_md_global_variables(md_param md_params )
{
    int no_md = md_params.no_md;// number of MD steps
    int gauge_scale = md_params.gauge_scale;   // Update fermions every gauge_scale gauge updates
    double t = md_params.t ;

    double epsilon_acc;
    d_complex ieps_acc,iepsh_acc;


    epsilon_acc = t/((double)(no_md));
    ieps_acc  = 0.0 + (epsilon_acc) * 1.0I;
    iepsh_acc = 0.0 + (epsilon_acc) * 0.5 * 1.0I;

    const double lambda=0.1931833275037836; // Omelyan Et Al.
    //const double lambda=0.1931833; // Omelyan Et Al.
    const double gs= 0.5/(double) gauge_scale;

    // used in the 'macro' step
    deltas_Omelyan[0]= - t/no_md * lambda;
    deltas_Omelyan[1]= - t/no_md * (1.0-2.0*lambda);
    deltas_Omelyan[2]= - t/no_md * 2.0*lambda;

    // used in the gauge substeps
    deltas_Omelyan[3]=   -t/(2*no_md*gauge_scale)*         lambda * BETA_BY_THREE;
    deltas_Omelyan[4]=0.5*t/(2*no_md*gauge_scale); // this is the only one which enters 
                                                   // directly the evolution of the links 
                                                   // and has the + sign,  and 
                                                   // no beta factor )
    deltas_Omelyan[5]=   -t/(2*no_md*gauge_scale)*(1.0-2.0*lambda)* BETA_BY_THREE;
    deltas_Omelyan[6]=   -t/(2*no_md*gauge_scale)*     2.0*lambda * BETA_BY_THREE;

    int iomelian;
    for(iomelian=0;iomelian<7;iomelian++)
     deltas_Omelyan_f[iomelian]= (float) deltas_Omelyan[iomelian];


#pragma acc enter data copyin(deltas_Omelyan[0:7])
#pragma acc enter data copyin(deltas_Omelyan_f[0:7])


}

int nMdInversionPerformed;// used to recycle inversion results
                              // after first force calculation are done
md_param md_parameters;

#endif
