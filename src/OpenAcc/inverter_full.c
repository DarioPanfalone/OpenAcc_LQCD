
#ifndef INVERTER_FULL_C_
#define INVERTER_FULL_C_

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#include "./inverter_full.h"
#include "../Include/inverter_tricks.h"
#include "../Mpi/multidev.h"


#ifndef __GNUC__
#include "openacc.h"
#endif

#define DEBUG_INVERTER_FULL_OPENACC


int ker_invert_openacc(__restrict const su3_soa * u, // non viene aggiornata mai qui dentro
        ferm_param *pars,
        __restrict vec3_soa * solution,
        __restrict const vec3_soa * in, // non viene aggiornato mai qui dentro
        double res,
        __restrict vec3_soa * loc_r,
        __restrict vec3_soa * loc_h,
        __restrict vec3_soa * loc_s,
        __restrict vec3_soa * loc_p,
        const int  max_cg,
        double shift  )
{

    int cg = 0;
    long int i;
    double delta, alpha, lambda, omega, gammag;

    double source_norm = l2norm2_global(in);

    do {

        fermion_matrix_multiplication_shifted(u,loc_s,solution,loc_h,pars,shift);

        combine_in1_minus_in2(in,loc_s,loc_r);
        assign_in_to_out(loc_r,loc_p);

        delta=l2norm2_global(loc_r);

        // loop over cg iterations
        if (verbosity_lv > 3 && 0==devinfo.myrank )
            printf("STARTING CG:\nCG\tR\n");

        int cg_restarted = 0 ;
        do {
            cg++;    
            cg_restarted++;
            // s=(M^dag M)p    alpha=(p,s)

            fermion_matrix_multiplication_shifted(u,loc_s,loc_p,loc_h,pars,shift);
            alpha = real_scal_prod_global(loc_p,loc_s);

            omega=delta/alpha;     
            // solution+=omega*p  r-=omega*s
            // lambda=(r,r);

            combine_in1xfactor_plus_in2(loc_p,omega,solution,solution);
            combine_in1xfactor_plus_in2(loc_s,-omega,loc_r,loc_r);

            lambda = l2norm2_global(loc_r);
            gammag=lambda/delta;
            delta=lambda;
            // p=r+gammag*p
            combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);

            if (verbosity_lv > 3 && cg%100==0 && 0==devinfo.myrank  ){

                printf("%d\t%1.1e\n",cg, sqrt(lambda/source_norm)/res);fflush(stdout);

            }

        } while( (sqrt(lambda/source_norm)>res) && 
                cg_restarted<inverter_tricks.restartingEvery);

    } while( (sqrt(lambda/source_norm)>res) && cg<max_cg);


    if (verbosity_lv > 3  && 0==devinfo.myrank ) printf("\n");
#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER_FULL_OPENACC))

    fermion_matrix_multiplication_shifted(u,loc_s,solution,loc_h,pars,shift);
    combine_in1_minus_in2(in,loc_s,loc_h); // r = s - y  
    double  giustoono=l2norm2_global(loc_h)/source_norm;
    if(verbosity_lv > 1 && 0==devinfo.myrank  ){
        printf("Terminated invert after   %d    iterations", cg);
        printf("[res/stop_res=  %e , stop_res=%e ]\n",
                sqrt(giustoono)/res,res);
    }
#endif
    if(cg==max_cg  && 0==devinfo.myrank )
    {
        printf("WARNING: maximum number of iterations reached in invert\n");
    }

    return cg;

}


#endif

