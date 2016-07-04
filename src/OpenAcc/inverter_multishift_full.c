#ifndef INVERTER_MULTISHIFT_FULL_C_
#define INVERTER_MULTISHIFT_FULL_C_

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "./fermion_matrix.h"
#include "./inverter_multishift_full.h"
#include "./inverter_full.h"

#include "../Mpi/multidev.h"

#define DEBUG_INVERTER_SHIFT_MULTI_FULL_OPENACC

extern int verbosity_lv;

int multishift_invert_iterations ; // global count of CG iterations

int multishift_invert(__restrict su3_soa * const u,
        __restrict ferm_param * pars,
        RationalApprox * approx,
        __restrict vec3_soa * out, // multi-fermion [nshifts]
        __restrict const vec3_soa * const in, // single ferm
        double residuo,
        __restrict vec3_soa * loc_r,
        __restrict vec3_soa * loc_h,
        __restrict vec3_soa * loc_s,
        __restrict vec3_soa * loc_p,
        __restrict vec3_soa * shiftferm, // multi-ferm [nshift]
        const int max_cg){
    /*********************
     * This function takes an input fermion 'in', a rational approximation
     * 'approx' and writes in 'out' a number of fermions, which are the
     * result of the inversions of all shifted matrices.
     * The result written in 'out' must then be summed to obtain the result of
     *
     * (M^\dagger M ) ^{fractional exponent} \psi 
     *
     * And this is done in the recombine_shifted_vec3_to_vec3_gl() 
     * function.
     ****************************/
    // AUXILIARY VARIABLES FOR THE INVERTER 
    if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF MULTISHIFT INVERTER\n");
    int  cg;
    //double *zeta_i,*zeta_ii,*zeta_iii,*omegas,*gammas;
    // DYNAMIC ALLOCATION OF THESE SMALL ARRAYS SEEMS TO FAIL.
    double zeta_i[MAX_APPROX_ORDER];
    double zeta_ii[MAX_APPROX_ORDER];
    double zeta_iii[MAX_APPROX_ORDER];
    double omegas[MAX_APPROX_ORDER];
    double gammas[MAX_APPROX_ORDER];
    int flag[MAX_APPROX_ORDER];

    //  printf("in pointer: %p\n", in);


    int iter;
    double alpha, delta, lambda, omega, omega_save, gammag, fact;
    alpha=0.0;
    // trial solution out = 0, set all flag to 1                                                                                                           
    for(iter=0; iter<(approx->approx_order); iter++){
        flag[iter]=1;
        set_vec3_soa_to_zero(&out[iter]);
    }

    // r=in, p=phi delta=(r,r)
    assign_in_to_out(in,loc_r);
    assign_in_to_out(loc_r,loc_p);
    delta=l2norm2_global(loc_r);
    double source_norm = delta;
    //printf("delta    %.18lf\n",delta);
    omega=1.0;

    //printf("Re in->c0[0]: %f\n"    ,creal(in->c0[0]) );
    //printf("Re loc_r->c0[0]: %f\n" ,creal(loc_r->c0[0]) );
    //printf("Re loc_p->c0[0]: %f\n" ,creal(loc_p->c0[0]) );

    for(iter=0; iter<(approx->approx_order); iter++){
        // ps_0=phi
        assign_in_to_out(in,&shiftferm[iter]);
        zeta_i[iter]=1.0;         // zeta_{-1}=1.0
        zeta_ii[iter]=1.0;        // zeta_{ 0}=1.0
        gammas[iter]=0.0;         // gammas_{-1}=0.0
    }
    gammag=0.0;
    cg=0;

    int maxiter; 
    //Find the max approx->approx_order where flag[approx->approx_order] == 1
    for(iter=0; iter<(approx->approx_order); iter++) {
        if(flag[iter]==1) {
            maxiter = iter+1;
        }
    }
    if (verbosity_lv > 0 && 0== devinfo.myrank ){
        printf("STARTING CG-M:\nCG\tR");
        for(iter=0; iter<(approx->approx_order); iter++)
            printf("\t%d",iter); printf("\n");
    } 
    do {      // loop over cg iterations
        cg++;

        // s=(M^dagM)p, alhpa=(p,s)=(p,Ap)
        //fermion_matrix_multiplication(su3_soa *u,vec3_soa *out,vec3_soa *in,vec3_soa *temp1, ferm_param *pars)
        fermion_matrix_multiplication(u,loc_s,loc_p,loc_h,pars);



        //      printf("component_loc_s[0]    %f\n",creal(loc_s->c0[0]));
        //      printf("component_loc_p[0]    %f\n",creal(loc_p->c0[0]));

        alpha = real_scal_prod_global(loc_p,loc_s);
        //      printf("alpha    %.18lf\n",alpha);
        //      printf("mass    %.18lf\n",pars->ferm_mass);

        omega_save=omega;   // omega_save=omega_(j-1)
        omega=-delta/alpha;  // omega = (r_j,r_j)/(p_j, Ap_j)               
        //      printf("omega    %.18lf\n",omega);

        // out-=omegas*ps
        for(iter=0; iter<maxiter; iter++){
            if(flag[iter]==1){
                zeta_iii[iter] = (zeta_i[iter]*zeta_ii[iter]*omega_save)/
                    ( omega*gammag*(zeta_i[iter]-zeta_ii[iter])+
                      zeta_i[iter]*omega_save*(1.0-(approx->RA_b[iter])*omega) );

                omegas[iter]=omega*zeta_iii[iter]/zeta_ii[iter];
            }
        }

        //multiple_combine_in1_minus_in2x_factor_back_into_in1(vec3_soa *out,vec3_soa *in,int maxiter,int *flag,double *omegas)
        multiple_combine_in1_minus_in2x_factor_back_into_in1(out,shiftferm,maxiter,flag,omegas);

        // r+=omega*s; lambda=(r,r)
        combine_add_factor_x_in2_to_in1(loc_r,loc_s,omega);
        lambda=l2norm2_global(loc_r);
        gammag=lambda/delta;

        // p=r+gammag*p
        combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);

        for(iter=0; iter<(approx->approx_order); iter++){
            if(flag[iter]==1){
                gammas[iter]=gammag*zeta_iii[iter]*omegas[iter]/(zeta_ii[iter]*omega);
            }
        }

        //multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1(vec3_soa *in1, int maxiter,int *flag,double *gammas,vec3_soa *in2,double *zeta_iii )
        multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1(shiftferm,maxiter,flag,gammas,loc_r,zeta_iii);

        maxiter = 0;
        for(iter=0; iter<(approx->approx_order); iter++){
            if(flag[iter]==1){
                fact=sqrt(delta*zeta_ii[iter]*zeta_ii[iter]/source_norm);
                if(fact<residuo) flag[iter]=0;
                else maxiter = iter+1;// modifying maxiter
                zeta_i[iter]=zeta_ii[iter];
                zeta_ii[iter]=zeta_iii[iter];
            }
        }
        delta=lambda;

        int print_condition = 0 ;
        if (verbosity_lv == 1  ) print_condition = (cg%500 == 0);
        if (verbosity_lv == 2  ) print_condition = (cg%200 == 0);
        if (verbosity_lv == 3  ) print_condition = (cg%100 == 0);
        if (verbosity_lv == 4  ) print_condition = (cg%50 == 0);
        if (verbosity_lv > 4  ) print_condition = (cg%10 == 0);


        if(print_condition && 0==devinfo.myrank  ){

            printf("%d\t%1.1e",cg, sqrt(lambda));
            for(iter=0; iter<(approx->approx_order); iter++){
                if(flag[iter]==0) printf("\t-");
                else printf("\t%1.1e", sqrt(delta*zeta_i[iter]*zeta_i[iter]) );
            }
            printf("\n");

        }
    } while(maxiter>0 && cg<max_cg); // end of cg iterations
    multishift_invert_iterations += cg ;  

    if(cg==max_cg && 0==devinfo.myrank )
    {
        printf("WARNING: maximum number of iterations reached in invert\n");
    }
    //      printf("\t CG count = %i \n",cg);


    if(verbosity_lv > 0 && 0==devinfo.myrank ) printf("Terminated multishift_invert ( target res = %1.1e,source_norm = %1.1e )\tCG count %d\n", residuo,source_norm,cg);
    // test 

    if(verbosity_lv > 2 && 0 == devinfo.myrank){
        printf("Shift:\t");
        for(iter=0; iter<approx->approx_order; iter++)printf("\t%d",iter);
        printf("\n");
    }

    printf("Relative Res:");
    for(iter=0; iter<approx->approx_order; iter++){
        if(verbosity_lv > 5 && 0 == devinfo.myrank)
            printf("Verifying result, shift %d\n", iter);
        assign_in_to_out(&out[iter],loc_p);
        fermion_matrix_multiplication_shifted(u,loc_s,loc_p,loc_h,pars,approx->RA_b[iter]);
        combine_in1_minus_in2(in,loc_s,loc_h); // r = s - y  
        double  giustoono=l2norm2_global(loc_h)/source_norm;

        if(verbosity_lv > 2 && 0 == devinfo.myrank){
            printf("\t%1.1e",sqrt(giustoono)/residuo);
        } 
        if(verbosity_lv > 5 && 0 == devinfo.myrank) printf("\n");
    }
    if(verbosity_lv > 2 && 0 == devinfo.myrank) printf("\n");

    return cg;

}

void recombine_shifted_vec3_to_vec3(const __restrict vec3_soa* const in_shifted /*multi-fermion*/, 
        const __restrict vec3_soa* const in, // [nshift]
        __restrict vec3_soa * const out, // [1] 
        const RationalApprox * const approx ){
    if(verbosity_lv > 3) printf("DOUBLE PRECISION VERSION OF RECOMBINE_SHIFTED_VEC3_TO_VEC3\n");
    int ih;
    int iter=0;
#pragma acc kernels present(out) present(in) present(in_shifted) present(approx)
#pragma acc loop independent
    for(ih=0; ih < sizeh; ih++){

        int ordine=approx->approx_order;
        out->c0[ih] =  in->c0[ih]*approx->RA_a0;
        out->c1[ih] =  in->c1[ih]*approx->RA_a0;
        out->c2[ih] =  in->c2[ih]*approx->RA_a0;

#pragma acc loop seq
        for(iter=0; iter<ordine; iter++){  // questo loop non lo vogliamo parallelizzare per forza ... forse puo andare bene cosi'
            out->c0[ih] +=  approx->RA_a[iter] * in_shifted[iter].c0[ih];
            out->c1[ih] +=  approx->RA_a[iter] * in_shifted[iter].c1[ih];
            out->c2[ih] +=  approx->RA_a[iter] * in_shifted[iter].c2[ih];
        }




    }
}



#endif

