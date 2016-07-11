#ifndef INVERTER_WRAPPERS_C_
#define INVERTER_WRAPPERS_C_

#include "../Include/fermion_parameters.h"
#include "../Include/inverter_tricks.h"
#include "../RationalApprox/rationalapprox.h"
#include "./alloc_vars.h"
#include "./geometry.h"
#include "./inverter_full.h"
#include "./inverter_mixedp.h"
#include "./inverter_multishift_full.h"
#include "./inverter_package.h"
#include "./inverter_wrappers.h"
#include "./sp_alloc_vars.h"
#include "./sp_inverter_full.h"
#include "./sp_inverter_multishift_full.h"
#include "./struct_c_def.h"
#include "./float_double_conv.h"



int multishift_invert_iterations ; // global count of CG iterations

int inverter_multishift_wrapper(inverter_package ip,
        ferm_param *pars,
        RationalApprox * approx,
        vec3_soa * out,
        const vec3_soa * in,
        double res,
        int max_cg)
{
    int total_iterations = 0;

    if(inverter_tricks.singlePInvAccelMultiInv){
        printf("Inverter multishift wrapper: using SINGLE precision multi shift inverter\n");

        convert_double_to_float_vec3_soa(in,aux1_f);
        // multishift inverter in single precision
        float singlePMultiInvTargetRes = 8e-7f*sqrtf(sizeh);
        if(singlePMultiInvTargetRes < res)
          singlePMultiInvTargetRes = res;
        //fferm_shiftmulti_acc_f,aux1_f: global variables
        total_iterations += multishift_invert_f(ip.u_f,pars,approx, 
                ferm_shiftmulti_acc_f,aux1_f, singlePMultiInvTargetRes, 
                ip.loc_r_f, ip.loc_h_f, ip.loc_s_f, ip.loc_p_f, ip.ferm_shift_temp_f, max_cg);
        int ishift;
        for(ishift=0;ishift< approx->approx_order;
                ishift++){

            double bshift = approx->RA_b[ishift];
            printf("Shift %d, %f\n", ishift,bshift);

            // this step may be inefficient (but it is needed)
            convert_float_to_double_vec3_soa(&(ferm_shiftmulti_acc_f[ishift]),
                    &(out[ishift]));

            total_iterations += inverter_wrapper(ip,pars,&(out[ishift]),
                        in,res,max_cg,bshift);

        }

    }
    else{ 
    
        printf("Inverter multishift wrapper: using DOUBLE precision multi shift inverter\n");

        total_iterations += multishift_invert(ip.u,pars,approx,out,in,res,
                ip.loc_r,ip.loc_h,ip.loc_s,ip.loc_p,ip.ferm_shift_temp,max_cg);

    }
    return total_iterations;

}

int inverter_wrapper(inverter_package ip,
        ferm_param *pars,
        vec3_soa * out,
        const vec3_soa * in,
        double res,
        int max_cg,
        double shift){

    int total_iterations = 0;
    float mixedPInvTargetRes = 8e-10f*sqrtf(sizeh);
    if(mixedPInvTargetRes<res)
        mixedPInvTargetRes = res;

    vec3_soa_f * out_f = ip.out_f;

    if(inverter_tricks.useMixedPrecision){
        convert_double_to_float_vec3_soa(out,out_f);// trial sol, hopefully close
        total_iterations += 
            inverter_mixed_precision(ip,pars,out_f,
                    in,mixedPInvTargetRes,max_cg,shift);

        convert_float_to_double_vec3_soa(out_f,out);// trial sol, hopefully close
    }
    total_iterations += ker_invert_openacc(ip.u, pars,
            out,in,res,ip.loc_r,ip.loc_h,ip.loc_s,ip.loc_p,max_cg,shift);


   return total_iterations; 
}



#endif


