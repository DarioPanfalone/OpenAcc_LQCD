#ifndef INVERTER_MULTISHIFT_WRAPPER_C_
#define INVERTER_MULTISHIFT_WRAPPER_C_

#include "./inverter_multishift_wrapper.h"
#include "./inverter_package.h"
#include "./inverter_multishift_full.h"
#include "./sp_inverter_multishift_full.h"
#include "./inverter_full.h"
#include "./sp_inverter_full.h"
#include "./inverter_mixedp.h"
#include "../Include/inverter_tricks.h"
#include "./alloc_vars.h"
#include "./sp_alloc_vars.h"
#include "../Include/fermion_parameters.h"
#include "../RationalApprox/rationalapprox.h"
#include "./struct_c_def.h"




void inverter_multishift_wrapper(inverter_package ip,
        ferm_param *pars,
        RationalApprox * approx,
        vec3_soa * out,
        vec3_soa * in,
        double res,
        int max_cg)
{


    if(inverter_tricks.singlePInvAccelMultiInv){

        convert_double_to_float_vec3_soa(in,aux1_f);
        // multishift inverter in single precision
        float singlePMultiInvTargetRes = 8e-7f*sqrtf(sizeh);
        if(singlePMultiInvTargetRes < res)
          singlePMultiInvTargetRes = res;
        float mixedPInvTargetRes = 8e-10f*sqrtf(sizeh);
        if(mixedPInvTargetRes<res)
          mixedPInvTargetRes = res;

        multishift_invert_f(ip.u_f,pars,approx, 
                ferm_shiftmulti_acc_f,aux1_f, singlePMultiInvTargetRes, 
                ip.loc_r_f, ip.loc_h_f, ip.loc_s_f, ip.loc_p_f, ip.ferm_shift_temp_f, max_cg);
        int ishift;
        for(ishift=0;ishift<tfermion_parameters[iflav].approx_md.approx_order;
                ishift++){

            double bshift = tfermion_parameters[iflav].approx_md.RA_b[ishift];
            printf("Shift %d, %f\n", ishift,bshift);

            if(inverter_tricks.useMixedPrecision)
                inverter_mixed_precision(ip,pars,&ferm_shiftmulti_acc_f[ishift],
                        in,mixedPInvTargetRes,max_cg,bshift);


            convert_float_to_double_vec3_soa(&ferm_shiftmulti_acc_f[ishift],
                    &out[ishift]);// trial sol, hopefully close

            ker_invert_openacc(conf_to_use, pars,
                    &out[ishift],in,res,
                    ip.loc_r,ip.loc_h,ip.loc_s,ip.loc_p,max_cg,bshift);


        }


    }
    else multishift_invert(ip.u,pars,approx,out,in,res,
            ip.loc_r,ip.loc_h,ip.loc_s,ip.loc_p,ip.ferm_shift_temp_f,max_cg);


}


#endif


