#ifndef INVERTER_MIXEDP_C_
#define INVERTER_MIXEDP_C_

#include "../Include/common_defines.h"
#include "../Include/inverter_tricks.h"
#include "./fermion_matrix.h"
#include "./fermionic_utilities.h"
#include "./float_double_conv.h"
#include "./sp_fermion_matrix.h"
#include "./sp_fermionic_utilities.h"
#include "./sp_struct_c_def.h"
#include "./struct_c_def.h"
#include "./inverter_package.h"

#ifndef __GNUC__
 #include "openacc.h"
 #define __restrict
#endif
#include "./sp_inverter_full.h"

#include "../Mpi/multidev.h"


void combine_add_in2_into_in1_mixed_precision(
        __restrict vec3_soa * in1,__restrict const  vec3_soa_f * in2){

#pragma acc kernels present(in1) present(in2)
#pragma acc loop independent
  for(int ih=LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih < LNH_SIZEH-LNH_VOL3/2*(D3_HALO-D3_FERMION_HALO); ih++) {
    in1->c0[ih] += (double)in2->c0[ih];
    in1->c1[ih] += (double)in2->c1[ih];
    in1->c2[ih] += (double)in2->c2[ih];
  }

}


// based on https://arxiv.org/pdf/0911.3191v2.pdf, pg.20


int inverter_mixed_precision(inverter_package ip,
			  ferm_param *pars,
			  __restrict vec3_soa * solution,// double precision output
			  __restrict const vec3_soa * in, // non viene aggiornato mai qui dentro
			  double res,
              const int  max_cg,
              double shift  )
{

  int cg;
  long int i;
  double delta, alpha, lambda, omega, gammag;
  double lastMaxResNorm=0;
  double source_norm = l2norm2_global(in);
  if (verbosity_lv > 3 && 0==devinfo.myrank )
  { 
      printf("source_norm: %e, ",source_norm);
      printf("res = %e, ",res);
      printf("target residue = %e\n", res*res*source_norm);
  }


  //// r = in  - M * solution
  // s = M * solution
  fermion_matrix_multiplication_shifted(ip.u,ip.loc_s,solution,ip.loc_h,pars,shift);
  // r = in - s
  combine_in1_minus_in2(in,ip.loc_s,ip.loc_r);
  
  // \hat{r} = (float)r
  convert_double_to_float_vec3_soa(ip.loc_r,ip.loc_r_f);
  // \hat{x} = 0 // 
  set_vec3_soa_to_zero_f(ip.out_f);

  // starting with CG, Single precision
  assign_in_to_out_f(ip.loc_r_f,ip.loc_p_f);
  delta=l2norm2_global_f(ip.loc_r_f);

  // loop over cg iterations
  cg=0;
    if (verbosity_lv > 3 && 0==devinfo.myrank )
      printf("STARTING CG:\nCG\tR\n");

  do {
    cg++;    
    // s=(M^dag M)p    alpha=(p,s)
    fermion_matrix_multiplication_shifted_f(ip.u_f,ip.loc_s_f,ip.loc_p_f,ip.loc_h_f,pars,shift);
    alpha = real_scal_prod_global_f(ip.loc_p_f,ip.loc_s_f);

    omega=delta/alpha;     
    // solution+=omega*p  r-=omega*s
    // lambda=(r,r);
    combine_in1xfactor_plus_in2_f(ip.loc_p_f,omega,ip.out_f,ip.out_f);
  
    // calculation of loc_r_f
    // EITHER THIS WAY 
    if(lastMaxResNorm < delta) lastMaxResNorm = delta;
    // the double precision magic touch
    if(delta < inverter_tricks.mixedPrecisionDelta * lastMaxResNorm){
        // adding the partial solution back to the full solution
        combine_add_in2_into_in1_mixed_precision(solution,ip.out_f);
        //// r = in  - M * solution
        // s = M * solution
        fermion_matrix_multiplication_shifted(ip.u,ip.loc_s,solution,ip.loc_h,pars,shift);
        // r = in - s
        combine_in1_minus_in2(in,ip.loc_s,ip.loc_r);

        // \hat{r} = (float)r
        //convert_double_to_float_vec3_soa(ip.loc_r,ip.loc_r_f);

        // \hat{x} = 0 // 
        set_vec3_soa_to_zero_f(ip.out_f);

        lastMaxResNorm=0;
        if(0==devinfo.myrank && verbosity_lv > 4){
            printf("Iteration %d, Inverter Mixed Precision: performed magic touch\n", cg);
            printf("Iteration %d, residue norm: %e (no mt)\n",cg,l2norm2_global_f(ip.loc_r_f));
            printf("Iteration %d, residue norm: %e (mt)\n",cg,l2norm2_global(ip.loc_r));
        }


    // OR THIS WAY
    }else combine_in1xfactor_plus_in2_f(ip.loc_s_f,-omega,ip.loc_r_f,ip.loc_r_f);

    
    lambda = l2norm2_global_f(ip.loc_r_f);
    gammag=lambda/delta;
    delta=lambda;
    if(0==devinfo.myrank && verbosity_lv > 4){
        printf("Iteration %d, residue norm (no magic touch): %e\n",cg,delta);
    }


    // p=r+gammag*p
    combine_in1xfactor_plus_in2_f(ip.loc_p_f,gammag,ip.loc_r_f,ip.loc_p_f);



    if (verbosity_lv > 3 && cg%100==0 && 0==devinfo.myrank  ){

        printf("%d\t%1.1e\n",cg, sqrt(lambda/source_norm)/res);fflush(stdout);
      
    }



  } while( (sqrt(lambda/source_norm)>res) && cg<max_cg);

  combine_add_in2_into_in1_mixed_precision(solution,ip.out_f);

  if (verbosity_lv > 3  && 0==devinfo.myrank ) printf("\n");
#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER_FULL_OPENACC))
  // s = M x
  fermion_matrix_multiplication_shifted(ip.u,ip.loc_s,solution,ip.loc_h,pars,shift);
  combine_in1_minus_in2(in,ip.loc_s,ip.loc_h); // r = s - y  
  double  giustoono=l2norm2_global(ip.loc_h)/source_norm;
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

