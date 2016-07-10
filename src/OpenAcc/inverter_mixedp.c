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
#endif
#include "./sp_inverter_full.h"

#include "../Mpi/multidev.h"


int inverter_mixed_precision(inverter_package ip,
			  ferm_param *pars,
			  __restrict vec3_soa_f * solution,// single precision output
			  __restrict const vec3_soa * in, // non viene aggiornato mai qui dentro
			  double res,
              const int  max_cg,
              double shift  )
{

  int cg;
  long int i;
  double delta, alpha, lambda, omega, gammag;

  convert_double_to_float_vec3_soa(in,ip.loc_p_f);


  fermion_matrix_multiplication_shifted_f(ip.u_f,ip.loc_s_f,solution,ip.loc_h_f,pars,shift);

  combine_in1_minus_in2_f(ip.loc_p_f,ip.loc_s_f,ip.loc_r_f);
  assign_in_to_out_f(ip.loc_r_f,ip.loc_p_f);

  delta=l2norm2_global_f(ip.loc_r_f);

  double source_norm = l2norm2_global(in);
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
    combine_in1xfactor_plus_in2_f(ip.loc_p_f,omega,solution,solution);

    if( cg % inverter_tricks.magicTouchEvery == 0 ){
        // calculation of r from "first principle" in double precision
        //here loc_r= tmp
        convert_float_to_double_vec3_soa(solution,ip.loc_h);
        fermion_matrix_multiplication_shifted(ip.u,ip.loc_s,ip.loc_h,ip.loc_r,pars,shift);
        combine_in1_minus_in2(in,ip.loc_s,ip.loc_r);
        convert_double_to_float_vec3_soa(ip.loc_r,ip.loc_r_f);

    }
    else combine_in1xfactor_plus_in2_f(ip.loc_s_f,-omega,ip.loc_r_f,ip.loc_r_f);

    lambda = l2norm2_global_f(ip.loc_r_f);
    gammag=lambda/delta;
    delta=lambda;
    // p=r+gammag*p
    combine_in1xfactor_plus_in2_f(ip.loc_p_f,gammag,ip.loc_r_f,ip.loc_p_f);


      if (verbosity_lv > 3 && cg%100==0 && 0==devinfo.myrank  ){

      printf("%d\t%1.1e\n",cg, sqrt(lambda/source_norm)/res);fflush(stdout);
      
      }



  } while( (sqrt(lambda/source_norm)>res) && cg<max_cg);


  if (verbosity_lv > 3  && 0==devinfo.myrank ) printf("\n");
#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER_FULL_OPENACC))
  // s = M x
  fermion_matrix_multiplication_shifted_f(ip.u_f,ip.loc_s_f,solution,ip.loc_h_f,pars,shift);
  convert_double_to_float_vec3_soa(in,ip.loc_p_f);  // y = in 
  combine_in1_minus_in2_f(ip.loc_p_f,ip.loc_s_f,ip.loc_h_f); // r = s - y  
  double  giustoono=l2norm2_global_f(ip.loc_h_f)/source_norm;
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

