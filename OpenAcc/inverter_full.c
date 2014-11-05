#include "./struct_c_def.c"
#include "./fermionic_utilities.c"
#include "./fermion_matrix.c"
#include "openacc.h"





void ker_invert_openacc(  const __restrict su3_soa * const u,
			  __restrict vec3_soa * const out,
			  const __restrict vec3_soa * const in,
			  double res,
			  const __restrict vec3_soa * const trialSolution,
			  __restrict vec3_soa * const loc_r,
			  __restrict vec3_soa * const loc_h,
			  __restrict vec3_soa * const loc_s,
			  __restrict vec3_soa * const loc_p
			  ){

  int cg;
  long int i;
  double delta, alpha, lambda, omega, gammag;
  d_complex aux_comp=0.0+0.0I;

  assign_in_to_out(trialSolution,out);

  acc_Doe(u,loc_h,out);
  acc_Deo(u,loc_s,loc_h);

  combine_in1xm2_minus_in2(out,loc_s,out);

  combine_in1_minus_in2(loc_r,in,loc_s);
  assign_in_to_out(loc_r,loc_p);

  assign_in_to_out(loc_r,loc_h);
  delta=   real_scal_prod_global(loc_h,loc_r);


  // loop over cg iterations
  cg=0;
  do {
    
    cg++;
    
    // s=(M^dag M)p    alpha=(p,s)
    acc_Doe(u,loc_h,loc_p);
    acc_Deo(u,loc_s,loc_h);
    combine_in1xm2_minus_in2(loc_p,loc_s,loc_s);
    alpha = real_scal_prod_global(loc_p,loc_s);
    omega=delta/alpha;     
    // out+=omega*p  r-=omega*s
    // lambda=(r,r);
    combine_in1xfactor_plus_in2(loc_p,omega,out,out);
    combine_in1xfactor_plus_in2(loc_s,-omega,loc_r,loc_r);

    assign_in_to_out(loc_r,loc_h);
    lambda = real_scal_prod_global(loc_h,loc_r);
    
    gammag=lambda/delta;
    delta=lambda;
    // p=r+gammag*p
    combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);
    
  } while( (sqrt(lambda)>res) && cg<max_cg);

#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER))
  cout << "\tterminated invert terminated in "<<cg<<" iterations [";

  acc_Doe(u,loc_h,out);
  acc_Deo(u,loc_s,loc_h);

  double giustoono;
  combine_in1xm2_minus_in2_minus_in3(out,loc_s,in,loc_p);
  assign_in_to_out(loc_p,loc_h);
  giustoono = real_scal_prod_global(loc_h,loc_p);
  // mettere il giusto print qui!!!
  //  cout << " res/stop_res="<< sqrt(giustoono)/res << " ,  stop_res="<< res << " ]"<<endl;
  /*
 if(cg==max_cg)
    {
    ofstream err_file;
    err_file.open(QUOTEME(ERROR_FILE), ios::app);   
    err_file << "WARNING: maximum number of iterations reached in invert\n";
    err_file.close();
    }

   */

  #endif


}




void invert_openacc_full( su3COM_soa  *conf,vec3COM_soa *out,vec3COM_soa *in,double res,vec3COM_soa *trialSol){

  su3_soa  * conf_acc;
  vec3_soa * ferm_in_acc;
  vec3_soa * ferm_out_acc;
  vec3_soa * ferm_try_acc;
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&ferm_in_acc , ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&ferm_out_acc, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&ferm_try_acc, ALIGN, sizeof(vec3_soa));

  // AUXILIARY FERMION FIELDS FOR THE INVERTER
  vec3_soa * kloc_r;
  vec3_soa * kloc_h;
  vec3_soa * kloc_s;
  vec3_soa * kloc_p;
  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa));


  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);
  convert_vec3COM_soa_to_vec3_soa(in,ferm_in_acc);
  convert_vec3COM_soa_to_vec3_soa(trialSol,ferm_try_acc);
  


#pragma acc data copyin(conf_acc[0:8]) copyin(ferm_in_acc[0:1]) copyin(ferm_try_acc[0:1]) copyout(ferm_out_acc[0:1]) create(kloc_r[0:1]) create(kloc_h[0:1]) create(kloc_s[0:1]) create(kloc_p[0:1])
  {
    ker_invert_openacc(conf_acc,ferm_out_acc,ferm_in_acc,res,ferm_try_acc,kloc_r,kloc_h,kloc_s,kloc_p);
  }

  for(dir=0;dir<8;dir++)  convert_su3_soa_to_su3COM_soa(&conf_acc[dir],&conf[dir]); // questa traduzione in linea di principio non serve --> toglierla ... ?
  convert_vec3_soa_to_vec3COM_soa(ferm_out_acc,out);

  free(conf_acc);
  free(ferm_in_acc);
  free(ferm_out_acc);
}
