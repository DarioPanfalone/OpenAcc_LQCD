#ifndef FIND_MIN_MAX_C_
#define FIND_MIN_MAX_C_


#include "./struct_c_def.h"
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"
#ifndef __GNUC__
 #include "openacc.h"
#endif
#include <stdio.h>
#include "./find_min_max.h"


//#define PRINT_EIGENVALUES


// find the maximum eigenvalue of the fermion matrix
// use loc_h, loc_p, loc_r
double ker_find_max_eigenvalue_openacc(  __restrict su3_soa * const u,
					 double_soa * const backfield,
					 ferm_param *pars,
					 __restrict vec3_soa * const loc_r,
					 __restrict vec3_soa * const loc_h,
					 __restrict vec3_soa * const loc_p
					 ){

  int loop_count;
  double norm, inorm, old_norm;
  double temp;
  // starting gauss vector p   (deve arrivargli gia' gaussiano)
  norm=sqrt(l2norm2_global(loc_p));
  printf("Norm:%lf\n",norm  );
  loop_count=0;
  // loop start
  do{
      // normalize  p
      // and r=p
      inorm=1.0/norm;
      multiply_fermion_x_doublefactor(loc_p,inorm);
      SETREQUESTED(loc_r);
      assign_in_to_out(loc_p,loc_r);

      old_norm=norm;
      // p=DeoDoe r
      SETREQUESTED(loc_h);
      acc_Doe(u,loc_h,loc_p,pars,backfield);
      SETFREE(loc_p);
      SETREQUESTED(loc_p);
      acc_Deo(u,loc_p,loc_h,pars,backfield);
      SETFREE(loc_h);

//#pragma acc update host(loc_p[0:1])
//      printf("%.18lf %.18lf\n",creal(loc_p->c0[0]),cimag(loc_p->c0[0]));



      // p=(M^dag M)r
      //      combine_in1xm2_minus_in2(loc_r,loc_p);
      combine_in1xferm_mass_minus_in2(loc_r,pars->ferm_mass*pars->ferm_mass,loc_p);
      SETFREE(loc_r);


      norm=sqrt(l2norm2_global(loc_p));
      old_norm=fabs(old_norm-norm);
      old_norm/=norm;
      loop_count++;
      if(loop_count %10 == 0)   printf("        iterations of max computat  = %d, norm = %lf, old_norm = %lf\n",loop_count, norm, old_norm);
  } while(old_norm>1.0e-5);    // loop end
  double max=norm;
  return max;
}

// find the minimum eigenvalue of the fermion matrix
// use loc_h, loc_p, loc_r
double ker_find_min_eigenvalue_openacc(  __restrict su3_soa * const u,
					 __restrict double_soa * const backfield,
					 __restrict ferm_param * const pars,
					 __restrict vec3_soa * const loc_r,
					 __restrict vec3_soa * const loc_h,
					 __restrict vec3_soa * const loc_p,
					 double max
					 ){
  int loop_count;
  double norm, inorm, old_norm, delta;
  delta=max-(pars->ferm_mass)*(pars->ferm_mass);
  norm=sqrt(l2norm2_global(loc_p));
  loop_count=0;
  // loop start
  do{
    // normalize  p 
    // and r=p 
    inorm=1.0/norm;
    multiply_fermion_x_doublefactor(loc_p,inorm);
    // r=p
    assign_in_to_out(loc_p,loc_r);
    old_norm=norm;
    // p=DeoDoe r 
    acc_Doe(u,loc_h,loc_p,pars,backfield);
    acc_Deo(u,loc_p,loc_h,pars,backfield);
    // p=max r - (M^dag M)r 
    combine_add_factor_x_in2_to_in1(loc_p,loc_r,delta);
    norm=sqrt(l2norm2_global(loc_p));
    old_norm=fabs(old_norm-norm);
    old_norm/=norm;
    loop_count++;
  }  while(old_norm>1.0e-5);   // loop end
  double  min=max-norm;
  return min;
}


void find_min_max_eigenvalue_soloopenacc(  __restrict su3_soa * const u,
					   double_soa * const backfield,
					   ferm_param *pars,
					   __restrict vec3_soa * const loc_r,
					   __restrict vec3_soa * const loc_h,
					   __restrict vec3_soa * const loc_p1,
					   __restrict vec3_soa * const loc_p2,
					   double *minmax
					   ){
  // minmax[0] --> minimo
  // minmax[1] --> massimo

  minmax[0] = pars->ferm_mass * pars->ferm_mass;
  minmax[1] = ker_find_max_eigenvalue_openacc(u,backfield,pars,loc_r,loc_h,loc_p1);

#ifdef PRINT_EIGENVALUES
  printf("    Computed the min and max eigs : [ min= %.18lf; max= %.18lf ]  \n",minmax[0],minmax[1]);
#endif
  //  ora il minimo e' messo a m*m, volendo lo si puo' calcolare con la routine seguente.
  //  minmax[0] = ker_find_min_eigenvalue_openacc(u,backfield,pars,loc_r,loc_h,loc_p2,minmax[1]); //--> si potrebbe mettere direttamente mass2

  SETFREE(loc_r);
  SETFREE(loc_h);
  SETFREE(loc_p1);
  SETFREE(loc_p2);
}



#endif

