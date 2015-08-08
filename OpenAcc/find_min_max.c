
#ifndef FIND_MIN_MAX_FULLACC_C_
#define FIND_MIN_MAX_FULLACC_C_


#ifndef INCLUDE_ACC
#define INCLUDE_ACC
#include "./struct_c_def.c"
#include "./fermionic_utilities.c"
#include "./fermion_matrix.c"
#include "openacc.h"
#endif
#include <stdio.h>


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
  loop_count=0;
  // loop start
  do{
      // normalize  p
      // and r=p
      inorm=1.0/norm;
      multiply_fermion_x_doublefactor(loc_p,inorm);
      assign_in_to_out(loc_p,loc_r);
      old_norm=norm;
      // p=DeoDoe r
      acc_Doe(u,loc_h,loc_p,pars,backfield);
      acc_Deo(u,loc_p,loc_h,pars,backfield);

      // p=(M^dag M)r
      combine_in1xm2_minus_in2(loc_r,loc_p);
      norm=sqrt(l2norm2_global(loc_p));
      old_norm=fabs(old_norm-norm);
      old_norm/=norm;
      loop_count++;
    } while(old_norm>1.0e-5);    // loop end
  double max=norm;
  return max;
}

// find the minimum eigenvalue of the fermion matrix
// use loc_h, loc_p, loc_r
double ker_find_min_eigenvalue_openacc(  __restrict su3_soa * const u,
					 double_soa * const backfield,
					 ferm_param *pars,
					 __restrict vec3_soa * const loc_r,
					 __restrict vec3_soa * const loc_h,
					 __restrict vec3_soa * const loc_p,
					 double max
					 ){
  int loop_count;
  double norm, inorm, old_norm, delta;
  delta=max-mass2;
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
					   int usestored,
					   double *minmax
					   ){
  // minmax[0] --> minimo
  // minmax[1] --> massimo
  if(usestored == 1){
    minmax[1] = ker_find_max_eigenvalue_openacc(u,backfield,pars,loc_r,loc_h,loc_p1);
    minmax[0] = ker_find_min_eigenvalue_openacc(u,backfield,pars,loc_r,loc_h,loc_p2,minmax[1]); //--> si potrebbe mettere direttamente mass2
  }
  // altrimenti se usestored == 0 allora lascia gli autocosi al valore che avevano gia' prima
}



#endif

