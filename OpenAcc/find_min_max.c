
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
					 __restrict vec3_soa * const loc_r,
					 __restrict vec3_soa * const loc_h,
					 __restrict vec3_soa * const loc_p
					 ){

  int loop_count;
  double norm, inorm, old_norm;
  double temp;
  // starting gauss vector p
  //   deve arrivargli gia' gaussiano

  norm=sqrt(l2norm2_global(loc_p));
  loop_count=0;
  //  printf("norm  out  %.18lf\n",norm);
  // loop start
  do{
      // normalize  p
      // and r=p
      inorm=1.0/norm;
      multiply_fermion_x_doublefactor(loc_p,inorm);
      assign_in_to_out(loc_p,loc_r);
      old_norm=norm;
      // p=DeoDoe r
      acc_Doe(u,loc_h,loc_p);
      acc_Deo(u,loc_p,loc_h);

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
    acc_Doe(u,loc_h,loc_p);
    acc_Deo(u,loc_p,loc_h);
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
					   __restrict vec3_soa * const loc_r,
					   __restrict vec3_soa * const loc_h,
					   __restrict vec3_soa * const loc_p1,
					   __restrict vec3_soa * const loc_p2,
					   int usestored,
					   double *minmax
					   ){

  if(usestored == 1){

    minmax[1] = ker_find_max_eigenvalue_openacc(u,loc_r,loc_h,loc_p1);
    minmax[0] = ker_find_min_eigenvalue_openacc(u,loc_r,loc_h,loc_p2,minmax[1]);


    //    minmax[0] =  mass2;
    //    minmax[1] =  ker_find_max_eigenvalue_openacc(u,loc_r,loc_h,loc_p);
  }
  // altrimenti se usestored == 0 allora lascia gli autocosi al valore che avevano gia' prima
}



void find_min_max_openacc(const  su3COM_soa  *conf, const vec3COM_soa *gaussian1, const vec3COM_soa *gaussian2,double *minmax){
  // minmax[0] --> minimo
  // minmax[1] --> massimo

  su3_soa  * conf_acc;
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));

  // AUXILIARY FERMION FIELDS FOR THE INVERTER
  vec3_soa * kloc_r;
  vec3_soa * kloc_h;
  vec3_soa * kloc_p1;
  vec3_soa * kloc_p2;
  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_p1, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_p2, ALIGN, sizeof(vec3_soa));


  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);
  convert_vec3COM_soa_to_vec3_soa(gaussian1,kloc_p1);
  convert_vec3COM_soa_to_vec3_soa(gaussian2,kloc_p2);

  struct timeval t0, t1;
  gettimeofday ( &t0, NULL );

#pragma acc data copyin(conf_acc[0:8]) create(kloc_r[0:1]) create(kloc_h[0:1]) copyin(kloc_p1[0:1]) copyin(kloc_p2[0:1])
  {
    minmax[1] = ker_find_max_eigenvalue_openacc(conf_acc,kloc_r,kloc_h,kloc_p1);
    minmax[0] = ker_find_min_eigenvalue_openacc(conf_acc,kloc_r,kloc_h,kloc_p2,minmax[1]);
  }

  printf("MIN eigen =  %f \n",minmax[0]);
  printf("MAX eigen =  %f \n",minmax[1]);

  gettimeofday ( &t1, NULL );
  double dt_tot = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  printf("FULL OPENACC FINDMINMAX times:        Tot time: %f sec  \n",dt_tot);

  free(conf_acc);
  free(kloc_r);
  free(kloc_h);
  free(kloc_p1);
  free(kloc_p2);
}

#endif

