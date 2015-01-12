#ifndef INVERTER_SHIFT_MULTI_FULL_C_
#define INVERTER_SHIFT_MULTI_FULL_C_

#ifndef INCLUDE_ACC_X_FERM_MATRIX
#define INCLUDE_ACC_X_FERM_MATRIX
#include "./struct_c_def.c"
#include "./fermionic_utilities.c"
#include "openacc.h"
#endif


#ifndef INCLUDE_ACC_FERM_MATRIX
#define INCLUDE_ACC_FERM_MATRIX
#include "./fermion_matrix.c"
#endif

#include "./inverter_full.c"
#include "./find_min_max.c"

#define DEBUG_INVERTER_SHIFT_MULTI_FULL_OPENACC

void ker_invert_openacc_shiftmulti(   const __restrict su3_soa * const u,
				     __restrict ACC_ShiftMultiFermion * const out,
				     const __restrict ACC_MultiFermion * const in,
				     double residuo,
				     const COM_RationalApprox * const approx,
				     __restrict vec3_soa * const loc_r,
				     __restrict vec3_soa * const loc_h,
				     __restrict vec3_soa * const loc_s,
				     __restrict vec3_soa * const loc_p,
				     __restrict ACC_ShiftFermion * const p_shiftferm
){
  // AUXILIARY VARIABLES FOR THE INVERTER 
  int  cg;
  int *cg_aux;
  posix_memalign((void **)&cg_aux, ALIGN, no_ps*sizeof(int));
  int *flag;
  posix_memalign((void **)&flag, ALIGN, max_approx_order*sizeof(int));

  double *zeta_i,*zeta_ii,*zeta_iii,*omegas,*gammas;
  posix_memalign((void **)&zeta_i,   ALIGN,  max_approx_order*sizeof(double));
  posix_memalign((void **)&zeta_ii,  ALIGN,  max_approx_order*sizeof(double));
  posix_memalign((void **)&zeta_iii, ALIGN,  max_approx_order*sizeof(double));
  posix_memalign((void **)&omegas,   ALIGN,  max_approx_order*sizeof(double));
  posix_memalign((void **)&gammas,   ALIGN,  max_approx_order*sizeof(double));
  printf("Allocated auxiliary variables \n");

  int pseudofermion,iter;
  double alpha, delta, lambda, omega, omega_save, gammag, fact;

  printf("Inside the kernel \n");
  printf("Ordine della approssimazione:   %i \n",approx[0].COM_approx_order);

  alpha=0.0;

  for(pseudofermion=0; pseudofermion<no_ps; pseudofermion++){

    // trial solution out = 0, set all flag to 1                                                                                                           
    for(iter=0; iter<(approx[0].COM_approx_order); iter++){
      flag[iter]=1;
      set_shiftmulti_to_zero(iter,pseudofermion,out);
    }

    // r=in, p=phi delta=(r,r)                                                                                                                             
    extract_pseudo_and_assign_to_fermion(in,pseudofermion,loc_r);
    assign_in_to_out(loc_r,loc_p);
    delta=l2norm2_global(loc_r);
    printf("delta    %.18lf\n",delta);
    omega=1.0;

    for(iter=0; iter<(approx[0].COM_approx_order); iter++){
      // ps_0=phi
      extract_pseudo_and_assign_to_shiftfermion(in,pseudofermion,p_shiftferm,iter);
      zeta_i[iter]=1.0;         // zeta_{-1}=1.0
      zeta_ii[iter]=1.0;        // zeta_{ 0}=1.0
      gammas[iter]=0.0;         // gammas_{-1}=0.0
    }
    gammag=0.0;
    cg=0;
    do {      // loop over cg iterations
      cg++;

      // s=(M^dagM)p, alhpa=(p,s)=(p,Ap)
      acc_Doe(u,loc_h,loc_p);
      acc_Deo(u,loc_s,loc_h);
      combine_in1xm2_minus_in2(loc_p,loc_s);

      printf("component_loc_s[0]    %.18e\n",creal(loc_s->c0[0]));
      printf("component_loc_p[0]    %.18e\n",creal(loc_p->c0[0]));

      alpha = real_scal_prod_global(loc_p,loc_s);
      printf("alpha    %.18lf\n",alpha);
      printf("mass2    %.18lf\n",mass2);

      omega_save=omega;   // omega_save=omega_(j-1)                                                                                                        
      omega=-delta/alpha;  // omega = (r_j,r_j)/(p_j, Ap_j)               
      printf("omega    %.18lf\n",omega);

      // out-=omegas*ps
      for(iter=0; iter<(approx[0].COM_approx_order); iter++){
	if(flag[iter]==1){
	  zeta_iii[iter] = (zeta_i[iter]*zeta_ii[iter]*omega_save)/
	    ( omega*gammag*(zeta_i[iter]-zeta_ii[iter])+
	      zeta_i[iter]*omega_save*(1.0-(approx[0].COM_RA_b[iter])*omega) );
	  omegas[iter]=omega*zeta_iii[iter]/zeta_ii[iter];
	  
	  combine_shiftmulti_minus_shiftfermion_x_factor_back_into_shiftmulti(out,p_shiftferm,pseudofermion,iter,omegas[iter]);
	}
      }

      // r+=omega*s; lambda=(r,r)
      combine_add_factor_x_in2_to_in1(loc_r,loc_s,omega);
      lambda=l2norm2_global(loc_r);
      gammag=lambda/delta;

      // p=r+gammag*p
      combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);

      for(iter=0; iter<(approx[0].COM_approx_order); iter++){
	if(flag[iter]==1){
	  gammas[iter]=gammag*zeta_iii[iter]*omegas[iter]/(zeta_ii[iter]*omega);
	  combine_shiftferm_x_fact1_plus_ferm_x_fact2_back_into_shiftferm(p_shiftferm,iter,gammas[iter],loc_r,zeta_iii[iter]);
	  
	  fact=sqrt(delta*zeta_ii[iter]*zeta_ii[iter]);
	  if(fact<residuo){
	    flag[iter]=0;
	  }
	  zeta_i[iter]=zeta_ii[iter];
	  zeta_ii[iter]=zeta_iii[iter];
	}
      }
      delta=lambda;
      printf("Pseudoferm: %i   Iteration: %i    --> residue = %e   (target = %e) \n", pseudofermion, cg, sqrt(lambda), residuo);

      printf("Inside multishift  ( step = %i )\n");
      printf("lambda   %.18lf\n",lambda );
      printf("omega    %.18lf\n",omega );
      printf("gammag   %.18lf\n",gammag );

      for(iter=0; iter<(approx[0].COM_approx_order); iter++){
	printf("zeta_i[%i] =  %.18lf\n",iter,zeta_i[iter]);
	printf("gammas[%i] =  %.18lf\n",iter,gammas[iter]);
	printf("omegas[%i] =  %.18lf\n",iter,omegas[iter]);
      }


    } while(sqrt(lambda)>residuo && cg<max_cg); // end of cg iterations 

    if(cg==max_cg)
      {
	printf("WARNING: maximum number of iterations reached in invert\n");
      }
    cg_aux[pseudofermion]=cg;



    
  } // end loop on pseudofermions

#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER_SHIFT_MULTI_FULL_OPENACC))
  printf("Terminated openacc_multips_shift_invert   ( target res = %f ) \n ", residuo);
  int i;
  for(i=0; i<no_ps; i++)
    {
      printf("\t CG count [pseudoferm n. %i ] = %i \n",i,cg_aux[i]);
    }
  // test 
  for(pseudofermion=0; pseudofermion<no_ps; pseudofermion++){
    for(iter=0; iter<(approx[0].COM_approx_order); iter++){
      extract_from_shiftmulti_and_assign_to_fermion(out,iter,pseudofermion,loc_p);
      acc_Doe(u,loc_h,loc_p);
      acc_Deo(u,loc_s,loc_h);
      combine_in1_x_fact_minus_in2_minus_multiin3_back_into_in1(loc_p,mass2+approx[0].COM_RA_b[iter],loc_s,in,pseudofermion);
      double  giustoono=l2norm2_global(loc_p);
      printf("\t\t pseudofermion= %i    iter_approx= %i      res/stop_res= %e          stop_res= %e \n",pseudofermion,iter,sqrt(giustoono)/residuo,residuo);
      printf("\t\t Shifted mass2  =  %f \n",approx[0].COM_RA_b[iter]);

    }
  }
#endif

  free(cg_aux);
  free(flag);
  free(zeta_i);
  free(zeta_ii);
  free(zeta_iii);
  free(omegas);
  free(gammas);

}

void ker_openacc_recombine_shiftmulti_to_multi( const __restrict ACC_ShiftMultiFermion * const in_shiftmulti,
						const __restrict ACC_MultiFermion * const in_multi,
						      __restrict ACC_MultiFermion * const out_multi,
						const COM_RationalApprox * const approx
						){
  int ips;
  int ih;
  int iter=0;
#pragma acc kernels present(out_multi) present(in_multi) present(in_shiftmulti)
#pragma acc loop independent
  for(ips=0; ips < no_ps; ips++){
#pragma acc loop independent
    for(ih=0; ih < sizeh; ih++){

      out_multi->multi[ips].c0[ih] =  (in_multi->multi[ips].c0[ih])*(approx[0].COM_RA_a0);
      out_multi->multi[ips].c1[ih] =  (in_multi->multi[ips].c1[ih])*(approx[0].COM_RA_a0);
      out_multi->multi[ips].c2[ih] =  (in_multi->multi[ips].c2[ih])*(approx[0].COM_RA_a0);

      for(iter=0; iter<(approx[0].COM_approx_order); iter++){  // questo loop non lo vogliamo parallelizzare per forza ... forse puo andare bene cosi'
	out_multi->multi[ips].c0[ih] +=  (approx[0].COM_RA_a[iter]) * (in_shiftmulti->shiftmulti[iter][ips].c0[ih]);
	out_multi->multi[ips].c1[ih] +=  (approx[0].COM_RA_a[iter]) * (in_shiftmulti->shiftmulti[iter][ips].c1[ih]);
	out_multi->multi[ips].c2[ih] +=  (approx[0].COM_RA_a[iter]) * (in_shiftmulti->shiftmulti[iter][ips].c2[ih]);
      }


    }
  }
}

void multips_invert_openacc_full(const su3COM_soa  *conf,COM_ShiftMultiFermion *out, COM_MultiFermion *in, double res, COM_RationalApprox *approx){
  
  printf("Residuo target = %f \n", res);

  su3_soa  * conf_acc;
  ACC_MultiFermion * ferm_in_acc;
  ACC_ShiftMultiFermion * ferm_out_acc;

  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&ferm_in_acc , ALIGN, sizeof(ACC_MultiFermion));
  posix_memalign((void **)&ferm_out_acc, ALIGN, sizeof(ACC_ShiftMultiFermion));

  printf("Allocated acc fermions \n");
  printf("Size of ACC_ShiftMultiFermion = %i \n",sizeof(ACC_ShiftMultiFermion));
  printf("Size of double = %i \n",sizeof(double));

  // AUXILIARY FERMION FIELDS FOR THE INVERTER (for the non-multishift inverter --> we need to change these variables)
  vec3_soa * kloc_r;
  vec3_soa * kloc_h;
  vec3_soa * kloc_s;
  vec3_soa * kloc_p;
  ACC_ShiftFermion *k_p_shiftferm;
  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&k_p_shiftferm, ALIGN, sizeof(ACC_ShiftFermion));
  printf("Allocated auxiliary fermions \n");


  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);
  convert_COM_MultiFermion_to_ACC_MultiFermion(in,ferm_in_acc);
  printf("Converted conf and multi fermion \n");


  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

  #pragma acc data copyin(conf_acc[0:8]) copyin(ferm_in_acc[0:1]) copyin(approx[0:1])  copyout(ferm_out_acc[0:1])  create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1])  create(k_p_shiftferm[0:1])
  {
  gettimeofday ( &t1, NULL );

    ker_invert_openacc_shiftmulti(conf_acc,ferm_out_acc,ferm_in_acc,res,approx,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
  gettimeofday ( &t2, NULL );

  }


  gettimeofday ( &t3, NULL );
  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("FULL OPENACC MULTISHIFT INVERSION times:        Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);
  
  
  
  convert_ACC_ShiftMultiFermion_to_COM_ShiftMultiFermion(ferm_out_acc,out);
  printf("Reconverted fermions \n");

  free(conf_acc);
  free(ferm_in_acc);
  free(ferm_out_acc);


  free(kloc_r);
  free(kloc_s);
  free(kloc_h);
  free(kloc_p);

  printf("Freed memory \n");

}


void first_inv_approx_calc_openacc(const su3COM_soa  *conf, COM_MultiFermion *out, const COM_MultiFermion *in, double res, const COM_RationalApprox *approx){
  
  printf("Residuo target = %f \n", res);

  su3_soa  * conf_acc;
  ACC_MultiFermion * ferm_in_acc;
  ACC_MultiFermion * ferm_out_acc;
  ACC_ShiftMultiFermion * ferm_shiftmulti_acc;

  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&ferm_in_acc  , ALIGN, sizeof(ACC_MultiFermion));
  posix_memalign((void **)&ferm_out_acc , ALIGN, sizeof(ACC_MultiFermion));
  posix_memalign((void **)&ferm_shiftmulti_acc, ALIGN, sizeof(ACC_ShiftMultiFermion));

  // AUXILIARY FERMION FIELDS FOR THE INVERTER (for the non-multishift inverter --> we need to change these variables)
  vec3_soa * kloc_r;
  vec3_soa * kloc_h;
  vec3_soa * kloc_s;
  vec3_soa * kloc_p;
  ACC_ShiftFermion *k_p_shiftferm;
  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&k_p_shiftferm, ALIGN, sizeof(ACC_ShiftFermion));
  printf("Allocated auxiliary fermions \n");


  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);
  convert_COM_MultiFermion_to_ACC_MultiFermion(in,ferm_in_acc);
  printf("Converted conf and multi fermion \n");


  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );

#pragma acc data copyin(conf_acc[0:8]) copyin(ferm_in_acc[0:1]) copyin(approx[0:1])  copyout(ferm_out_acc[0:1])  create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1])  create(k_p_shiftferm[0:1]) create(ferm_shiftmulti_acc[0:1])
  {
  gettimeofday ( &t1, NULL );

  ker_invert_openacc_shiftmulti(conf_acc,ferm_shiftmulti_acc,ferm_in_acc,res,approx,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);

  ker_openacc_recombine_shiftmulti_to_multi(ferm_shiftmulti_acc,ferm_in_acc,ferm_out_acc,approx);

  gettimeofday ( &t2, NULL );

  }


  int ips = 1;
  int ish = 108986;

  printf("%e    %e \n",creal(ferm_out_acc->multi[ips].c0[ish]),cimag(ferm_out_acc->multi[ips].c0[ish]));
  printf("%e    %e \n",creal(ferm_out_acc->multi[ips].c1[ish]),cimag(ferm_out_acc->multi[ips].c1[ish]));
  printf("%e    %e \n",creal(ferm_out_acc->multi[ips].c2[ish]),cimag(ferm_out_acc->multi[ips].c2[ish]));
  convert_ACC_MultiFermion_to_COM_MultiFermion(ferm_out_acc,out);


  gettimeofday ( &t3, NULL );
  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("FULL FIRST INV APPROX CALC                      Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);
  
  
  
  printf("Reconverted fermions \n");

  free(conf_acc);
  free(ferm_in_acc);
  free(ferm_out_acc);
  free(ferm_shiftmulti_acc);

  free(kloc_r);
  free(kloc_s);
  free(kloc_h);
  free(kloc_p);

  printf("Freed memory \n");

}

#endif

