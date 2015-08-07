#ifndef INVERTER_SHIFT_MULTI_FULL_C_
#define INVERTER_SHIFT_MULTI_FULL_C_

#define DEBUG_INVERTER_SHIFT_MULTI_FULL_OPENACC


void ker_invert_openacc_shiftmulti(   __restrict su3_soa * const u, // non viene mai aggiornata qui dentro
				      __restrict ACC_ShiftMultiFermion * const out,
				      __restrict ACC_MultiFermion * const in, // non viene mai aggiornato qui dentro
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

  int pseudofermion, iter, maxiter;
  double alpha, delta, lambda, omega, omega_save, gammag, fact;

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

    //Find the max approx[0].COM_approx_order where flag[approx[0].COM_approx_order] == 1
    for(iter=0; iter<(approx[0].COM_approx_order); iter++) {
      if(flag[iter]==1) {
        maxiter = iter+1;
      }
    }

    do {      // loop over cg iterations
      cg++;

      // s=(M^dagM)p, alhpa=(p,s)=(p,Ap)
      acc_Doe(u,loc_h,loc_p);
      acc_Deo(u,loc_s,loc_h);
      combine_in1xm2_minus_in2(loc_p,loc_s);

      alpha = real_scal_prod_global(loc_p,loc_s);

      omega_save=omega;   // omega_save=omega_(j-1) 
      omega=-delta/alpha;  // omega = (r_j,r_j)/(p_j, Ap_j)               

      // out-=omegas*ps
      //for (iter=0; iter<(approx[0].COM_approx_order); iter++) {
      for (iter=0; iter<maxiter; iter++) {
     	  if (flag[iter]==1) {
	        zeta_iii[iter] = (zeta_i[iter]*zeta_ii[iter]*omega_save)/
	                         ( omega*gammag*(zeta_i[iter]-zeta_ii[iter])+
	                         zeta_i[iter]*omega_save*(1.0-(approx[0].COM_RA_b[iter])*omega) );
	      
            omegas[iter] = omega*zeta_iii[iter]/zeta_ii[iter];
	      }
      }

	    combine_shiftmulti_minus_shiftfermion_x_factor_back_into_shiftmulti_all(out,p_shiftferm,pseudofermion,maxiter,flag,omegas);
  
//	    combine_shiftmulti_minus_shiftfermion_x_factor_back_into_shiftmulti(out,p_shiftferm,pseudofermion,iter,omegas[iter]);

//	      }

//      }

      // r+=omega*s; lambda=(r,r)
      combine_add_factor_x_in2_to_in1(loc_r,loc_s,omega);
      lambda=l2norm2_global(loc_r);
      gammag=lambda/delta;

      // p=r+gammag*p
      combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);

//ORIGINALE
/*
      for(iter=0; iter<(approx[0].COM_approx_order); iter++) {
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
*/
     for(iter=0; iter<maxiter; iter++) {
        if(flag[iter]==1){
          gammas[iter]=gammag*zeta_iii[iter]*omegas[iter]/(zeta_ii[iter]*omega);
        }
     }

     combine_shiftferm_x_fact1_plus_ferm_x_fact2_back_into_shiftferm_all(p_shiftferm,maxiter,flag,gammas,loc_r,zeta_iii);

     for(iter=0; iter<maxiter; iter++) {
        if(flag[iter]==1){

          fact=sqrt(delta*zeta_ii[iter]*zeta_ii[iter]);

          if(fact<residuo){
            flag[iter]=0;
          }

          zeta_i[iter]=zeta_ii[iter];
          zeta_ii[iter]=zeta_iii[iter];
        }
     }


//




     delta=lambda;

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

int multishift_invert(const __restrict su3_soa * const u,
                     ferm_param * pars,
                     const double_soa * backfield,
                     const __restrict vec3_soa * const out, // multi-fermion
				     const __restrict vec3_soa * const in, // single ferm
				     double residuo,
				     __restrict vec3_soa * const loc_r,
				     __restrict vec3_soa * const loc_h,
				     __restrict vec3_soa * const loc_s,
				     __restrict vec3_soa * const loc_p,
				     __restrict vec3_soa * const shiftferm // multi-ferm
){
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
  RationalApprox *approx = pars->approx,
  int  cg;
  double *zeta_i,*zeta_ii,*zeta_iii,*omegas,*gammas;
  int *flag;
  posix_memalign((void **)&flag,    ALIGN,approx->approx_order*sizeof(int));
  posix_memalign((void **)&zeta_i,  ALIGN,approx->approx_order*sizeof(double));
  posix_memalign((void **)&zeta_ii, ALIGN,approx->approx_order*sizeof(double));
  posix_memalign((void **)&zeta_iii,ALIGN,approx->approx_order*sizeof(double));
  posix_memalign((void **)&omegas,  ALIGN,approx->approx_order*sizeof(double));
  posix_memalign((void **)&gammas,  ALIGN,approx->approx_order*sizeof(double));
  //  printf("Allocated auxiliary variables \n");
  int iter;
  double alpha, delta, lambda, omega, omega_save, gammag, fact;
  //  printf("Inside the kernel \n");
  //  printf("Ordine della approssimazione:   %i \n",approx[0].COM_approx_order);
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
    //printf("delta    %.18lf\n",delta);
    omega=1.0;

    for(iter=0; iter<(approx->approx_order); iter++){
      // ps_0=phi
      assign_in_to_out(in,&shiftferm[iter]);
      zeta_i[iter]=1.0;         // zeta_{-1}=1.0
      zeta_ii[iter]=1.0;        // zeta_{ 0}=1.0
      gammas[iter]=0.0;         // gammas_{-1}=0.0
    }
    gammag=0.0;
    cg=0;
    do {      // loop over cg iterations
      cg++;

      // s=(M^dagM)p, alhpa=(p,s)=(p,Ap)
      fermion_matrix_multiplication(u,loc_s,loc_p,loc_h,pars,backfield);
      
//      printf("component_loc_s[0]    %f\n",creal(loc_s->c0[0]));
//      printf("component_loc_p[0]    %f\n",creal(loc_p->c0[0]));

      alpha = real_scal_prod_global(loc_p,loc_s);
//      printf("alpha    %.18lf\n",alpha);
//      printf("mass2    %.18lf\n",mass2);

      omega_save=omega;   // omega_save=omega_(j-1)
      omega=-delta/alpha;  // omega = (r_j,r_j)/(p_j, Ap_j)               
//      printf("omega    %.18lf\n",omega);

      // out-=omegas*ps
      for(iter=0; iter<(approx->approx_order); iter++){
          if(flag[iter]==1){
              zeta_iii[iter] = (zeta_i[iter]*zeta_ii[iter]*omega_save)/
                  ( omega*gammag*(zeta_i[iter]-zeta_ii[iter])+
                    zeta_i[iter]*omega_save*(1.0-(approx->RA_b[iter])*omega) );
              omegas[iter]=omega*zeta_iii[iter]/zeta_ii[iter];
          }
      }
      multiple_combine_in1_minus_in2x_factor_back_into_in1(out,shiftferm,approx->approx_order,flag,omegas);

      // r+=omega*s; lambda=(r,r)
      combine_add_factor_x_in2_to_in1(loc_r,loc_s,omega);
      lambda=l2norm2_global(loc_r);
      gammag=lambda/delta;

      // p=r+gammag*p
      combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);

      for(iter=0; iter<(approx->approx_order); iter++)
          if(flag[iter]==1)
              gammas[iter]=gammag*zeta_iii[iter]*omegas[iter]/(zeta_ii[iter]*omega);
      multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1(shiftferm,approx->approx_order,flag,gammas,loc_r,zeta_iii);


      for(iter=0; iter<(approx->approx_order); iter++)
          if(flag[iter]==1){
              fact=sqrt(delta*zeta_ii[iter]*zeta_ii[iter]);
              if(fact<residuo){
                  flag[iter]=0;
              }
              zeta_i[iter]=zeta_ii[iter];
              zeta_ii[iter]=zeta_iii[iter];
          }
      delta=lambda;
      printf("Iteration: %i    --> residue = %e   (target = %e) \n", cg, sqrt(lambda), residuo);
//      printf("Inside multishift  ( step = %i )\n");
//      printf("lambda   %.18lf\n",lambda );
//      printf("omega    %.18lf\n",omega );
//      printf("gammag   %.18lf\n",gammag );
//      for(iter=0; iter<(approx->approx_order); iter++){
//		printf("zeta_i[%i] =  %.18lf\n",iter,zeta_i[iter]);
//		printf("gammas[%i] =  %.18lf\n",iter,gammas[iter]);
//		printf("omegas[%i] =  %.18lf\n",iter,omegas[iter]);
//      }//DEBUG


    } while(sqrt(lambda)>residuo && cg<max_cg); // end of cg iterations 

    if(cg==max_cg)
      {
	printf("WARNING: maximum number of iterations reached in invert\n");
      }
    

#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER_SHIFT_MULTI_FULL_OPENACC))
  printf("Terminated multishift_invert_gl ( target res = %f ) \n ", residuo);
  int i;
      printf("\t CG count = %i \n",cg);
  // test 
    for(iter=0; iter<approx->approx_order; iter++){
      assign_in_to_out(&out[iter],loc_p);
      fermion_matrix_multiplication_shifted(u,loc_s,loc_p,loc_h,backfield,approx->RA_b[iter]);
      combine_in1_minus_in2(in,loc_s,loc_h); // r = s - y  
      double  giustoono=l2norm2_global(loc_h);
      printf("\t\titer_approx= %i      res/stop_res= %e        stop_res= %e \n",iter,sqrt(giustoono)/residuo,residuo);
      printf("\t\t Shifted mass2  =  %f \n",approx->RA_b[iter]);
    }
#endif

  free(flag);
  free(zeta_i);
  free(zeta_ii);
  free(zeta_iii);
  free(omegas);
  free(gammas);

  return cg;

}
void recombine_shifted_vec3_to_vec3(const __restrict vec3_soa* const in_shifted /*multi-fermion*/, const __restrict vec3_soa* const in /*single fermion*/, __restrict vec3_soa * const out /*multi fermion*/, const RationalApprox * const approx ){
  int ih;
  int iter=0;
#pragma acc kernels present(out) present(in) present(in_shifted)
#pragma acc loop independent
    for(ih=0; ih < SIZEH; ih++){

      out->c0[ih] =  in->c0[ih]*approx->RA_a0;
      out->c1[ih] =  in->c1[ih]*approx->RA_a0;
      out->c2[ih] =  in->c2[ih]*approx->RA_a0;

      for(iter=0; iter<approx->approx_order; iter++){  // questo loop non lo vogliamo parallelizzare per forza ... forse puo andare bene cosi'
	out->c0[ih] +=  approx->RA_a[iter] * in_shifted[iter].c0[ih];
	out->c1[ih] +=  approx->RA_a[iter] * in_shifted[iter].c1[ih];
	out->c2[ih] +=  approx->RA_a[iter] * in_shifted[iter].c2[ih];
      }
    }
}

static inline void vec1_directprod_conj_vec2_into_mat1( __restrict su3_soa * const aux_u,
							int idxh,
							__restrict vec3_soa  * const fer_l, // questo fermione e' costante e non viene modificato qui dentro
							int idl,
							__restrict vec3_soa  * const fer_r, // questo fermione e' costante e non viene modificato qui dentro
							int idr,
							double factor){
  // forse si possono risparmiare un po di conti andando a prendere le sole parti reale e immaginarie 
  //  e scrivendo i prodotti in modo molto piu' sbrodolato
  // in particolare per gli elementi lungo la diagonale
  d_complex r0 = factor * conj(fer_r->c0[idr]);
  d_complex r1 = factor * conj(fer_r->c1[idr]);
  d_complex r2 = factor * conj(fer_r->c2[idr]);
  d_complex l0 = fer_l->c0[idl];
  d_complex l1 = fer_l->c1[idl];
  d_complex l2 = fer_l->c2[idl];

  aux_u->r0.c0[idxh] += l0*r0;
  aux_u->r0.c1[idxh] += l0*r1;
  aux_u->r0.c2[idxh] += l0*r2;
  aux_u->r1.c0[idxh] += l1*r0;
  aux_u->r1.c1[idxh] += l1*r1;
  aux_u->r1.c2[idxh] += l1*r2;
  aux_u->r2.c0[idxh] += l2*r0;
  aux_u->r2.c1[idxh] += l2*r1;
  aux_u->r2.c2[idxh] += l2*r2;

}


#ifdef BACKFIELD
static inline  void mat1_times_auxmat_into_tamat(  __restrict su3_soa * const mat1, // e' costante e non viene modificato
						   const  int idx,
						   const  int eta,
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   const  int idx_aux,
						   __restrict tamat_soa * const ipdot,
						   const  int idipdot,
						   d_complex phase){
  d_complex mat1_00 = mat1->r0.c0[idx] * phase;
  d_complex mat1_01 = mat1->r0.c1[idx] * phase;
  d_complex mat1_02 = mat1->r0.c2[idx] * phase;
  d_complex mat1_10 = mat1->r1.c0[idx] * phase;
  d_complex mat1_11 = mat1->r1.c1[idx] * phase;
  d_complex mat1_12 = mat1->r1.c2[idx] * phase;
  //Compute 3rd matrix row from the first two
  d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
  //Multiply 3rd row by eta
  mat1_20 = (mat1_20)*eta * phase;
  mat1_21 = (mat1_21)*eta * phase;
  mat1_22 = (mat1_22)*eta * phase;
#else
static inline  void mat1_times_auxmat_into_tamat(  __restrict su3_soa * const mat1, // e' costante e non viene modificato
						   const  int idx,
						   const  int eta,
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   const  int idx_aux,
						   __restrict tamat_soa * const ipdot,
						   const  int idipdot){
  d_complex mat1_00 = mat1->r0.c0[idx];
  d_complex mat1_01 = mat1->r0.c1[idx];
  d_complex mat1_02 = mat1->r0.c2[idx];
  d_complex mat1_10 = mat1->r1.c0[idx];
  d_complex mat1_11 = mat1->r1.c1[idx];
  d_complex mat1_12 = mat1->r1.c2[idx];
  //Compute 3rd matrix row from the first two
  d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
  //Multiply 3rd row by eta
  mat1_20 = (mat1_20)*eta;
  mat1_21 = (mat1_21)*eta;
  mat1_22 = (mat1_22)*eta;
#endif

  d_complex auxmat_00 = auxmat->r0.c0[idx_aux];
  d_complex auxmat_01 = auxmat->r0.c1[idx_aux];
  d_complex auxmat_02 = auxmat->r0.c2[idx_aux];
  d_complex auxmat_10 = auxmat->r1.c0[idx_aux];
  d_complex auxmat_11 = auxmat->r1.c1[idx_aux];
  d_complex auxmat_12 = auxmat->r1.c2[idx_aux];
  d_complex auxmat_20 = auxmat->r2.c0[idx_aux];
  d_complex auxmat_21 = auxmat->r2.c1[idx_aux];
  d_complex auxmat_22 = auxmat->r2.c2[idx_aux];

  // product;
  d_complex a00 = mat1_00 * auxmat_00 + mat1_01 * auxmat_10 + mat1_02 * auxmat_20;
  d_complex a01 = mat1_00 * auxmat_01 + mat1_01 * auxmat_11 + mat1_02 * auxmat_21;
  d_complex a02 = mat1_00 * auxmat_02 + mat1_01 * auxmat_12 + mat1_02 * auxmat_22;

  mat1_00 = mat1_10 * auxmat_00 + mat1_11 * auxmat_10 + mat1_12 * auxmat_20;
  mat1_01 = mat1_10 * auxmat_01 + mat1_11 * auxmat_11 + mat1_12 * auxmat_21;
  mat1_02 = mat1_10 * auxmat_02 + mat1_11 * auxmat_12 + mat1_12 * auxmat_22;

  mat1_10 = mat1_20 * auxmat_00 + mat1_21 * auxmat_10 + mat1_22 * auxmat_20;
  mat1_11 = mat1_20 * auxmat_01 + mat1_21 * auxmat_11 + mat1_22 * auxmat_21; 
  mat1_12 = mat1_20 * auxmat_02 + mat1_21 * auxmat_12 + mat1_22 * auxmat_22;

  ipdot->c01[idipdot]  -= 0.5*((a01) - conj(mat1_00));
  ipdot->c02[idipdot]  -= 0.5*((a02) - conj(mat1_10));
  ipdot->c12[idipdot]  -= 0.5*((mat1_02) - conj(mat1_11));
  ipdot->rc00[idipdot] -= cimag(a00)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));
  ipdot->rc11[idipdot] -= cimag(mat1_01)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));

#ifdef BACKFIELD
}
#else
}
#endif


void direct_product_of_fermions_into_auxmat(__restrict vec3_soa  * const loc_s, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict vec3_soa  * const loc_h, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict su3_soa * const aux_u,
					    const RationalApprox * const approx,
					    int iter){
  
  //   ////////////////////////////////////////////////////////////////////////   //
  //    Riflettere se conviene tenere i loop sui siti pari e su quelli dispari    //
  //    separati come sono adesso o se invece conviene fare un unico kernel!!!    //
  //   ////////////////////////////////////////////////////////////////////////   //

  //LOOP SUI SITI PARI
  int xh, y, z, t;
#pragma acc kernels present(loc_s) present(loc_h)  present(approx) present(aux_u)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(xh=0; xh < nxh; xh++) {
          int idxh,idxpmu,x;
          int parity;
	  int dir_mu;
	  int mu;
	  x = 2*xh + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);  // r
	  //  parity = (x+y+z+t) % 2;
	  parity = 0; // la fisso cosi' perche' sto prendendo il sito pari

	  for(mu=0;mu<4;mu++){
	    idxpmu = nnp_openacc[idxh][mu][parity];// r+mu        
	    dir_mu = 2*mu +  parity;
	    vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_h,idxpmu,loc_s,idxh,approx->RA_a[iter]);
	  }//mu

        }  // x     
      }  // y       
    }  // z         
  }  // t

  //LOOP SUI SITI DISPARI
#pragma acc kernels present(loc_s) present(loc_h)  present(approx) present(aux_u)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(xh=0; xh < nxh; xh++) {
          int idxh,idxpmu,x;
          int parity;
	  int dir_mu;
	  int mu;
	  x = 2*xh + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);  // r
	  //  parity = (x+y+z+t) % 2;
	  parity = 1; // la fisso cosi' perche' sto prendendo il sito dispari
#pragma acc loop independent
	  for(mu=0;mu<4;mu++){
	    idxpmu = nnp_openacc[idxh][mu][parity];// r+mu        
	    dir_mu = 2*mu +  parity;
	    vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_s,idxpmu,loc_h,idxh,-approx->RA_a[iter]);
	  }//mu
        }  // x     
      }  // y       
    }  // z         
  }  // t
}// closes routine     

static inline void assign_zero_to_tamat_soa_component(__restrict tamat_soa * const matrix_comp,
						      int idx){
  matrix_comp->c01[idx]=0.0+I*0.0;
  matrix_comp->c02[idx]=0.0+I*0.0;
  matrix_comp->c12[idx]=0.0+I*0.0;
  matrix_comp->rc00[idx]=0.0;
  matrix_comp->rc11[idx]=0.0;
}

void set_tamat_soa_to_zero( __restrict tamat_soa * const matrix){
  int hx, y, z, t;
  int mu;
#pragma acc kernels present(matrix)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
	for(hx=0; hx < nxh; hx++) {
	  int x,idxh;
	  x = 2*hx + ((y+z+t) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  for(mu=0; mu<8; mu++) {
	    assign_zero_to_tamat_soa_component(&matrix[mu],idxh);
	  }
	}  // x
      }  // y
    }  // z
  }  // t
}



void multiply_conf_times_force_and_take_ta_even(__restrict su3_soa * const u, // la conf e' costante e non viene modificata
						__restrict ferm_param * const tpars,
						__restrict double_soa * const backfield,
						__restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
						__restrict tamat_soa * const ipdot){
  int hx,y,z,t,idxh;
#ifdef BACKFIELD
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars) present(backfield)
#else
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars)
#endif
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

          double arg;
          d_complex phase;
#ifdef IMCHEMPOT
          double imchempot = tpars->ferm_im_chem_pot/((double)(nt));
#endif

          //even sites
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);

          // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  eta = 1;
#ifdef BACKFIELD
	  arg = backfield[0].d[idxh] * tpars->ferm_charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh);
#endif


          // dir  2  =  y even
          eta = 1 - ( 2*(x & 0x1) );
#ifdef BACKFIELD
	  arg = backfield[2].d[idxh] * tpars->ferm_charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh);
#endif


          // dir  4  =  z even
          eta = 1 - ( 2*((x+y) & 0x1) );
#ifdef BACKFIELD
	  arg = backfield[4].d[idxh] * tpars->ferm_charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh);
#endif

          // dir  6  =  t even
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  arg = 0;
#ifdef BACKFIELD
	  arg += backfield[6].d[idxh] * tpars->ferm_charge;
#endif
#ifdef IMCHEMPOT
          arg += imchempot;
#endif
#ifdef PHASE_MAT_VEC_MULT
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh);
#endif

        }
      }
    }
  }

}



void multiply_conf_times_force_and_take_ta_odd(  __restrict su3_soa * const u, // e' costante e non viene modificata
						 __restrict ferm_param * const tpars,
						 __restrict double_soa * const backfield,
					         __restrict su3_soa * const auxmat, // e' costante e non viene modificata
					         __restrict tamat_soa * const ipdot){
  int hx,y,z,t,idxh;
#ifdef BACKFIELD
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars) present(backfield)
#else
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars)
#endif
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

          //odd sites
          x = 2*hx + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);
          // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  eta = 1;
#ifdef BACKFIELD
          arg = backfield[1].d[idxh] * tpars->ferm_charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh);
#endif

          // dir  3  =  y odd
          eta = 1 - ( 2*(x & 0x1) );
#ifdef BACKFIELD
          arg = backfield[3].d[idxh] * tpars->ferm_charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh);
#endif

          // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
#ifdef BACKFIELD
          arg = backfield[5].d[idxh] * tpars->ferm_charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh);
#endif

          // dir  7  =  t odd
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
          arg = 0;
#ifdef BACKFIELD
          arg += backfield[7].d[idxh] * tpars->ferm_charge;
#endif
#ifdef IMCHEMPOT
          arg += imchempot;
#endif
#ifdef PHASE_MAT_VEC_MULT
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh);
#endif


        }
      }
    }
  }

}




void ker_openacc_compute_fermion_force( __restrict su3_soa * const u, // e' costante e non viene mai modificato qui dentro
					__restrict double_soa * const backfield,
					__restrict su3_soa * const aux_u,
					__restrict vec3_soa * const in_shiftmulti,  // e' costante e non viene mai modificato qui dentro
					__restrict vec3_soa  * const loc_s,
					__restrict vec3_soa  * const loc_h,

					const RationalApprox * const approx  // e' costante e non viene mai modificato qui dentro
					){
  int ih;
  int iter=0;

  for(iter=0; iter<approx->approx_order; iter++){
      assign_in_to_out(&in_shiftmulti[iter],loc_s);
      acc_Doe(u,loc_h,loc_s,pars,backfield);
      direct_product_of_fermions_into_auxmat(loc_s,loc_h,aux_u,approx,iter);
    }
  }

 }


void fermion_force_soloopenacc(__restrict su3_soa  *conf_acc, // la configurazione qui dentro e' costante e non viene modificata
                   __restrict double_soa *backfield,
			       __restrict tamat_soa *ipdot_acc,
                   __restrict ferm_param *tpseudofermion_parameters,
                   int no_tot_ps,
			       __restrict vec3_soa *ferm_in_acc, // questo multifermione e' costante e non viene modificato
			       //__restrict ACC_MultiFermion *ferm_in_acc, // questo multifermione e' costante e non viene modificato
			       double res,
			       //const COM_RationalApprox *approx,
			       __restrict vec3_soa * ferm_out_acc,
			       //__restrict ACC_MultiFermion * ferm_out_acc,
			       __restrict su3_soa  * aux_conf_acc,
			       __restrict vec3_soa ** ferm_shiftmulti_acc,//parking variable
			       //__restrict ACC_ShiftMultiFermion * ferm_shiftmulti_acc,
			       __restrict vec3_soa * kloc_r,
			       __restrict vec3_soa * kloc_h,
			       __restrict vec3_soa * kloc_s,
			       __restrict vec3_soa * kloc_p,
			       __restrict vec3_soa *k_p_shiftferm//parking variable
                   ){
			       //__restrict ACC_ShiftFermion *k_p_shiftferm){

  printf("############################################ \n");
  printf("#### Inside fermion force soloopenacc ###### \n");
  printf("############################################ \n");

  struct timeval t1,t2;
  gettimeofday ( &t1, NULL );
 
  for(int ips = 0; ips < no_tot_ps; ips++) 
      multishift_invert(conf_acc,&tpseudofermion_parameters[ips],backfield,ferm_shiftmulti_acc[ips],&ferm_in_acc[ips],res,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);

  set_su3_soa_to_zero(aux_u);

  for(int ips = 0; ips < no_tot_ps; ips++) 
  ker_openacc_compute_fermion_force(conf_acc,backfield,aux_conf_acc,ferm_shiftmulti_acc[ips],kloc_s,kloc_h,&(tpseudofermion_parameters[approx]));
  set_tamat_soa_to_zero(ipdot_acc);
  multiply_conf_times_force_and_take_ta_even(conf_acc,aux_conf_acc,ipdot_acc);
  multiply_conf_times_force_and_take_ta_odd(conf_acc,aux_conf_acc,ipdot_acc);
  gettimeofday ( &t2, NULL );
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  printf("FULL FERMION FORCE COMPUTATION                  PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("########################################### \n");
  printf("#### Completed fermion force openacc ###### \n");
  printf("########################################### \n");

}


#endif

