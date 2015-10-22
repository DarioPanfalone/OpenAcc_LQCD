#ifndef INVERTER_SHIFT_MULTI_FULL_C_
#define INVERTER_SHIFT_MULTI_FULL_C_

#include "./struct_c_def.c"

//#define DEBUG_INVERTER_SHIFT_MULTI_FULL_OPENACC

int multishift_invert(__restrict su3_soa * const u,
		      ferm_param * pars,
		      RationalApprox * approx,
		      double_soa * backfield,
		      __restrict vec3_soa * const out, // multi-fermion [nshifts]
		      __restrict vec3_soa * const in, // single ferm
		      double residuo,
		      __restrict vec3_soa * const loc_r,
		      __restrict vec3_soa * const loc_h,
		      __restrict vec3_soa * const loc_s,
		      __restrict vec3_soa * const loc_p,
		      __restrict vec3_soa * const shiftferm // multi-ferm [nshift]
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
//  printf("INSIDE MULTISHIFT_INVERT \n");
  int  cg;
  //double *zeta_i,*zeta_ii,*zeta_iii,*omegas,*gammas;
  // DYNAMIC ALLOCATION OF THESE SMALL ARRAYS SEEMS TO FAIL.
  double zeta_i[max_approx_order];
  double zeta_ii[max_approx_order];
  double zeta_iii[max_approx_order];
  double omegas[max_approx_order];
  double gammas[max_approx_order];
  int flag[max_approx_order];

//  printf("in pointer: %p\n", in);


  int iter;
  double alpha, delta, lambda, omega, omega_save, gammag, fact;
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

    //printf("Re in->c0[0]: %f\n"    ,creal(in->c0[0]) );
    //printf("Re loc_r->c0[0]: %f\n" ,creal(loc_r->c0[0]) );
    //printf("Re loc_p->c0[0]: %f\n" ,creal(loc_p->c0[0]) );

    for(iter=0; iter<(approx->approx_order); iter++){
      // ps_0=phi
      assign_in_to_out(in,&shiftferm[iter]);
      zeta_i[iter]=1.0;         // zeta_{-1}=1.0
      zeta_ii[iter]=1.0;        // zeta_{ 0}=1.0
      gammas[iter]=0.0;         // gammas_{-1}=0.0
    }
    gammag=0.0;
    cg=0;
 
    int maxiter; 
//Find the max approx->approx_order where flag[approx->approx_order] == 1
    for(iter=0; iter<(approx->approx_order); iter++) {
      if(flag[iter]==1) {
        maxiter = iter+1;
      }
    }
  
    do {      // loop over cg iterations
      cg++;

      // s=(M^dagM)p, alhpa=(p,s)=(p,Ap)
//fermion_matrix_multiplication(su3_soa *u,vec3_soa *out,vec3_soa *in,vec3_soa *temp1, ferm_param *pars,double_soa * backfield)
      fermion_matrix_multiplication(u,loc_s,loc_p,loc_h,pars,backfield);
      
      
      
//      printf("component_loc_s[0]    %f\n",creal(loc_s->c0[0]));
//      printf("component_loc_p[0]    %f\n",creal(loc_p->c0[0]));

      alpha = real_scal_prod_global(loc_p,loc_s);
//      printf("alpha    %.18lf\n",alpha);
//      printf("mass    %.18lf\n",pars->ferm_mass);

      omega_save=omega;   // omega_save=omega_(j-1)
      omega=-delta/alpha;  // omega = (r_j,r_j)/(p_j, Ap_j)               
//      printf("omega    %.18lf\n",omega);

      // out-=omegas*ps
      for(iter=0; iter<maxiter; iter++){
          if(flag[iter]==1){
              zeta_iii[iter] = (zeta_i[iter]*zeta_ii[iter]*omega_save)/
                  ( omega*gammag*(zeta_i[iter]-zeta_ii[iter])+
                    zeta_i[iter]*omega_save*(1.0-(approx->RA_b[iter])*omega) );

              omegas[iter]=omega*zeta_iii[iter]/zeta_ii[iter];
          }
      }

//multiple_combine_in1_minus_in2x_factor_back_into_in1(vec3_soa *out,vec3_soa *in,int maxiter,int *flag,double *omegas)
      multiple_combine_in1_minus_in2x_factor_back_into_in1(out,shiftferm,maxiter,flag,omegas);

      // r+=omega*s; lambda=(r,r)
      combine_add_factor_x_in2_to_in1(loc_r,loc_s,omega);
      lambda=l2norm2_global(loc_r);
      gammag=lambda/delta;

      // p=r+gammag*p
      combine_in1xfactor_plus_in2(loc_p,gammag,loc_r,loc_p);

      for(iter=0; iter<(approx->approx_order); iter++){
          if(flag[iter]==1){
              gammas[iter]=gammag*zeta_iii[iter]*omegas[iter]/(zeta_ii[iter]*omega);
          }
      }
      
//multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1(vec3_soa *in1, int maxiter,int *flag,double *gammas,vec3_soa *in2,double *zeta_iii )
      multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1(shiftferm,maxiter,flag,gammas,loc_r,zeta_iii);

      for(iter=0; iter<(approx->approx_order); iter++){
          if(flag[iter]==1){
              fact=sqrt(delta*zeta_ii[iter]*zeta_ii[iter]);
              if(fact<residuo) flag[iter]=0;
              else maxiter = iter+1;// modifying maxiter
              zeta_i[iter]=zeta_ii[iter];
              zeta_ii[iter]=zeta_iii[iter];
          }
      }
      delta=lambda;
//      printf("Iteration: %i    --> residue = %e   (target = %e) \n", cg, sqrt(lambda), residuo);

    } while(sqrt(lambda)>residuo && cg<max_cg); // end of cg iterations 

    if(cg==max_cg)
      {
	printf("WARNING: maximum number of iterations reached in invert\n");
      }
      printf("\t CG count = %i \n",cg);
    

#if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER_SHIFT_MULTI_FULL_OPENACC))
  printf("Terminated multishift_invert ( target res = %f ) \n ", residuo);
  int i;
      printf("\t CG count = %i \n",cg);
  // test 


    printf("\t\tnshift\tres/stop_res\tstop_res\tshift\n");
    for(iter=0; iter<approx->approx_order; iter++){
      assign_in_to_out(&out[iter],loc_p);
      fermion_matrix_multiplication_shifted(u,loc_s,loc_p,loc_h,pars,backfield,approx->RA_b[iter]);
      combine_in1_minus_in2(in,loc_s,loc_h); // r = s - y  
      double  giustoono=l2norm2_global(loc_h);
      printf("\t\t%i\t%e\t%1.1e\t\t%e\n",iter,sqrt(giustoono)/residuo,residuo,approx->RA_b[iter]);
    }
#endif

  return cg;

}

void recombine_shifted_vec3_to_vec3(const __restrict vec3_soa* const in_shifted /*multi-fermion*/, 
				    const __restrict vec3_soa* const in, // [nshift]
				    __restrict vec3_soa * const out, // [1] 
				    const RationalApprox * const approx ){
  int ih;
  int iter=0;
#pragma acc kernels present(out) present(in) present(in_shifted) present(approx)
#pragma acc loop independent
    for(ih=0; ih < sizeh; ih++){

      int ordine=approx->approx_order;
      out->c0[ih] =  in->c0[ih]*approx->RA_a0;
      out->c1[ih] =  in->c1[ih]*approx->RA_a0;
      out->c2[ih] =  in->c2[ih]*approx->RA_a0;

#pragma acc loop seq
      for(iter=0; iter<ordine; iter++){  // questo loop non lo vogliamo parallelizzare per forza ... forse puo andare bene cosi'
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


static inline  void mat1_times_auxmat_into_tamat(  __restrict su3_soa * const mat1, // e' costante e non viene modificato
						   const  int idx,
						   const  int eta,
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   const  int idx_aux,
						   __restrict tamat_soa * const ipdot,
						   const  int idipdot
#if defined(BACKFIELD) || defined(IMCHEMPOT)
						   ,d_complex phase
#endif
                           ){
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
  mat1_20 *= eta;
  mat1_21 *= eta;
  mat1_22 *= eta;

#if defined(BACKFIELD) || defined(IMCHEMPOT)
  mat1_00 *= phase;
  mat1_01 *= phase;
  mat1_02 *= phase;
  mat1_10 *= phase;
  mat1_11 *= phase;
  mat1_12 *= phase;
  mat1_20 *= phase;
  mat1_21 *= phase;
  mat1_22 *= phase;
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

}


#ifdef STOUT_FERMIONS
static inline  void mat1_times_auxmat_into_tamat_nophase(  __restrict su3_soa * const mat1, // e' costante e non viene modificato
						   const  int idx,
						   const  int eta,
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   const  int idx_aux,
						   __restrict tamat_soa * const ipdot,
						   const  int idipdot
                           ){
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
  mat1_20 *= eta;
  mat1_21 *= eta;
  mat1_22 *= eta;

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

}
#endif // ifdef STOUT_FERMIONS


#if defined(BACKFIELD) || defined(IMCHEMPOT)
static inline  void phase_times_auxmat_into_auxmat(
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   __restrict su3_soa * const pseudo_ipdot,
						   const  int idx,
                           d_complex phase
                           ){

    pseudo_ipdot->r0.c0[idx] += phase * auxmat->r0.c0[idx];
    pseudo_ipdot->r0.c1[idx] += phase * auxmat->r0.c1[idx];
    pseudo_ipdot->r0.c2[idx] += phase * auxmat->r0.c2[idx];
    pseudo_ipdot->r1.c0[idx] += phase * auxmat->r1.c0[idx];
    pseudo_ipdot->r1.c1[idx] += phase * auxmat->r1.c1[idx];
    pseudo_ipdot->r1.c2[idx] += phase * auxmat->r1.c2[idx];
    pseudo_ipdot->r2.c0[idx] += phase * auxmat->r2.c0[idx];
    pseudo_ipdot->r2.c1[idx] += phase * auxmat->r2.c1[idx];
    pseudo_ipdot->r2.c2[idx] += phase * auxmat->r2.c2[idx];


}

#else  //if defined(BACKFIELD) || defined(IMCHEMPOT)
static inline void accumulate_auxmat1_into_auxmat2(
						   __restrict su3_soa * const auxmat1,  // e' costante e non viene modificato
						   __restrict su3_soa * const auxmat2,
						   const  int idx
                           ){

    auxmat2->r0.c0[idx] += auxmat1->r0.c0[idx];
    auxmat2->r0.c1[idx] += auxmat1->r0.c1[idx];
    auxmat2->r0.c2[idx] += auxmat1->r0.c2[idx];
    auxmat2->r1.c0[idx] += auxmat1->r1.c0[idx];
    auxmat2->r1.c1[idx] += auxmat1->r1.c1[idx];
    auxmat2->r1.c2[idx] += auxmat1->r1.c2[idx];
    auxmat2->r2.c0[idx] += auxmat1->r2.c0[idx];
    auxmat2->r2.c1[idx] += auxmat1->r2.c1[idx];
    auxmat2->r2.c2[idx] += auxmat1->r2.c2[idx];

}

#endif  //if defined(BACKFIELD) || defined(IMCHEMPOT)




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
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

          double arg;
          d_complex phase;
#ifdef IMCHEMPOT
          double imchempot = tpars->ferm_im_chem_pot/((double)(nt));
#endif

#ifdef BACKFIELD
          double charge = (double)(tpars->ferm_charge);
#endif

          //even sites
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);

          // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  eta = 1;
#ifdef BACKFIELD
	  arg = backfield[0].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh);
#endif


          // dir  2  =  y even
          eta = 1 - ( 2*(x & 0x1) );
#ifdef BACKFIELD
	  arg = backfield[2].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh);
#endif


          // dir  4  =  z even
          eta = 1 - ( 2*((x+y) & 0x1) );
#ifdef BACKFIELD
	  arg = backfield[4].d[idxh] * charge;
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
	  arg += backfield[6].d[idxh] * charge;
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
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;
	  double arg;
	  d_complex phase;
#ifdef IMCHEMPOT
          double imchempot = tpars->ferm_im_chem_pot/((double)(nt));
#endif
#ifdef BACKFIELD
          double charge = (double)(tpars->ferm_charge);
#endif

          //odd sites
          x = 2*hx + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);
          // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  eta = 1;
#ifdef BACKFIELD
          arg = backfield[1].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh);
#endif

          // dir  3  =  y odd
          eta = 1 - ( 2*(x & 0x1) );
#ifdef BACKFIELD
          arg = backfield[3].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh);
#endif

          // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
#ifdef BACKFIELD
          arg = backfield[5].d[idxh] * charge;
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
          arg += backfield[7].d[idxh] * charge;
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


        } //hx
      } //y
    } // z
  } // t
} // end  multiply_conf_times_force_and_take_ta_odd()


#ifdef STOUT_FERMIONS
void multiply_conf_times_force_and_take_ta_even_nophase(__restrict su3_soa * const u, // la conf e' costante e non viene modificata
							__restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante; non viene modificata
							__restrict tamat_soa * const ipdot){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u) present(auxmat) present(ipdot) 
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
          for(hx=0; hx < nxh; hx++) {
              int x,eta;

              //even sites
              x = 2*hx + ((y+z+t) & 0x1);
              idxh = snum_acc(x,y,z,t);

              // dir  0  =  x even   --> eta = 1 , no multiplication needed
              eta = 1;
              mat1_times_auxmat_into_tamat_nophase(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh); 
              // dir  2  =  y even
              eta = 1 - ( 2*(x & 0x1) );
              mat1_times_auxmat_into_tamat_nophase(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh);
              // dir  4  =  z even
              eta = 1 - ( 2*((x+y) & 0x1) );
              mat1_times_auxmat_into_tamat_nophase(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh);

              // dir  6  =  t even
              eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
              eta *= (1- 2*(int)(t/(nt-1)));
#endif
              mat1_times_auxmat_into_tamat_nophase(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh);

          } // hx
      } // y
    } // z
  } // t
}


void multiply_conf_times_force_and_take_ta_odd_nophase(  __restrict su3_soa * const u, // e' costante e non viene modificata
							 __restrict su3_soa * const auxmat, // e' costante e non viene modificata
							 __restrict tamat_soa * const ipdot){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u) present(auxmat) present(ipdot)
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;
          //odd sites
          x = 2*hx + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);
          // dir  1  =  x odd    --> eta = 1 , no multiplication needed
          eta = 1;
          mat1_times_auxmat_into_tamat_nophase(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh);

          // dir  3  =  y odd
          eta = 1 - ( 2*(x & 0x1) );
          mat1_times_auxmat_into_tamat_nophase(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh);

          // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_auxmat_into_tamat_nophase(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh);

          // dir  7  =  t odd
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
          mat1_times_auxmat_into_tamat_nophase(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh);

        } //hx
      } //y
    } // z
  } // t
} // end  multiply_conf_times_force_and_take_ta_odd()
#endif //ifdef STOUT_FERMIONS


#if defined(IMCHEMPOT) || defined(BACKFIELD)
void multiply_backfield_times_force(__restrict ferm_param * const tpars,
        __restrict double_soa * const backfield,
        __restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
        __restrict su3_soa * const pseudo_ipdot){

  double arg;
  d_complex phase;
#ifdef IMCHEMPOT
  double imchempot = tpars->ferm_im_chem_pot/((double)(nt));
#endif
  double charge = (double)(tpars->ferm_charge);
  int idxh;
#pragma acc data present(backfield) present(auxmat) present(pseudo_ipdot)
#pragma acc loop independent
  for(int dirindex = 0 ; dirindex < 8 ; dirindex++){
#pragma acc loop independent
    for( idxh = 0 ; idxh < sizeh; idxh++){
      
#ifdef BACKFIELD
      arg = backfield[dirindex].d[idxh] * charge;
#else
      arg = 0.0;
#endif

#ifdef IMCHEMPOT
      if(dirindex>6) arg += imchempot;
#endif
      phase = cos(arg) + I * sin(arg);
      phase_times_auxmat_into_auxmat(&auxmat[dirindex],&pseudo_ipdot[dirindex],idxh,phase);
    }
  }
} // end multiply_backfield_times_force()

#else 
void accumulate_gl3soa_into_gl3soa(
        __restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
        __restrict su3_soa * const pseudo_ipdot){

    int idxh;
#pragma acc data present(auxmat) present(pseudo_ipdot)
#pragma acc loop independent 
    for(int dirindex = 0 ; dirindex < 8 ; dirindex++){
#pragma acc loop independent
        for( idxh = 0 ; idxh < sizeh; idxh++){
            accumulate_auxmat1_into_auxmat2(&auxmat[dirindex],&pseudo_ipdot[dirindex],idxh);
        }
    }
}

#endif

void ker_openacc_compute_fermion_force( __restrict su3_soa * const u, // e' costante e non viene mai modificato qui dentro
					double_soa * backfield,
					__restrict su3_soa * const aux_u,
					__restrict vec3_soa * const in_shiftmulti,  // e' costante e non viene mai modificato qui dentro
					__restrict vec3_soa  * const loc_s,
					__restrict vec3_soa  * const loc_h,
					ferm_param  *  tpars
					){
  int ih;
  int iter=0;

  for(iter=0; iter<tpars->approx_md.approx_order; iter++){
    assign_in_to_out(&in_shiftmulti[iter],loc_s);
    acc_Doe(u,loc_h,loc_s,tpars,backfield);
    direct_product_of_fermions_into_auxmat(loc_s,loc_h,aux_u,&(tpars->approx_md),iter);
  }
}

#endif

