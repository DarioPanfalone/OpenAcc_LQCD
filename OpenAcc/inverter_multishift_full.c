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
  //  printf("Allocated auxiliary variables \n");

  int pseudofermion,iter;
  double alpha, delta, lambda, omega, omega_save, gammag, fact;

  //  printf("Inside the kernel \n");
  //  printf("Ordine della approssimazione:   %i \n",approx[0].COM_approx_order);

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
    //    printf("delta    %.18lf\n",delta);
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

      //      printf("component_loc_s[0]    %.18e\n",creal(loc_s->c0[0]));
      //      printf("component_loc_p[0]    %.18e\n",creal(loc_p->c0[0]));

      alpha = real_scal_prod_global(loc_p,loc_s);
      //      printf("alpha    %.18lf\n",alpha);
      //      printf("mass2    %.18lf\n",mass2);

      omega_save=omega;   // omega_save=omega_(j-1)                                                                                                        
      omega=-delta/alpha;  // omega = (r_j,r_j)/(p_j, Ap_j)               
      //      printf("omega    %.18lf\n",omega);

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
      //      printf("Pseudoferm: %i   Iteration: %i    --> residue = %e   (target = %e) \n", pseudofermion, cg, sqrt(lambda), residuo);

      //      printf("Inside multishift  ( step = %i )\n");
      //      printf("lambda   %.18lf\n",lambda );
      //      printf("omega    %.18lf\n",omega );
      //      printf("gammag   %.18lf\n",gammag );

      //      for(iter=0; iter<(approx[0].COM_approx_order); iter++){
	//	printf("zeta_i[%i] =  %.18lf\n",iter,zeta_i[iter]);
	//	printf("gammas[%i] =  %.18lf\n",iter,gammas[iter]);
	//	printf("omegas[%i] =  %.18lf\n",iter,omegas[iter]);
      //      }


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


void direct_product_of_fermions_into_auxmat(__restrict vec3_soa  * const loc_s, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict vec3_soa  * const loc_h, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict su3_soa * const aux_u,
					    const COM_RationalApprox * const approx,
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
	    vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_h,idxpmu,loc_s,idxh,approx[0].COM_RA_a[iter]);
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
	    vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_s,idxpmu,loc_h,idxh,-approx[0].COM_RA_a[iter]);
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


/*
// Questa routine fa tutto quello che dovrebbero fare 
//    - multiply_conf_times_force_and_take_ta_even
//    - multiply_conf_times_force_and_take_ta_odd
// pero' NON COMPILA perche' l'istruzione "mat1_times_auxmat_into_tamat"
// e' ripetuta troppe volte (8). Invece quando la ripeto solo (6) volte compila... mah!!

// multiply the whole configuration for the staggered phases field
void multiply_conf_times_force_and_take_ta(const __restrict su3_soa * const u,
					   const __restrict su3_soa * const auxmat,
					   __restrict tamat_soa * const ipdot){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u) present(auxmat) present(ipdot)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

          //even sites
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);

          // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  eta = 1;
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh);

          // dir  2  =  y even
          eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh);

          // dir  4  =  z even
          eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh);

          // dir  6  =  t even
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_auxmat_into_tamat(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh);


          //odd sites
          x = 2*hx + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);
          // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  eta = 1;
	  mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh);

          // dir  3  =  y odd
          eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh);

          // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh);

          // dir  7  =  t odd
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_auxmat_into_tamat(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh);


        }
      }
    }
  }

}


*/


void multiply_conf_times_force_and_take_ta_even(__restrict su3_soa * const u, // la conf e' costante e non viene modificata
						__restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
						  __restrict tamat_soa * const ipdot){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u) present(auxmat) present(ipdot)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

          //even sites
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);

          // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  eta = 1;
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh);

          // dir  2  =  y even
          eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh);

          // dir  4  =  z even
          eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh);

          // dir  6  =  t even
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_auxmat_into_tamat(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh);

        }
      }
    }
  }

}



void multiply_conf_times_force_and_take_ta_odd(  __restrict su3_soa * const u, // e' costante e non viene modificata
					         __restrict su3_soa * const auxmat, // e' costante e non viene modificata
					         __restrict tamat_soa * const ipdot){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u) present(auxmat) present(ipdot)
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
	  mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh);

          // dir  3  =  y odd
          eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh);

          // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_auxmat_into_tamat(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh);

          // dir  7  =  t odd
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_auxmat_into_tamat(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh);


        }
      }
    }
  }

}




void ker_openacc_compute_fermion_force( __restrict su3_soa * const u, // e' costante e non viene mai modificato qui dentro
					__restrict su3_soa * const aux_u,
					__restrict ACC_ShiftMultiFermion * const in_shiftmulti,  // e' costante e non viene mai modificato qui dentro
					__restrict vec3_soa  * const loc_s,
					__restrict vec3_soa  * const loc_h,
					const COM_RationalApprox * const approx  // e' costante e non viene mai modificato qui dentro
					){
  int ips;
  int ih;
  int iter=0;
  //set_tamat_soa_to_zero(ipdot);
  set_su3_soa_to_zero(aux_u);

  for(ips=0; ips < no_ps; ips++){
    for(iter=0; iter<(approx[0].COM_approx_order); iter++){
      extract_from_shiftmulti_and_assign_to_fermion(in_shiftmulti,iter,ips,loc_s);
      acc_Doe(u,loc_h,loc_s);
      direct_product_of_fermions_into_auxmat(loc_s,loc_h,aux_u,approx,iter);
    }
  }
  //  multiply_conf_times_force_and_take_ta(u,aux_u,ipdot);

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
  int ish = 0;

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


void fermion_force_soloopenacc(__restrict su3_soa  *conf_acc, // la configurazione qui dentro e' costante e non viene modificata
			       __restrict tamat_soa *ipdot_acc,
			       __restrict ACC_MultiFermion *ferm_in_acc, // questo multifermione e' costante e non viene modificato
			       double res,
			       const COM_RationalApprox *approx,
			       __restrict ACC_MultiFermion * ferm_out_acc,
			       __restrict su3_soa  * aux_conf_acc,
			       __restrict ACC_ShiftMultiFermion * ferm_shiftmulti_acc,
			       __restrict vec3_soa * kloc_r,
			       __restrict vec3_soa * kloc_h,
			       __restrict vec3_soa * kloc_s,
			       __restrict vec3_soa * kloc_p,
			       __restrict ACC_ShiftFermion *k_p_shiftferm){
  printf("############################################ \n");
  printf("#### Inside fermion force soloopenacc ###### \n");
  printf("############################################ \n");

  struct timeval t1,t2;
  gettimeofday ( &t1, NULL );
  ker_invert_openacc_shiftmulti(conf_acc,ferm_shiftmulti_acc,ferm_in_acc,res,approx,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);
  ker_openacc_compute_fermion_force(conf_acc,aux_conf_acc,ferm_shiftmulti_acc,kloc_s,kloc_h,approx);
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


void fermion_force_openacc(const su3COM_soa  *conf, tamatCOM_soa *out, const COM_MultiFermion *in, double res, const COM_RationalApprox *approx){
  printf("######################################## \n");
  printf("#### Inside fermion force openacc ###### \n");
  printf("######################################## \n");
  printf("Residuo target = %f \n", res);


  tamat_soa * ipdot_acc;
  su3_soa  * conf_acc;
  su3_soa  * aux_conf_acc;
  ACC_MultiFermion * ferm_in_acc;
  ACC_MultiFermion * ferm_out_acc;
  ACC_ShiftMultiFermion * ferm_shiftmulti_acc;

  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&ipdot_acc, ALIGN, 8*sizeof(tamat_soa));
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
  int index = 19;

  struct timeval t0, t1,t2,t3;
  gettimeofday ( &t0, NULL );


  int dir;
  for(dir=0;dir<8;dir++)  convert_su3COM_soa_to_su3_soa(&conf[dir],&conf_acc[dir]);
  convert_COM_MultiFermion_to_ACC_MultiFermion(in,ferm_in_acc);
  printf("Converted conf and multi fermion \n");




#pragma acc data copyin(conf_acc[0:8]) copyout(aux_conf_acc[0:8]) copyin(ferm_in_acc[0:1]) copyin(approx[0:1])  copyout(ferm_out_acc[0:1])  create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1])  create(k_p_shiftferm[0:1]) create(ferm_shiftmulti_acc[0:1]) copy(ipdot_acc[0:8])
  {
  gettimeofday ( &t1, NULL );

    ker_invert_openacc_shiftmulti(conf_acc,ferm_shiftmulti_acc,ferm_in_acc,res,approx,kloc_r,kloc_h,kloc_s,kloc_p,k_p_shiftferm);

   // I think that for the force reconstruction this is not necessary
  //  ker_openacc_recombine_shiftmulti_to_multi(ferm_shiftmulti_acc,ferm_in_acc,ferm_out_acc,approx); 
   ker_openacc_compute_fermion_force(conf_acc,aux_conf_acc,ferm_shiftmulti_acc,kloc_s,kloc_h,approx);

   //   multiply_conf_times_force_and_take_ta(conf_acc,aux_conf_acc,ipdot_acc);
   multiply_conf_times_force_and_take_ta_even(conf_acc,aux_conf_acc,ipdot_acc);
   multiply_conf_times_force_and_take_ta_odd(conf_acc,aux_conf_acc,ipdot_acc);
  gettimeofday ( &t2, NULL );

  }



  printf("Ferm - Ipdot  OPENACC \n");
  printf("Ipdot 00 = ( %.18lf )\n", ipdot_acc[0].rc00[index]);
  printf("Ipdot 11 = ( %.18lf )\n", ipdot_acc[0].rc11[index]);
  printf("Ipdot 01 = ( %.18lf , %.18lf )\n", creal(ipdot_acc[0].c01[index]) , cimag(ipdot_acc[0].c01[index]));
  printf("Ipdot 02 = ( %.18lf , %.18lf )\n", creal(ipdot_acc[0].c02[index]) , cimag(ipdot_acc[0].c02[index]));
  printf("Ipdot 12 = ( %.18lf , %.18lf )\n", creal(ipdot_acc[0].c12[index]) , cimag(ipdot_acc[0].c12[index]));


  for(dir=0;dir<8;dir++)  convert_tamat_soa_to_tamatCOM_soa(&ipdot_acc[dir],&out[dir]);

  printf("Ipdot 00 = ( %.18lf )\n", out[0].rc00[index]);     
  printf("Ipdot 11 = ( %.18lf )\n", out[0].rc11[index]);     
  printf("Ipdot 01 = ( %.18lf , %.18lf )\n", out[0].c01[index].Re , out[0].c01[index].Im);     
  printf("Ipdot 02 = ( %.18lf , %.18lf )\n", out[0].c02[index].Re , out[0].c02[index].Im);     
  printf("Ipdot 12 = ( %.18lf , %.18lf )\n", out[0].c12[index].Re , out[0].c12[index].Im);     


  gettimeofday ( &t3, NULL );
  double dt_tot = (double)(t3.tv_sec - t0.tv_sec) + ((double)(t3.tv_usec - t0.tv_usec)/1.0e6);
  double dt_pretrans_to_preker = (double)(t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
  double dt_preker_to_postker = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec)/1.0e6);
  double dt_postker_to_posttrans = (double)(t3.tv_sec - t2.tv_sec) + ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);

  printf("FULL FERMION FORCE COMPUTATION                  Tot time          : %f sec  \n",dt_tot);
  printf("                                                PreTrans->Preker  : %f sec  \n",dt_pretrans_to_preker);
  printf("                                                PreKer->PostKer   : %f sec  \n",dt_preker_to_postker);
  printf("                                                PostKer->PostTrans: %f sec  \n",dt_postker_to_posttrans);
  
  
  
  printf("Reconverted fermions \n");


  free(conf_acc);
  free(aux_conf_acc);
  free(ipdot_acc);
  free(ferm_in_acc);
  free(ferm_out_acc);
  free(ferm_shiftmulti_acc);


  free(kloc_r);
  free(kloc_s);
  free(kloc_h);
  free(kloc_p);

  printf("Freed memory \n");

  printf("########################################### \n");
  printf("#### Completed fermion force openacc ###### \n");
  printf("########################################### \n");

}

#endif

