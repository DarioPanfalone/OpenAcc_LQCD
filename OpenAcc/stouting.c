#ifndef STOUTING_C
#define STOUTING_C

#include "./cayley_hamilton.c"

#pragma acc routine seq
static inline void su3_soa_to_single_su3(__restrict su3_soa * const in,
					int idx,
					single_su3 * out){
  out->comp[0][0] = in->r0.c0[idx];
  out->comp[0][1] = in->r0.c1[idx];
  out->comp[0][2] = in->r0.c2[idx];
  out->comp[1][0] = in->r1.c0[idx];
  out->comp[1][1] = in->r1.c1[idx];
  out->comp[1][2] = in->r1.c2[idx];
}

#pragma acc routine seq
static inline void tamat_soa_to_single_tamat( __restrict tamat_soa * const in,
					      int idx,
					      single_tamat * out){
  out->c01  = in->c01[idx];
  out->c02  = in->c02[idx];
  out->c12  = in->c12[idx];
  out->rc00 = in->rc00[idx];
  out->rc11 = in->rc11[idx];

}


static inline void conf_left_exp_multiply_to_su3_soa(__restrict su3_soa * const cnf,
						     const int idx,
						     __restrict su3_soa * const  EXP,
						     __restrict su3_soa * const cnf_out){

  single_su3 AUX;
  //Multiply: U_new = EXP * U_old
  AUX.comp[0][0] = cnf->r0.c0[idx];
  AUX.comp[0][1] = cnf->r0.c1[idx];
  AUX.comp[0][2] = cnf->r0.c2[idx];
  AUX.comp[1][0] = cnf->r1.c0[idx];
  AUX.comp[1][1] = cnf->r1.c1[idx];
  AUX.comp[1][2] = cnf->r1.c2[idx];
  // ricostruisco la terza
  AUX.comp[2][0] = conj(AUX.comp[0][1] * AUX.comp[1][2] - AUX.comp[0][2] * AUX.comp[1][1]);
  AUX.comp[2][1] = conj(AUX.comp[0][2] * AUX.comp[1][0] - AUX.comp[0][0] * AUX.comp[1][2]);
  AUX.comp[2][2] = conj(AUX.comp[0][0] * AUX.comp[1][1] - AUX.comp[0][1] * AUX.comp[1][0]);

  // moltiplica
  cnf_out->r0.c0[idx] = EXP->r0.c0[idx] * AUX.comp[0][0] + EXP->r0.c1[idx] * AUX.comp[1][0] + EXP->r0.c2[idx] * AUX.comp[2][0];
  cnf_out->r0.c1[idx] = EXP->r0.c0[idx] * AUX.comp[0][1] + EXP->r0.c1[idx] * AUX.comp[1][1] + EXP->r0.c2[idx] * AUX.comp[2][1];
  cnf_out->r0.c2[idx] = EXP->r0.c0[idx] * AUX.comp[0][2] + EXP->r0.c1[idx] * AUX.comp[1][2] + EXP->r0.c2[idx] * AUX.comp[2][2];

  cnf_out->r1.c0[idx] = EXP->r1.c0[idx] * AUX.comp[0][0] + EXP->r1.c1[idx] * AUX.comp[1][0] + EXP->r1.c2[idx] * AUX.comp[2][0];
  cnf_out->r1.c1[idx] = EXP->r1.c0[idx] * AUX.comp[0][1] + EXP->r1.c1[idx] * AUX.comp[1][1] + EXP->r1.c2[idx] * AUX.comp[2][1];
  cnf_out->r1.c2[idx] = EXP->r1.c0[idx] * AUX.comp[0][2] + EXP->r1.c1[idx] * AUX.comp[1][2] + EXP->r1.c2[idx] * AUX.comp[2][2];

}


void exp_minus_QA_times_conf(__restrict su3_soa * const tu,
			     __restrict tamat_soa * const QA,
			     __restrict su3_soa * const tu_out,
			     __restrict su3_soa * const exp_aux){
  
  int x, y, z, t;
#pragma acc kernels present(tu) present(QA) present(tu_out) present(exp_aux)
#pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector
        for(x=0; x < nx; x++) {
	  const  int idxh   = snum_acc(x,y,z,t);  // r
          const  int parity = (x+y+z+t) % 2;
          int dir_link;
          int mu;

          for(mu=0;mu<4;mu++){

            dir_link = 2*mu + parity;
	    //	    single_su3 tempU;
	    //	    single_su3 EXP;
	    //	    single_tamat tempQA;
	    //	    su3_soa_to_single_su3(&tu[dir_link],idxh,&tempU);
	    //	    tamat_soa_to_single_tamat(&QA[dir_link],idxh,&tempQA);

	    CH_exponential_antihermitian_soa(&exp_aux[dir_link],&QA[dir_link],idxh); // qui la variabile tu_out viene usata come output temporaneo
	    //	    CH_exponential_antihermitian(&EXP,&tempQA);
	    conf_left_exp_multiply_to_su3_soa(&tu[dir_link],idxh,&exp_aux[dir_link],&tu_out[dir_link]);
          }

        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine



void stout_isotropic( __restrict su3_soa * const u,
		      __restrict su3_soa * const uprime,
		      __restrict su3_soa * const local_staples,
		      __restrict su3_soa * const auxiliary,
		      __restrict tamat_soa * const tipdot){

  set_su3_soa_to_zero(local_staples);
  mult_conf_times_stag_phases(u);

  calc_loc_staples_removing_stag_phases_nnptrick_all(u,local_staples);

  RHO_times_conf_times_staples_ta_part(u,local_staples,tipdot);

  exp_minus_QA_times_conf(u,tipdot,uprime,auxiliary);

  mult_conf_times_stag_phases(u);
  mult_conf_times_stag_phases(uprime);

}




static inline d_complex  b1(double denom,
			    double u,
			    double w,
			    d_complex r_1,
			    d_complex r_2,
			    d_complex f){
  return  0.5*denom*denom*(2.0*u*r_1 + (3.0*u*u-w*w)*r_2 -2.0*(15.0*u*u+w*w)*f); // (57)
}

static inline d_complex  b2(double denom,
			    double u,
			    d_complex r_1,
			    d_complex r_2,
			    d_complex f){
  return   0.5*denom*denom*(r_1 - 3.0*u*r_2 -24.0*u*f); // (58)
}


//calcolo di lambda
static inline void compute_Lambda(__restrict thmat_soa * const L, // la Lambda --> ouput
				  __restrict tamat_soa * const SP, // Sigma primo --> input
				  __restrict su3_soa   * const U,    // la configurazione di gauge --> input
				  __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input
				  __restrict su3_soa   * const TMP,  // variabile di parcheggio
				  int idx
				  ){
  
  double c0 = det_i_times_QA_soa(QA,idx); //(14)
  double c1  = 0.5 * Tr_i_times_QA_sq_soa(QA,idx); // (15)
  //  double c0max = 2*pow(c1/3,1.5); // (17)
  double theta = homebrew_acos(c0/(2*pow(c1/3,1.5)));//(25)

  double u = sqrt(c1/3) * cos(theta/3) ;//(23)
  double w = sqrt(c1) * sin(theta/3) ;//(23)
  //  double xi0_A = 1 - w*w/6*(1-w*w/20*(1-w*w/42)); // (33 e seguenti)
  //  double xi0_B = sin(w)/w; // (33 e seguenti)

  double xi0 = (1 - w*w/6*(1-w*w/20*(1-w*w/42))) * (((int) (400*w*w-1) >> 31) & 0x1) +
    (sin(w)/w) * (((int) (1-400*w*w) >> 31) & 0x1) ; // (33 e seguenti)
  double xi1 = (cos(w)-xi0)/(w*w); // (67)

  d_complex expmiu = cos(u) - sin(u)*I;
  d_complex exp2iu = cos(2.0*u) + sin(2.0*u)*I;

  double denom = 1.0/(9*u*u - w*w);

  d_complex f0 =   denom * ((u*u - w*w) * exp2iu + expmiu*( // (30)
		  8*u*u *cos(w) + 2 * u * (3*u*u+ w*w) * xi0 * I ));
  d_complex f1 = denom * (2*u*exp2iu - expmiu* (  // (31)
		    2*u*cos(w) - (3*u*u-w*w)* xi0 * I )) ;
  d_complex f2 = denom * (exp2iu - expmiu* (cos(w)+ 3*u*xi0*I)); // (32)

  d_complex r0_1 = 2.0*(u+(u*u-w*w)*I)*exp2iu+2.0*expmiu*(4.0*u*(2.0-u*I)*cos(w) + (9*u*u+w*w-(3.0*u*u+w*w)*u*I)*xi0*I); //(60)
  d_complex r1_1 = 2.0*(1.0+2.0*u*I)*exp2iu+expmiu*(-2.0*(1.0-u*I)*cos(w) + (6.0*u+(w*w-3.0*u*u)*I)*xi0*I) ; // (61)
  d_complex r2_1 = (2.0*I)*exp2iu+(1.0*I)*expmiu*(cos(w)-3.0*(1.0-u*I)*xi0); // (62)

  d_complex r0_2 = -2.0*exp2iu+(2.0*I)*u*expmiu*(cos(w)+(1.0+4.0*u*I)*xi0+3.0*u*u*xi1); // (63)
  d_complex r1_2 = (-1.0*I)*expmiu*(cos(w)+(1.0+2.0*u*I)*xi0-3.0*u*u*xi1); // (64)
  d_complex r2_2 = expmiu*(xi0-3.0*u*xi1*I); // (65)

  d_complex b10 = 0.5*denom*denom*(2.0*u*r0_1 + (3.0*u*u-w*w)*r0_2 -2.0*(15.0*u*u+w*w)*f0); // (57)
  d_complex b11 = 0.5*denom*denom*(2.0*u*r1_1 + (3.0*u*u-w*w)*r1_2 -2.0*(15.0*u*u+w*w)*f1); // (57)
  d_complex b12 = 0.5*denom*denom*(2.0*u*r2_1 + (3.0*u*u-w*w)*r2_2 -2.0*(15.0*u*u+w*w)*f2); // (57)
  //  b10 = b1(denom,u,w,r0_1,r0_2,f0);
  //  b11 = b1(denom,u,w,r1_1,r1_2,f1);
  //  b12 = b1(denom,u,w,r2_1,r2_2,f2);

  //  d_complex b20 = 0.5*denom*denom*(r0_1 - 3.0*u*r0_2 -24.0*u*f0); // (58)
  //  d_complex b21 = 0.5*denom*denom*(r1_1 - 3.0*u*r1_2 -24.0*u*f1); // (58)
  //  d_complex b22 = 0.5*denom*denom*(r2_1 - 3.0*u*r2_2 -24.0*u*f2); // (58)
  //  b20 = b1(denom,u,r0_1,r0_2,f0);
  //  b21 = b1(denom,u,r1_1,r1_2,f1);
  //  b22 = b1(denom,u,r2_1,r2_2,f2);
  
  ////////////CALCOLO DI B1   (EQ 69) ////////////////////////////////////////
  TMP->r0.c0[idx] =b10 -   b11*QA->rc00[idx]                + b12*(QA->rc00[idx]*QA->rc00[idx]
								  + QA->c01[idx] *conj(QA->c01[idx])+ QA->c02[idx] *conj(QA->c02[idx]));
  TMP->r0.c1[idx] =      ( b11*I) *  QA->c01[idx]           + b12*(QA->c02[idx]*conj(QA->c12[idx])+(-1.0*I)*QA->c01[idx]*(QA->rc00[idx]+QA->rc11[idx]));
  TMP->r0.c2[idx] =      ( b11*I) * QA->c02[idx]            + b12*(-QA->c01[idx] * QA->c12[idx] + ( 1.0*I)* QA->c02[idx] * QA->rc11[idx]);
  /////////
  TMP->r1.c0[idx] =      (-b11*I) * conj(QA->c01[idx])     + b12*(QA->c12[idx]*conj(QA->c02[idx])+(1.0*I)*conj(QA->c01[idx])*(QA->rc00[idx]+QA->rc11[idx]));
  TMP->r1.c1[idx] =b10 -   b11*QA->rc11[idx]                + b12*( QA->rc11[idx]*QA->rc11[idx]
								   +QA->c01[idx]*conj(QA->c01[idx]) + QA->c12[idx] * conj(QA->c12[idx]));
  TMP->r1.c2[idx] =      ( b11*I)*QA->c12[idx]              + b12*((1.0*I)*QA->rc00[idx] * QA->c12[idx] + QA->c02[idx] * conj(QA->c01[idx]));
  /////////
  TMP->r2.c0[idx] =      (-b11*I)*conj(QA->c02[idx])        + b12*((-1.0*I)*QA->rc11[idx]*conj(QA->c02[idx]) - conj(QA->c01[idx]*QA->c12[idx]));
  TMP->r2.c1[idx] =      (-b11*I)*conj(QA->c12[idx])        + b12*((-1.0*I)*QA->rc00[idx]*conj(QA->c12[idx]) + QA->c01[idx]*conj(QA->c02[idx]));
  TMP->r2.c2[idx] =b10 +   b11*(QA->rc00[idx]+QA->rc11[idx])+ b12*((QA->rc00[idx]+QA->rc11[idx])*(QA->rc00[idx]+QA->rc11[idx])
								  + QA->c02[idx] * conj(QA->c02[idx])+ QA->c12[idx] * conj(QA->c12[idx]));
  //////////////////////////////////////////////////////////////////////////////

  //ricostruisco la terza riga della conf
  d_complex U20 = conj(U->r0.c1[idx] * U->r1.c2[idx] - U->r0.c2[idx] * U->r1.c1[idx]);
  d_complex U21 = conj(U->r0.c2[idx] * U->r1.c0[idx] - U->r0.c0[idx] * U->r1.c2[idx]);
  d_complex U22 = conj(U->r0.c0[idx] * U->r1.c1[idx] - U->r0.c1[idx] * U->r1.c0[idx]);

  /////////////////////////////////
  // CALCOLO DELLA TRACCIA => tr1 = Tr(Sigma' * B1 * U)
  d_complex tr1 =       SP->rc00[idx]   *  (TMP->r0.c0[idx]*U->r0.c0[idx] + TMP->r0.c1[idx]*U->r1.c0[idx] + TMP->r0.c2[idx]*U20)
                 +      SP->c01[idx]    *  (TMP->r1.c0[idx]*U->r0.c0[idx] + TMP->r1.c1[idx]*U->r1.c0[idx] + TMP->r1.c2[idx]*U20)
                 +      SP->c02[idx]    *  (TMP->r2.c0[idx]*U->r0.c0[idx] + TMP->r2.c1[idx]*U->r1.c0[idx] + TMP->r2.c2[idx]*U20)
                 - conj(SP->c01[idx])   *  (TMP->r0.c0[idx]*U->r0.c1[idx] + TMP->r0.c1[idx]*U->r1.c1[idx] + TMP->r0.c2[idx]*U21)
                 +      SP->rc11[idx]   *  (TMP->r1.c0[idx]*U->r0.c1[idx] + TMP->r1.c1[idx]*U->r1.c1[idx] + TMP->r1.c2[idx]*U21)
                 +      SP->c12[idx]    *  (TMP->r2.c0[idx]*U->r0.c1[idx] + TMP->r2.c1[idx]*U->r1.c1[idx] + TMP->r2.c2[idx]*U21)
                 - conj(SP->c02[idx])   *  (TMP->r0.c0[idx]*U->r0.c2[idx] + TMP->r0.c1[idx]*U->r1.c2[idx] + TMP->r0.c2[idx]*U22)
                 - conj(SP->c12[idx])   *  (TMP->r1.c0[idx]*U->r0.c2[idx] + TMP->r1.c1[idx]*U->r1.c2[idx] + TMP->r1.c2[idx]*U22)
       -(SP->rc00[idx]+SP->rc11[idx])   *  (TMP->r2.c0[idx]*U->r0.c2[idx] + TMP->r2.c1[idx]*U->r1.c2[idx] + TMP->r2.c2[idx]*U22);
  /////////////////////////////////


  /// questi sono i coefficienti b2j (58)
  b10 = 0.5*denom*denom*(r0_1 - 3.0*u*r0_2 -24.0*u*f0); // (58)
  b11 = 0.5*denom*denom*(r1_1 - 3.0*u*r1_2 -24.0*u*f1); // (58)
  b12 = 0.5*denom*denom*(r2_1 - 3.0*u*r2_2 -24.0*u*f2); // (58)
  ////////////CALCOLO DI B2   (EQ 69) ////////////////////////////////////////
  TMP->r0.c0[idx] =b10 -   b11*QA->rc00[idx]                + b12*(QA->rc00[idx]*QA->rc00[idx]
								  + QA->c01[idx] *conj(QA->c01[idx])+ QA->c02[idx] *conj(QA->c02[idx]));
  TMP->r0.c1[idx] =      ( b11*I) *  QA->c01[idx]           + b12*(QA->c02[idx]*conj(QA->c12[idx])+(-1.0*I)*QA->c01[idx]*(QA->rc00[idx]+QA->rc11[idx]));
  TMP->r0.c2[idx] =      ( b11*I) * QA->c02[idx]            + b12*(-QA->c01[idx] * QA->c12[idx] + ( 1.0*I)* QA->c02[idx] * QA->rc11[idx]);
  /////////
  TMP->r1.c0[idx] =      (-b11*I) * conj(QA->c01[idx])     + b12*(QA->c12[idx]*conj(QA->c02[idx])+(1.0*I)*conj(QA->c01[idx])*(QA->rc00[idx]+QA->rc11[idx]));
  TMP->r1.c1[idx] =b10 -   b11*QA->rc11[idx]                + b12*( QA->rc11[idx]*QA->rc11[idx]
								   +QA->c01[idx]*conj(QA->c01[idx]) + QA->c12[idx] * conj(QA->c12[idx]));
  TMP->r1.c2[idx] =      ( b11*I)*QA->c12[idx]              + b12*((1.0*I)*QA->rc00[idx] * QA->c12[idx] + QA->c02[idx] * conj(QA->c01[idx]));
  /////////
  TMP->r2.c0[idx] =      (-b11*I)*conj(QA->c02[idx])        + b12*((-1.0*I)*QA->rc11[idx]*conj(QA->c02[idx]) - conj(QA->c01[idx]*QA->c12[idx]));
  TMP->r2.c1[idx] =      (-b11*I)*conj(QA->c12[idx])        + b12*((-1.0*I)*QA->rc00[idx]*conj(QA->c12[idx]) + QA->c01[idx]*conj(QA->c02[idx]));
  TMP->r2.c2[idx] =b10 +   b11*(QA->rc00[idx]+QA->rc11[idx])+ b12*((QA->rc00[idx]+QA->rc11[idx])*(QA->rc00[idx]+QA->rc11[idx])
								  + QA->c02[idx] * conj(QA->c02[idx])+ QA->c12[idx] * conj(QA->c12[idx]));
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////
  // CALCOLO DELLA TRACCIA => tr2 = Tr(Sigma' * B2 * U)
  d_complex tr2 =       SP->rc00[idx]   *  (TMP->r0.c0[idx]*U->r0.c0[idx] + TMP->r0.c1[idx]*U->r1.c0[idx] + TMP->r0.c2[idx]*U20)
                 +      SP->c01[idx]    *  (TMP->r1.c0[idx]*U->r0.c0[idx] + TMP->r1.c1[idx]*U->r1.c0[idx] + TMP->r1.c2[idx]*U20)
                 +      SP->c02[idx]    *  (TMP->r2.c0[idx]*U->r0.c0[idx] + TMP->r2.c1[idx]*U->r1.c0[idx] + TMP->r2.c2[idx]*U20)
                 - conj(SP->c01[idx])   *  (TMP->r0.c0[idx]*U->r0.c1[idx] + TMP->r0.c1[idx]*U->r1.c1[idx] + TMP->r0.c2[idx]*U21)
                 +      SP->rc11[idx]   *  (TMP->r1.c0[idx]*U->r0.c1[idx] + TMP->r1.c1[idx]*U->r1.c1[idx] + TMP->r1.c2[idx]*U21)
                 +      SP->c12[idx]    *  (TMP->r2.c0[idx]*U->r0.c1[idx] + TMP->r2.c1[idx]*U->r1.c1[idx] + TMP->r2.c2[idx]*U21)
                 - conj(SP->c02[idx])   *  (TMP->r0.c0[idx]*U->r0.c2[idx] + TMP->r0.c1[idx]*U->r1.c2[idx] + TMP->r0.c2[idx]*U22)
                 - conj(SP->c12[idx])   *  (TMP->r1.c0[idx]*U->r0.c2[idx] + TMP->r1.c1[idx]*U->r1.c2[idx] + TMP->r1.c2[idx]*U22)
       -(SP->rc00[idx]+SP->rc11[idx])   *  (TMP->r2.c0[idx]*U->r0.c2[idx] + TMP->r2.c1[idx]*U->r1.c2[idx] + TMP->r2.c2[idx]*U22);
  /////////////////////////////////

  ////////////////CALCOLO U*Sigma' e lo metto in TMP che non serve piu'
  TMP->r0.c0[idx] = U->r0.c0[idx] * SP->rc00[idx]  - U->r0.c1[idx]*conj(SP->c01[idx]) - U->r0.c2[idx]*conj(SP->c02[idx]);
  TMP->r0.c1[idx] = U->r0.c0[idx] * SP->c01[idx]   + U->r0.c1[idx]*     SP->rc11[idx] - U->r0.c2[idx]*conj(SP->c12[idx]);
  TMP->r0.c2[idx] = U->r0.c0[idx] * SP->c02[idx]   + U->r0.c1[idx]*     SP->c12[idx]  - U->r0.c2[idx]*(SP->rc00[idx]+SP->rc11[idx]);

  TMP->r1.c0[idx] = U->r1.c0[idx] * SP->rc00[idx]  - U->r1.c1[idx]*conj(SP->c01[idx]) - U->r1.c2[idx]*conj(SP->c02[idx]);
  TMP->r1.c1[idx] = U->r1.c0[idx] * SP->c01[idx]   + U->r1.c1[idx]*     SP->rc11[idx] - U->r1.c2[idx]*conj(SP->c12[idx]);
  TMP->r1.c2[idx] = U->r1.c0[idx] * SP->c02[idx]   + U->r1.c1[idx]*     SP->c12[idx]  - U->r1.c2[idx]*(SP->rc00[idx]+SP->rc11[idx]);

  TMP->r2.c0[idx] = U20           * SP->rc00[idx]  - U21          *conj(SP->c01[idx]) - U22          *conj(SP->c02[idx]);
  TMP->r2.c1[idx] = U20           * SP->c01[idx]   + U21          *     SP->rc11[idx] - U22          *conj(SP->c12[idx]);
  TMP->r1.c2[idx] = U20           * SP->c02[idx]   + U21          *     SP->c12[idx]  - U22          *(SP->rc00[idx]+SP->rc11[idx]);
  /////////////////////////////////////////////////////////////////////


  ///////////////// CALCOLO DI GAMMA = tr1*Q + tr2*Q^2 + f1* (U * Sigma') + f2 * (Q * U * Sigma' + U * Sigma' * Q ) = 
  /////////////////                  = tr1*Q + tr2*Q^2 + f1* TMP + f2 * (Q * TMP + TMP * Q)
  /////////////////  GAMMA_00 = r0_1
  /////////////////  GAMMA_01 = r1_1
  /////////////////  GAMMA_02 = r2_1
  /////////////////  GAMMA_10 = r0_2
  /////////////////  GAMMA_11 = r1_2
  /////////////////  GAMMA_12 = r2_2
  /////////////////  GAMMA_20 = U20
  /////////////////  GAMMA_21 = U21
  /////////////////  GAMMA_22 = U22

  r0_1 = -tr1 * QA->rc00[idx]   + tr2 * (QA->rc00[idx]*QA->rc00[idx] + QA->c01[idx] *conj(QA->c01[idx])+ QA->c02[idx] *conj(QA->c02[idx]))
    + f1 * TMP->r0.c0[idx] + f2  * (-2.0*TMP->r0.c0[idx]*QA->rc00[idx]+(1.0*I)*(QA->c01[idx]*TMP->r1.c0[idx]-conj(QA->c01[idx])*TMP->r0.c1[idx]
										+QA->c02[idx]*TMP->r2.c0[idx]-conj(QA->c02[idx])*TMP->r0.c2[idx]));
  r1_1 =(tr1*I)*QA->c01[idx]    + tr2 * (QA->c02[idx]*conj(QA->c12[idx])+(-1.0*I)*QA->c01[idx]*(QA->rc00[idx]+QA->rc11[idx]))
    + f1 * TMP->r0.c1[idx] + f2  * (-TMP->r0.c1[idx]/*..........*/);






  ///////////////// INFINE CALCOLO DI LAMBDA = 0.5*(GAMMA + GAMMA^CROCE) - (1/6)* Id * Tr(GAMMA + GAMMA^CROCE)





}







#endif
