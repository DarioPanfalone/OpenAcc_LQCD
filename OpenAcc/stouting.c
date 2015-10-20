#ifndef STOUTING_C
#define STOUTING_C

#include "./cayley_hamilton.c"
#include "./struct_c_def.c"

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
	    CH_exponential_antihermitian_soa_nissalike(&exp_aux[dir_link],&QA[dir_link],idxh);
	    conf_left_exp_multiply_to_su3_soa(&tu[dir_link],idxh,&exp_aux[dir_link],&tu_out[dir_link]);
          }

        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine



void stout_isotropic( __restrict su3_soa * const u,               // --> input conf
		      __restrict su3_soa * const uprime,          // --> output conf [stouted]
		      __restrict su3_soa * const local_staples,   // --> parking variable
		      __restrict su3_soa * const auxiliary,       // --> parking variable
		      __restrict tamat_soa * const tipdot){       // --> parking variable

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
#pragma acc routine seq
static inline void compute_loc_Lambda(__restrict thmat_soa * const L, // la Lambda --> ouput
				      __restrict su3_soa   * const SP, // Sigma primo --> input
				      __restrict su3_soa   * const U,    // la configurazione di gauge --> input
				      __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input
				      __restrict su3_soa   * const TMP,  // variabile di parcheggio
				      int idx
				      ){
  
  double c0 = det_i_times_QA_soa(QA,idx);
  double c1  = 0.5 * Tr_i_times_QA_sq_soa(QA,idx);
  double c0max = 2*pow(c1/3,1.5);
  d_complex f0,f1,f2;
  d_complex b10,b11,b12,b20,b21,b22;
  d_complex  r0_1,r1_1,r2_1,r0_2,r1_2,r2_2;
  if(c1<4e-3)
    {
      f0 = (1-c0*c0/720) + (1.0*I)*(-c0*(1-c1*(1-c1/42)/20)/6);
      f1 = (c0*(1-c1*(1-3*c1/112)/15)/24) + (1.0*I)*(1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040);
      f2 = (0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160)) + (1.0*I)*(0.5*(c0*(1-c1*(1-c1/48)/21)/60));

      //from nissa
      /* 
      b[0][0][RE]=0;
      b[0][0][IM]=c0/120*(1-c1/21);
      b[0][1][RE]=-c0/360*(1-3*c1/56);
      b[0][1][IM]=-1.0/6*(1-c1/10*(1.0-c1/28));
      b[0][2][RE]=0.5*(1.0/12*(1-2*c1/30*(1-3*c1/112)));
      b[0][2][IM]=0.5*(-c0/1260*(1-c1/24));

      b[1][0][RE]=-c0/360;
      b[1][0][IM]=-1.0/6*(1-c1/20*(1-c1/42));
      b[1][1][RE]=1.0/24*(1-c1/15*(1-3*c1/112));
      b[1][1][IM]=-c0/2520;
      b[1][2][RE]=0.5*c0/10080;
      b[1][2][IM]=0.5*(1.0/60*(1-c1/21*(1-c1/48)));
      */

      b10 = 0.0 + (1.0*I)*(c0/120*(1-c1/21));
      b11 = (-c0/360*(1-3*c1/56)) + (1.0*I)*(-1.0/6*(1-c1/10*(1.0-c1/28)));
      b12 = (0.5*(1.0/12*(1-2*c1/30*(1-3*c1/112)))) + (1.0*I)*(0.5*(-c0/1260*(1-c1/24)));

      b20 = (-c0/360) + (1.0*I)*(-1.0/6*(1-c1/20*(1-c1/42)));
      b21 = (1.0/24*(1-c1/15*(1-3*c1/112))) + (1.0*I)*(-c0/2520);
      b22 = (0.5*c0/10080) + (1.0*I)*(0.5*(1.0/60*(1-c1/21*(1-c1/48))));

    }
  else
    {
      int segno=1;
      if(c0<0)
        {
          segno=-1;
          c0=-c0;
        }

      double eps=(c0max-c0)/c0max;
      double theta;
      if(eps<0) theta=0.0;
      else
        if(eps<1e-3) theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
        else theta=acos(c0/c0max);
      double u = sqrt(c1/3) * cos(theta/3) ;
      double w = sqrt(c1) * sin(theta/3) ;
      double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
      double cu=cos(u),c2u=cos(2*u);
      double su=sin(u),s2u=sin(2*u);
      double cw=cos(w);
      double xi0w;
      if(fabs(w)<0.05)
        {
          double temp0=w*w,temp1=1-temp0/42,temp2=1.0-temp0/20*temp1;
          xi0w=1-temp0/6*temp2;
        }
      else xi0w=sin(w)/w;
      double xi1w; //eq. (67)
      if(fabs(w)<0.05) xi1w=-(1-w2*(1-w2*(1-w2/54)/28)/10)/3;
      else xi1w=cw/w2-sin(w)/(w2*w);
      
      double denom = 1/(9*u*u - w*w);
      
      f0=(u2mw2*c2u+ //(u2-w2)*cos(2u)
	  cu*8*u2*cw+ //cos(u)*8*u2*cos(w)
	  2*su*u*w2p3u2*xi0w) //sin(u)*2*mu*(3*u2+w2)*xi0(w)
        +(1.0*I)*(u2mw2*s2u+ //(u2-w2)*sin(2u)
		  -su*8*u2*cw+ //-sin(u)*8*u2*cos(w)
		  cu*2*u*w2p3u2*xi0w); //cos(u)*2*u*(3*u2+w2)*xi0(w)
      f0 *= denom;
      f1=(2*u*c2u+ //2*u*cos(2u)
	  -cu*2*u*cw+ //cos(u)*2*u*cos(w)
	  -su*w2m3u2*xi0w) //sin(u)*(u2-3*w2)*xi0(w)
        +(1.0*I)*(2*u*s2u+ //2*u*sin(2u)
		  su*2*u*cos(w)+ //sin(u)*2*u*cos(w)
                  -cu*w2m3u2*xi0w);//cos(u)*(3*u2-w2)*xi0(w)
      f1 *= denom;
      f2=(c2u+ //cos(2u)
	  -cu*cw+ //-cos(u)*cos(w)
          -3*su*u*xi0w) //-3*sin(u)*u*xi0(w)
        +(1.0*I)*(s2u+ //sin(2u)
		  su*cw+ //sin(w)*cos(w)
                  -cu*3*u*xi0w);//-cos(u)*3*u*xi0(w)
      f2 *= denom;

      //from nissa
      /*
       complex r[2][3]=
	 {{{2*c2u*u+s2u*(-2*u2+2*w2)+2*cu*u*(8*cw+3*u2*xi0w+w2*xi0w)+su*(-8*cw*u2+18*u2*xi0w+2*w2*xi0w),
	    -8*cw*(2*su*u+cu*u2)+2*(s2u*u+c2u*u2-c2u*w2)+2*(9*cu*u2-3*su*u*u2+cu*w2 -su*u*w2)*xi0w},
	   {2*c2u-4*s2u*u+su*(2*cw*u+6*u*xi0w)+cu*(-2*cw+3*u2*xi0w-w2*xi0w),
	    2*s2u+4*c2u*u+2*cw*(su+cu*u)+(6*cu*u-3*su*u2+su*w2)*xi0w},
	   {-2*s2u+cw*su-3*(su+cu*u)*xi0w,
	    2*c2u+cu*cw+(-3*cu+3*su*u)*xi0w}},

	  {{-2*c2u+2*cw*su*u+2*su*u*xi0w-8*cu*u2*xi0w+6*su*u*u2*xi1w,
	    2*(-s2u+4*su*u2*xi0w+cu*u*(cw+xi0w+3*u2*xi1w))},
	   {2*cu*u*xi0w+su*(-cw-xi0w+3*u2*xi1w),
	    -2*su*u*xi0w-cu*(cw+xi0w-3*u2*xi1w)},
	   {cu*xi0w-3*su*u*xi1w,
	    -(su*xi0w)-3*cu*u*xi1w}}};
      */

      r0_1 = (2*c2u*u+s2u*(-2*u2+2*w2)+2*cu*u*(8*cw+3*u2*xi0w+w2*xi0w)+su*(-8*cw*u2+18*u2*xi0w+2*w2*xi0w))
	+ (1.0*I)*(-8*cw*(2*su*u+cu*u2)+2*(s2u*u+c2u*u2-c2u*w2)+2*(9*cu*u2-3*su*u*u2+cu*w2 -su*u*w2)*xi0w);
      r1_1 = (2*c2u-4*s2u*u+su*(2*cw*u+6*u*xi0w)+cu*(-2*cw+3*u2*xi0w-w2*xi0w))
	+ (1.0*I)*(2*s2u+4*c2u*u+2*cw*(su+cu*u)+(6*cu*u-3*su*u2+su*w2)*xi0w);
      r2_1 = (-2*s2u+cw*su-3*(su+cu*u)*xi0w)
	+ (1.0*I)*(2*c2u+cu*cw+(-3*cu+3*su*u)*xi0w);

      r0_2 = (-2*c2u+2*cw*su*u+2*su*u*xi0w-8*cu*u2*xi0w+6*su*u*u2*xi1w)
	+ (1.0*I)*(2*(-s2u+4*su*u2*xi0w+cu*u*(cw+xi0w+3*u2*xi1w)));
      r1_2 = (2*cu*u*xi0w+su*(-cw-xi0w+3*u2*xi1w))
	+ (1.0*I)*(-2*su*u*xi0w-cu*(cw+xi0w-3*u2*xi1w));
      r2_2 = (cu*xi0w-3*su*u*xi1w)
	+ (1.0*I)*(-(su*xi0w)-3*cu*u*xi1w);

      b10 = 0.5*denom*denom*(2.0*u*r0_1 + (3.0*u*u-w*w)*r0_2 -2.0*(15.0*u*u+w*w)*f0); // (57) 
      b11 = 0.5*denom*denom*(2.0*u*r1_1 + (3.0*u*u-w*w)*r1_2 -2.0*(15.0*u*u+w*w)*f1); // (57)
      b12 = 0.5*denom*denom*(2.0*u*r2_1 + (3.0*u*u-w*w)*r2_2 -2.0*(15.0*u*u+w*w)*f2); // (57)

      b20 = 0.5*denom*denom*(r0_1 - 3.0*u*r0_2 -24.0*u*f0); // (58)
      b21 = 0.5*denom*denom*(r1_1 - 3.0*u*r1_2 -24.0*u*f1); // (58)
      b22 = 0.5*denom*denom*(r2_1 - 3.0*u*r2_2 -24.0*u*f2); // (58)

      if(segno==-1){
	b10 = +conj(b10);
	b11 = -conj(b11);
	b12 = +conj(b12);
	b20 = -conj(b20);
	b21 = +conj(b21);
	b22 = -conj(b22);
	
	
        f0 =  conj(f0);
        f1 = -conj(f1);
        f2 =  conj(f2);
      }
    }

  //ricostruisco la terza riga della conf
  d_complex U20 = conj(U->r0.c1[idx] * U->r1.c2[idx] - U->r0.c2[idx] * U->r1.c1[idx]);
  d_complex U21 = conj(U->r0.c2[idx] * U->r1.c0[idx] - U->r0.c0[idx] * U->r1.c2[idx]);
  d_complex U22 = conj(U->r0.c0[idx] * U->r1.c1[idx] - U->r0.c1[idx] * U->r1.c0[idx]);

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

  /////////////////////////////////
  // CALCOLO DELLA TRACCIA => tr1 = Tr(Sigma' * B1 * U)
  d_complex tr1 =    SP->r0.c0[idx]   *  (TMP->r0.c0[idx]*U->r0.c0[idx] + TMP->r0.c1[idx]*U->r1.c0[idx] + TMP->r0.c2[idx]*U20)
                 +   SP->r0.c1[idx]   *  (TMP->r1.c0[idx]*U->r0.c0[idx] + TMP->r1.c1[idx]*U->r1.c0[idx] + TMP->r1.c2[idx]*U20)
                 +   SP->r0.c2[idx]   *  (TMP->r2.c0[idx]*U->r0.c0[idx] + TMP->r2.c1[idx]*U->r1.c0[idx] + TMP->r2.c2[idx]*U20)
                 +   SP->r1.c0[idx]   *  (TMP->r0.c0[idx]*U->r0.c1[idx] + TMP->r0.c1[idx]*U->r1.c1[idx] + TMP->r0.c2[idx]*U21)
                 +   SP->r1.c1[idx]   *  (TMP->r1.c0[idx]*U->r0.c1[idx] + TMP->r1.c1[idx]*U->r1.c1[idx] + TMP->r1.c2[idx]*U21)
                 +   SP->r1.c2[idx]   *  (TMP->r2.c0[idx]*U->r0.c1[idx] + TMP->r2.c1[idx]*U->r1.c1[idx] + TMP->r2.c2[idx]*U21)
                 +   SP->r2.c0[idx]   *  (TMP->r0.c0[idx]*U->r0.c2[idx] + TMP->r0.c1[idx]*U->r1.c2[idx] + TMP->r0.c2[idx]*U22)
                 +   SP->r2.c1[idx]   *  (TMP->r1.c0[idx]*U->r0.c2[idx] + TMP->r1.c1[idx]*U->r1.c2[idx] + TMP->r1.c2[idx]*U22)
                 +   SP->r2.c2[idx]   *  (TMP->r2.c0[idx]*U->r0.c2[idx] + TMP->r2.c1[idx]*U->r1.c2[idx] + TMP->r2.c2[idx]*U22);
  /////////////////////////////////


  ////////////CALCOLO DI B2   (EQ 69) ////////////////////////////////////////
  TMP->r0.c0[idx] =b20 -   b21*QA->rc00[idx]                + b22*(QA->rc00[idx]*QA->rc00[idx]
								  + QA->c01[idx] *conj(QA->c01[idx])+ QA->c02[idx] *conj(QA->c02[idx]));
  TMP->r0.c1[idx] =      ( b21*I) *  QA->c01[idx]           + b22*(QA->c02[idx]*conj(QA->c12[idx])+(-1.0*I)*QA->c01[idx]*(QA->rc00[idx]+QA->rc11[idx]));
  TMP->r0.c2[idx] =      ( b21*I) * QA->c02[idx]            + b22*(-QA->c01[idx] * QA->c12[idx] + ( 1.0*I)* QA->c02[idx] * QA->rc11[idx]);
  /////////
  TMP->r1.c0[idx] =      (-b21*I) * conj(QA->c01[idx])     + b22*(QA->c12[idx]*conj(QA->c02[idx])+(1.0*I)*conj(QA->c01[idx])*(QA->rc00[idx]+QA->rc11[idx]));
  TMP->r1.c1[idx] =b20 -   b21*QA->rc11[idx]                + b22*( QA->rc11[idx]*QA->rc11[idx]
								   +QA->c01[idx]*conj(QA->c01[idx]) + QA->c12[idx] * conj(QA->c12[idx]));
  TMP->r1.c2[idx] =      ( b21*I)*QA->c12[idx]              + b22*((1.0*I)*QA->rc00[idx] * QA->c12[idx] + QA->c02[idx] * conj(QA->c01[idx]));
  /////////
  TMP->r2.c0[idx] =      (-b21*I)*conj(QA->c02[idx])        + b22*((-1.0*I)*QA->rc11[idx]*conj(QA->c02[idx]) - conj(QA->c01[idx]*QA->c12[idx]));
  TMP->r2.c1[idx] =      (-b21*I)*conj(QA->c12[idx])        + b22*((-1.0*I)*QA->rc00[idx]*conj(QA->c12[idx]) + QA->c01[idx]*conj(QA->c02[idx]));
  TMP->r2.c2[idx] =b20 +   b21*(QA->rc00[idx]+QA->rc11[idx])+ b22*((QA->rc00[idx]+QA->rc11[idx])*(QA->rc00[idx]+QA->rc11[idx])
								  + QA->c02[idx] * conj(QA->c02[idx])+ QA->c12[idx] * conj(QA->c12[idx]));
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////
  // CALCOLO DELLA TRACCIA => tr2 = Tr(Sigma' * B2 * U)
  d_complex tr2 =    SP->r0.c0[idx]   *  (TMP->r0.c0[idx]*U->r0.c0[idx] + TMP->r0.c1[idx]*U->r1.c0[idx] + TMP->r0.c2[idx]*U20)
                 +   SP->r0.c1[idx]   *  (TMP->r1.c0[idx]*U->r0.c0[idx] + TMP->r1.c1[idx]*U->r1.c0[idx] + TMP->r1.c2[idx]*U20)
                 +   SP->r0.c2[idx]   *  (TMP->r2.c0[idx]*U->r0.c0[idx] + TMP->r2.c1[idx]*U->r1.c0[idx] + TMP->r2.c2[idx]*U20)
                 +   SP->r1.c0[idx]   *  (TMP->r0.c0[idx]*U->r0.c1[idx] + TMP->r0.c1[idx]*U->r1.c1[idx] + TMP->r0.c2[idx]*U21)
                 +   SP->r1.c1[idx]   *  (TMP->r1.c0[idx]*U->r0.c1[idx] + TMP->r1.c1[idx]*U->r1.c1[idx] + TMP->r1.c2[idx]*U21)
                 +   SP->r1.c2[idx]   *  (TMP->r2.c0[idx]*U->r0.c1[idx] + TMP->r2.c1[idx]*U->r1.c1[idx] + TMP->r2.c2[idx]*U21)
                 +   SP->r2.c0[idx]   *  (TMP->r0.c0[idx]*U->r0.c2[idx] + TMP->r0.c1[idx]*U->r1.c2[idx] + TMP->r0.c2[idx]*U22)
                 +   SP->r2.c1[idx]   *  (TMP->r1.c0[idx]*U->r0.c2[idx] + TMP->r1.c1[idx]*U->r1.c2[idx] + TMP->r1.c2[idx]*U22)
                 +   SP->r2.c2[idx]   *  (TMP->r2.c0[idx]*U->r0.c2[idx] + TMP->r2.c1[idx]*U->r1.c2[idx] + TMP->r2.c2[idx]*U22);
  /////////////////////////////////

  ////////////////CALCOLO U*Sigma' e lo metto in TMP che non serve piu'
  TMP->r0.c0[idx] = U->r0.c0[idx] * SP->r0.c0[idx]  - U->r0.c1[idx] * SP->r1.c0[idx] - U->r0.c2[idx] * SP->r2.c0[idx];
  TMP->r0.c1[idx] = U->r0.c0[idx] * SP->r0.c1[idx]  - U->r0.c1[idx] * SP->r1.c1[idx] - U->r0.c2[idx] * SP->r2.c1[idx];
  TMP->r0.c2[idx] = U->r0.c0[idx] * SP->r0.c2[idx]  - U->r0.c1[idx] * SP->r1.c2[idx] - U->r0.c2[idx] * SP->r2.c2[idx];

  TMP->r1.c0[idx] = U->r1.c0[idx] * SP->r0.c0[idx]  - U->r1.c1[idx] * SP->r1.c0[idx] - U->r1.c2[idx] * SP->r2.c0[idx];
  TMP->r1.c1[idx] = U->r1.c0[idx] * SP->r0.c1[idx]  - U->r1.c1[idx] * SP->r1.c1[idx] - U->r1.c2[idx] * SP->r2.c1[idx];
  TMP->r1.c2[idx] = U->r1.c0[idx] * SP->r0.c2[idx]  - U->r1.c1[idx] * SP->r1.c2[idx] - U->r1.c2[idx] * SP->r2.c2[idx];

  TMP->r2.c0[idx] = U->r2.c0[idx] * SP->r0.c0[idx]  - U->r2.c1[idx] * SP->r1.c0[idx] - U->r2.c2[idx] * SP->r2.c0[idx];
  TMP->r2.c1[idx] = U->r2.c0[idx] * SP->r0.c1[idx]  - U->r2.c1[idx] * SP->r1.c1[idx] - U->r2.c2[idx] * SP->r2.c1[idx];
  TMP->r2.c2[idx] = U->r2.c0[idx] * SP->r0.c2[idx]  - U->r2.c1[idx] * SP->r1.c2[idx] - U->r2.c2[idx] * SP->r2.c2[idx];
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
  // prima riga
  r0_1 = -tr1 * QA->rc00[idx]   + tr2 * (QA->rc00[idx]*QA->rc00[idx] + QA->c01[idx] *conj(QA->c01[idx])+ QA->c02[idx] *conj(QA->c02[idx]))
    + f1 * TMP->r0.c0[idx] + f2  * (-2.0*TMP->r0.c0[idx]*QA->rc00[idx]+(1.0*I)*(QA->c01[idx]*TMP->r1.c0[idx]-conj(QA->c01[idx])*TMP->r0.c1[idx]
										+QA->c02[idx]*TMP->r2.c0[idx]-conj(QA->c02[idx])*TMP->r0.c2[idx]));
  r1_1 =(tr1*I)*QA->c01[idx]    + tr2 * (QA->c02[idx]*conj(QA->c12[idx])+(-1.0*I)*QA->c01[idx]*(QA->rc00[idx]+QA->rc11[idx]))
    + f1 * TMP->r0.c1[idx] + f2  * (-TMP->r0.c1[idx]*(QA->rc11[idx]+QA->rc00[idx])+(1.0*I)*(QA->c01[idx]*(TMP->r0.c0[idx]+TMP->r1.c1[idx])
									    +QA->c02[idx]*TMP->r2.c1[idx]-conj(QA->c12[idx])*TMP->r0.c2[idx]));
  r2_1 =(tr1*I)*QA->c02[idx]    + tr2 * (-QA->c01[idx] * QA->c12[idx] + ( 1.0*I)* QA->c02[idx] * QA->rc11[idx])
    + f1 * TMP->r0.c2[idx] + f2  * (QA->rc11[idx]*TMP->r0.c2[idx]+(1.0*I)*(QA->c02[idx]*(TMP->r0.c0[idx]+TMP->r2.c2[idx])+
									   QA->c01[idx]*TMP->r1.c2[idx]+QA->c12[idx]*TMP->r0.c1[idx]));
  // seconda riga
  r0_2 = (-tr1*I) * conj(QA->c01[idx])  + tr2*(QA->c12[idx]*conj(QA->c02[idx])+(1.0*I)*conj(QA->c01[idx])*(QA->rc00[idx]+QA->rc11[idx]))
    + f1 * TMP->r1.c0[idx] + f2 *(-TMP->r1.c0[idx]*(QA->rc00[idx]+QA->rc11[idx])+(1.0*I)*(-conj(QA->c01[idx])*(TMP->r0.c0[idx]+TMP->r1.c1[idx])  
									       + QA->c12[idx]*TMP->r2.c0[idx] - conj(QA->c02[idx])*TMP->r1.c2[idx]));

  r1_2 = -tr1*QA->rc11[idx]   + tr2*( QA->rc11[idx]*QA->rc11[idx]+QA->c01[idx]*conj(QA->c01[idx]) + QA->c12[idx] * conj(QA->c12[idx]))
    + f1 * TMP->r1.c1[idx] + f2*(-2.0*TMP->r1.c1[idx]*QA->rc11[idx]+(1.0*I)*(QA->c12[idx]*TMP->r2.c1[idx]-conj(QA->c01[idx])*TMP->r0.c1[idx]+
									     QA->c01[idx]*TMP->r1.c0[idx]-conj(QA->c12[idx])*TMP->r1.c2[idx]));
  r2_2 = (tr1*I)*QA->c12[idx] + tr2*((1.0*I)*QA->rc00[idx] * QA->c12[idx] + QA->c02[idx] * conj(QA->c01[idx]))
    + f1 * TMP->r1.c2[idx] + f2*(QA->rc00[idx]*TMP->r1.c2[idx]+(1.0*I)*(QA->c12[idx]*(TMP->r1.c1[idx]+TMP->r2.c2[idx]) 
									+QA->c02[idx]*TMP->r1.c0[idx]-conj(QA->c01[idx])*TMP->r0.c2[idx]));
  // terza riga
  U20  = (-tr1*I)*conj(QA->c02[idx]) + tr2*((-1.0*I)*QA->rc11[idx]*conj(QA->c02[idx]) - conj(QA->c01[idx]*QA->c12[idx]))
    + f1 * TMP->r2.c0[idx] + f2 *(QA->rc11[idx]*TMP->r2.c0[idx]+(-1.0*I)*(conj(QA->c02[idx])*(TMP->r0.c0[idx]+TMP->r2.c2[idx])+
									  conj(QA->c01[idx])*TMP->r2.c1[idx]+conj(QA->c12[idx])*TMP->r1.c0[idx]));
  U21  = (-tr1*I)*conj(QA->c12[idx]) + tr2*((-1.0*I)*QA->rc00[idx]*conj(QA->c12[idx]) + QA->c01[idx]*conj(QA->c02[idx]))
    + f1 * TMP->r2.c1[idx] + f2 *(QA->rc00[idx]*TMP->r2.c1[idx]+(1.0*I)*(QA->c01[idx]*TMP->r2.c0[idx]-conj(QA->c02[idx])*TMP->r0.c1[idx]
									 -conj(QA->c12[idx])*(TMP->r1.c1[idx]+TMP->r2.c2[idx])));
  U22  = tr1*(QA->rc00[idx]+QA->rc11[idx])+ tr2*((QA->rc00[idx]+QA->rc11[idx])*(QA->rc00[idx]+QA->rc11[idx])
						 + QA->c02[idx] * conj(QA->c02[idx])+ QA->c12[idx] * conj(QA->c12[idx]))
    + f1 * TMP->r2.c2[idx] + f2 *(2.0*(QA->rc00[idx]+QA->rc11[idx])*TMP->r2.c2[idx]+(1.0*I)*(QA->c02[idx]*TMP->r2.c0[idx]+QA->c12[idx]*TMP->r2.c1[idx]
											     -conj(QA->c02[idx])*TMP->r0.c2[idx]
											     -conj(QA->c12[idx])*TMP->r1.c2[idx]));

  /////////////////////////////////////////
  ///////////////// INFINE CALCOLO DI LAMBDA = 0.5*(GAMMA + GAMMA^CROCE) - (1/6)* Id * Tr(GAMMA + GAMMA^CROCE)
  ///// LAMBDA_00 = (2*re(G00)-re(G11)-re(G22))/3
  ///// LAMBDA_11 = (2*re(G11)-re(G00)-re(G22))/3
  ///// LAMBDA_01 = (G01+conj(G10))/2
  ///// LAMBDA_02 = (G02+conj(G20))/2
  ///// LAMBDA_12 = (G12+conj(G21))/2

  L->rc00[idx] = (2*creal(r0_1)-creal(r1_2)-creal(U22))*ONE_BY_THREE;
  L->rc11[idx] = (2*creal(r1_2)-creal(r0_1)-creal(U22))*ONE_BY_THREE;
  L->c01[idx]  = (r1_1+conj(r0_2))*0.5;
  L->c02[idx]  = (r2_1+conj(U20))*0.5;
  L->c12[idx]  = (r2_2+conj(U21))*0.5;

}



void compute_lambda(__restrict thmat_soa * const L, // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
		    __restrict su3_soa   * const SP, // Sigma primo --> input (forza fermionica del passo precedente)
		    __restrict su3_soa   * const U,    // la configurazione di gauge --> input
		    __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
		    __restrict su3_soa   * const TMP  // variabile di parcheggio
		    ){


  int x, y, z, t;
#pragma acc kernels present(L)  present(SP)  present(U)  present(QA)  present(TMP)
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
#pragma acc loop seq
          for(mu=0;mu<4;mu++){
            dir_link = 2*mu + parity;
	    compute_loc_Lambda(&L[dir_link],&SP[dir_link],&U[dir_link],&QA[dir_link],&TMP[dir_link],idxh);
          }
        }  // x
      }  // y
    }  // z
  }  // t
}


#ifdef STOT_FERMIONS
inline void stout_wrapper(su3_soa * tconf_acc, su3_soa * tstout_conf_acc_arr){


    stout_isotropic(tconf_acc, tstout_conf_acc_arr, auxbis_conf_acc, glocal_staples, gtipdot );
    for(int stoutlevel=1;stoutlevel < STOUT_STEPS; stoutlevel++)
        stout_isotropic(&(tstout_conf_acc_arr[8*(stoutlevel-1)]),&(tstout_conf_acc_arr[8*stoutlevel]),auxbis_conf_acc, glocal_staples,  gtipdot );

}
#endif





#pragma acc routine seq
void compute_sigma_local_PEZZO1(__restrict thmat_soa * const L,  // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
				__restrict su3_soa   * const U,  // la configurazione di gauge --> input
				__restrict su3_soa   * const SP,  // entra Sigma primo (input: fermforce del passo precedente) ED esce Sigma --> sia input che ouput
				__restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
				__restrict su3_soa   * const TMP, // variabile di parcheggio
				int idx
				){


  double c0 = det_i_times_QA_soa(QA,idx);
  double c1  = 0.5 * Tr_i_times_QA_sq_soa(QA,idx);
  double c0max = 2*pow(c1/3,1.5);
  d_complex f0,f1,f2;
  if(c1<4e-3)
    {
      f0 = (1-c0*c0/720) + (1.0*I)*(-c0*(1-c1*(1-c1/42)/20)/6);
      f1 = (c0*(1-c1*(1-3*c1/112)/15)/24) + (1.0*I)*(1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040);
      f2 = (0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160)) + (1.0*I)*(0.5*(c0*(1-c1*(1-c1/48)/21)/60));
    }
  else
    {
      int segno=1;
      if(c0<0)
        {
          segno=-1;
          c0=-c0;
        }

      double eps=(c0max-c0)/c0max;
      double theta;
      if(eps<0) theta=0.0;
      else
        if(eps<1e-3) theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
        else theta=acos(c0/c0max);
      double u = sqrt(c1/3) * cos(theta/3) ;
      double w = sqrt(c1) * sin(theta/3) ;
      double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
      double cu=cos(u),c2u=cos(2*u);
      double su=sin(u),s2u=sin(2*u);
      double cw=cos(w);
      double xi0w;
      if(fabs(w)<0.05)
        {
          double temp0=w*w,temp1=1-temp0/42,temp2=1.0-temp0/20*temp1;
          xi0w=1-temp0/6*temp2;
        }
      else xi0w=sin(w)/w;
      
      double denom = 1/(9*u*u - w*w);
      
      f0=(u2mw2*c2u+ //(u2-w2)*cos(2u)
	  cu*8*u2*cw+ //cos(u)*8*u2*cos(w)
	  2*su*u*w2p3u2*xi0w) //sin(u)*2*mu*(3*u2+w2)*xi0(w)
        +(1.0*I)*(u2mw2*s2u+ //(u2-w2)*sin(2u)
		  -su*8*u2*cw+ //-sin(u)*8*u2*cos(w)
		  cu*2*u*w2p3u2*xi0w); //cos(u)*2*u*(3*u2+w2)*xi0(w)
      f0 *= denom;
      f1=(2*u*c2u+ //2*u*cos(2u)
	  -cu*2*u*cw+ //cos(u)*2*u*cos(w)
	  -su*w2m3u2*xi0w) //sin(u)*(u2-3*w2)*xi0(w)
        +(1.0*I)*(2*u*s2u+ //2*u*sin(2u)
		  su*2*u*cos(w)+ //sin(u)*2*u*cos(w)
                  -cu*w2m3u2*xi0w);//cos(u)*(3*u2-w2)*xi0(w)
      f1 *= denom;
      f2=(c2u+ //cos(2u)
	  -cu*cw+ //-cos(u)*cos(w)
          -3*su*u*xi0w) //-3*sin(u)*u*xi0(w)
        +(1.0*I)*(s2u+ //sin(2u)
		  su*cw+ //sin(w)*cos(w)
                  -cu*3*u*xi0w);//-cos(u)*3*u*xi0(w)
      f2 *= denom;

      if(segno==-1){
        f0 =  conj(f0);
        f1 = -conj(f1);
        f2 =  conj(f2);
      }
    }

    ///   SIGMA = SIGMA'*EXP(iQ) + i * RHO * STAPLE(complete, non parte TA)*LAMBDA
    ///           - i * RHO * (  U Udag Udag LAMBDA + Udag Udag LAMBDA U + Udag LAMBDA Udag U
    ///                        - Udag Udag LAMBDA U - LAMBDA U Udag Udag - U Udag LAMBDA Udag ) =
    ///            PEZZO1 + PEZZO2 + PEZZO3

    //Calcolo il PEZZO 1 e lo metto in RES

    //////  CALCOLO TMP = EXP(iQ)
    TMP->r0.c0[idx] =f0-f1  *QA->rc00[idx]+f2*(QA->rc00[idx]*QA->rc00[idx]+QA->c01[idx]*conj(QA->c01[idx]) + QA->c02[idx]*conj(QA->c02[idx]));
    TMP->r0.c1[idx] =( f1*I)*QA->c01[idx] +f2*(QA->c02[idx]*conj(QA->c12[idx])+(-1.0*I)* QA->c01[idx] * ( QA->rc00[idx] + QA->rc11[idx]));
    TMP->r0.c2[idx] =( f1*I)*QA->c02[idx] +f2*(-QA->c01[idx] * QA->c12[idx]+( 1.0*I)* QA->c02[idx] * QA->rc11[idx]);
    TMP->r1.c0[idx] =(-f1*I)*conj(QA->c01[idx])+ f2*(QA->c12[idx]*conj(QA->c02[idx])+(1.0*I) * conj(QA->c01[idx]) * ( QA->rc00[idx] + QA->rc11[idx]));
    TMP->r1.c1[idx] =f0-f1*QA->rc11[idx]  +f2*(QA->rc11[idx] * QA->rc11[idx]+ QA->c01[idx] * conj(QA->c01[idx])+ QA->c12[idx] * conj(QA->c12[idx]));
    TMP->r1.c2[idx] =( f1*I)*QA->c12[idx] + f2*((1.0*I)*QA->rc00[idx] * QA->c12[idx]+ QA->c02[idx] * conj(QA->c01[idx]));
    TMP->r2.c0[idx] = conj(TMP->r0.c1[idx] * TMP->r1.c2[idx]  - TMP->r0.c2[idx]  * TMP->r1.c1[idx]);
    TMP->r2.c1[idx] = conj(TMP->r0.c2[idx] * TMP->r1.c0[idx]  - TMP->r0.c0[idx]  * TMP->r1.c2[idx]);
    TMP->r2.c2[idx] = conj(TMP->r0.c0[idx] * TMP->r1.c1[idx]  - TMP->r0.c1[idx]  * TMP->r1.c0[idx]);

    // PEZZO2 e PEZZO3 sono gia' stati calcolati da un'altra routine (a meno del fatto di moltiplicare per RHO)
    // Adesso quindi faccio RES = RES * RHO + PEZZO1
    //////////  CALCOLO PEZZO1 ==>  SIGMA' * EXP(IQ) = SIGMA' * TMP
    d_complex s00 = SP->r0.c0[idx] * TMP->r0.c0[idx] + SP->r0.c1[idx] * TMP->r1.c0[idx] + SP->r0.c2[idx] * TMP->r2.c0[idx];
    d_complex s01 = SP->r0.c0[idx] * TMP->r0.c1[idx] + SP->r0.c1[idx] * TMP->r1.c1[idx] + SP->r0.c2[idx] * TMP->r2.c1[idx];
    d_complex s02 = SP->r0.c0[idx] * TMP->r0.c2[idx] + SP->r0.c1[idx] * TMP->r1.c2[idx] + SP->r0.c2[idx] * TMP->r2.c2[idx];
    SP->r0.c0[idx] = s00;
    SP->r0.c1[idx] = s01;
    SP->r0.c2[idx] = s02;


    s00 = SP->r1.c0[idx] * TMP->r0.c0[idx] + SP->r1.c1[idx] * TMP->r1.c0[idx] + SP->r1.c2[idx] * TMP->r2.c0[idx];
    s01 = SP->r1.c0[idx] * TMP->r0.c1[idx] + SP->r1.c1[idx] * TMP->r1.c1[idx] + SP->r1.c2[idx] * TMP->r2.c1[idx];
    s02 = SP->r1.c0[idx] * TMP->r0.c2[idx] + SP->r1.c1[idx] * TMP->r1.c2[idx] + SP->r1.c2[idx] * TMP->r2.c2[idx];
    SP->r1.c0[idx] = s00;
    SP->r1.c1[idx] = s01;
    SP->r1.c2[idx] = s02;

    s00 = SP->r2.c0[idx] * TMP->r0.c0[idx] + SP->r2.c1[idx] * TMP->r1.c0[idx] + SP->r2.c2[idx] * TMP->r2.c0[idx];
    s01 = SP->r2.c0[idx] * TMP->r0.c1[idx] + SP->r2.c1[idx] * TMP->r1.c1[idx] + SP->r2.c2[idx] * TMP->r2.c1[idx];
    s02 = SP->r2.c0[idx] * TMP->r0.c2[idx] + SP->r2.c1[idx] * TMP->r1.c2[idx] + SP->r2.c2[idx] * TMP->r2.c2[idx];
    SP->r2.c0[idx] = s00;
    SP->r2.c1[idx] = s01;
    SP->r2.c2[idx] = s02;
}



#pragma acc routine seq
static inline void RIGHT_iABC_times_DminusE_absent_stag_phases(  __restrict su3_soa * const UA,
								 const int idxA,
								 __restrict su3_soa * const UB,
								 const int idxB,
								 __restrict su3_soa * const UC,
								 const int idxC,
								 __restrict thmat_soa * const LD,
								 const int idxD,
								 __restrict thmat_soa * const LE,
								 const int idxE,
								 __restrict su3_soa * const RES,
								 const int idxRES){
  // Cosa calcoliamo in questa routine:
  //  RES += UA * dag(UB) * dag(UC) * ((RHO*I)*(LD - LE))
  d_complex matA_00 = UA->r0.c0[idxA];
  d_complex matA_01 = UA->r0.c1[idxA];
  d_complex matA_02 = UA->r0.c2[idxA];
  d_complex matA_10 = UA->r1.c0[idxA];
  d_complex matA_11 = UA->r1.c1[idxA];
  d_complex matA_12 = UA->r1.c2[idxA];

  // construct (into the variables matB_ij) the hermitian conjugate of the UB matrix
  d_complex matB_00 = conj( UB->r0.c0[idxB] ) ;
  d_complex matB_10 = conj( UB->r0.c1[idxB] ) ;
  d_complex matB_20 = conj( UB->r0.c2[idxB] ) ;
  d_complex matB_01 = conj( UB->r1.c0[idxB] ) ;
  d_complex matB_11 = conj( UB->r1.c1[idxB] ) ;
  d_complex matB_21 = conj( UB->r1.c2[idxB] ) ;
  //Compute 3rd matB column from the first two
  d_complex matB_02 = conj( ( matB_10 * matB_21 ) - ( matB_20 * matB_11) ) ;
  d_complex matB_12 = conj( ( matB_20 * matB_01 ) - ( matB_00 * matB_21) ) ;
  d_complex matB_22 = conj( ( matB_00 * matB_11 ) - ( matB_10 * matB_01) ) ;

  //Compute the first two rows of the result of UA * ~UB and assign to matT
  d_complex matT_00 = matA_00 * matB_00 + matA_01 * matB_10 + matA_02 * matB_20 ;
  d_complex matT_01 = matA_00 * matB_01 + matA_01 * matB_11 + matA_02 * matB_21 ;
  d_complex matT_02 = matA_00 * matB_02 + matA_01 * matB_12 + matA_02 * matB_22 ;
  d_complex matT_10 = matA_10 * matB_00 + matA_11 * matB_10 + matA_12 * matB_20 ;
  d_complex matT_11 = matA_10 * matB_01 + matA_11 * matB_11 + matA_12 * matB_21 ;
  d_complex matT_12 = matA_10 * matB_02 + matA_11 * matB_12 + matA_12 * matB_22 ;

  // construct (into the variables matB_ij) the hermitian conjugate of the matC matrix
  matB_00 = conj( UC->r0.c0[idxC] ) ;
  matB_10 = conj( UC->r0.c1[idxC] ) ;
  matB_20 = conj( UC->r0.c2[idxC] ) ;
  matB_01 = conj( UC->r1.c0[idxC] ) ;
  matB_11 = conj( UC->r1.c1[idxC] ) ;
  matB_21 = conj( UC->r1.c2[idxC] ) ;
  //Compute 3rd UC column from the first two
  matB_02 = conj( ( matB_10 * matB_21 ) - ( matB_20 * matB_11) ) ;
  matB_12 = conj( ( matB_20 * matB_01 ) - ( matB_00 * matB_21) ) ;
  matB_22 = conj( ( matB_00 * matB_11 ) - ( matB_10 * matB_01) ) ;


  //Compute the first two rows of the result of matT * ~UC and assign to matA
  matA_00 = matT_00 * matB_00 + matT_01 * matB_10 + matT_02 * matB_20 ;
  matA_01 = matT_00 * matB_01 + matT_01 * matB_11 + matT_02 * matB_21 ;
  matA_02 = matT_00 * matB_02 + matT_01 * matB_12 + matT_02 * matB_22 ;

  matA_10 = matT_10 * matB_00 + matT_11 * matB_10 + matT_12 * matB_20 ;
  matA_11 = matT_10 * matB_01 + matT_11 * matB_11 + matT_12 * matB_21 ;
  matA_12 = matT_10 * matB_02 + matT_11 * matB_12 + matT_12 * matB_22 ;
  // la terza riga la metto in matT_0j
  matT_00 = conj( ( matA_01 * matA_12 ) - ( matA_02 * matA_11) );
  matT_01 = conj( ( matA_02 * matA_10 ) - ( matA_00 * matA_12) ) ;
  matT_02 = conj( ( matA_00 * matA_11 ) - ( matA_10 * matA_01) ) ;

  // write into RES the product ABC * ((RHO*I)*(LD-LE))
  //  RES += matA
  RES->r0.c0[idxRES]+=(RHO*I)*(matA_00*(LD->rc00[idxD]-LE->rc00[idxE])+matA_01*conj(LD->c01[idxD]-LE->c01[idxE])+matA_02*conj(LD->c02[idxD]-LE->c02[idxE]));
  RES->r1.c0[idxRES]+=(RHO*I)*(matA_10*(LD->rc00[idxD]-LE->rc00[idxE])+matA_11*conj(LD->c01[idxD]-LE->c01[idxE])+matA_12*conj(LD->c02[idxD]-LE->c02[idxE]));
  RES->r2.c0[idxRES]+=(RHO*I)*(matT_00*(LD->rc00[idxD]-LE->rc00[idxE])+matT_01*conj(LD->c01[idxD]-LE->c01[idxE])+matT_02*conj(LD->c02[idxD]-LE->c02[idxE]));

  RES->r0.c1[idxRES]+=(RHO*I)*(matA_00*(LD->c01[idxD]-LE->c01[idxE])+matA_01*(LD->rc11[idxD]-LE->rc11[idxE])+matA_02*conj(LD->c12[idxD]-LE->c12[idxE]));
  RES->r1.c1[idxRES]+=(RHO*I)*(matA_10*(LD->c01[idxD]-LE->c01[idxE])+matA_11*(LD->rc11[idxD]-LE->rc11[idxE])+matA_12*conj(LD->c12[idxD]-LE->c12[idxE]));
  RES->r2.c1[idxRES]+=(RHO*I)*(matT_00*(LD->c01[idxD]-LE->c01[idxE])+matT_01*(LD->rc11[idxD]-LE->rc11[idxE])+matT_02*conj(LD->c12[idxD]-LE->c12[idxE]));

  RES->r0.c2[idxRES]+=(RHO*I)*(matA_00*(LD->c02[idxD]-LE->c02[idxE])+matA_01*(LD->c12[idxD]-LE->c12[idxE])-matA_02*(LD->rc00[idxD]-LE->rc00[idxE]+LD->rc11[idxD]-LE->rc11[idxE]));
  RES->r1.c2[idxRES]+=(RHO*I)*(matA_10*(LD->c02[idxD]-LE->c02[idxE])+matA_11*(LD->c12[idxD]-LE->c12[idxE])-matA_12*(LD->rc00[idxD]-LE->rc00[idxE]+LD->rc11[idxD]-LE->rc11[idxE]));
  RES->r2.c2[idxRES]+=(RHO*I)*(matT_00*(LD->c02[idxD]-LE->c02[idxE])+matT_01*(LD->c12[idxD]-LE->c12[idxE])-matT_02*(LD->rc00[idxD]-LE->rc00[idxE]+LD->rc11[idxD]-LE->rc11[idxE]));
}


// Questa e' di quelle della categoria RIGHT 
#pragma acc routine seq
static inline void RIGHT_iFABC_absent_stag_phases(  __restrict su3_soa * const UA,
						    const int idxA,
						    __restrict su3_soa * const UB,
						    const int idxB,
						    __restrict su3_soa * const UC,
						    const int idxC,
						    __restrict thmat_soa * const LF,
						    const int idxF,
						    __restrict su3_soa * const RES,
						    const int idxRES){
  // Cosa calcoliamo in questa routine:
  //  RES += ((RHO*I)*LF) * UA * dag(UB) * dag(UC)

  d_complex matA_00 = UA->r0.c0[idxA];
  d_complex matA_01 = UA->r0.c1[idxA];
  d_complex matA_02 = UA->r0.c2[idxA];
  d_complex matA_10 = UA->r1.c0[idxA];
  d_complex matA_11 = UA->r1.c1[idxA];
  d_complex matA_12 = UA->r1.c2[idxA];
  //Compute 3rd matA row from the first two
  d_complex matA_20 = conj( ( matA_01 * matA_12 ) - ( matA_02 * matA_11) ) ;
  d_complex matA_21 = conj( ( matA_02 * matA_10 ) - ( matA_00 * matA_12) ) ;
  d_complex matA_22 = conj( ( matA_00 * matA_11 ) - ( matA_01 * matA_10) ) ;

  //Compute the first two rows of the result of LF * UA and assign to matT
  d_complex matT_00 = (RHO*I)*(LF->rc00[idxF]     *matA_00 + LF->c01[idxF]      *matA_10 + LF->c02[idxF]                  *matA_20);
  d_complex matT_01 = (RHO*I)*(LF->rc00[idxF]     *matA_01 + LF->c01[idxF]      *matA_11 + LF->c02[idxF]                  *matA_21);
  d_complex matT_02 = (RHO*I)*(LF->rc00[idxF]     *matA_02 + LF->c01[idxF]      *matA_12 + LF->c02[idxF]                  *matA_22);
  d_complex matT_10 = (RHO*I)*(conj(LF->c01[idxF])*matA_00 + LF->rc11[idxF]     *matA_10 + LF->c12[idxF]                  *matA_20);
  d_complex matT_11 = (RHO*I)*(conj(LF->c01[idxF])*matA_01 + LF->rc11[idxF]     *matA_11 + LF->c12[idxF]                  *matA_21);
  d_complex matT_12 = (RHO*I)*(conj(LF->c01[idxF])*matA_02 + LF->rc11[idxF]     *matA_12 + LF->c12[idxF]                  *matA_22);
  d_complex matT_20 = (RHO*I)*(conj(LF->c02[idxF])*matA_00 + conj(LF->c12[idxF])*matA_10 - (LF->rc00[idxF]+LF->rc11[idxF])*matA_20);
  d_complex matT_21 = (RHO*I)*(conj(LF->c02[idxF])*matA_01 + conj(LF->c12[idxF])*matA_11 - (LF->rc00[idxF]+LF->rc11[idxF])*matA_21);
  d_complex matT_22 = (RHO*I)*(conj(LF->c02[idxF])*matA_02 + conj(LF->c12[idxF])*matA_12 - (LF->rc00[idxF]+LF->rc11[idxF])*matA_22);
  // construct (into the variables matB_ij) the hermitian conjugate of the UB matrix
  // --> non le salvo in temporanei [che lascio commentati], ma le scrivo direttamente nel prodotto
  //  matB_00 = conj( UB->r0.c0[idxB] ) ;
  //  matB_10 = conj( UB->r0.c1[idxB] ) ;
  //  matB_20 = conj( UB->r0.c2[idxB] ) ;
  //  matB_01 = conj( UB->r1.c0[idxB] ) ;
  //  matB_11 = conj( UB->r1.c1[idxB] ) ;
  //  matB_21 = conj( UB->r1.c2[idxB] ) ;
  //Compute 3rd matB column from the first two
  //  matB_02 = ( UB->r0.c1[idxB] * UB->r1.c2[idxB]  - UB->r0.c2[idxB] * UB->r1.c1[idxB] ) ;
  //  matB_12 = ( UB->r0.c2[idxB] * UB->r1.c0[idxB]  - UB->r0.c0[idxB] * UB->r1.c2[idxB] ) ;
  //  matB_22 = ( UB->r0.c0[idxB] * UB->r1.c1[idxB]  - UB->r0.c1[idxB] * UB->r1.c0[idxB] ) ;

  //Compute the first two rows of the result of (RHO*I)*LF * UA * ~UB = T * ~UB and assign to matA
  matA_00 = matT_00 * conj( UB->r0.c0[idxB] ) + matT_01 * conj( UB->r0.c1[idxB] ) + matT_02 * conj( UB->r0.c2[idxB] ) ;
  matA_10 = matT_10 * conj( UB->r0.c0[idxB] ) + matT_11 * conj( UB->r0.c1[idxB] ) + matT_12 * conj( UB->r0.c2[idxB] ) ;
  matA_20 = matT_20 * conj( UB->r0.c0[idxB] ) + matT_21 * conj( UB->r0.c1[idxB] ) + matT_22 * conj( UB->r0.c2[idxB] ) ;

  matA_01 = matT_00 * conj( UB->r1.c0[idxB] ) + matT_01 * conj( UB->r1.c1[idxB] ) + matT_02 * conj( UB->r1.c2[idxB] ) ;
  matA_11 = matT_10 * conj( UB->r1.c0[idxB] ) + matT_11 * conj( UB->r1.c1[idxB] ) + matT_12 * conj( UB->r1.c2[idxB] ) ;
  matA_21 = matT_20 * conj( UB->r1.c0[idxB] ) + matT_21 * conj( UB->r1.c1[idxB] ) + matT_22 * conj( UB->r1.c2[idxB] ) ;

  matA_02 = matT_00 * ( UB->r0.c1[idxB] * UB->r1.c2[idxB]  - UB->r0.c2[idxB] * UB->r1.c1[idxB] ) + matT_01 * ( UB->r0.c2[idxB] * UB->r1.c0[idxB]  - UB->r0.c0[idxB] * UB->r1.c2[idxB] ) + matT_02 * ( UB->r0.c0[idxB] * UB->r1.c1[idxB]  - UB->r0.c1[idxB] * UB->r1.c0[idxB] ) ;
  matA_12 = matT_10 * ( UB->r0.c1[idxB] * UB->r1.c2[idxB]  - UB->r0.c2[idxB] * UB->r1.c1[idxB] ) + matT_11 * ( UB->r0.c2[idxB] * UB->r1.c0[idxB]  - UB->r0.c0[idxB] * UB->r1.c2[idxB] ) + matT_12 * ( UB->r0.c0[idxB] * UB->r1.c1[idxB]  - UB->r0.c1[idxB] * UB->r1.c0[idxB] ) ;
  matA_22 = matT_20 * ( UB->r0.c1[idxB] * UB->r1.c2[idxB]  - UB->r0.c2[idxB] * UB->r1.c1[idxB] ) + matT_21 * ( UB->r0.c2[idxB] * UB->r1.c0[idxB]  - UB->r0.c0[idxB] * UB->r1.c2[idxB] ) + matT_22 * ( UB->r0.c0[idxB] * UB->r1.c1[idxB]  - UB->r0.c1[idxB] * UB->r1.c0[idxB] ) ;
  
  // construct (into the variables matC_ij) the hermitian conjugate of the matC matrix
  // --> non le salvo in temporanei [che lascio commentati], ma le scrivo direttamente nel prodotto
  //  matC_00 = conj( UC->r0.c0[idxC] ) ;
  //  matC_10 = conj( UC->r0.c1[idxC] ) ;
  //  matC_20 = conj( UC->r0.c2[idxC] ) ;
  //  matC_01 = conj( UC->r1.c0[idxC] ) ;
  //  matC_11 = conj( UC->r1.c1[idxC] ) ;
  //  matC_21 = conj( UC->r1.c2[idxC] ) ;
  //Compute 3rd UC column from the first two
  //  matC_02 = ( UC->r0.c1[idxC] * UC->r1.c2[idxC]  - UC->r0.c2[idxC] * UC->r1.c1[idxC] ) ;
  //  matC_12 = ( UC->r0.c2[idxC] * UC->r1.c0[idxC]  - UC->r0.c0[idxC] * UC->r1.c2[idxC] ) ;
  //  matC_22 = ( UC->r0.c0[idxC] * UC->r1.c1[idxC]  - UC->r0.c1[idxC] * UC->r1.c0[idxC] ) ;
  // add to RES the product I * F * ABC 
  //  RES += matA * ~UC
  RES->r0.c0[idxRES]+= matA_00 * conj( UC->r0.c0[idxB] ) + matA_01 * conj( UC->r0.c1[idxB] ) + matA_02 * conj( UC->r0.c2[idxB] ) ;
  RES->r1.c0[idxRES]+= matA_10 * conj( UC->r0.c0[idxB] ) + matA_11 * conj( UC->r0.c1[idxB] ) + matA_12 * conj( UC->r0.c2[idxB] ) ;
  RES->r2.c0[idxRES]+= matA_20 * conj( UC->r0.c0[idxB] ) + matA_21 * conj( UC->r0.c1[idxB] ) + matA_22 * conj( UC->r0.c2[idxB] ) ;

  RES->r0.c1[idxRES]+= matA_00 * conj( UC->r1.c0[idxB] ) + matA_01 * conj( UC->r1.c1[idxB] ) + matA_02 * conj( UC->r1.c2[idxB] ) ;
  RES->r1.c1[idxRES]+= matA_10 * conj( UC->r1.c0[idxB] ) + matA_11 * conj( UC->r1.c1[idxB] ) + matA_12 * conj( UC->r1.c2[idxB] ) ;
  RES->r2.c1[idxRES]+= matA_20 * conj( UC->r1.c0[idxB] ) + matA_21 * conj( UC->r1.c1[idxB] ) + matA_22 * conj( UC->r1.c2[idxB] ) ;

  RES->r0.c2[idxRES]+= matA_00 * ( UC->r0.c1[idxB] * UC->r1.c2[idxB]  - UC->r0.c2[idxB] * UC->r1.c1[idxB] ) + matA_01 * ( UC->r0.c2[idxB] * UC->r1.c0[idxB]  - UC->r0.c0[idxB] * UC->r1.c2[idxB] ) + matA_02 * ( UC->r0.c0[idxB] * UC->r1.c1[idxB]  - UC->r0.c1[idxB] * UC->r1.c0[idxB] ) ;
  RES->r1.c2[idxRES]+=matA_10 * ( UC->r0.c1[idxB] * UC->r1.c2[idxB]  - UC->r0.c2[idxB] * UC->r1.c1[idxB] )  + matA_11 * ( UC->r0.c2[idxB] * UC->r1.c0[idxB]  - UC->r0.c0[idxB] * UC->r1.c2[idxB] ) + matA_12 * ( UC->r0.c0[idxB] * UC->r1.c1[idxB]  - UC->r0.c1[idxB] * UC->r1.c0[idxB] ) ;
  RES->r2.c2[idxRES]+= matA_20 * ( UC->r0.c1[idxB] * UC->r1.c2[idxB]  - UC->r0.c2[idxB] * UC->r1.c1[idxB] ) + matA_21 * ( UC->r0.c2[idxB] * UC->r1.c0[idxB]  - UC->r0.c0[idxB] * UC->r1.c2[idxB] ) + matA_22 * ( UC->r0.c0[idxB] * UC->r1.c1[idxB]  - UC->r0.c1[idxB] * UC->r1.c0[idxB] ) ;
}



// Questa e' di quelle della categoria RIGHT 
#pragma acc routine seq
static inline void RIGHT_miABGC_absent_stag_phases(  __restrict su3_soa * const UA,
						     const int idxA,
						     __restrict su3_soa * const UB,
						     const int idxB,
						     __restrict su3_soa * const UC,
						     const int idxC,
						     __restrict thmat_soa * const LG,
						     const int idxG,
						     __restrict su3_soa * const RES,
						     const int idxRES){
  // Cosa calcoliamo in questa routine:
  //  RES +=  UA * dag(UB) * ((-RHO*I)*LG) * dag(UC)

  d_complex matA_00 = UA->r0.c0[idxA];
  d_complex matA_01 = UA->r0.c1[idxA];
  d_complex matA_02 = UA->r0.c2[idxA];
  d_complex matA_10 = UA->r1.c0[idxA];
  d_complex matA_11 = UA->r1.c1[idxA];
  d_complex matA_12 = UA->r1.c2[idxA];

  // construct (into the variables matB_ij) the hermitian conjugate of the UB matrix
  d_complex matB_00 = conj( UB->r0.c0[idxB] ) ;
  d_complex matB_10 = conj( UB->r0.c1[idxB] ) ;
  d_complex matB_20 = conj( UB->r0.c2[idxB] ) ;
  d_complex matB_01 = conj( UB->r1.c0[idxB] ) ;
  d_complex matB_11 = conj( UB->r1.c1[idxB] ) ;
  d_complex matB_21 = conj( UB->r1.c2[idxB] ) ;
  //Compute 3rd matB column from the first two
  d_complex matB_02 = conj( ( matB_10 * matB_21 ) - ( matB_20 * matB_11) ) ;
  d_complex matB_12 = conj( ( matB_20 * matB_01 ) - ( matB_00 * matB_21) ) ;
  d_complex matB_22 = conj( ( matB_00 * matB_11 ) - ( matB_10 * matB_01) ) ;

  //Compute the first two rows of the result of UA * ~UB and assign to matT and matA
  //   the result is stored in the variables:
  // T00 T01 T02
  // A00 A01 A02
  // A10 A11 A12
  d_complex matT_00 = matA_00 * matB_00 + matA_01 * matB_10 + matA_02 * matB_20 ;
  d_complex matT_01 = matA_00 * matB_01 + matA_01 * matB_11 + matA_02 * matB_21 ;
  d_complex matT_02 = matA_00 * matB_02 + matA_01 * matB_12 + matA_02 * matB_22 ;
  /*
  d_complex matT_10 = matA_10 * matB_00 + matA_11 * matB_10 + matA_12 * matB_20 ; // --> matA_00
  d_complex matT_11 = matA_10 * matB_01 + matA_11 * matB_11 + matA_12 * matB_21 ; // --> matA_01
  d_complex matT_12 = matA_10 * matB_02 + matA_11 * matB_12 + matA_12 * matB_22 ; // --> matA_02
  d_complex matT_20 = matA_20 * matB_00 + matA_21 * matB_10 + matA_22 * matB_20 ; // --> matA_10
  d_complex matT_21 = matA_20 * matB_01 + matA_21 * matB_11 + matA_22 * matB_21 ; // --> matA_11
  d_complex matT_22 = matA_20 * matB_02 + matA_21 * matB_12 + matA_22 * matB_22 ; // --> matA_12
  */
  matA_00 = matA_10 * matB_00 + matA_11 * matB_10 + matA_12 * matB_20 ;
  matA_01 = matA_10 * matB_01 + matA_11 * matB_11 + matA_12 * matB_21 ;
  matA_02 = matA_10 * matB_02 + matA_11 * matB_12 + matA_12 * matB_22 ;
  // prodotto vettore per ricostruire la terza dalle prime due righe
  matA_10 = conj(matT_01*matA_02-matT_02*matA_01);
  matA_11 = conj(matT_02*matA_00-matT_00*matA_02);
  matA_12 = conj(matT_00*matA_01-matT_01*matA_00);


  matB_00=(-RHO*I)*(matT_00*(LG->rc00[idxG])+matT_01*conj(LG->c01[idxG])+matT_02*conj(LG->c02[idxG]));
  matB_10=(-RHO*I)*(matA_00*(LG->rc00[idxG])+matA_01*conj(LG->c01[idxG])+matA_02*conj(LG->c02[idxG]));
  matB_20=(-RHO*I)*(matA_10*(LG->rc00[idxG])+matA_11*conj(LG->c01[idxG])+matA_12*conj(LG->c02[idxG]));
  matB_01=(-RHO*I)*(matT_00*(LG->c01[idxG])+matT_01*(LG->rc11[idxG])+matT_02*conj(LG->c12[idxG]));
  matB_11=(-RHO*I)*(matA_00*(LG->c01[idxG])+matA_01*(LG->rc11[idxG])+matA_02*conj(LG->c12[idxG]));
  matB_21=(-RHO*I)*(matA_10*(LG->c01[idxG])+matA_11*(LG->rc11[idxG])+matA_12*conj(LG->c12[idxG]));
  matB_02=(-RHO*I)*(matT_00*(LG->c02[idxG])+matT_01*(LG->c12[idxG])-matT_02*(LG->rc00[idxG]+LG->rc11[idxG]));
  matB_12=(-RHO*I)*(matA_00*(LG->c02[idxG])+matA_01*(LG->c12[idxG])-matA_02*(LG->rc00[idxG]+LG->rc11[idxG]));
  matB_22=(-RHO*I)*(matA_10*(LG->c02[idxG])+matA_11*(LG->c12[idxG])-matA_12*(LG->rc00[idxG]+LG->rc11[idxG]));

  //  RES += matB * ~UC
  RES->r0.c0[idxRES]+= matB_00 * conj( UC->r0.c0[idxC] ) + matB_01 * conj( UC->r0.c1[idxC] ) + matB_02 * conj( UC->r0.c2[idxC] ) ;
  RES->r1.c0[idxRES]+= matB_10 * conj( UC->r0.c0[idxC] ) + matB_11 * conj( UC->r0.c1[idxC] ) + matB_12 * conj( UC->r0.c2[idxC] ) ;
  RES->r2.c0[idxRES]+= matB_20 * conj( UC->r0.c0[idxC] ) + matB_21 * conj( UC->r0.c1[idxC] ) + matB_22 * conj( UC->r0.c2[idxC] ) ;

  RES->r0.c1[idxRES]+= matB_00 * conj( UC->r1.c0[idxC] ) + matB_01 * conj( UC->r1.c1[idxC] ) + matB_02 * conj( UC->r1.c2[idxC] ) ;
  RES->r1.c1[idxRES]+= matB_10 * conj( UC->r1.c0[idxC] ) + matB_11 * conj( UC->r1.c1[idxC] ) + matB_12 * conj( UC->r1.c2[idxC] ) ;
  RES->r2.c1[idxRES]+= matB_20 * conj( UC->r1.c0[idxC] ) + matB_21 * conj( UC->r1.c1[idxC] ) + matB_22 * conj( UC->r1.c2[idxC] ) ;

  RES->r0.c2[idxRES]+= matB_00 * ( UC->r0.c1[idxC] * UC->r1.c2[idxC]  - UC->r0.c2[idxC] * UC->r1.c1[idxC] ) + matB_01 * ( UC->r0.c2[idxC] * UC->r1.c0[idxC]- UC->r0.c0[idxC] * UC->r1.c2[idxC] ) + matB_02 * ( UC->r0.c0[idxC] * UC->r1.c1[idxC]  - UC->r0.c1[idxC] * UC->r1.c0[idxC] ) ;
  RES->r1.c2[idxRES]+=matB_10 * ( UC->r0.c1[idxC] * UC->r1.c2[idxC]  - UC->r0.c2[idxC] * UC->r1.c1[idxC] )  + matB_11 * ( UC->r0.c2[idxC] * UC->r1.c0[idxC]- UC->r0.c0[idxC] * UC->r1.c2[idxC] ) + matB_12 * ( UC->r0.c0[idxC] * UC->r1.c1[idxC]  - UC->r0.c1[idxC] * UC->r1.c0[idxC] ) ;
  RES->r2.c2[idxRES]+= matB_20 * ( UC->r0.c1[idxC] * UC->r1.c2[idxC]  - UC->r0.c2[idxC] * UC->r1.c1[idxC] ) + matB_21 * ( UC->r0.c2[idxC] * UC->r1.c0[idxC]- UC->r0.c0[idxC] * UC->r1.c2[idxC] ) + matB_22 * ( UC->r0.c0[idxC] * UC->r1.c1[idxC]  - UC->r0.c1[idxC] * UC->r1.c0[idxC] ) ;
}



#pragma acc routine seq
static inline void LEFT_iAB_times_GminusE_times_C_absent_stag_phases(  __restrict su3_soa * const UA,
								       const int idxA,
								       __restrict su3_soa * const UB,
								       const int idxB,
								       __restrict su3_soa * const UC,
								       const int idxC,
								       __restrict thmat_soa * const LG,
								       const int idxG,
								       __restrict thmat_soa * const LE,
								       const int idxE,
								       __restrict su3_soa * const RES,
								       const int idxRES){
  // Cosa calcoliamo in questa routine:
  //  RES += dag(UA) * dag(UB) * ((RHO*I)*(LG - LE)) * UC
  // construct (into the variables matA_ij) the hermitian conjugate of the UA matrix
  d_complex matA_00 = conj( UA->r0.c0[idxA] ) ;
  d_complex matA_10 = conj( UA->r0.c1[idxA] ) ;
  d_complex matA_20 = conj( UA->r0.c2[idxA] ) ;
  d_complex matA_01 = conj( UA->r1.c0[idxA] ) ;
  d_complex matA_11 = conj( UA->r1.c1[idxA] ) ;
  d_complex matA_21 = conj( UA->r1.c2[idxA] ) ;
  //Compute 3rd matA column from the first two
  d_complex matA_02 = conj( ( matA_10 * matA_21 ) - ( matA_20 * matA_11) ) ;
  d_complex matA_12 = conj( ( matA_20 * matA_01 ) - ( matA_00 * matA_21) ) ;
  d_complex matA_22 = conj( ( matA_00 * matA_11 ) - ( matA_10 * matA_01) ) ;//questo nel passo dopo non serve

  // construct (into the variables matB_ij) the hermitian conjugate of the UB matrix
  d_complex matB_00 = conj( UB->r0.c0[idxB] ) ;
  d_complex matB_10 = conj( UB->r0.c1[idxB] ) ;
  d_complex matB_20 = conj( UB->r0.c2[idxB] ) ;
  d_complex matB_01 = conj( UB->r1.c0[idxB] ) ;
  d_complex matB_11 = conj( UB->r1.c1[idxB] ) ;
  d_complex matB_21 = conj( UB->r1.c2[idxB] ) ;
  //Compute 3rd matB column from the first two
  d_complex matB_02 = conj( ( matB_10 * matB_21 ) - ( matB_20 * matB_11) ) ;
  d_complex matB_12 = conj( ( matB_20 * matB_01 ) - ( matB_00 * matB_21) ) ;
  d_complex matB_22 = conj( ( matB_00 * matB_11 ) - ( matB_10 * matB_01) ) ;

  //Compute the first two rows of the result of ~UA * ~UB and assign to matT
  d_complex matT_00 = matA_00 * matB_00 + matA_01 * matB_10 + matA_02 * matB_20 ;
  d_complex matT_01 = matA_00 * matB_01 + matA_01 * matB_11 + matA_02 * matB_21 ;
  d_complex matT_02 = matA_00 * matB_02 + matA_01 * matB_12 + matA_02 * matB_22 ;
  d_complex matT_10 = matA_10 * matB_00 + matA_11 * matB_10 + matA_12 * matB_20 ;
  d_complex matT_11 = matA_10 * matB_01 + matA_11 * matB_11 + matA_12 * matB_21 ;
  d_complex matT_12 = matA_10 * matB_02 + matA_11 * matB_12 + matA_12 * matB_22 ;
  //matT_20 = (conj(matT_01*matT_12-matT_02*matT_11));
  //matT_21 = (conj(matT_02*matT_10-matT_00*matT_12));
  //matT_22 = (conj(matT_00*matT_11-matT_01*matT_10));

  // write into RES the product (~A*~B) * ((RHO*I)*(LG-LE)) = matT * ((RHO*I)*(LG-LE))
  //  RES += matA
  matA_00 = (RHO*I)*(matT_00*(LG->rc00[idxG]-LE->rc00[idxE])+matT_01*conj(LG->c01[idxG]-LE->c01[idxE])+matT_02*conj(LG->c02[idxG]-LE->c02[idxE]));
  matA_10 = (RHO*I)*(matT_10*(LG->rc00[idxG]-LE->rc00[idxG])+matT_11*conj(LG->c01[idxG]-LE->c01[idxE])+matT_12*conj(LG->c02[idxG]-LE->c02[idxE]));
  matA_20 = (RHO*I)*((conj(matT_01*matT_12-matT_02*matT_11))*(LG->rc00[idxG]-LE->rc00[idxE])+(conj(matT_02*matT_10-matT_00*matT_12))*conj(LG->c01[idxG]-LE->c01[idxE])+(conj(matT_00*matT_11-matT_01*matT_10))*conj(LG->c02[idxG]-LE->c02[idxE]));

  matA_01 = (RHO*I)*(matT_00*(LG->c01[idxG]-LE->c01[idxE])+matT_01*(LG->rc11[idxG]-LE->rc11[idxE])+matT_02*conj(LG->c12[idxG]-LE->c12[idxE]));
  matA_11 = (RHO*I)*(matT_10*(LG->c01[idxG]-LE->c01[idxE])+matT_11*(LG->rc11[idxG]-LE->rc11[idxE])+matT_12*conj(LG->c12[idxG]-LE->c12[idxE]));
  matA_21 = (RHO*I)*((conj(matT_01*matT_12-matT_02*matT_11))*(LG->c01[idxG]-LE->c01[idxE])+(conj(matT_02*matT_10-matT_00*matT_12))*(LG->rc11[idxG]-LE->rc11[idxE])+(conj(matT_00*matT_11-matT_01*matT_10))*conj(LG->c12[idxG]-LE->c12[idxE]));

  matA_02 = (RHO*I)*(matT_00*(LG->c02[idxG]-LE->c02[idxE])+matT_01*(LG->c12[idxG]-LE->c12[idxE])-matT_02*(LG->rc00[idxG]-LE->rc00[idxE]+LG->rc11[idxG]-LE->rc11[idxE]));
  matA_12 = (RHO*I)*(matT_10*(LG->c02[idxG]-LE->c02[idxE])+matT_11*(LG->c12[idxG]-LE->c12[idxE])-matT_12*(LG->rc00[idxG]-LE->rc00[idxE]+LG->rc11[idxG]-LE->rc11[idxE]));
  matA_22 = (RHO*I)*((conj(matT_01*matT_12-matT_02*matT_11))*(LG->c02[idxG]-LE->c02[idxE])+(conj(matT_02*matT_10-matT_00*matT_12))*(LG->c12[idxG]-LE->c12[idxE])-(conj(matT_00*matT_11-matT_01*matT_10))*(LG->rc00[idxG]-LE->rc00[idxE]+LG->rc11[idxG]-LE->rc11[idxE]));


  /////
  //  RES += matA * UC                                                                 /////3rd row!! /////////
  RES->r0.c0[idxRES]+= matA_00 *  UC->r0.c0[idxC]+ matA_01 *  UC->r1.c0[idxC]+ matA_02 *conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );
  RES->r1.c0[idxRES]+= matA_10 *  UC->r0.c0[idxC]+ matA_11 *  UC->r1.c0[idxC]+ matA_12 *conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );
  RES->r2.c0[idxRES]+= matA_20 *  UC->r0.c0[idxC]+ matA_21 *  UC->r1.c0[idxC]+ matA_22 *conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );

  RES->r0.c1[idxRES]+= matA_00 *  UC->r0.c1[idxC]+ matA_01 *  UC->r1.c1[idxC]+ matA_02 *conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) );
  RES->r1.c1[idxRES]+= matA_10 *  UC->r0.c1[idxC]+ matA_11 *  UC->r1.c1[idxC]+ matA_12 *conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) );
  RES->r2.c1[idxRES]+= matA_20 *  UC->r0.c1[idxC]+ matA_21 *  UC->r1.c1[idxC]+ matA_22 *conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) );

  RES->r0.c2[idxRES]+= matA_00 *  UC->r0.c2[idxC]+ matA_01 *  UC->r1.c2[idxC]+ matA_02 *conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );
  RES->r1.c2[idxRES]+= matA_10 *  UC->r0.c2[idxC]+ matA_11 *  UC->r1.c2[idxC]+ matA_12 *conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );
  RES->r2.c2[idxRES]+= matA_20 *  UC->r0.c2[idxC]+ matA_21 *  UC->r1.c2[idxC]+ matA_22 *conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );

}

#pragma acc routine seq
static inline void LEFT_iABCD_absent_stag_phases(  __restrict su3_soa * const UA,
						   const int idxA,
						   __restrict su3_soa * const UB,
						   const int idxB,
						   __restrict su3_soa * const UC,
						   const int idxC,
						   __restrict thmat_soa * const LD,
						   const int idxD,
						   __restrict su3_soa * const RES,
						   const int idxRES){
  // Cosa calcoliamo in questa routine:
  //  RES += dag(UA) * dag(UB) * UC * ((RHO*I)*LD
  // construct (into the variables matA_ij) the hermitian conjugate of the UA matrix
  d_complex matA_00 = conj( UA->r0.c0[idxA] ) ;
  d_complex matA_10 = conj( UA->r0.c1[idxA] ) ;
  d_complex matA_20 = conj( UA->r0.c2[idxA] ) ;
  d_complex matA_01 = conj( UA->r1.c0[idxA] ) ;
  d_complex matA_11 = conj( UA->r1.c1[idxA] ) ;
  d_complex matA_21 = conj( UA->r1.c2[idxA] ) ;
  //Compute 3rd matA column from the first two
  d_complex matA_02 = conj( ( matA_10 * matA_21 ) - ( matA_20 * matA_11) ) ;
  d_complex matA_12 = conj( ( matA_20 * matA_01 ) - ( matA_00 * matA_21) ) ;
  d_complex matA_22 = conj( ( matA_00 * matA_11 ) - ( matA_10 * matA_01) ) ;

  // construct (into the variables matB_ij) the hermitian conjugate of the UB matrix
  d_complex matB_00 = conj( UB->r0.c0[idxB] ) ;
  d_complex matB_10 = conj( UB->r0.c1[idxB] ) ;
  d_complex matB_20 = conj( UB->r0.c2[idxB] ) ;
  d_complex matB_01 = conj( UB->r1.c0[idxB] ) ;
  d_complex matB_11 = conj( UB->r1.c1[idxB] ) ;
  d_complex matB_21 = conj( UB->r1.c2[idxB] ) ;
  //Compute 3rd matB column from the first two
  d_complex matB_02 = conj( ( matB_10 * matB_21 ) - ( matB_20 * matB_11) ) ;
  d_complex matB_12 = conj( ( matB_20 * matB_01 ) - ( matB_00 * matB_21) ) ;
  d_complex matB_22 = conj( ( matB_00 * matB_11 ) - ( matB_10 * matB_01) ) ;
  
  //Compute the first two rows of the result of ~UA * ~UB and assign to matT
  d_complex matT_00 = matA_00 * matB_00 + matA_01 * matB_10 + matA_02 * matB_20 ;
  d_complex matT_01 = matA_00 * matB_01 + matA_01 * matB_11 + matA_02 * matB_21 ;
  d_complex matT_02 = matA_00 * matB_02 + matA_01 * matB_12 + matA_02 * matB_22 ;
  d_complex matT_10 = matA_10 * matB_00 + matA_11 * matB_10 + matA_12 * matB_20 ;
  d_complex matT_11 = matA_10 * matB_01 + matA_11 * matB_11 + matA_12 * matB_21 ;
  d_complex matT_12 = matA_10 * matB_02 + matA_11 * matB_12 + matA_12 * matB_22 ;
  //matT_20 = (conj(matT_01*matT_12-matT_02*matT_11));
  //matT_21 = (conj(matT_02*matT_10-matT_00*matT_12));
  //matT_22 = (conj(matT_00*matT_11-matT_01*matT_10));

  // compute  A = (~A*~B) * C = matT * UC [only the first two rows]
  matA_00 = matT_00 * UC->r0.c0[idxC] + matT_01 * UC->r1.c0[idxC] + matT_02 * conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );
  matA_10 = matT_10 * UC->r0.c0[idxC] + matT_11 * UC->r1.c0[idxC] + matT_12 * conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );

  matA_01 = matT_00 * UC->r0.c1[idxC] + matT_01 * UC->r1.c1[idxC] + matT_02 * conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) );
  matA_11 = matT_10 * UC->r0.c1[idxC] + matT_11 * UC->r1.c1[idxC] + matT_12 * conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) ); 

  matA_02 = matT_00 * UC->r0.c2[idxC] + matT_01 * UC->r1.c2[idxC] + matT_02 * conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );
  matA_12 = matT_10 * UC->r0.c2[idxC] + matT_11 * UC->r1.c2[idxC] + matT_12 * conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );

  matT_00 = conj( ( matA_01 * matA_12 ) - ( matA_02 * matA_11) );
  matT_01 = conj( ( matA_02 * matA_10 ) - ( matA_00 * matA_12) ) ;
  matT_02 = conj( ( matA_00 * matA_11 ) - ( matA_10 * matA_01) ) ;

  /////
  //  RES += matA * ((RHO*I)*D)
  RES->r0.c0[idxRES]+=(RHO*I)*(matA_00*(LD->rc00[idxD])+matA_01*conj(LD->c01[idxD])+matA_02*conj(LD->c02[idxD]));
  RES->r1.c0[idxRES]+=(RHO*I)*(matA_10*(LD->rc00[idxD])+matA_11*conj(LD->c01[idxD])+matA_12*conj(LD->c02[idxD]));
  RES->r2.c0[idxRES]+=(RHO*I)*(matT_00*(LD->rc00[idxD])+matT_01*conj(LD->c01[idxD])+matT_02*conj(LD->c02[idxD]));

  RES->r0.c1[idxRES]+=(RHO*I)*(matA_00*(LD->c01[idxD])+matA_01*(LD->rc11[idxD])+matA_02*conj(LD->c12[idxD]));
  RES->r1.c1[idxRES]+=(RHO*I)*(matA_10*(LD->c01[idxD])+matA_11*(LD->rc11[idxD])+matA_12*conj(LD->c12[idxD]));
  RES->r2.c1[idxRES]+=(RHO*I)*(matT_00*(LD->c01[idxD])+matT_01*(LD->rc11[idxD])+matT_02*conj(LD->c12[idxD]));

  RES->r0.c2[idxRES]+=(RHO*I)*(matA_00*(LD->c02[idxD])+matA_01*(LD->c12[idxD])-matA_02*(LD->rc00[idxD]+LD->rc11[idxD]));
  RES->r1.c2[idxRES]+=(RHO*I)*(matA_10*(LD->c02[idxD])+matA_11*(LD->c12[idxD])-matA_12*(LD->rc00[idxD]+LD->rc11[idxD]));
  RES->r2.c2[idxRES]+=(RHO*I)*(matT_00*(LD->c02[idxD])+matT_01*(LD->c12[idxD])-matT_02*(LD->rc00[idxD]+LD->rc11[idxD]));

}



#pragma acc routine seq
static inline void LEFT_miAFBC_absent_stag_phases(  __restrict su3_soa * const UA,
						    const int idxA,
						    __restrict su3_soa * const UB,
						    const int idxB,
						    __restrict su3_soa * const UC,
						    const int idxC,
						    __restrict thmat_soa * const LF,
						    const int idxF,
						    __restrict su3_soa * const RES,
						    const int idxRES){
  // Cosa calcoliamo in questa routine:
  //  RES += dag(UA) * (-RHO*I)*LF * dag(UB) * UC
  // construct (into the variables matA_ij) the hermitian conjugate of the UA matrix
  d_complex matA_00 = conj( UA->r0.c0[idxA] ) ;
  d_complex matA_10 = conj( UA->r0.c1[idxA] ) ;
  d_complex matA_20 = conj( UA->r0.c2[idxA] ) ;
  d_complex matA_01 = conj( UA->r1.c0[idxA] ) ;
  d_complex matA_11 = conj( UA->r1.c1[idxA] ) ;
  d_complex matA_21 = conj( UA->r1.c2[idxA] ) ;
  //Compute 3rd matA column from the first two
  d_complex matA_02 = conj( ( matA_10 * matA_21 ) - ( matA_20 * matA_11) ) ;
  d_complex matA_12 = conj( ( matA_20 * matA_01 ) - ( matA_00 * matA_21) ) ;
  d_complex matA_22 = conj( ( matA_00 * matA_11 ) - ( matA_10 * matA_01) ) ;

  // LF00  = (-RHO*I)*LF->rc00[idxF]
  // LF01  = (-RHO*I)*LF->c01[idxF]
  // LF02  = (-RHO*I)*LF->c02[idxF]
  // LF10  = (-RHO*I)*conj(LF->c01[idxF])
  // LF11  = (-RHO*I)*LF->rc11[idxF]
  // LF12  = (-RHO*I)*LF->c12[idxF]
  // LF20  = (-RHO*I)*conj(LF->c02[idxF])
  // LF21  = (-RHO*I)*conj(LF->c12[idxF])
  // LF22  = (RHO*I)*(LF->rc00[idxF]+LF->rc11[idxF])

  //Compute the first two rows of the result of ~UA * (-RHO*I)*LF and assign to matT
  d_complex matT_00 = matA_00 * (-RHO*I)*LF->rc00[idxF] + matA_01 * (-RHO*I)*conj(LF->c01[idxF]) + matA_02 * (-RHO*I)*conj(LF->c02[idxF]) ;
  d_complex matT_01 = matA_00 * (-RHO*I)*LF->c01[idxF]  + matA_01 * (-RHO*I)*LF->rc11[idxF]      + matA_02 * (-RHO*I)*conj(LF->c12[idxF]) ;
  d_complex matT_02 = matA_00 * (-RHO*I)*LF->c02[idxF]  + matA_01 * (-RHO*I)*LF->c12[idxF]       + matA_02 * (RHO*I)*(LF->rc00[idxF]+LF->rc11[idxF]) ;
  d_complex matT_10 = matA_10 * (-RHO*I)*LF->rc00[idxF] + matA_11 * (-RHO*I)*conj(LF->c01[idxF]) + matA_12 * (-RHO*I)*conj(LF->c02[idxF]) ;
  d_complex matT_11 = matA_10 * (-RHO*I)*LF->c01[idxF]  + matA_11 * (-RHO*I)*LF->rc11[idxF]      + matA_12 * (-RHO*I)*conj(LF->c12[idxF]) ;
  d_complex matT_12 = matA_10 * (-RHO*I)*LF->c02[idxF]  + matA_11 * (-RHO*I)*LF->c12[idxF]       + matA_12 * (RHO*I)*(LF->rc00[idxF]+LF->rc11[idxF]) ;
  d_complex matT_20 = matA_20 * (-RHO*I)*LF->rc00[idxF] + matA_21 * (-RHO*I)*conj(LF->c01[idxF]) + matA_22 * (-RHO*I)*conj(LF->c02[idxF]) ;
  d_complex matT_21 = matA_20 * (-RHO*I)*LF->c01[idxF]  + matA_21 * (-RHO*I)*LF->rc11[idxF]      + matA_22 * (-RHO*I)*conj(LF->c12[idxF]) ;
  d_complex matT_22 = matA_20 * (-RHO*I)*LF->c02[idxF]  + matA_21 * (-RHO*I)*LF->c12[idxF]       + matA_22 * (RHO*I)*(LF->rc00[idxF]+LF->rc11[idxF]) ;


  // construct (into the variables matB_ij) the hermitian conjugate of the UB matrix
  d_complex matB_00 = conj( UB->r0.c0[idxB] ) ;
  d_complex matB_10 = conj( UB->r0.c1[idxB] ) ;
  d_complex matB_20 = conj( UB->r0.c2[idxB] ) ;
  d_complex matB_01 = conj( UB->r1.c0[idxB] ) ;
  d_complex matB_11 = conj( UB->r1.c1[idxB] ) ;
  d_complex matB_21 = conj( UB->r1.c2[idxB] ) ;
  //Compute 3rd matB column from the first two
  d_complex matB_02 = conj( ( matB_10 * matB_21 ) - ( matB_20 * matB_11) ) ;
  d_complex matB_12 = conj( ( matB_20 * matB_01 ) - ( matB_00 * matB_21) ) ;
  d_complex matB_22 = conj( ( matB_00 * matB_11 ) - ( matB_10 * matB_01) ) ;  

  // compute ~UA * (-RHO*I)*LF * ~UB = T * ~UB into matA
  matA_00 = matT_00 * matB_00 + matT_01 * matB_10 + matT_02 * matB_20;
  matA_01 = matT_00 * matB_01 + matT_01 * matB_11 + matT_02 * matB_21;
  matA_02 = matT_00 * matB_02 + matT_01 * matB_12 + matT_02 * matB_22;
  matA_10 = matT_10 * matB_00 + matT_11 * matB_10 + matT_12 * matB_20;
  matA_11 = matT_10 * matB_01 + matT_11 * matB_11 + matT_12 * matB_21;
  matA_12 = matT_10 * matB_02 + matT_11 * matB_12 + matT_12 * matB_22;
  matA_20 = matT_20 * matB_00 + matT_21 * matB_10 + matT_22 * matB_20;
  matA_21 = matT_20 * matB_01 + matT_21 * matB_11 + matT_22 * matB_21;
  matA_22 = matT_20 * matB_02 + matT_21 * matB_12 + matT_22 * matB_22;

  //  RES += matA * UC                                                                 /////3rd row!! /////////
  RES->r0.c0[idxRES]+= matA_00 *  UC->r0.c0[idxC]+ matA_01 *  UC->r1.c0[idxC]+ matA_02 *conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );
  RES->r1.c0[idxRES]+= matA_10 *  UC->r0.c0[idxC]+ matA_11 *  UC->r1.c0[idxC]+ matA_12 *conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );
  RES->r2.c0[idxRES]+= matA_20 *  UC->r0.c0[idxC]+ matA_21 *  UC->r1.c0[idxC]+ matA_22 *conj( (UC->r0.c1[idxC] * UC->r1.c2[idxC] ) - ( UC->r0.c2[idxC] * UC->r1.c1[idxC]) );

  RES->r0.c1[idxRES]+= matA_00 *  UC->r0.c1[idxC]+ matA_01 *  UC->r1.c1[idxC]+ matA_02 *conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) );
  RES->r1.c1[idxRES]+= matA_10 *  UC->r0.c1[idxC]+ matA_11 *  UC->r1.c1[idxC]+ matA_12 *conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) );
  RES->r2.c1[idxRES]+= matA_20 *  UC->r0.c1[idxC]+ matA_21 *  UC->r1.c1[idxC]+ matA_22 *conj( (UC->r0.c2[idxC] * UC->r1.c0[idxC] ) - ( UC->r0.c0[idxC] * UC->r1.c2[idxC]) );

  RES->r0.c2[idxRES]+= matA_00 *  UC->r0.c2[idxC]+ matA_01 *  UC->r1.c2[idxC]+ matA_02 *conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );
  RES->r1.c2[idxRES]+= matA_10 *  UC->r0.c2[idxC]+ matA_11 *  UC->r1.c2[idxC]+ matA_12 *conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );
  RES->r2.c2[idxRES]+= matA_20 *  UC->r0.c2[idxC]+ matA_21 *  UC->r1.c2[idxC]+ matA_22 *conj( (UC->r0.c0[idxC] * UC->r1.c1[idxC] ) - ( UC->r0.c1[idxC] * UC->r1.c0[idxC]) );

}


void compute_sigma(__restrict thmat_soa * const L,  // la Lambda --> ouput  (una cosa che serve per calcolare la forza fermionica successiva)
		   __restrict su3_soa   * const U,  // la configurazione di gauge --> input
		   __restrict su3_soa   * const S,  // entra Sigma primo (input: fermforce del passo precedente) ED esce Sigma --> sia input che ouput
		   __restrict tamat_soa * const QA, // gli stessi Q che arrivano a Cayley hamilton --> input (sostanzialmente sono rho*ta(staples))
		   __restrict su3_soa   * const TMP // variabile di parcheggio
		   ){
  
  int x, y, z, t, mu, iter;

#pragma acc kernels present(L) present(U) present(nnp_openacc) present(nnm_openacc) present(S) present(QA) present(TMP)
 #pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32)
        for(x=0; x < nx; x++) {
     #pragma acc loop seq
          for(mu=0; mu<4; mu++){

	    const int indice = snum_acc(x,y,z,t); // anche questa e' calcolata pure dopo ...
	    const int parita = (x+y+z+t) % 2;
	    const int dir_mu = 2*mu + parita; // definisco una variabile di tipo dirlink pure qui (ce ne e' anche una dentro il loop dopo)
	    
	    // in questa routine faccio RES = PEZZO1
	    compute_sigma_local_PEZZO1(&L[dir_mu],&U[dir_mu],&S[dir_mu],&QA[dir_mu],&TMP[dir_mu],indice);
	    
      #pragma acc loop seq
            for(iter=0; iter<3; iter++){
              int nu;
              if (mu==0) { nu = iter + 1; }
              else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
              else if (mu==2) { nu = iter + (iter >> 1); }
              else if (mu==3) { nu = iter; }
              else { //error
              }
	      
              const int idxh = snum_acc(x,y,z,t);  // r
              const int parity = (x+y+z+t) % 2;
#pragma acc cache (nnp_openacc[idxh:8])

              const int dir_link = 2*mu + parity;
              const int dir_mu_2R = 2*mu + !parity;
              const int dir_mu_2L = 2*mu + !parity;
              const int idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])
              const int dir_nu_1R = 2*nu + !parity;
              const int dir_nu_3R = 2*nu +  parity;
              const int dir_nu_1L = 2*nu +  parity;
              const int dir_nu_3L = 2*nu + !parity;

	      const int idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu

              //  computation of the >>>Right<<< part of the staple with LAMBDAs
	      //  su3
	      //  A = U_nu(x+mu)
	      //  B = dagger U_mu(x+nu)
	      //  C = dagger U_nu(x)
	      //  thmat
	      //  D = LAMBDA_mu(x)
	      //  E = LAMBDA_nu(x)
	      //  F = LAMBDA_nu(x+mu)
	      //  G = LAMBDA_mu(x+nu)
	      //  iABCD and -iABCE
              RIGHT_iABC_times_DminusE_absent_stag_phases(&U[dir_nu_1R],       idx_pmu, // A
							  &U[dir_mu_2R],       idx_pnu, // B
							  &U[dir_nu_3R],       idxh,    // C
							  &L[dir_link],        idxh,    // D
							  &L[dir_nu_3R],       idxh,    // E
							  &S[dir_link],        idxh);
	      //  iFABC
              RIGHT_iFABC_absent_stag_phases(&U[dir_nu_1R],       idx_pmu, // A
					     &U[dir_mu_2R],       idx_pnu, // B
					     &U[dir_nu_3R],       idxh,    // C
					     &L[dir_nu_1R],       idx_pmu, // F
					     &S[dir_link],        idxh);
	      // -iABGC
	      RIGHT_miABGC_absent_stag_phases(&U[dir_nu_1R],       idx_pmu, // A
					      &U[dir_mu_2R],       idx_pnu, // B
					      &U[dir_nu_3R],       idxh,    // C
					      &L[dir_mu_2R],       idx_pnu, // G
					      &S[dir_link],        idxh);
	      
              const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
              const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu

              //  computation of the >>>Left<<< part of the staple with LAMBDAs
	      //  su3
	      //  A = dagger U_nu(x+mu-nu)
	      //  B = dagger U_mu(x-nu)
	      //  C = U_nu(x-nu)
	      //  thmat
	      //  D = LAMBDA_mu(x)
	      //  E = LAMBDA_mu(x-nu)
	      //  F = LAMBDA_nu(x+mu-nu)
	      //  G = LAMBDA_nu(x-nu)
	      //  iABGC and -iABEC

	      LEFT_iAB_times_GminusE_times_C_absent_stag_phases(&U[dir_nu_1L],       idx_pmu_mnu, // A
							   &U[dir_mu_2L],       idx_mnu,     // B
							   &U[dir_nu_3L],       idx_mnu,     // C
							   &L[dir_mu_2L],       idx_mnu,     // E
							   &L[dir_nu_3L],       idx_mnu,     // G
							   &S[dir_link],        idxh);

	      //  iABCD

	      LEFT_iABCD_absent_stag_phases(&U[dir_nu_1L],       idx_pmu_mnu, // A
					    &U[dir_mu_2L],       idx_mnu,     // B
					    &U[dir_nu_3L],       idx_mnu,     // C
					    &L[dir_link],        idxh,        // D
					    &S[dir_link],        idxh);
	      // -iAFBC
	      LEFT_miAFBC_absent_stag_phases(&U[dir_nu_1L],       idx_pmu_mnu, // A
					     &U[dir_mu_2L],       idx_mnu,     // B
					     &U[dir_nu_3L],       idx_mnu,     // C
					     &L[dir_nu_1L],       idx_pmu_mnu, // F
					     &S[dir_link],        idxh);
	      
            }  // iter
	    

          } // mu

        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine



#endif

