#ifndef CAYLEY_HAMILTON_C
#define CAYLEY_HAMILTON_C


// This is based on http://arXiv.org/abs/hep-lat/0311018v1,
// "Analytic Smearing of SU(3) Link Variables in Lattice QCD",
//  Morningstar & Peardon (2008)

#include <accelmath.h>
#include <complex.h>
#include "./struct_c_def.c"
#include "./single_types.c"
#include "../DbgTools/debug_macros_glvarcheck.c"


inline void Itamat_2ndDeg_poly_no3rdrow(d_complex f0, d_complex f1, d_complex f2, __restrict single_tamat * const QA, single_su3 * const out){
// this is a degree 2 polynomial in Q= iQA, that is
// P = f0 + f1 * Q + f2 * Q^2 = f0 + i * f1 * QA - f2 * QA^2


    out->comp[0][0] = f0 -   f1    * QA->rc00
        + f2 * (  QA->rc00 *      QA->rc00
                + QA->c01  * conj(QA->c01)
                + QA->c02  * conj(QA->c02));

    out->comp[0][1] =      ( f1*I) *  QA->c01
        + f2    * (QA->c02 * conj(QA->c12)
                + (-1.0*I)* QA->c01 * ( QA->rc00 + QA->rc11));


    out->comp[0][2] =      ( f1*I) * QA->c02        
        + f2    * (-QA->c01 * QA->c12
                + ( 1.0*I)* QA->c02 * QA->rc11);


    out->comp[1][0] =      (-f1*I) * conj(QA->c01)  
        + f2    * ( QA->c12 * conj(QA->c02)
                + ( 1.0*I) * conj(QA->c01) * ( QA->rc00 + QA->rc11));


    out->comp[1][1] = f0 -   f1    * QA->rc11       
        + f2    * ( QA->rc11 * QA->rc11
                + QA->c01 * conj(QA->c01)
                + QA->c12 * conj(QA->c12));

    out->comp[1][2] =      ( f1*I) * QA->c12        
        + f2 * ( (1.0*I)*QA->rc00 * QA->c12 
                + QA->c02 * conj(QA->c01));

}
inline void Itamat_2ndDeg_poly(d_complex f0, d_complex f1, d_complex f2, single_tamat * QA, single_su3 * out){
   // See comments in Itamat_2ndDeg_poly_no3rdrow(f0,f1,f2,QA,out);
   
    Itamat_2ndDeg_poly_no3rdrow(f0,f1,f2,QA,out);

    out->comp[2][0] = (f1*I) * (-conj(QA->c02)) 
        - f2 * conj(QA->c01 * QA->c12 - QA->c02 * QA->rc11 * I );

    out->comp[2][1] = (f1*I) * (-conj(QA->c12))
        - f2 * (-conj(QA->c02) * QA->c01 + conj(QA->c12)* QA->rc00 * I );


    out->comp[2][2] = f0 + f1*(QA->rc00 + QA->rc11)
        +  f2 * (QA->c02 * conj(QA->c02) + QA->c12 * conj(QA->c12)
               +  (QA->rc00 + QA->rc11 ) *  (QA->rc00 + QA->rc11 ) );


}

#pragma acc routine seq
static inline double det_i_times_QA_soa( __restrict tamat_soa * const QA,const int idx){

    double rc22 = -QA->rc00[idx]-QA->rc11[idx] ;

    return -creal(          QA->rc00[idx] *      QA->rc11[idx] * rc22
		  + 2*cimag(QA->c01[idx]  *      QA->c12[idx]  * conj(QA->c02[idx]))
		  -         QA->rc00[idx] *      QA->c12[idx]  * conj(QA->c12[idx])
		  -         QA->rc11[idx] *      QA->c02[idx]  * conj(QA->c02[idx])
      	          -         QA->c01[idx]  * conj(QA->c01[idx]) * rc22  );

}
#pragma acc routine seq
static inline double Tr_i_times_QA_sq_soa( __restrict tamat_soa * const QA,const int idx){
  // computes Tr( (i*QA)^2 )
    return 2 * creal( QA->rc00[idx] * QA->rc00[idx] +
		      QA->rc11[idx] * QA->rc11[idx] +
		      QA->rc00[idx] * QA->rc11[idx] +
		      QA->c01[idx]  * conj(QA->c01[idx]) +
		      QA->c02[idx]  * conj(QA->c02[idx]) +
		      QA->c12[idx]  * conj(QA->c12[idx]) );
}

#pragma acc routine seq
static inline void CH_exponential_antihermitian_soa_nissalike(__restrict su3_soa * const exp_out,
							      __restrict tamat_soa * const QA,const int idx){
  // exp( - QA) , where Q=i*QA ==> exp(-QA) = exp(i*Q)
  //        ~~>  QA is anti-hermitian
  //        ~~>  hence Q=i*QA is hermitian
  //   based on Sez. III of http://arXiv.org/abs/hep-lat/0311018v1

  double c0 = det_i_times_QA_soa(QA,idx); //(14) // cosi calcolo direttamente il det(Q)
  double c1  = 0.5 * Tr_i_times_QA_sq_soa(QA,idx); // (15)
  double c0max = 2*pow(c1/3,1.5); // (17) // forse e' meglio mettere (c1/3)*sqrt(c1/3) ?!?!

  d_complex f0,f1,f2;
  //consider the case in which c1<4*10^-3 apart, as done in MILC [TAKEN FROM NISSA]
  if(c1<4e-3)
    {
      f0 = (1-c0*c0/720) + (1.0*I)*(-c0*(1-c1*(1-c1/42)/20)/6);
      f1 = (c0*(1-c1*(1-3*c1/112)/15)/24) + (1.0*I)*(1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040);
      f2 = (0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160)) + (1.0*I)*(0.5*(c0*(1-c1*(1-c1/48)/21)/60));
    }
  else
    {
      //take c0 module and write separately its sign (see note after eq. 34)
      int segno=1;
      if(c0<0)
	{
	  segno=-1;
	  c0=-c0;
	}
      
      //check rounding error
      double eps=(c0max-c0)/c0max;
      //(eqs. 23-24)
      double theta;
      if(eps<0) theta=0.0; //only possible as an effect of rounding error when c0/c0_max=1
      else
	if(eps<1e-3) theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
	else theta=acos(c0/c0max);
      double u = sqrt(c1/3) * cos(theta/3) ;//(23)
      double w = sqrt(c1) * sin(theta/3) ;//(23)
      
      //auxiliary variables for the computation of h0, h1, h2 
      double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
      double cu=cos(u),c2u=cos(2*u);
      double su=sin(u),s2u=sin(2*u);
      double cw=cos(w);

      //compute xi function defined after (eq. 33)
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
  // first term in eq. 19
  exp_out->r0.c0[idx] = f0 -   f1    * QA->rc00[idx]
                      + f2 * (  QA->rc00[idx] *      QA->rc00[idx]
			      + QA->c01[idx]  * conj(QA->c01[idx])
		              + QA->c02[idx]  * conj(QA->c02[idx]));

  exp_out->r0.c1[idx] =      ( f1*I) *  QA->c01[idx]
                             + f2    * (QA->c02[idx] * conj(QA->c12[idx])
                  	                + (-1.0*I)* QA->c01[idx] * ( QA->rc00[idx] + QA->rc11[idx]));


  exp_out->r0.c2[idx] =      ( f1*I) * QA->c02[idx]        
                             + f2    * (-QA->c01[idx] * QA->c12[idx]
					+ ( 1.0*I)* QA->c02[idx] * QA->rc11[idx]);


  exp_out->r1.c0[idx] =      (-f1*I) * conj(QA->c01[idx])  
                             + f2    * ( QA->c12[idx] * conj(QA->c02[idx])
				        + ( 1.0*I) * conj(QA->c01[idx]) * ( QA->rc00[idx] + QA->rc11[idx]));


  exp_out->r1.c1[idx] = f0 -   f1    * QA->rc11[idx]       
                             + f2    * ( QA->rc11[idx] * QA->rc11[idx]
                                       + QA->c01[idx] * conj(QA->c01[idx])
					 + QA->c12[idx] * conj(QA->c12[idx]));

  exp_out->r1.c2[idx] =      ( f1*I) * QA->c12[idx]        
                                        + f2 * ( (1.0*I)*QA->rc00[idx] * QA->c12[idx] 
					+ QA->c02[idx] * conj(QA->c01[idx]));


  /*
  if(idx==0){
    printf("c0 = %.18lf\n",c0);
    printf("c1 = %.18lf\n",c1);
    printf("f0 = (%.18lf) + (%.18lf)*I\n",creal(f0),cimag(f0));
    printf("f1 = (%.18lf) + (%.18lf)*I\n",creal(f1),cimag(f1));
    printf("f2 = (%.18lf) + (%.18lf)*I\n",creal(f2),cimag(f2));
    printf("\n");

  }
  */

}

#pragma acc routine seq
static inline void CH_exponential_antihermitian_nissalike(single_su3 * const exp_out,
							      __restrict single_tamat * const QA){
  // exp( - QA) , where Q=i*QA ==> exp(-QA) = exp(i*Q)
  //        ~~>  QA is anti-hermitian
  //        ~~>  hence Q=i*QA is hermitian
  //   based on Sez. III of http://arXiv.org/abs/hep-lat/0311018v1



  double c0 = det_i_times_QA(QA); //(14) // cosi calcolo direttamente il det(Q)
  double c1  = 0.5 * Tr_i_times_QA_sq(QA); // (15)
  double c0max = 2*pow(c1/3,1.5); // (17) // forse e' meglio mettere (c1/3)*sqrt(c1/3) ?!?!

  d_complex f0,f1,f2;
  //consider the case in which c1<4*10^-3 apart, as done in MILC [TAKEN FROM NISSA]
  if(c1<4e-3)
    {
      f0 = (1-c0*c0/720) + (1.0*I)*(-c0*(1-c1*(1-c1/42)/20)/6);
      f1 = (c0*(1-c1*(1-3*c1/112)/15)/24) + (1.0*I)*(1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040);
      f2 = (0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160)) + (1.0*I)*(0.5*(c0*(1-c1*(1-c1/48)/21)/60));
    }
  else
    {
      //take c0 module and write separately its sign (see note after eq. 34)
      int segno=1;
      if(c0<0)
	{
	  segno=-1;
	  c0=-c0;
	}
      
      //check rounding error
      double eps=(c0max-c0)/c0max;
      //(eqs. 23-24)
      double theta;
      if(eps<0) theta=0.0; //only possible as an effect of rounding error when c0/c0_max=1
      else
	if(eps<1e-3) theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
	else theta=acos(c0/c0max);
      double u = sqrt(c1/3) * cos(theta/3) ;//(23)
      double w = sqrt(c1) * sin(theta/3) ;//(23)
      
      //auxiliary variables for the computation of h0, h1, h2 
      double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
      double cu=cos(u),c2u=cos(2*u);
      double su=sin(u),s2u=sin(2*u);
      double cw=cos(w);

      //compute xi function defined after (eq. 33)
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
  // first term in eq. 19

  
  Itamat_2ndDeg_poly_no3rdrow(f0,f1,f2,QA,exp_out);

  /*
    printf("c0 = %.18lf\n",c0);
    printf("c1 = %.18lf\n",c1);
    printf("f0 = (%.18lf) + (%.18lf)*I\n",creal(f0),cimag(f0));
    printf("f1 = (%.18lf) + (%.18lf)*I\n",creal(f1),cimag(f1));
    printf("f2 = (%.18lf) + (%.18lf)*I\n",creal(f2),cimag(f2));
    printf("\n");

  */


}
#endif



