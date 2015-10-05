#ifndef CAYLEY_HAMILTON_C
#define CAYLEY_HAMILTON_C


// This is based on http://arXiv.org/abs/hep-lat/0311018v1,
// "Analytic Smearing of SU(3) Link Variables in Lattice QCD",
//  Morningstar & Peardon (2008)

//#include <math.h>
#include <complex.h>
#include "../OpenAcc/struct_c_def.c"

static inline void thmat_to_su3(single_su3 * out, single_thmat * Q){

    out->comp[0][0] =      Q->rc00 ; 
    out->comp[0][1] =      Q->c01 ;
    out->comp[0][2] =      Q->c02 ;
    out->comp[1][0] = conj(Q->c01) ;
    out->comp[1][1] =      Q->rc11 ; 
    out->comp[1][2] =      Q->c12 ;
    out->comp[2][0] = conj(Q->c02);
    out->comp[2][1] = conj(Q->c12);
    out->comp[2][2] =    - Q->rc00 - Q->rc11 ; 


}
static inline void i_times_tamat_to_su3(single_su3 * out, single_tamat * QA){
  //I*tamat is a thmat
    out->comp[0][0] =           - QA->rc00 ; 
    out->comp[0][1] =   (1.0*I) * QA->c01  ;
    out->comp[0][2] =   (1.0*I) * QA->c02  ;

    out->comp[1][0] =   (-1.0*I) * conj(QA->c01) ;
    out->comp[1][1] =           - QA->rc11 ; 
    out->comp[1][2] =   (1.0*I) * QA->c12  ;

    out->comp[2][0] =   (-1.0*I) * conj(QA->c02) ;
    out->comp[2][1] =   (-1.0*I) * conj(QA->c12) ;
    out->comp[2][2] =   QA->rc00 + QA->rc11; 


}

static inline void i_times_tamat_soa_to_su3(single_su3 * out, __restrict tamat_soa * const QA,const int idx){
  //I*tamat is a thmat
    out->comp[0][0] =           - QA->rc00[idx] ; 
    out->comp[0][1] =   (1.0*I) * QA->c01[idx]  ;
    out->comp[0][2] =   (1.0*I) * QA->c02[idx]  ;

    out->comp[1][0] =   (-1.0*I) * conj(QA->c01[idx]) ;
    out->comp[1][1] =           - QA->rc11[idx] ; 
    out->comp[1][2] =   (1.0*I) * QA->c12[idx]  ;

    out->comp[2][0] =   (-1.0*I) * conj(QA->c02[idx]) ;
    out->comp[2][1] =   (-1.0*I) * conj(QA->c12[idx]) ;
    out->comp[2][2] =   QA->rc00[idx] + QA->rc11[idx]; 


}

static inline double detQ(single_thmat *Q){

    double rc22 = -Q->rc00-Q->rc11 ; 
    return creal(Q->rc00*Q->rc11*rc22 + 2*creal(Q->c01*Q->c12*conj(Q->c02)) -
        Q->rc00 * Q->c12 * conj(Q->c12) - rc22 * Q->c01 * conj(Q->c01) - 
        Q->rc11 * Q->c02 * conj(Q->c02));

}

static inline double det_i_times_QA(single_tamat *QA){

    double rc22 = -QA->rc00-QA->rc11 ; 

    return -creal(  QA->rc00 * QA->rc11 * rc22
		  + 2*cimag(QA->c01 * QA->c12 * conj(QA->c02))
		  - QA->rc00 * QA->c12 * conj(QA->c12)
		  - QA->rc11 * QA->c02 * conj(QA->c02)
      	          - rc22 * QA->c01 * conj(QA->c01) );

}

static inline double det_i_times_QA_soa( __restrict tamat_soa * const QA,const int idx){

    double rc22 = -QA->rc00[idx]-QA->rc11[idx] ;

    return -creal(          QA->rc00[idx] *      QA->rc11[idx] * rc22
		  + 2*cimag(QA->c01[idx]  *      QA->c12[idx]  * conj(QA->c02[idx]))
		  -         QA->rc00[idx] *      QA->c12[idx]  * conj(QA->c12[idx])
		  -         QA->rc11[idx] *      QA->c02[idx]  * conj(QA->c02[idx])
      	          -         QA->c01[idx]  * conj(QA->c01[idx]) * rc22  );

}

static inline d_complex detSu3(single_su3 *m){


    return  
    m->comp[0][0]* m->comp[1][1] * m->comp[2][2] +
    m->comp[0][1]* m->comp[1][2] * m->comp[2][0] +
    m->comp[0][2]* m->comp[1][0] * m->comp[2][1] -
    m->comp[0][0]* m->comp[1][2] * m->comp[2][1] -
    m->comp[0][1]* m->comp[1][0] * m->comp[2][2] -
    m->comp[0][2]* m->comp[1][1] * m->comp[2][0] ;

}
static inline double TrQsq(single_thmat *Q){

  return 2 * creal( Q->rc00 * Q->rc00 +
                    Q->rc11 * Q->rc11 +
                    Q->c01 * conj(Q->c01) +
                    Q->c02 * conj(Q->c02) +
                    Q->c12 * conj(Q->c12) +
	            Q->rc00 * Q->rc11 );

}
static inline double Tr_i_times_QA_sq(single_tamat *QA){
  // computes Tr( (i*QA)^2 )
    return 2 * creal( QA->rc00 * QA->rc00 +
		      QA->rc11 * QA->rc11 +
		      QA->rc00 * QA->rc11 +
		      QA->c01  * conj(QA->c01) +
		      QA->c02  * conj(QA->c02) +
		      QA->c12  * conj(QA->c12) );
}

static inline double Tr_i_times_QA_sq_soa( __restrict tamat_soa * const QA,const int idx){
  // computes Tr( (i*QA)^2 )
    return 2 * creal( QA->rc00[idx] * QA->rc00[idx] +
		      QA->rc11[idx] * QA->rc11[idx] +
		      QA->rc00[idx] * QA->rc11[idx] +
		      QA->c01[idx]  * conj(QA->c01[idx]) +
		      QA->c02[idx]  * conj(QA->c02[idx]) +
		      QA->c12[idx]  * conj(QA->c12[idx]) );
}




static inline void single_su3_times_scalar(single_su3 * m , d_complex scalar){

   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++)
     m->comp[r][c] *= scalar;

}
static inline void single_su3xsu3(single_su3 * out , single_su3 *m1, single_su3 *m2){

   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++){
        out->comp[r][c] = 0;
        for(int d=0;d<3;d++) out->comp[r][c] += m1->comp[r][d] * m2->comp[d][c] ;

    }
}
void print_su3_stdout(single_su3 *m){

    printf("\n");
   for(int r=0;r<3;r++){
    for(int c=0;c<3;c++) printf("%.18lf + %.18lf I   ",  creal(m->comp[r][c]), cimag(m->comp[r][c]));
    printf("\n");

   }
}
static inline void single_su3add(single_su3 * out , single_su3 *m){

   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++)
        out->comp[r][c] += m->comp[r][c];

}

static inline void single_su3add_no3rdrow(single_su3 * out , single_su3 *m){

   for(int r=0;r<2;r++)
    for(int c=0;c<3;c++)
        out->comp[r][c] += m->comp[r][c];

}


void CH_exponential(single_su3 * out, single_thmat * Q ){ // exp(iQ)
    // based on Sez. III of http://arXiv.org/abs/hep-lat/0311018v1

    double c0 = detQ(Q); //(14)
    double c1  = 0.5 * TrQsq(Q); // (15)
    double c0max = 2*pow(c1/3,1.5); // (17)

    double theta = acos(c0/c0max);//(25)
    double u = sqrt(c1/3) * cos(theta/3) ;//(23)
    double w = sqrt(c1) * sin(theta/3) ;//(23)
    //    printf("c0:%f c1:%f \n" ,c0,c1);
    //    printf("u:%f w:%f \n" ,u,w);

    // see (23), 
    double xi1 = 1 - w*w/6*(1-w*w/20*(1-w*w/42));
    double xi2 = sin(w)/w;
    double xi0 = xi1 * (((int) (20*w-1) >> 31) & 0x1) +
                xi2 * (((int) (1-20*w) >> 31) & 0x1) ;

    d_complex h0 = (u*u - w*w) * cexp(2*u*I) + cexp(-u*I)*( // (30)
        8*u*u *cos(w) + 2 * u * (3*u*u+ w*w) * xi0 * I );

    d_complex h1 = 2*u*cexp(2*u*I) - cexp(-u*I)* (  // (31)
            2*u*cos(w) - (3*u*u-w*w)* xi0 * I ) ; 
    d_complex h2 = cexp(2*u*I) - cexp(-u*I)* (cos(w)+ 3*u*xi0*I); // (32)

    double denom = (9*u*u - w*w);
    
    single_su3 Q1, Q2;
    thmat_to_su3(&Q1,Q);
    single_su3xsu3(&Q2,&Q1,&Q1);

    single_su3_times_scalar(&Q1,h1/denom);
    single_su3_times_scalar(&Q2,h2/denom);

    d_complex f0 = h0/denom ;

    // first term in eq. 19
    out->comp[0][0] = f0;
    out->comp[0][1] = 0;
    out->comp[0][2] = 0;
    out->comp[1][0] = 0;
    out->comp[1][1] = f0;
    out->comp[1][2] = 0;
    out->comp[2][0] = 0;
    out->comp[2][1] = 0;
    out->comp[2][2] = f0;

//    print_su3_stdout(out);
//    print_su3_stdout(&Q1);
//    print_su3_stdout(&Q2);

    single_su3add(out, &Q1);// second term in (19)
    single_su3add(out, &Q2);// third

}

#pragma acc routine seq
static inline void CH_exponential_antihermitian(single_su3 * out, single_tamat * QA){
  // exp( - QA) , where Q=i*QA ==> exp(-QA) = exp(i*Q)
  //        ~~>  QA is anti-hermitian
  //        ~~>  hence Q=i*QA is hermitian
  //   based on Sez. III of http://arXiv.org/abs/hep-lat/0311018v1
  
  double c0 = det_i_times_QA(QA); //(14) // cosi calcolo direttamente il det(Q)
  double c1  = 0.5 * Tr_i_times_QA_sq(QA); // (15)
  double c0max = 2*pow(c1/3,1.5); // (17) // forse e' meglio mettere (c1/3)*sqrt(c1/3) ?!?!
  double theta = homebrew_acos(c0/c0max);//(25)


  double u = sqrt(c1/3) * cos(theta/3) ;//(23)
  double w = sqrt(c1) * sin(theta/3) ;//(23)

  //  printf("c0:%f c1:%f \n" ,c0,c1);
  //  printf("u:%f w:%f \n" ,u,w);
  
  // see (23), 
  double xi1 = 1 - w*w/6*(1-w*w/20*(1-w*w/42));
  double xi2 = sin(w)/w;
  double xi0 = xi1 * (((int) (20*w-1) >> 31) & 0x1) +
    xi2 * (((int) (1-20*w) >> 31) & 0x1) ;
  
  d_complex expmiu = cos(u) - sin(u)*I;
  d_complex exp2iu = cos(2.0*u) + sin(2.0*u)*I;


  d_complex h0 =   (u*u - w*w) * exp2iu + expmiu*( // (30)
		    8*u*u *cos(w) + 2 * u * (3*u*u+ w*w) * xi0 * I );
  
  d_complex h1 = 2*u*exp2iu - expmiu* (  // (31)
		 2*u*cos(w) - (3*u*u-w*w)* xi0 * I ) ; 
  d_complex h2 = exp2iu - expmiu* (cos(w)+ 3*u*xi0*I); // (32)
  
  double denom = (9*u*u - w*w);
  
  single_su3 Q1, Q2;
  i_times_tamat_to_su3(&Q1,QA);
  single_su3xsu3(&Q2,&Q1,&Q1);
  
  single_su3_times_scalar(&Q1,h1/denom);
  single_su3_times_scalar(&Q2,h2/denom);

  
  d_complex f0 = h0/denom ;

  // first term in eq. 19
  out->comp[0][0] = f0;
  out->comp[0][1] = 0;
  out->comp[0][2] = 0;
  out->comp[1][0] = 0;
  out->comp[1][1] = f0;
  out->comp[1][2] = 0;
  //  out->comp[2][0] = 0;
  //  out->comp[2][1] = 0;
  //  out->comp[2][2] = f0;
  
  //    print_su3_stdout(out);
  //    print_su3_stdout(&Q1);
  //    print_su3_stdout(&Q2);
  

  single_su3add_no3rdrow(out, &Q1);// second term in (19)
  single_su3add_no3rdrow(out, &Q2);// third

}



#pragma acc routine seq
static inline void CH_exponential_antihermitian_soa(__restrict su3_soa * const exp_out,  __restrict tamat_soa * const QA,const int idx){
  // exp( - QA) , where Q=i*QA ==> exp(-QA) = exp(i*Q)
  //        ~~>  QA is anti-hermitian
  //        ~~>  hence Q=i*QA is hermitian
  //   based on Sez. III of http://arXiv.org/abs/hep-lat/0311018v1

  double c0 = det_i_times_QA_soa(QA,idx); //(14) // cosi calcolo direttamente il det(Q)


  double c1  = 0.5 * Tr_i_times_QA_sq_soa(QA,idx); // (15)
  double c0max = 2*pow(c1/3,1.5); // (17) // forse e' meglio mettere (c1/3)*sqrt(c1/3) ?!?!
  double theta = homebrew_acos(c0/c0max);//(25)



  double u = sqrt(c1/3) * cos(theta/3) ;//(23)
  double w = sqrt(c1) * sin(theta/3) ;//(23)

  //  printf("c0:%f c1:%f \n" ,c0,c1);
  //  printf("u:%f w:%f \n" ,u,w);
  
  // see (23), 
  double xi1 = 1 - w*w/6*(1-w*w/20*(1-w*w/42));
  double xi2 = sin(w)/w;
  double xi0 = xi1 * (((int) (20*w-1) >> 31) & 0x1) +
    xi2 * (((int) (1-20*w) >> 31) & 0x1) ;
  
  d_complex expmiu = cos(u) - sin(u)*I;
  d_complex exp2iu = cos(2.0*u) + sin(2.0*u)*I;


  double denom = 1.0/(9*u*u - w*w);

  d_complex f0 =   denom * ((u*u - w*w) * exp2iu + expmiu*( // (30)
		   8*u*u *cos(w) + 2 * u * (3*u*u+ w*w) * xi0 * I ));
  
  d_complex f1 = denom * (2*u*exp2iu - expmiu* (  // (31)
		 2*u*cos(w) - (3*u*u-w*w)* xi0 * I )) ; 
  d_complex f2 = denom * (exp2iu - expmiu* (cos(w)+ 3*u*xi0*I)); // (32)
  
  

  //  single_su3 Q1, Q2;
  //  Q1.comp[0][0] = - QA->rc00[idx];
  /*
  i_times_tamat_soa_to_su3(&Q1,&QA,idx);

  single_su3xsu3(&Q2,&Q1,&Q1);
  */
  
  //  single_su3_times_scalar(&Q1,h1/denom);
  //  single_su3_times_scalar(&Q2,h2/denom);

  

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
                             + f2    * ( QA->rc11[idx] * QA->rc11[idx])
                                       + QA->c01[idx] * conj(QA->c01[idx])
                                       + QA->c12[idx] * conj(QA->c12[idx]);

  exp_out->r1.c2[idx] =      ( f1*I) * QA->c12[idx]        
                             + f2    * ( QA->rc00[idx] * QA->c12[idx] 
					 + QA->c02[idx] * conj(QA->c01[idx]));



  //  out->comp[2][0] = 0;
  //  out->comp[2][1] = 0;
  //  out->comp[2][2] = f0;
  
  //    print_su3_stdout(out);
  //    print_su3_stdout(&Q1);
  //    print_su3_stdout(&Q2);

  
  //  single_su3add_no3rdrow(out, &Q1);// second term in (19)
  //  single_su3add_no3rdrow(out, &Q2);// third



}





#endif
