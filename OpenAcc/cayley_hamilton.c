// This is based on http://arXiv.org/abs/hep-lat/0311018v1,
// "Analytic Smearing of SU(3) Link Variables in Lattice QCD",
//  Morningstar & Peardon (2008)

#include <math.h>
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
static inline double detQ(single_thmat *Q){

    double rc22 = -Q->rc00-Q->rc11 ; 
    return creal(Q->rc00*Q->rc11*rc22 + 2*creal(Q->c01*Q->c12*conj(Q->c02)) -
        Q->rc00 * Q->c12 * conj(Q->c12) - rc22 * Q->c01 * conj(Q->c01) - 
        Q->rc11 * Q->c02 * conj(Q->c02));

}
static inline double detSu3(single_su3 *m){


    return  
    m->comp[0][0]* m->comp[1][1] * m->comp[2][2] +
    m->comp[0][1]* m->comp[1][2] * m->comp[2][0] +
    m->comp[0][2]* m->comp[1][0] * m->comp[2][1] -
    m->comp[0][0]* m->comp[1][2] * m->comp[2][1] -
    m->comp[0][1]* m->comp[1][0] * m->comp[2][2] -
    m->comp[0][2]* m->comp[1][1] * m->comp[2][0] ;

}
static inline double TrQsq(single_thmat *Q){

    return 2* ( Q->rc00 * Q->rc00 +
         Q->rc11 * Q->rc11 +
         Q->c01 * conj(Q->c01) +
         Q->c02 * conj(Q->c02) +
         Q->c12 * conj(Q->c12) +
         Q->rc00 * Q->rc11 ); 

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
    for(int c=0;c<3;c++) printf("%f + %f I   ",  creal(m->comp[r][c]), cimag(m->comp[r][c]));
    printf("\n");

   }
}
static inline void single_su3add(single_su3 * out , single_su3 *m){

   for(int r=0;r<3;r++)
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
    printf("c0:%f c1:%f \n" ,c0,c1);
    printf("u:%f w:%f \n" ,u,w);

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
