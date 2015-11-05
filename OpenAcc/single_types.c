#ifndef SINGLE_TYPES_H
#define SINGLE_TYPES_H

#include <complex.h>
#include "./struct_c_def.c"

/****************************************************
 * (Perhaps) Check if adding                        *
 *                                                  *
 * #pragma acc routine seq                          *
 *                                                  *
 *before each function can improve the performances.*
 ****************************************************/

//Type definitions
typedef struct single_su3_t {
  d_complex comp[3][3];
} single_su3;
typedef struct single_tamat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double rc00;   // Im(comp_00)
  double rc11;   // Im(comp_11)
} single_tamat;
typedef struct single_thmat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double rc00;   // Re(comp_00)
  double rc11;   // Re(comp_11)
} single_thmat;

// 3rd row reconstruction on site
void static inline rebuild3row(single_su3 *AUX){

  AUX->comp[2][0] = conj(AUX->comp[0][1] * AUX->comp[1][2] - AUX->comp[0][2] * AUX->comp[1][1]);
  AUX->comp[2][1] = conj(AUX->comp[0][2] * AUX->comp[1][0] - AUX->comp[0][0] * AUX->comp[1][2]);
  AUX->comp[2][2] = conj(AUX->comp[0][0] * AUX->comp[1][1] - AUX->comp[0][1] * AUX->comp[1][0]);

}
//Conversion functions
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
static inline void tamat_to_su3(single_su3 * out, single_thmat * QA){

    out->comp[0][0] =     I* QA->rc00 ; 
    out->comp[0][1] =      QA->c01 ;
    out->comp[0][2] =      QA->c02 ;
    out->comp[1][0] = -conj(QA->c01) ;
    out->comp[1][1] =     I* QA->rc11 ; 
    out->comp[1][2] =      QA->c12 ;
    out->comp[2][0] = -conj(QA->c02);
    out->comp[2][1] = -conj(QA->c12);
    out->comp[2][2] =    (- QA->rc00 - QA->rc11)*I ; 


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
static inline void gl3_to_tamat(single_su3 * in , single_tamat * out){

     double im_trace_third = (cimag(in->comp[0][0]) +  cimag(in->comp[1][1]) + cimag(in->comp[2][2]))/3;

     out->c01 = 0.5*(in->comp[0][1] - conj(in->comp[1][0]));
     out->c02 = 0.5*(in->comp[0][2] - conj(in->comp[2][0]));
     out->c12 = 0.5*(in->comp[1][2] - conj(in->comp[2][1]));

     out->rc00 = cimag(in->comp[0][0])-im_trace_third;
     out->rc11 = cimag(in->comp[1][1])-im_trace_third;

}
static inline void gl3_to_thmat(single_su3 * in , single_thmat * out){

     double re_trace_third = (creal(in->comp[0][0]) +  creal(in->comp[1][1]) + creal(in->comp[2][2]))/3;

     out->c01 = 0.5*(in->comp[0][1]+ conj(in->comp[1][0]));
     out->c02 = 0.5*(in->comp[0][2]+ conj(in->comp[2][0]));
     out->c12 = 0.5*(in->comp[1][2]+ conj(in->comp[2][1]));

     out->rc00 = creal(in->comp[0][0])-re_trace_third;
     out->rc11 = creal(in->comp[1][1])-re_trace_third;

}

//extraction from soas
static inline void single_su3_from_su3_soa( __restrict su3_soa * const mat, const int idx_mat,						  single_su3 * Omat){
  Omat->comp[0][0] = mat->r0.c0[idx_mat];
  Omat->comp[0][1] = mat->r0.c1[idx_mat];
  Omat->comp[0][2] = mat->r0.c2[idx_mat];
  Omat->comp[1][0] = mat->r1.c0[idx_mat];
  Omat->comp[1][1] = mat->r1.c1[idx_mat];
  Omat->comp[1][2] = mat->r1.c2[idx_mat];

}
static inline void single_gl3_from_su3_soa( __restrict su3_soa * const mat, const int idx_mat, single_su3 * Omat){
    single_su3_from_su3_soa(mat,idx_mat,Omat);
    Omat->comp[2][0] = mat->r2.c0[idx_mat];
    Omat->comp[2][1] = mat->r2.c1[idx_mat];
    Omat->comp[2][2] = mat->r2.c2[idx_mat];

}
static inline void single_tamat_from_tamat_soa(__restrict tamat_soa * in, int idx, single_tamat * out){

    out->rc00 = in->rc00[idx];
    out->rc11 = in->rc11[idx];
    out->c01 = in->c01[idx];
    out->c02 = in->c02[idx];
    out->c12 = in->c12[idx];

}
static inline void single_thmat_from_thmat_soa(__restrict thmat_soa * in, int idx, single_thmat * out){

    out->rc00 = in->rc00[idx];
    out->rc11 = in->rc11[idx];
    out->c01 = in->c01[idx];
    out->c02 = in->c02[idx];
    out->c12 = in->c12[idx];

}



//insertion into soas
static inline void single_su3_into_su3_soa( __restrict su3_soa * const mat, const int idx_mat,  single_su3 * Imat){
  mat->r0.c0[idx_mat] = Imat->comp[0][0];
  mat->r0.c1[idx_mat] = Imat->comp[0][1];
  mat->r0.c2[idx_mat] = Imat->comp[0][2];
  mat->r1.c0[idx_mat] = Imat->comp[1][0];
  mat->r1.c1[idx_mat] = Imat->comp[1][1];
  mat->r1.c2[idx_mat] = Imat->comp[1][2];
}
static inline void single_su3_addinto_su3_soa( __restrict su3_soa * const mat, const int idx_mat,  single_su3 * Imat){
  mat->r0.c0[idx_mat] += Imat->comp[0][0];
  mat->r0.c1[idx_mat] += Imat->comp[0][1];
  mat->r0.c2[idx_mat] += Imat->comp[0][2];
  mat->r1.c0[idx_mat] += Imat->comp[1][0];
  mat->r1.c1[idx_mat] += Imat->comp[1][1];
  mat->r1.c2[idx_mat] += Imat->comp[1][2];
}
static inline void single_gl3_into_su3_soa( __restrict su3_soa * const mat, const int idx_mat,						   single_su3 * Imat){

 single_su3_into_su3_soa(mat,idx_mat,Imat);
  mat->r2.c0[idx_mat] = Imat->comp[2][0];
  mat->r2.c1[idx_mat] = Imat->comp[2][1];
  mat->r2.c2[idx_mat] = Imat->comp[2][2];
}
static inline void single_gl3_addinto_su3_soa( __restrict su3_soa * const mat, const int idx_mat,						   single_su3 * Imat){

 single_su3_addinto_su3_soa(mat,idx_mat,Imat);
  mat->r2.c0[idx_mat] += Imat->comp[2][0];
  mat->r2.c1[idx_mat] += Imat->comp[2][1];
  mat->r2.c2[idx_mat] += Imat->comp[2][2];
}
static inline void single_tamat_into_tamat_soa(__restrict tamat_soa * out, int idx, single_tamat * in){

    out->rc00[idx] = in->rc00;
    out->rc11[idx] = in->rc11;
    out->c01[idx]  = in->c01  ;
    out->c02[idx]  = in->c02  ;
    out->c12[idx]  = in->c12  ;

}
static inline void single_thmat_into_thmat_soa(__restrict thmat_soa * out, int idx, single_thmat * in){

    out->rc00[idx] = in->rc00;
    out->rc11[idx] = in->rc11;
    out->c01[idx]  = in->c01  ;
    out->c02[idx]  = in->c02  ;
    out->c12[idx]  = in->c12  ;

}


//misc. extraction and conversion
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


// functions
static inline void herm_conj_gl3(single_su3 * in, single_su3 *out){

    out->comp[0][0] = conj(in->comp[0][0]);
    out->comp[0][1] = conj(in->comp[1][0]);
    out->comp[0][2] = conj(in->comp[2][0]);
    out->comp[1][0] = conj(in->comp[0][1]);
    out->comp[1][1] = conj(in->comp[1][1]);
    out->comp[1][2] = conj(in->comp[2][1]);
    out->comp[2][0] = conj(in->comp[0][2]);
    out->comp[2][1] = conj(in->comp[1][2]);
    out->comp[2][2] = conj(in->comp[2][2]);

}
static inline double detQ(single_thmat *Q){

    double rc22 = -Q->rc00-Q->rc11 ; 
    return creal(Q->rc00*Q->rc11*rc22 + 2*creal(Q->c01*Q->c12*conj(Q->c02)) -
        Q->rc00 * Q->c12 * conj(Q->c12) - rc22 * Q->c01 * conj(Q->c01) - 
        Q->rc11 * Q->c02 * conj(Q->c02));

}
static inline double det_i_times_QA(__restrict single_tamat * const QA){

    double rc22 = -QA->rc00-QA->rc11 ; 

    return -creal(  QA->rc00 * QA->rc11 * rc22
		  + 2*cimag(QA->c01 * QA->c12 * conj(QA->c02))
		  - QA->rc00 * QA->c12 * conj(QA->c12)
		  - QA->rc11 * QA->c02 * conj(QA->c02)
      	          - rc22 * QA->c01 * conj(QA->c01) );

}
void print_su3_stdout(single_su3 *m){

    printf("\n");
   for(int r=0;r<3;r++){
    for(int c=0;c<3;c++) printf("%.18lf + %.18lf I   ",  creal(m->comp[r][c]), cimag(m->comp[r][c]));
    printf("\n");

   }
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
static inline void single_su3xsu3_add_to_out(single_su3 * out , single_su3 *m1, single_su3 *m2){

   for(int r=0;r<3;r++)
    for(int c=0;c<3;c++){
        //   out->comp[r][c] = 0; // add to out!
        for(int d=0;d<3;d++) out->comp[r][c] += m1->comp[r][d] * m2->comp[d][c] ;

    }
}
static inline void single_su3add(single_su3 * out , single_su3 *m){

   for(int r=0;r<3;r++)// Magari fino alla seconda riga?
   //for(int r=0;r<2;r++) //??
    for(int c=0;c<3;c++)
        out->comp[r][c] += m->comp[r][c];

}
static inline double TrQsq(single_thmat *Q){

  return 2 * creal( Q->rc00 * Q->rc00 +
                    Q->rc11 * Q->rc11 +
                    Q->c01 * conj(Q->c01) +
                    Q->c02 * conj(Q->c02) +
                    Q->c12 * conj(Q->c12) +
	            Q->rc00 * Q->rc11 );

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
static inline double Tr_i_times_QA_sq(__restrict single_tamat * const QA){
  // computes Tr( (i*QA)^2 )
    return 2 * creal( QA->rc00 * QA->rc00 +
		      QA->rc11 * QA->rc11 +
		      QA->rc00 * QA->rc11 +
		      QA->c01  * conj(QA->c01) +
		      QA->c02  * conj(QA->c02) +
		      QA->c12  * conj(QA->c12) );
}
static inline void single_thmatAeqAmB(single_thmat *A, single_thmat* B){

    A->rc00 -= B->rc00;
    A->rc11 -= B->rc11;
    A->c01 -= B->c01;
    A->c02 -= B->c02;
    A->c12 -= B->c12;

}

static inline void gl3_dagger(single_su3 * inout){

    inout->comp[0][0] = conj(inout->comp[0][0]);
    inout->comp[1][1] = conj(inout->comp[1][1]);
    inout->comp[2][2] = conj(inout->comp[2][2]);

    d_complex temp;

    temp = inout->comp[0][1]; inout->comp[0][1] = conj(inout->comp[1][0]);
    inout->comp[1][0] = conj(temp);

    temp = inout->comp[0][2]; inout->comp[0][2] = conj(inout->comp[2][0]);
    inout->comp[2][0] = conj(temp);

    temp = inout->comp[2][1]; inout->comp[2][1] = conj(inout->comp[1][2]);
    inout->comp[1][2] = conj(temp);


}
static inline void single_su3_copy(const single_su3 * const in, single_su3 * out){

    for(int i =0;i<3;i++)    
        for(int j =0;j<3;j++)    
            out->comp[i][j] = in->comp[i][j];

}


//no 3rd row
static inline void single_su3xsu3_no3rdrow(single_su3 * out , single_su3 *m1, single_su3 *m2){

   for(int r=0;r<2;r++)
    for(int c=0;c<3;c++){
        out->comp[r][c] = 0;
        for(int d=0;d<3;d++) out->comp[r][c] += m1->comp[r][d] * m2->comp[d][c] ;

    }
}
static inline void single_su3_times_scalar_no3rdrow(single_su3 * m , d_complex scalar){

   for(int r=0;r<2;r++)
    for(int c=0;c<3;c++)
     m->comp[r][c] *= scalar;

}
static inline void single_su3add_no3rdrow(single_su3 * out , single_su3 *m){

   for(int r=0;r<2;r++)
    for(int c=0;c<3;c++)
        out->comp[r][c] += m->comp[r][c];

}
#endif
