#ifndef FERMION_MATRIX_C
#define FERMION_MATRIX_C

#include "../Include/common_defines.h"
#include "./geometry.h"
#include "./struct_c_def.h"
#ifndef __GNUC__
#include "openacc.h"
#endif
#include "./fermionic_utilities.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "./fermion_matrix.h"




static inline vec3 mat_vec_mul( __restrict su3_soa * const matrix,
        const int idx_mat,
        __restrict vec3_soa * const in_vect,
        const int idx_vect,
        d_complex phase){
    vec3 out_vect;

    d_complex vec0 = (in_vect->c0[idx_vect])*phase;
    d_complex vec1 = (in_vect->c1[idx_vect])*phase;
    d_complex vec2 = (in_vect->c2[idx_vect])*phase;

    d_complex mat00 = matrix->r0.c0[idx_mat];
    d_complex mat01 = matrix->r0.c1[idx_mat];
    d_complex mat02 = matrix->r0.c2[idx_mat];

    d_complex mat10 = matrix->r1.c0[idx_mat];
    d_complex mat11 = matrix->r1.c1[idx_mat];
    d_complex mat12 = matrix->r1.c2[idx_mat];

    // Load 3rd matrix row from global memory
    //  d_complex mat20 = matrix->r2.c0[idx_mat];
    //  d_complex mat21 = matrix->r2.c1[idx_mat];
    //  d_complex mat22 = matrix->r2.c2[idx_mat];

    //Compute 3rd matrix row from the first two  
    d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) ) ;
    d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) ) ;
    d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) ) ;

    out_vect.c0 = ( mat00 * vec0 ) + ( mat01 * vec1 ) + ( mat02 * vec2 );
    out_vect.c1 = ( mat10 * vec0 ) + ( mat11 * vec1 ) + ( mat12 * vec2 );
    out_vect.c2 = ( mat20 * vec0 ) + ( mat21 * vec1 ) + ( mat22 * vec2 );

    return out_vect;

}

static inline vec3 conjmat_vec_mul( __restrict su3_soa * const matrix,
        const int idx_mat,
        __restrict vec3_soa * const in_vect,
        const int idx_vect,
        d_complex phase) {
    vec3 out_vect;


    d_complex vec0 = in_vect->c0[idx_vect]*conj(phase);
    d_complex vec1 = in_vect->c1[idx_vect]*conj(phase);
    d_complex vec2 = in_vect->c2[idx_vect]*conj(phase);


    d_complex mat00 = matrix->r0.c0[idx_mat];
    d_complex mat01 = matrix->r0.c1[idx_mat];
    d_complex mat02 = matrix->r0.c2[idx_mat];

    d_complex mat10 = matrix->r1.c0[idx_mat];
    d_complex mat11 = matrix->r1.c1[idx_mat];
    d_complex mat12 = matrix->r1.c2[idx_mat];

    // Load 3rd matrix row from global memory
    //  d_complex mat20 = matrix->r2.c0[idx_mat];
    //  d_complex mat21 = matrix->r2.c1[idx_mat];
    //  d_complex mat22 = matrix->r2.c2[idx_mat];

    //Compute 3rd matrix row from the first two  
    d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) );
    d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) );
    d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) );

    out_vect.c0 = ( conj(mat00) * vec0 ) + ( conj(mat10) * vec1 ) + ( conj(mat20) * vec2 );
    out_vect.c1 = ( conj(mat01) * vec0 ) + ( conj(mat11) * vec1 ) + ( conj(mat21) * vec2 );
    out_vect.c2 = ( conj(mat02) * vec0 ) + ( conj(mat12) * vec1 ) + ( conj(mat22) * vec2 );

    return out_vect;

}

static inline vec3 mat_vec_mul_arg( __restrict su3_soa * const matrix,
        const int idx_mat,
        __restrict vec3_soa * const in_vect,
        const int idx_vect,
        __restrict double_soa * arrarg){
    vec3 out_vect;
    double arg = arrarg->d[idx_mat];
    d_complex phase = cos(arg) + I * sin(arg);

    d_complex vec0 = (in_vect->c0[idx_vect])*phase;
    d_complex vec1 = (in_vect->c1[idx_vect])*phase;
    d_complex vec2 = (in_vect->c2[idx_vect])*phase;

    d_complex mat00 = matrix->r0.c0[idx_mat];
    d_complex mat01 = matrix->r0.c1[idx_mat];
    d_complex mat02 = matrix->r0.c2[idx_mat];

    d_complex mat10 = matrix->r1.c0[idx_mat];
    d_complex mat11 = matrix->r1.c1[idx_mat];
    d_complex mat12 = matrix->r1.c2[idx_mat];

    // Load 3rd matrix row from global memory
    //  d_complex mat20 = matrix->r2.c0[idx_mat];
    //  d_complex mat21 = matrix->r2.c1[idx_mat];
    //  d_complex mat22 = matrix->r2.c2[idx_mat];

    //Compute 3rd matrix row from the first two  
    d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) ) ;
    d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) ) ;
    d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) ) ;

    out_vect.c0 = ( mat00 * vec0 ) + ( mat01 * vec1 ) + ( mat02 * vec2 );
    out_vect.c1 = ( mat10 * vec0 ) + ( mat11 * vec1 ) + ( mat12 * vec2 );
    out_vect.c2 = ( mat20 * vec0 ) + ( mat21 * vec1 ) + ( mat22 * vec2 );

    return out_vect;

}

static inline vec3 conjmat_vec_mul_arg( __restrict su3_soa * const matrix,
        const int idx_mat,
        __restrict vec3_soa * const in_vect,
        const int idx_vect,
        __restrict double_soa* arrarg  ) {
    vec3 out_vect;
    double arg = arrarg->d[idx_mat];
    d_complex phase = cos(arg) + I * sin(arg);


    d_complex vec0 = in_vect->c0[idx_vect]*conj(phase);
    d_complex vec1 = in_vect->c1[idx_vect]*conj(phase);
    d_complex vec2 = in_vect->c2[idx_vect]*conj(phase);


    d_complex mat00 = matrix->r0.c0[idx_mat];
    d_complex mat01 = matrix->r0.c1[idx_mat];
    d_complex mat02 = matrix->r0.c2[idx_mat];

    d_complex mat10 = matrix->r1.c0[idx_mat];
    d_complex mat11 = matrix->r1.c1[idx_mat];
    d_complex mat12 = matrix->r1.c2[idx_mat];

    // Load 3rd matrix row from global memory
    //  d_complex mat20 = matrix->r2.c0[idx_mat];
    //  d_complex mat21 = matrix->r2.c1[idx_mat];
    //  d_complex mat22 = matrix->r2.c2[idx_mat];

    //Compute 3rd matrix row from the first two  
    d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) );
    d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) );
    d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) );

    out_vect.c0 = ( conj(mat00) * vec0 ) + ( conj(mat10) * vec1 ) + ( conj(mat20) * vec2 );
    out_vect.c1 = ( conj(mat01) * vec0 ) + ( conj(mat11) * vec1 ) + ( conj(mat21) * vec2 );
    out_vect.c2 = ( conj(mat02) * vec0 ) + ( conj(mat12) * vec1 ) + ( conj(mat22) * vec2 );

    return out_vect;

}


static inline vec3 sumResult ( vec3 aux, vec3 aux_tmp) {

    aux.c0 += aux_tmp.c0;
    aux.c1 += aux_tmp.c1;
    aux.c2 += aux_tmp.c2;

    return aux;

}

static inline vec3 subResult ( vec3 aux, vec3 aux_tmp) {

    aux.c0 -= aux_tmp.c0;
    aux.c1 -= aux_tmp.c1;
    aux.c2 -= aux_tmp.c2;

    return aux;

}

void acc_Deo( __restrict su3_soa * const u, 
        __restrict vec3_soa * const out, 
        __restrict vec3_soa * const in,double_soa * backfield) {
    SETINUSE(out);
    int hd0, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent gang(nd3)
    for(d3=0; d3<nd3; d3++) {
#pragma acc loop independent gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector(DIM_BLOCK_X) 
                for(hd0=0; hd0 < nd0h; hd0++) {

                    int d0, d0m, d1m, d2m, d3m, d0p, d1p, d2p, d3p, idxh, matdir,dirindex;
                    vec3 aux=(vec3) {0,0,0};

                    // (d0+d1+d2+d3) even
                    d0 = 2*hd0 + ((d1+d2+d3) & 0x1);

                    d0m = d0 - 1;
                    d0m = d0m + (((d0m >> 31) & 0x1) * nd0);
                    d1m = d1 - 1;
                    d1m = d1m + (((d1m >> 31) & 0x1) * nd1);
                    d2m = d2 - 1;
                    d2m = d2m + (((d2m >> 31) & 0x1) * nd2);
                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d0p = d0 + 1;
                    d0p *= (((d0p-nd0) >> 31) & 0x1);
                    d1p = d1 + 1;
                    d1p *= (((d1p-nd1) >> 31) & 0x1);
                    d2p = d2 + 1;
                    d2p *= (((d2p-nd2) >> 31) & 0x1);
                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

#define SUB_RESULT_DEO(matdir, indexm) \
                    aux = subResult(aux,conjmat_vec_mul_arg( &u[matdir],indexm,\
                                in,indexm,&backfield[matdir]));  

                    SUB_RESULT_DEO(1,snum_acc(d0m,d1,d2,d3));
                    SUB_RESULT_DEO(3,snum_acc(d0,d1m,d2,d3));
                    SUB_RESULT_DEO(5,snum_acc(d0,d1,d2m,d3));
                    SUB_RESULT_DEO(7,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DEO

                    //      for(dirindex=0;dirindex<4;dirindex++) {   
                    //          matdir = 2*dirindex+1;
                    //          aux = subResult(aux,
                    //             conjmat_vec_mul_arg( &u[matdir],idxhm[dirindex],
                    //                   in,idxhm[dirindex],&backfield[matdir]));
                    //      }                                


                    //////////////////////////////////////////////////////////////
                    idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DEO(matdir, indexp) \
                    aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,\
                                in,indexp,&backfield[matdir])); 

                    SUM_RESULT_DEO(0,snum_acc(d0p,d1,d2,d3));
                    SUM_RESULT_DEO(2,snum_acc(d0,d1p,d2,d3));
                    SUM_RESULT_DEO(4,snum_acc(d0,d1,d2p,d3));
                    SUM_RESULT_DEO(6,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DEO

                    //      for(dirindex=0;dirindex<4;dirindex++) {   
                    //	      matdir = 2*dirindex;
                    //         aux   = sumResult(aux,
                    //            mat_vec_mul_arg(&u[matdir],idxh,
                    //                in,idxhp[dirindex],&backfield[matdir])); 
                    //      }	  

                    ///////////////////////////////////////////////////////////// 

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1   
        } // Loop over nd2    
    } // Loop over nd3      
}
void acc_Doe( __restrict su3_soa * const u,
        __restrict vec3_soa * const out,
        __restrict vec3_soa * const in,double_soa * backfield) {

    SETINUSE(out);
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent gang(nd3)
    for(d3=0; d3<nd3; d3++) {
#pragma acc loop independent gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d0m, d1m, d2m, d3m, d0p, d1p, d2p, d3p, idxh, matdir,dirindex;
                    vec3 aux=(vec3) {0,0,0};

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    d0m = d0 - 1;
                    d0m = d0m + (((d0m >> 31) & 0x1) * nd0);
                    d1m = d1 - 1;
                    d1m = d1m + (((d1m >> 31) & 0x1) * nd1);
                    d2m = d2 - 1;
                    d2m = d2m + (((d2m >> 31) & 0x1) * nd2);
                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d0p = d0 + 1;
                    d0p *= (((d0p-nd0) >> 31) & 0x1);
                    d1p = d1 + 1;
                    d1p *= (((d1p-nd1) >> 31) & 0x1);
                    d2p = d2 + 1;
                    d2p *= (((d2p-nd2) >> 31) & 0x1);
                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

#define SUB_RESULT_DOE(matdir,index)\
                    aux = subResult(aux, conjmat_vec_mul_arg( &u[matdir],index,\
                                in,index,&backfield[matdir]));

                    SUB_RESULT_DOE(0,snum_acc(d0m,d1,d2,d3));
                    SUB_RESULT_DOE(2,snum_acc(d0,d1m,d2,d3));
                    SUB_RESULT_DOE(4,snum_acc(d0,d1,d2m,d3));
                    SUB_RESULT_DOE(6,snum_acc(d0,d1,d2,d3m));
#undef SUB_RESULT_DOE        



                    //        for(dirindex=0;dirindex<4;dirindex++) {   
                    //            matdir = 2*dirindex;
                    //            aux = subResult(aux,
                    //                 conjmat_vec_mul_arg( &u[matdir],idxhm[dirindex],
                    //                        in,idxhm[dirindex],&backfield[matdir]));
                    //        }                                

                    ///////////////////////////////////////////////////////////////////

                    idxh = snum_acc(d0,d1,d2,d3);

#define SUM_RESULT_DOE(matdir,index)\
                    aux   = sumResult(aux, mat_vec_mul_arg(&u[matdir],idxh,\
                                in,index,&backfield[matdir]));

                    SUM_RESULT_DOE(1,snum_acc(d0p,d1,d2,d3));
                    SUM_RESULT_DOE(3,snum_acc(d0,d1p,d2,d3));
                    SUM_RESULT_DOE(5,snum_acc(d0,d1,d2p,d3));
                    SUM_RESULT_DOE(7,snum_acc(d0,d1,d2,d3p));
#undef SUM_RESULT_DOE


                    //        for(dirindex=0;dirindex<4;dirindex++) {
                    //            matdir = 2*dirindex+1;
                    //            aux   = sumResult(aux,
                    //                    mat_vec_mul_arg(&u[matdir],idxh,
                    //                        in,idxhp[dirindex],&backfield[matdir]));
                    //        }	  


                    /////////////////////////////////////////////////////////////////

                    out->c0[idxh] = aux.c0*0.5;
                    out->c1[idxh] = aux.c1*0.5;
                    out->c2[idxh] = aux.c2*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
}


void dM_dmu_eo( __restrict su3_soa * const u,
        vec3_soa * const out,
        vec3_soa * const in,
        double_soa * backfield)
{
    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent gang(nd3)
    for(d3=0; d3<nd3; d3++) {
#pragma acc loop independent gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
                for(d0h=0; d0h < nd0h; d0h++) {


                    int d0, d3m, d3p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) even
                    d0 = 2*d0h + ((d1+d2+d3) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

                    //direzione tempo positiva,
                    aux = mat_vec_mul_arg( &u[geom_par.tmap*2], idxh, in, snum_acc(d0,d1,d2,d3p),
                            &backfield[geom_par.tmap*2] );

                    // direzione tempo negativa
                    aux = sumResult(aux, conjmat_vec_mul_arg( &u[geom_par.tmap*2+1],
                                snum_acc(d0,d1,d2,d3m),in,snum_acc(d0,d1,d2,d3m),&backfield[geom_par.tmap*2+1]));

                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
}

void dM_dmu_oe( __restrict su3_soa * const u,
        vec3_soa * const out,
        vec3_soa * const in,
        double_soa * backfield)
{

    int d0h, d1, d2, d3;

#pragma acc kernels present(u) present(out) present(in) present(backfield)
#pragma acc loop independent gang(nd3)
    for(d3=0; d3<nd3; d3++) {
#pragma acc loop independent gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
                for(d0h=0; d0h < nd0h; d0h++) {

                    int d0, d3m, d3p, idxh;
                    vec3 aux;

                    // (d0+d1+d2+d3) odd
                    d0 = 2*d0h + ((d1+d2+d3+1) & 0x1);

                    idxh = snum_acc(d0,d1,d2,d3);

                    d3m = d3 - 1;
                    d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

                    d3p = d3 + 1;
                    d3p *= (((d3p-nd3) >> 31) & 0x1);

                    //direzione tempo negativa
                    aux = conjmat_vec_mul_arg( &u[geom_par.tmap*2],snum_acc(d0,d1,d2,d3m), in,
                            snum_acc(d0,d1,d2,d3m), &backfield[geom_par.tmap*2]);

                    // direzione tempo positiva
                    aux = sumResult(aux, mat_vec_mul_arg( &u[geom_par.tmap*2+1],idxh,
                                in,snum_acc(d0,d1,d2,d3p),&backfield[geom_par.tmap*2+1]));

                    /////////////////////////////////////////////////////////////////////     

                    out->c0[idxh] = (aux.c0)*0.5;
                    out->c1[idxh] = (aux.c1)*0.5;
                    out->c2[idxh] = (aux.c2)*0.5;

                } // Loop over nd0h
            } // Loop over nd1
        } // Loop over nd2
    } // Loop over nd3
}


inline void fermion_matrix_multiplication( __restrict su3_soa * const u, 
        __restrict vec3_soa * const out,  __restrict vec3_soa * const in, 
        __restrict vec3_soa * const temp1, ferm_param *pars){
    SETREQUESTED(temp1);
    SETREQUESTED(out);
    acc_Doe(u,temp1,in,pars->phases);
    acc_Deo(u,out,temp1,pars->phases);
    combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass,out);// Nuova funzione in OpenAcc/fermionic_utilities.c
    SETFREE(temp1);
}
inline void fermion_matrix_multiplication_shifted( __restrict su3_soa * const u, 
        __restrict vec3_soa * const out,  __restrict vec3_soa * const in, 
        __restrict vec3_soa * const temp1, ferm_param *pars, double shift){
    SETREQUESTED(temp1);
    SETREQUESTED(out);
    acc_Doe(u,temp1,in,pars->phases);
    acc_Deo(u,out,temp1,pars->phases);
    combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass+shift,out);// Nuova funzione in OpenAcc/fermionic_utilities.c
    SETFREE(temp1);


}









#endif

