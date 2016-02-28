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
#include "./matvecmul.h"


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

