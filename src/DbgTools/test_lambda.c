/***********************
 * Test per la routine 
 *
 * (DEOTT)_compute_loc_Lambda()
 *
 ************************/




#include "openacc.h"
#include <complex.h>
#include "../OpenAcc/struct_c_def.c"
#include "../OpenAcc/dbgtools.c"
#include "../OpenAcc/single_types.c"
#include "../OpenAcc/cayley_hamilton.c"
#include "../OpenAcc/stouting_deottimizzato.c"
#include "../OpenAcc/stouting.c"



int main(){

    
    // Lambda 
    single_su3 gl3_temp0;
    single_su3 sSP;
    single_tamat sQA;
    su3_soa * U =   (su3_soa * ) malloc(sizeof(su3_soa)); 
    su3_soa * SP =  (su3_soa * ) malloc(sizeof(su3_soa));
    su3_soa * tmp = (su3_soa * ) malloc(sizeof(su3_soa));
    thmat_soa * L =  (thmat_soa * ) malloc(sizeof(thmat_soa));
    tamat_soa * QA = (tamat_soa * ) malloc(sizeof(tamat_soa));
    int idx;


    // FOR Lambda testing
    idx = 0;

    sQA.rc00 = -0.5;
    sQA.rc11 = 0.5+1.2*I;
    sQA.c01 = -1.0*I;
    sQA.c02 = 1.2+2*I;
    sQA.c12 = 2.3*I;

    CH_exponential_antihermitian_nissalike(&gl3_temp0,&sQA);
    // sSP for lambda
    sSP.comp[0][0] =  0.4-4.4*I;
    sSP.comp[0][1] =  1.5+5.1*I;
    sSP.comp[0][2] = -0.2+0.5*I;
    sSP.comp[1][0] =  0.7-1.7*I;
    sSP.comp[1][1] = -1.8+0.8*I;
    sSP.comp[1][2] =  1.9+1.5*I;
    sSP.comp[2][0] =  0.6-0.4*I;
    sSP.comp[2][1] = -1.5+1.1*I;
    sSP.comp[2][2] =  2.4+5.2*I;
    single_gl3_into_su3_soa(SP,idx,&sSP);
    // U for lambda
    single_su3_into_su3_soa(U,idx,&gl3_temp0);
    // QA for lambda
    single_tamat_into_tamat_soa(QA,idx,&sQA);
    // temp should not necessitate any initialization, but...
    // ad cazzum initialization
    single_gl3_into_su3_soa(tmp,idx,&sSP);


//#define DEOTT

#ifdef DEOTT

    DEOTT_compute_loc_Lambda(
#else
         compute_loc_Lambda(
#endif
            L, 
            SP,
            U,
            QA,
            tmp,
            idx );

 print_1thmat_soa(L,"Lambda.thmat");

return 0;

}


