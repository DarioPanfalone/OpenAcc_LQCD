/*****************
Test per le routines

 void (DEOTT)_RIGHT_iABC_times_DminusE_absent_stag_phases()
 void (DEOTT)_RIGHT_iFABC_absent_stag_phases()
 void (DEOTT)_RIGHT_miABGC_absent_stag_phases()
 void (DEOTT)_LEFT_iAB_times_GminusE_times_C_absent_stag_phases()
 void (DEOTT)_LEFT_iABCD_absent_stag_phases()
 void (DEOTT)_LEFT_miAFBC_absent_stag_phases()
***********************/

double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);


#include "openacc.h"
#include <complex.h>
#include "../OpenAcc/struct_c_def.c"
#include "../OpenAcc/dbgtools.c"
#include "../OpenAcc/single_types.c"
#include "../OpenAcc/cayley_hamilton.c"
#include "../OpenAcc/stouting_deottimizzato.c"
#include "../OpenAcc/stouting.c"



int main(){

    su3_soa * gauge_matrix = (su3_soa * ) malloc(sizeof(su3_soa));
    su3_soa * result = (su3_soa * ) malloc(sizeof(su3_soa));
    single_su3 gl3_temp0;
    //just one array
    thmat_soa * thmats = (thmat_soa * ) malloc(sizeof(thmat_soa));
    singel_tamat QA;
    single_thmat thmat_temp0;

    const double a = 1/sqrt(2.0);


    int idxA,idxB,idxC,idxLD,idxLE,idxRES;

    // Insertion of the right matrices in the right places in 
    // gauge_matrix and thmats.


    //A
    idxA = 0 ; 
    QA.rc00 = 1; QA.rc11 = I; QA.rc01 = 3.2*I; QA.rc02 = 5; QA.rc12 = 1;
    CH_exponential_antihermitian_nissalike(&gl3_temp0,&QA);
    single_su3_into_su3_soa(gauge_matrix,idxA,&gl3_temp0);
    //B
    idxB = 1 ; 
    QA.rc00 = 0.1; QA.rc11 = 3+I; QA.rc01 = 2+I; QA.rc02 =0.1*I; QA.rc12 = -0.02;
    CH_exponential_antihermitian_nissalike(&gl3_temp0,&QA);
    single_su3_into_su3_soa(gauge_matrix,idxB,&gl3_temp0);
    //C
    idxC = 2 ; 
    QA.rc00 = -0.5; QA.rc11 = 0.5+1.2*I; QA.rc01 = -1.0*I; QA.rc02 = 1.2+2*I; QA.rc12 = 2.3*I;
    CH_exponential_antihermitian_nissalike(&gl3_temp0,&QA);
    single_su3_into_su3_soa(gauge_matrix,idxC,&gl3_temp0);
    //LD
    idxLD = 0 ; 
    thmat_temp0.rc00 = 1 ; 
    thmat_temp0.rc11 = 2.1 ; 
    thmat_temp0.c01 = 1.5*I ; 
    thmat_temp0.c02 = -1.2 ; 
    thmat_temp0.c12 = 0.1 ; 
    single_thmat_into_thmat_soa(thmats,idxLD,&thmat_temp0);
    //LE
    idxLE = 1 ; 
    thmat_temp0.rc00 = .02 ; 
    thmat_temp0.rc11 = 1.2 ; 
    thmat_temp0.c01 = -1.2*I ; 
    thmat_temp0.c02 = 0.8 ; 
    thmat_temp0.c12 = -I*1.5 ; 
    single_thmat_into_thmat_soa(thmats,idxLE,&thmat_temp0);

    //res
    idxRES = 0 ; 
    gl3_temp0.comp[0][0] = 0 ; 
    gl3_temp0.comp[0][1] = 0 ; 
    gl3_temp0.comp[0][2] = 0 ; 
    gl3_temp0.comp[1][0] = 0 ; 
    gl3_temp0.comp[1][1] = 0 ; 
    gl3_temp0.comp[1][2] = 0 ; 
    gl3_temp0.comp[2][0] = 0 ; 
    gl3_temp0.comp[2][1] = 0 ; 
    gl3_temp0.comp[2][2] = 0 ; 
    single_gl3_into_su3_soa(result,idxRES,&gl3_temp0);

//#define DEOTT
#define NFUNC 6

#ifdef DEOTT

  #if NFUNC == 1
  DEOTT_RIGHT_iABC_times_DminusE_absent_stag_phases(//UA*dag(UB)*dag(UC)*RHO*I*(LD - LE))
  #elif NFUNC == 2
  DEOTT_RIGHT_iFABC_absent_stag_phases(//((RHO*I)*LF)*UA*dag(UB)*dag(UC)
  #elif NFUNC == 3
  DEOTT_RIGHT_miABGC_absent_stag_phases(//-UA*dag(UB)*RHO*I*LG*dag(UC)
  #elif NFUNC == 4
  DEOTT_LEFT_iAB_times_GminusE_times_C_absent_stag_phases(//dag(UA)*dag(UB)*RHO*I*(LG-LE)*UC
  #elif NFUNC == 5
  DEOTT_LEFT_iABCD_absent_stag_phases(//dag(UA)*dag(UB)*UC*RHO*I*LD
  #elif NFUNC == 6
  DEOTT_LEFT_miAFBC_absent_stag_phases(//-dag(UA)*RHO*I*LF*dag(UB)*UC
  #endif

#else

  #if NFUNC == 1
  RIGHT_iABC_times_DminusE_absent_stag_phases(//UA*dag(UB)*dag(UC)*RHO*I*(LD - LE))
  #elif NFUNC == 2
  RIGHT_iFABC_absent_stag_phases(//((RHO*I)*LF)*UA*dag(UB)*dag(UC)
  #elif NFUNC == 3
  RIGHT_miABGC_absent_stag_phases(//-UA*dag(UB)*RHO*I*LG*dag(UC)
  #elif NFUNC == 4
  LEFT_iAB_times_GminusE_times_C_absent_stag_phases(//dag(UA)*dag(UB)*RHO*I*(LG-LE)*UC
  #elif NFUNC == 5
  LEFT_iABCD_absent_stag_phases(//dag(UA)*dag(UB)*UC*RHO*I*LD
  #elif NFUNC == 6
  LEFT_miAFBC_absent_stag_phases(//-dag(UA)*RHO*I*LF*dag(UB)*UC
  #endif
  
#endif
            gauge_matrix,
            idxA,
            gauge_matrix,
            idxB,
            gauge_matrix,
            idxC,
            thmats,
            idxLD,
#if (NFUNC == 1 ) || (NFUNC == 4)
            thmats,
            idxLE,
#endif
            result,
            idxRES);

 print_1su3_soa(gauge_matrix,"GaugeMatric.su3");
 print_1su3_soa(result,"result.gl3");
 print_1thmat_soa(thmats,"thmat");


 





return 0;

}


