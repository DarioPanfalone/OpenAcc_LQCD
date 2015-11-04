/*****************
Test per le routines

 void DEOTT_RIGHT_iABC_times_DminusE_absent_stag_phases()
 void DEOTT_RIGHT_iFABC_absent_stag_phases()
 void DEOTT_RIGHT_miABGC_absent_stag_phases()
 void DEOTT_LEFT_iAB_times_GminusE_times_C_absent_stag_phases()
 void DEOTT_LEFT_iABCD_absent_stag_phases()
 void DEOTT_LEFT_miAFBC_absent_stag_phases()
***********************/

double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);


#include "openacc.h"
#include "../OpenAcc/struct_c_def.c"
#include "../OpenAcc/dbgtools.c"
#include "../OpenAcc/single_types.c"
#include "../OpenAcc/cayley_hamilton.c"
#include "../OpenAcc/stouting_deottimizzato.c"



int main(){

    su3_soa * gauge_matrix = (su3_soa * ) malloc(sizeof(su3_soa));
    su3_soa * result = (su3_soa * ) malloc(sizeof(su3_soa));
    single_su3 gl3_temp0;
    //just one array
    thmat_soa * thmats = (thmat_soa * ) malloc(sizeof(thmat_soa));
    single_thmat thmat_temp0;


    int idxA,idxB,idxC,idxLD,idxLE,idxRES;

    // Insertion of the right matrices in the right places in 
    // gauge_matrix and thmats.

    //A
    idxA = 11 ; 
    gl3_temp0.comp[0][0] = 1 ; 
    gl3_temp0.comp[0][1] = 0 ; 
    gl3_temp0.comp[0][2] = 0 ; 
    gl3_temp0.comp[1][0] = 0 ; 
    gl3_temp0.comp[1][1] = 1 ; 
    gl3_temp0.comp[1][2] = 0 ; 
    single_su3_into_su3_soa(gauge_matrix,idxA,&gl3_temp0);
    //B
    idxB = 15 ; 
    gl3_temp0.comp[0][0] = 1 ; 
    gl3_temp0.comp[0][1] = 0 ; 
    gl3_temp0.comp[0][2] = 0 ; 
    gl3_temp0.comp[1][0] = 0 ; 
    gl3_temp0.comp[1][1] = 1 ; 
    gl3_temp0.comp[1][2] = 0 ; 
    single_su3_into_su3_soa(gauge_matrix,idxB,&gl3_temp0);
    //C
    idxC = 100 ; 
    gl3_temp0.comp[0][0] = 1 ; 
    gl3_temp0.comp[0][1] = 0 ; 
    gl3_temp0.comp[0][2] = 0 ; 
    gl3_temp0.comp[1][0] = 0 ; 
    gl3_temp0.comp[1][1] = 1 ; 
    gl3_temp0.comp[1][2] = 0 ; 
    single_su3_into_su3_soa(gauge_matrix,idxC,&gl3_temp0);
    //LD
    idxLD = 100 ; 
    thmat_temp0.rc00 = 1 ; 
    thmat_temp0.rc11 = 0 ; 
    thmat_temp0.c01 = 0 ; 
    thmat_temp0.c02 = 0 ; 
    thmat_temp0.c12 = 1 ; 
    single_thmat_into_thmat_soa(thmats,idxLD,&thmat_temp0);
    //LE
    idxLE = 23 ; 
    thmat_temp0.rc00 = 1 ; 
    thmat_temp0.rc11 = 0 ; 
    thmat_temp0.c01 = 0 ; 
    thmat_temp0.c02 = 0 ; 
    thmat_temp0.c12 = 1 ; 
    single_thmat_into_thmat_soa(thmats,idxLE,&thmat_temp0);

    //res
    idxRES = 23 ; 
    gl3_temp0.comp[0][0] = 0 ; 
    gl3_temp0.comp[0][1] = 0 ; 
    gl3_temp0.comp[0][2] = 0 ; 
    gl3_temp0.comp[1][0] = 0 ; 
    gl3_temp0.comp[1][1] = 0 ; 
    gl3_temp0.comp[1][2] = 0 ; 
    gl3_temp0.comp[2][0] = 0 ; 
    gl3_temp0.comp[2][1] = 0 ; 
    gl3_temp0.comp[2][2] = 0 ; 
    single_gl3_into_su3_soa(gauge_matrix,idxRES,&gl3_temp0);








//DEOTT_RIGHT_iABC_times_DminusE_absent_stag_phases(//UA*dag(UB)*dag(UC)*RHO*I*(LD - LE))
DEOTT_RIGHT_iFABC_absent_stag_phases(//((RHO*I)*LF)*UA*dag(UB)*dag(UC)
//DEOTT_RIGHT_miABGC_absent_stag_phases(//-UA*dag(UB)*RHO*I*LG*dag(UC)
//DEOTT_LEFT_iAB_times_GminusE_times_C_absent_stag_phases(//dag(UA)*dag(UB)*RHO*I*(LG-LE)*UC
//DEOTT_LEFT_iABCD_absent_stag_phases(//dag(UA)*dag(UB)*UC*RHO*I*LD
//DEOTT_LEFT_miAFBC_absent_stag_phases(//-dag(UA)*RHO*I*LF*dag(UB)*UC
            gauge_matrix,
            idxA,
            gauge_matrix,
            idxB,
            gauge_matrix,
            idxC,
            thmats,
            idxLD,
//          thmats,
//          idxLE,
            result,
            idxRES);

 print_1su3_soa(gauge_matrix,"GaugeMatric.su3");
 print_1su3_soa(result,"result.gl3");
 print_1thmat_soa(thmats,"thmat");



return 0;

}


