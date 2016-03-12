#ifndef PLAQUETTES_H_
#define PLAQUETTES_H_


#include "../OpenAcc/struct_c_def.h"


// routine for the computation of the average of the plaquettes computed on the plane mu-nu
// 1) all the plaquettes on the plane mu-nu are computed and saved locally
// 2) finally the reduction of the traces is performed
// routine to compute the staples for each site on a given plane mu-nu and sum the result to the local stored staples
void calc_loc_staples_removing_stag_phases_nnptrick( 
        __restrict su3_soa * const u, __restrict su3_soa * const loc_stap,
        const int mu,const int nu);
void calc_loc_staples_removing_stag_phases_nnptrick_all(  
        __restrict su3_soa * const u, __restrict su3_soa * const loc_stap );
void calc_loc_staples_removing_stag_phases_nnptrick_all_only_even( 
        __restrict su3_soa * const u,  __restrict su3_soa * const loc_stap );
void calc_loc_staples_removing_stag_phases_nnptrick_all_only_odd(
        __restrict su3_soa * const u,__restrict su3_soa * const loc_stap); 
double calc_loc_plaquettes_removing_stag_phases_nnptrick(
        __restrict su3_soa * const u,
        __restrict su3_soa * const loc_plaq,
        dcomplex_soa * const tr_local_plaqs,
        const int mu, const int nu);




#endif
