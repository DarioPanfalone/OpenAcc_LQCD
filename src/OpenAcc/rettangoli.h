#ifndef RETTANGOLI_H
#define RETTANGOLI_H

#include "./struct_c_def.h"
#include "./rettangoli.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif



double calc_loc_rectangles_2x1_removing_stag_phases_nnptrick( 
        __restrict su3_soa * const u,__restrict su3_soa * const loc_plaq, 
        dcomplex_soa * const tr_local_plaqs,const int mu,	const int nu);



double calc_loc_rectangles_2x1_removing_stag_phases_nnptrick(  
        __restrict su3_soa * const u,__restrict su3_soa * const loc_plaq,
        dcomplex_soa * const tr_local_plaqs,const int mu,const int nu);



#pragma acc routine seq
void    PPMMM_5mat_prod_addto_mat6_absent_stag_phases(  
        __restrict su3_soa * const mat1, const int idx_mat1,
        __restrict su3_soa * const mat2, const int idx_mat2,
        __restrict su3_soa * const mat3, const int idx_mat3,
        __restrict su3_soa * const mat4, const int idx_mat4,
        __restrict su3_soa * const mat5, const int idx_mat5,
        __restrict su3_soa * const mat6, const int idx_mat6);

#pragma acc routine seq
void    PMMMP_5mat_prod_addto_mat6_absent_stag_phases( 
        __restrict su3_soa * const mat1, const int idx_mat1,
        __restrict su3_soa * const mat2, const int idx_mat2,
        __restrict su3_soa * const mat3, const int idx_mat3,
        __restrict su3_soa * const mat4, const int idx_mat4,
        __restrict su3_soa * const mat5, const int idx_mat5,
        __restrict su3_soa * const mat6, const int idx_mat6);

void    MMMPP_5mat_prod_addto_mat6_absent_stag_phases( 
        __restrict su3_soa * const mat1, const int idx_mat1,
        __restrict su3_soa * const mat2, const int idx_mat2,
        __restrict su3_soa * const mat3, const int idx_mat3,
        __restrict su3_soa * const mat4, const int idx_mat4,
        __restrict su3_soa * const mat5, const int idx_mat5,
        __restrict su3_soa * const mat6, const int idx_mat6);

void calc_loc_improved_staples_typeA_removing_stag_phases_nnptrick_all( 
        __restrict su3_soa * const u, __restrict su3_soa * const loc_stap );

void calc_loc_improved_staples_typeB_removing_stag_phases_nnptrick_all( 
        __restrict su3_soa * const u, __restrict su3_soa * const loc_stap );

void calc_loc_improved_staples_typeC_removing_stag_phases_nnptrick_all( 
        __restrict su3_soa * const u, __restrict su3_soa * const loc_stap );

void calc_loc_improved_staples_typeABC_removing_stag_phases_nnptrick_all( 
        __restrict su3_soa * const u,__restrict su3_soa * const loc_stap );

double  calc_rettangolo_soloopenacc( __restrict  su3_soa * const tconf_acc,
        __restrict su3_soa * const local_plaqs, 
        dcomplex_soa * const tr_local_plaqs);

#endif
