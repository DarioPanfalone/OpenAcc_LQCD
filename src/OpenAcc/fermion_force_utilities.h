#ifndef FERMION_FORCE_UTILITIES_H
#define FERMION_FORCE_UTILITIES_H

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif


#pragma acc routine seq
static inline void vec1_directprod_conj_vec2_into_mat1( __restrict su3_soa * const aux_u,
							int idxh,
							__restrict vec3_soa  * const fer_l, // questo fermione e' costante e non viene modificato qui dentro
							int idl,
							__restrict vec3_soa  * const fer_r, // questo fermione e' costante e non viene modificato qui dentro
							int idr,
							double factor){
  // forse si possono risparmiare un po di conti andando a prendere le sole parti reale e immaginarie 
  //  e scrivendo i prodotti in modo molto piu' sbrodolato
  // in particolare per gli elementi lungo la diagonale
  d_complex r0 = factor * conj(fer_r->c0[idr]);
  d_complex r1 = factor * conj(fer_r->c1[idr]);
  d_complex r2 = factor * conj(fer_r->c2[idr]);
  d_complex l0 = fer_l->c0[idl];
  d_complex l1 = fer_l->c1[idl];
  d_complex l2 = fer_l->c2[idl];

  aux_u->r0.c0[idxh] += l0*r0;
  aux_u->r0.c1[idxh] += l0*r1;
  aux_u->r0.c2[idxh] += l0*r2;
  aux_u->r1.c0[idxh] += l1*r0;
  aux_u->r1.c1[idxh] += l1*r1;
  aux_u->r1.c2[idxh] += l1*r2;
  aux_u->r2.c0[idxh] += l2*r0;
  aux_u->r2.c1[idxh] += l2*r1;
  aux_u->r2.c2[idxh] += l2*r2;

}


#pragma acc routine seq
static inline  void mat1_times_auxmat_into_tamat(  __restrict su3_soa * const mat1, // e' costante e non viene modificato
						   const  int idx,
						   const  int eta,
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   const  int idx_aux,
						   __restrict tamat_soa * const ipdot,
						   const  int idipdot
#if defined(BACKFIELD) || defined(IMCHEMPOT)
						   ,d_complex phase
#endif
                           ){
  d_complex mat1_00 = mat1->r0.c0[idx];
  d_complex mat1_01 = mat1->r0.c1[idx];
  d_complex mat1_02 = mat1->r0.c2[idx];
  d_complex mat1_10 = mat1->r1.c0[idx];
  d_complex mat1_11 = mat1->r1.c1[idx];
  d_complex mat1_12 = mat1->r1.c2[idx];
  //Compute 3rd matrix row from the first two
  d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
  //Multiply 3rd row by eta
  mat1_20 *= eta;
  mat1_21 *= eta;
  mat1_22 *= eta;

#if defined(BACKFIELD) || defined(IMCHEMPOT)
  mat1_00 *= phase;
  mat1_01 *= phase;
  mat1_02 *= phase;
  mat1_10 *= phase;
  mat1_11 *= phase;
  mat1_12 *= phase;
  mat1_20 *= phase;
  mat1_21 *= phase;
  mat1_22 *= phase;
#endif

  d_complex auxmat_00 = auxmat->r0.c0[idx_aux];
  d_complex auxmat_01 = auxmat->r0.c1[idx_aux];
  d_complex auxmat_02 = auxmat->r0.c2[idx_aux];
  d_complex auxmat_10 = auxmat->r1.c0[idx_aux];
  d_complex auxmat_11 = auxmat->r1.c1[idx_aux];
  d_complex auxmat_12 = auxmat->r1.c2[idx_aux];
  d_complex auxmat_20 = auxmat->r2.c0[idx_aux];
  d_complex auxmat_21 = auxmat->r2.c1[idx_aux];
  d_complex auxmat_22 = auxmat->r2.c2[idx_aux];

  // product;
  d_complex a00 = mat1_00 * auxmat_00 + mat1_01 * auxmat_10 + mat1_02 * auxmat_20;
  d_complex a01 = mat1_00 * auxmat_01 + mat1_01 * auxmat_11 + mat1_02 * auxmat_21;
  d_complex a02 = mat1_00 * auxmat_02 + mat1_01 * auxmat_12 + mat1_02 * auxmat_22;

  mat1_00 = mat1_10 * auxmat_00 + mat1_11 * auxmat_10 + mat1_12 * auxmat_20;
  mat1_01 = mat1_10 * auxmat_01 + mat1_11 * auxmat_11 + mat1_12 * auxmat_21;
  mat1_02 = mat1_10 * auxmat_02 + mat1_11 * auxmat_12 + mat1_12 * auxmat_22;

  mat1_10 = mat1_20 * auxmat_00 + mat1_21 * auxmat_10 + mat1_22 * auxmat_20;
  mat1_11 = mat1_20 * auxmat_01 + mat1_21 * auxmat_11 + mat1_22 * auxmat_21; 
  mat1_12 = mat1_20 * auxmat_02 + mat1_21 * auxmat_12 + mat1_22 * auxmat_22;

  ipdot->c01[idipdot]  -= 0.5*((a01) - conj(mat1_00));
  ipdot->c02[idipdot]  -= 0.5*((a02) - conj(mat1_10));
  ipdot->c12[idipdot]  -= 0.5*((mat1_02) - conj(mat1_11));
  ipdot->rc00[idipdot] -= cimag(a00)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));
  ipdot->rc11[idipdot] -= cimag(mat1_01)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));

}


#ifdef STOUT_FERMIONS
#pragma acc routine seq
static inline  void mat1_times_auxmat_into_tamat_nophase(  __restrict su3_soa * const mat1, // e' costante e non viene modificato
						   const  int idx,
						   const  int eta,
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   const  int idx_aux,
						   __restrict tamat_soa * const ipdot,
						   const  int idipdot
                           ){
  d_complex mat1_00 = mat1->r0.c0[idx];
  d_complex mat1_01 = mat1->r0.c1[idx];
  d_complex mat1_02 = mat1->r0.c2[idx];
  d_complex mat1_10 = mat1->r1.c0[idx];
  d_complex mat1_11 = mat1->r1.c1[idx];
  d_complex mat1_12 = mat1->r1.c2[idx];
  //Compute 3rd matrix row from the first two
  d_complex mat1_20 = conj( ( mat1_01 * mat1_12 ) - ( mat1_02 * mat1_11) ) ;
  d_complex mat1_21 = conj( ( mat1_02 * mat1_10 ) - ( mat1_00 * mat1_12) ) ;
  d_complex mat1_22 = conj( ( mat1_00 * mat1_11 ) - ( mat1_01 * mat1_10) ) ;
  //Multiply 3rd row by eta
  mat1_20 *= eta;
  mat1_21 *= eta;
  mat1_22 *= eta;

  d_complex auxmat_00 = auxmat->r0.c0[idx_aux];
  d_complex auxmat_01 = auxmat->r0.c1[idx_aux];
  d_complex auxmat_02 = auxmat->r0.c2[idx_aux];
  d_complex auxmat_10 = auxmat->r1.c0[idx_aux];
  d_complex auxmat_11 = auxmat->r1.c1[idx_aux];
  d_complex auxmat_12 = auxmat->r1.c2[idx_aux];
  d_complex auxmat_20 = auxmat->r2.c0[idx_aux];
  d_complex auxmat_21 = auxmat->r2.c1[idx_aux];
  d_complex auxmat_22 = auxmat->r2.c2[idx_aux];

  // product;
  d_complex a00 = mat1_00 * auxmat_00 + mat1_01 * auxmat_10 + mat1_02 * auxmat_20;
  d_complex a01 = mat1_00 * auxmat_01 + mat1_01 * auxmat_11 + mat1_02 * auxmat_21;
  d_complex a02 = mat1_00 * auxmat_02 + mat1_01 * auxmat_12 + mat1_02 * auxmat_22;

  mat1_00 = mat1_10 * auxmat_00 + mat1_11 * auxmat_10 + mat1_12 * auxmat_20;
  mat1_01 = mat1_10 * auxmat_01 + mat1_11 * auxmat_11 + mat1_12 * auxmat_21;
  mat1_02 = mat1_10 * auxmat_02 + mat1_11 * auxmat_12 + mat1_12 * auxmat_22;

  mat1_10 = mat1_20 * auxmat_00 + mat1_21 * auxmat_10 + mat1_22 * auxmat_20;
  mat1_11 = mat1_20 * auxmat_01 + mat1_21 * auxmat_11 + mat1_22 * auxmat_21; 
  mat1_12 = mat1_20 * auxmat_02 + mat1_21 * auxmat_12 + mat1_22 * auxmat_22;

  ipdot->c01[idipdot]  -= 0.5*((a01) - conj(mat1_00));
  ipdot->c02[idipdot]  -= 0.5*((a02) - conj(mat1_10));
  ipdot->c12[idipdot]  -= 0.5*((mat1_02) - conj(mat1_11));
  ipdot->rc00[idipdot] -= cimag(a00)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));
  ipdot->rc11[idipdot] -= cimag(mat1_01)-ONE_BY_THREE*(cimag(a00)+cimag(mat1_01)+cimag(mat1_12));

}
#endif // ifdef STOUT_FERMIONS


#if defined(BACKFIELD) || defined(IMCHEMPOT)
#pragma acc routine seq
static inline  void phase_times_auxmat_into_auxmat(
						   __restrict su3_soa * const auxmat,  // e' costante e non viene modificato
						   __restrict su3_soa * const pseudo_ipdot,
						   const  int idx,
                           d_complex phase
                           ){

    pseudo_ipdot->r0.c0[idx] += phase * auxmat->r0.c0[idx];
    pseudo_ipdot->r0.c1[idx] += phase * auxmat->r0.c1[idx];
    pseudo_ipdot->r0.c2[idx] += phase * auxmat->r0.c2[idx];
    pseudo_ipdot->r1.c0[idx] += phase * auxmat->r1.c0[idx];
    pseudo_ipdot->r1.c1[idx] += phase * auxmat->r1.c1[idx];
    pseudo_ipdot->r1.c2[idx] += phase * auxmat->r1.c2[idx];
    pseudo_ipdot->r2.c0[idx] += phase * auxmat->r2.c0[idx];
    pseudo_ipdot->r2.c1[idx] += phase * auxmat->r2.c1[idx];
    pseudo_ipdot->r2.c2[idx] += phase * auxmat->r2.c2[idx];


}

#else  //if defined(BACKFIELD) || defined(IMCHEMPOT)
//#pragma acc routine seq
static inline void accumulate_auxmat1_into_auxmat2(
						   __restrict su3_soa * const auxmat1,  // e' costante e non viene modificato
						   __restrict su3_soa * const auxmat2,
						   const  int idx
                           ){

    auxmat2->r0.c0[idx] += auxmat1->r0.c0[idx];
    auxmat2->r0.c1[idx] += auxmat1->r0.c1[idx];
    auxmat2->r0.c2[idx] += auxmat1->r0.c2[idx];
    auxmat2->r1.c0[idx] += auxmat1->r1.c0[idx];
    auxmat2->r1.c1[idx] += auxmat1->r1.c1[idx];
    auxmat2->r1.c2[idx] += auxmat1->r1.c2[idx];
    auxmat2->r2.c0[idx] += auxmat1->r2.c0[idx];
    auxmat2->r2.c1[idx] += auxmat1->r2.c1[idx];
    auxmat2->r2.c2[idx] += auxmat1->r2.c2[idx];

}

#endif  //if defined(BACKFIELD) || defined(IMCHEMPOT)






#pragma acc routine seq
static inline void assign_zero_to_tamat_soa_component(__restrict tamat_soa * const matrix_comp,
						      int idx){
  matrix_comp->c01[idx]=0.0+I*0.0;
  matrix_comp->c02[idx]=0.0+I*0.0;
  matrix_comp->c12[idx]=0.0+I*0.0;
  matrix_comp->rc00[idx]=0.0;
  matrix_comp->rc11[idx]=0.0;
}

void set_tamat_soa_to_zero( __restrict tamat_soa * const matrix);



void direct_product_of_fermions_into_auxmat(__restrict vec3_soa  * const loc_s, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict vec3_soa  * const loc_h, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict su3_soa * const aux_u,
					    const RationalApprox * const approx,
					    int iter);



void multiply_conf_times_force_and_take_ta_even(__restrict su3_soa * const u, // la conf e' costante e non viene modificata
						__restrict ferm_param * const tpars,
						__restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
						__restrict tamat_soa * const ipdot);



void multiply_conf_times_force_and_take_ta_odd(  __restrict su3_soa * const u, // e' costante e non viene modificata
						 __restrict ferm_param * const tpars,
					         __restrict su3_soa * const auxmat, // e' costante e non viene modificata
					         __restrict tamat_soa * const ipdot);


#ifdef STOUT_FERMIONS
void multiply_conf_times_force_and_take_ta_even_nophase(__restrict su3_soa * const u, // la conf e' costante e non viene modificata
							__restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante; non viene modificata
							__restrict tamat_soa * const ipdot);


void multiply_conf_times_force_and_take_ta_odd_nophase(  __restrict su3_soa * const u, // e' costante e non viene modificata
							 __restrict su3_soa * const auxmat, // e' costante e non viene modificata
							 __restrict tamat_soa * const ipdot);

#endif //ifdef STOUT_FERMIONS


#if defined(IMCHEMPOT) || defined(BACKFIELD)
void multiply_backfield_times_force(__restrict ferm_param * const tpars,
        __restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
        __restrict su3_soa * const pseudo_ipdot)i;

#else 
void accumulate_gl3soa_into_gl3soa(
        __restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
        __restrict su3_soa * const pseudo_ipdot);

#endif

void ker_openacc_compute_fermion_force( __restrict su3_soa * const u, // e' costante e non viene mai modificato qui dentro
					__restrict su3_soa * const aux_u,
					__restrict vec3_soa * const in_shiftmulti,  // e' costante e non viene mai modificato qui dentro
					__restrict vec3_soa  * const loc_s,
					__restrict vec3_soa  * const loc_h,
					__restrict ferm_param  *  tpars
					);

#endif
