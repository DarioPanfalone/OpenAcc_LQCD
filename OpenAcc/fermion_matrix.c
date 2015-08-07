
#ifndef FERMION_MATRIX_C
#define FERMION_MATRIX_C

#include "./struct_c_def.c"
#include "openacc.h"
#include "./fermionic_utilities.c"

static inline vec3 mat_vec_mul( __restrict su3_soa * const matrix,
                                const int idx_mat,
                                const int eta,
                                __restrict vec3_soa * const in_vect,
                                const int idx_vect) {

  vec3 out_vect;

  d_complex vec0 = in_vect->c0[idx_vect];
  d_complex vec1 = in_vect->c1[idx_vect];
  d_complex vec2 = in_vect->c2[idx_vect];

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
  //Multiply 3rd row by eta 

  mat20 = (mat20)*eta;
  mat21 = (mat21)*eta;
  mat22 = (mat22)*eta;

  out_vect.c0 = ( mat00 * vec0 ) + ( mat01 * vec1 ) + ( mat02 * vec2 );
  out_vect.c1 = ( mat10 * vec0 ) + ( mat11 * vec1 ) + ( mat12 * vec2 );
  out_vect.c2 = ( mat20 * vec0 ) + ( mat21 * vec1 ) + ( mat22 * vec2 );

  return out_vect;

}


static inline vec3 conjmat_vec_mul( __restrict su3_soa * const matrix,
                                    const int idx_mat,
                                    const int eta,
                                    __restrict vec3_soa * const in_vect,
                                    const int idx_vect) {

  vec3 out_vect;

  d_complex vec0 = in_vect->c0[idx_vect];
  d_complex vec1 = in_vect->c1[idx_vect];
  d_complex vec2 = in_vect->c2[idx_vect];

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
  //Multiply 3rd row by eta
  mat20 = (mat20)*eta;
  mat21 = (mat21)*eta;
  mat22 = (mat22)*eta;

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


//void acc_Deo( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in, ferm_param pars) {
void acc_Deo( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in,double_soa * backfield) {

  int hx, y, z, t;
#pragma acc kernels present(u) present(out) present(in)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {

          int x, xm, ym, zm, tm, xp, yp, zp, tp, idxh, eta, matdir;
          vec3 aux;

          x = 2*hx + ((y+z+t) & 0x1);

          idxh = snum_acc(x,y,z,t);

          xm = x - 1;
          xm = xm + (((xm >> 31) & 0x1) * nx);
          ym = y - 1;
          ym = ym + (((ym >> 31) & 0x1) * ny);
          zm = z - 1;
          zm = zm + (((zm >> 31) & 0x1) * nz);
          tm = t - 1;
          tm = tm + (((tm >> 31) & 0x1) * nt);

          xp = x + 1;
          xp *= (((xp-nx) >> 31) & 0x1);
          yp = y + 1;
          yp *= (((yp-ny) >> 31) & 0x1);
          zp = z + 1;
          zp *= (((zp-nz) >> 31) & 0x1);
          tp = t + 1;
          tp *= (((tp-nt) >> 31) & 0x1);

          matdir = 0;
          eta = 1;
	  // mat_vec_mul( &(u_work[snum_acc(x,y,z,t)       ]), &(in[snum_acc(xp,y,z,t)]), &aux_tmp );
          aux = mat_vec_mul( &u[matdir], idxh, eta, in, snum_acc(xp,y,z,t) );

          matdir = 2;
          eta = 1 - ( 2*(x & 0x1) ); // if (x % 2 = 0) eta = 1 else -1                       
	  // mat_vec_mul( &(u_work[snum_acc(x,y,z,t) + size ]), &(in[snum_acc(x,yp,z,t)]), &aux_tmp );
          aux = sumResult(aux, mat_vec_mul( &u[matdir], idxh, eta, in, snum_acc(x,yp,z,t) ) );

          matdir = 4;
          eta = 1 - ( 2*((x+y) & 0x1) );
	  // mat_vec_mul( &(u_work[snum_acc(x,y,z,t) + size2]), &(in[snum_acc(x,y,zp,t)]), &aux_tmp );
          aux = sumResult(aux, mat_vec_mul( &u[matdir], idxh, eta, in, snum_acc(x,y,zp,t) ) );

          matdir = 6;
          eta = 1 - ( 2*((x+y+z) & 0x1) );

#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  // mat_vec_mul( &(u_work[snum_acc(x,y,z,t) + size3]), &(in[snum_acc(x,y,z,tp)]), &aux_tmp );
          aux = sumResult(aux, mat_vec_mul( &u[matdir], idxh, eta, in, snum_acc(x,y,z,tp) ) );

	  //////////////////////////////////////////////////////////////////////////////////////////////

          matdir = 1;
          eta = 1;
	  // conjmat_vec_mul( &(u_work[sizeh + snum_acc(xm,y,z,t)      ]), &(in[ snum_acc(xm,y,z,t) ]), &aux_tmp );
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(xm,y,z,t), eta, in, snum_acc(xm,y,z,t) ) );

          matdir = 3;
          eta = 1 - ( 2*(x & 0x1) );
	  // conjmat_vec_mul( &(u_work[sizeh + snum_acc(x,ym,z,t) + size ]), &(in[ snum_acc(x,ym,z,t) ]), &aux_tmp );
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(x,ym,z,t), eta, in, snum_acc(x,ym,z,t) ) );

          matdir = 5;
          eta = 1 - ( 2*((x+y) & 0x1) );
	  // conjmat_vec_mul( &(u_work[sizeh + snum_acc(x,y,zm,t) + size2]), &(in[ snum_acc(x,y,zm,t) ]), &aux_tmp );
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(x,y,zm,t), eta, in, snum_acc(x,y,zm,t) ) );

          matdir = 7;
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(tm/(nt-1)));
#endif
	  // conjmat_vec_mul( &(u_work[sizeh + snum_acc(x,y,z,tm) + size3]), &(in[ snum_acc(x,y,z,tm) ]), &aux_tmp );
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(x,y,z,tm), eta, in, snum_acc(x,y,z,tm) ) );

	  //////////////////////////////////////////////////////////////////////////////////////////////     

          out->c0[idxh] = (aux.c0)*0.5;
          out->c1[idxh] = (aux.c1)*0.5;
          out->c2[idxh] = (aux.c2)*0.5;

        } // Loop over nxh
      } // Loop over ny   
    } // Loop over nz    
  } // Loop over nt      
}


//void acc_Doe(__restrict su3_soa * const u, __restrict vec3_soa * const out, __restrict vec3_soa * const in, ferm_param *pars, int ferm_id) {
void acc_Doe(__restrict su3_soa * const u, __restrict vec3_soa * const out, __restrict vec3_soa * const in) {

  int hx, y, z, t;

#pragma acc kernels present(u) present(out) present(in)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {

          int x, xm, ym, zm, tm, xp, yp, zp, tp, idxh, eta, matdir;
          vec3 aux;

          x = 2*hx + ((y+z+t+1) & 0x1);

          idxh = snum_acc(x,y,z,t);


          xm = x - 1;
          xm = xm + (((xm >> 31) & 0x1) * nx);
          ym = y - 1;
          ym = ym + (((ym >> 31) & 0x1) * ny);
          zm = z - 1;
          zm = zm + (((zm >> 31) & 0x1) * nz);
          tm = t - 1;
          tm = tm + (((tm >> 31) & 0x1) * nt);

          xp = x + 1;
          xp *= (((xp-nx) >> 31) & 0x1);
          yp = y + 1;
          yp *= (((yp-ny) >> 31) & 0x1);
          zp = z + 1;
          zp *= (((zp-nz) >> 31) & 0x1);
          tp = t + 1;
          tp *= (((tp-nt) >> 31) & 0x1);


          matdir = 1;
          eta = 1;
          aux = mat_vec_mul( &u[1], idxh, eta, in, snum_acc(xp,y,z,t));

          matdir = 3;
          eta = 1 - ( 2*(x & 0x1) );
          aux = sumResult(aux, mat_vec_mul( &u[matdir], idxh, eta, in, snum_acc(x,yp,z,t)) );

          matdir = 5;
          eta = 1 - ( 2*((x+y) & 0x1) );
          aux = sumResult(aux, mat_vec_mul( &u[matdir], idxh, eta, in, snum_acc(x,y,zp,t)) );

          matdir = 7;
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
          aux = sumResult(aux, mat_vec_mul( &u[matdir], idxh, eta, in, snum_acc(x,y,z,tp)) );

	  //////////////////////////////////////////////////////////////////////////////////////////////

          matdir = 0;
          eta = 1;
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(xm,y,z,t), eta, in, snum_acc(xm,y,z,t)) );

          matdir = 2;
          eta = 1 - ( 2*(x & 0x1) );
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(x,ym,z,t), eta, in, snum_acc(x,ym,z,t)) );

          matdir = 4;
          eta = 1 - ( 2*((x+y) & 0x1) );
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(x,y,zm,t), eta, in, snum_acc(x,y,zm,t)) );

          matdir = 6;
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(tm/(nt-1)));
#endif
          aux = subResult(aux, conjmat_vec_mul( &u[matdir], snum_acc(x,y,z,tm), eta, in, snum_acc(x,y,z,tm)) );

	  //////////////////////////////////////////////////////////////////////////////////////////////

          out->c0[idxh] = aux.c0*0.5;
          out->c1[idxh] = aux.c1*0.5;
          out->c2[idxh] = aux.c2*0.5;

        } // Loop over nxh
      } // Loop over ny
    } // Loop over nz
  } // Loop over nt
}


inline void fermion_matrix_multiplication( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in, __restrict vec3_soa * const temp1, ferm_param *pars,double_soa * backfield){
    acc_Doe(u,temp1,in,backfield);
    acc_Deo(u,out,temp1,backfield);
    combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass,out);// Nuova funzione in OpenAcc/fermionic_utilities.c

}
inline void fermion_matrix_multiplication_shifted( __restrict su3_soa * const u, __restrict vec3_soa * const out,  __restrict vec3_soa * const in, __restrict vec3_soa * const temp1, ferm_param *pars,double_soa * backfield, const double shift){
    acc_Doe(u,temp1,in,backfield);
    acc_Deo(u,out,temp1,backfield);
    combine_in1xferm_mass_minus_in2(in,pars->ferm_mass*pars->ferm_mass+shift,out);// Nuova funzione in OpenAcc/fermionic_utilities.c

}




#endif

