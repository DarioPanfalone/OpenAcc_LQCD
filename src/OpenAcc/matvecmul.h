#ifndef MATVECMUL_H_
#define MATVECMUL_H_

#include "./struct_c_def.h"
#include "math.h"

#pragma acc routine seq
static inline vec3 mat_vec_mul(__restrict const su3_soa * const matrix,
															 const int idx_mat,
															 __restrict const vec3_soa * const in_vect,
															 const int idx_vect,
															 const d_complex phase)
{
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

	// load 3rd matrix row from global memory
	// d_complex mat20 = matrix->r2.c0[idx_mat];
	// d_complex mat21 = matrix->r2.c1[idx_mat];
	// d_complex mat22 = matrix->r2.c2[idx_mat];

	// compute 3rd matrix row from the first two  
	d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) ) ;
	d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) ) ;
	d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) ) ;

	out_vect.c0 = ( mat00 * vec0 ) + ( mat01 * vec1 ) + ( mat02 * vec2 );
	out_vect.c1 = ( mat10 * vec0 ) + ( mat11 * vec1 ) + ( mat12 * vec2 );
	out_vect.c2 = ( mat20 * vec0 ) + ( mat21 * vec1 ) + ( mat22 * vec2 );

	return out_vect;

}

#pragma acc routine seq
static inline vec3 conjmat_vec_mul(__restrict const su3_soa *const matrix,
																	 const int idx_mat,
																	 __restrict const vec3_soa * const in_vect,
																	 const int idx_vect,
																	 const d_complex phase)
{
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

	// load 3rd matrix row from global memory
	// d_complex mat20 = matrix->r2.c0[idx_mat];
	// d_complex mat21 = matrix->r2.c1[idx_mat];
	// d_complex mat22 = matrix->r2.c2[idx_mat];

	// compute 3rd matrix row from the first two  
	d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) );
	d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) );
	d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) );

	out_vect.c0 = ( conj(mat00) * vec0 ) + ( conj(mat10) * vec1 ) + ( conj(mat20) * vec2 );
	out_vect.c1 = ( conj(mat01) * vec0 ) + ( conj(mat11) * vec1 ) + ( conj(mat21) * vec2 );
	out_vect.c2 = ( conj(mat02) * vec0 ) + ( conj(mat12) * vec1 ) + ( conj(mat22) * vec2 );

	return out_vect;

}

#pragma acc routine seq
static inline vec3 mat_vec_mul_arg(__restrict const su3_soa *const matrix,
																	 const int idx_mat,
																	 __restrict const vec3_soa * const in_vect,
																	 const int idx_vect,
																	 __restrict const double_soa * arrarg)
{
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

	// load 3rd matrix row from global memory
	// d_complex mat20 = matrix->r2.c0[idx_mat];
	// d_complex mat21 = matrix->r2.c1[idx_mat];
	// d_complex mat22 = matrix->r2.c2[idx_mat];

	// compute 3rd matrix row from the first two  
	d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) ) ;
	d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) ) ;
	d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) ) ;

	out_vect.c0 = ( mat00 * vec0 ) + ( mat01 * vec1 ) + ( mat02 * vec2 );
	out_vect.c1 = ( mat10 * vec0 ) + ( mat11 * vec1 ) + ( mat12 * vec2 );
	out_vect.c2 = ( mat20 * vec0 ) + ( mat21 * vec1 ) + ( mat22 * vec2 );

	return out_vect;

}

#pragma acc routine seq
static inline vec3 conjmat_vec_mul_arg(__restrict const su3_soa * const matrix,
																			 const int idx_mat,
																			 __restrict const vec3_soa * const in_vect,
																			 const int idx_vect,
																			 __restrict const double_soa * const arrarg  )
{

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

	// load 3rd matrix row from global memory
	// d_complex mat20 = matrix->r2.c0[idx_mat];
	// d_complex mat21 = matrix->r2.c1[idx_mat];
	// d_complex mat22 = matrix->r2.c2[idx_mat];

	// compute 3rd matrix row from the first two  
	d_complex cmat20 = ( mat01 * mat12 ) - ( mat02 * mat11) ;
	d_complex cmat21 = ( mat02 * mat10 ) - ( mat00 * mat12) ;
	d_complex cmat22 = ( mat00 * mat11 ) - ( mat01 * mat10) ;

	out_vect.c0 = ( conj(mat00) * vec0 ) + ( conj(mat10) * vec1 ) + ( cmat20 * vec2 );
	out_vect.c1 = ( conj(mat01) * vec0 ) + ( conj(mat11) * vec1 ) + ( cmat21 * vec2 );
	out_vect.c2 = ( conj(mat02) * vec0 ) + ( conj(mat12) * vec1 ) + ( cmat22 * vec2 );

	return out_vect;

}

// 'wf = with a field (for magnetic susceptibility)'
#pragma acc routine seq
static inline vec3 mat_vec_mul_arg_wf(__restrict const su3_soa *const matrix,
																			const int idx_mat,
																			__restrict const vec3_soa * const in_vect,
																			const int idx_vect,
																			__restrict const double_soa * arrarg,
																			__restrict const double_soa* field_re, // e.g. e^(i phi_1 ) - e^(i phi_2),
																			__restrict const double_soa* field_im) // or dphi/dbz
{
	vec3 out_vect;
	double arg = arrarg->d[idx_mat];
	d_complex phase = cos(arg) + I * sin(arg);
	phase *= (field_re->d[idx_mat] + I * field_im->d[idx_mat]); // not in U(1) anymore


	d_complex vec0 = (in_vect->c0[idx_vect])*phase;
	d_complex vec1 = (in_vect->c1[idx_vect])*phase;
	d_complex vec2 = (in_vect->c2[idx_vect])*phase;

	d_complex mat00 = matrix->r0.c0[idx_mat];
	d_complex mat01 = matrix->r0.c1[idx_mat];
	d_complex mat02 = matrix->r0.c2[idx_mat];

	d_complex mat10 = matrix->r1.c0[idx_mat];
	d_complex mat11 = matrix->r1.c1[idx_mat];
	d_complex mat12 = matrix->r1.c2[idx_mat];

	// load 3rd matrix row from global memory
	// d_complex mat20 = matrix->r2.c0[idx_mat];
	// d_complex mat21 = matrix->r2.c1[idx_mat];
	// d_complex mat22 = matrix->r2.c2[idx_mat];

	// compute 3rd matrix row from the first two  
	d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) ) ;
	d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) ) ;
	d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) ) ;

	out_vect.c0 = ( mat00 * vec0 ) + ( mat01 * vec1 ) + ( mat02 * vec2 );
	out_vect.c1 = ( mat10 * vec0 ) + ( mat11 * vec1 ) + ( mat12 * vec2 );
	out_vect.c2 = ( mat20 * vec0 ) + ( mat21 * vec1 ) + ( mat22 * vec2 );

	return out_vect;

}

// 'wf = with a field (for magnetic susceptibility)'
#pragma acc routine seq
static inline vec3 conjmat_vec_mul_arg_wf(__restrict const su3_soa * const matrix,
																					const int idx_mat,
																					__restrict const vec3_soa * const in_vect,
																					const int idx_vect,
																					__restrict const double_soa* arrarg, 
																					__restrict const double_soa* field_re,   // e.g. e^(i phi_1 ) - e^(i phi_2),
																					__restrict const double_soa* field_im)   // or dphi/dbz
{
	vec3 out_vect;
	double arg = arrarg->d[idx_mat];
	d_complex phase = cos(arg) + I * sin(arg);
	phase *= (field_re->d[idx_mat] + I * field_im->d[idx_mat]); // not in U(1) anymore


	d_complex vec0 = in_vect->c0[idx_vect]*conj(phase);
	d_complex vec1 = in_vect->c1[idx_vect]*conj(phase);
	d_complex vec2 = in_vect->c2[idx_vect]*conj(phase);


	d_complex mat00 = matrix->r0.c0[idx_mat];
	d_complex mat01 = matrix->r0.c1[idx_mat];
	d_complex mat02 = matrix->r0.c2[idx_mat];

	d_complex mat10 = matrix->r1.c0[idx_mat];
	d_complex mat11 = matrix->r1.c1[idx_mat];
	d_complex mat12 = matrix->r1.c2[idx_mat];

	// load 3rd matrix row from global memory
	// d_complex mat20 = matrix->r2.c0[idx_mat];
	// d_complex mat21 = matrix->r2.c1[idx_mat];
	// d_complex mat22 = matrix->r2.c2[idx_mat];

	// compute 3rd matrix row from the first two  
	d_complex mat20 = conj( ( mat01 * mat12 ) - ( mat02 * mat11) );
	d_complex mat21 = conj( ( mat02 * mat10 ) - ( mat00 * mat12) );
	d_complex mat22 = conj( ( mat00 * mat11 ) - ( mat01 * mat10) );

	out_vect.c0 = ( conj(mat00) * vec0 ) + ( conj(mat10) * vec1 ) + ( conj(mat20) * vec2 );
	out_vect.c1 = ( conj(mat01) * vec0 ) + ( conj(mat11) * vec1 ) + ( conj(mat21) * vec2 );
	out_vect.c2 = ( conj(mat02) * vec0 ) + ( conj(mat12) * vec1 ) + ( conj(mat22) * vec2 );

	return out_vect;

}




#pragma acc routine seq
static inline vec3 sumResult ( vec3 aux, vec3 aux_tmp)
{
	aux.c0 += aux_tmp.c0;
	aux.c1 += aux_tmp.c1;
	aux.c2 += aux_tmp.c2;

	return aux;
}

#pragma acc routine seq
static inline vec3 subResult ( vec3 aux, vec3 aux_tmp)
{
	aux.c0 -= aux_tmp.c0;
	aux.c1 -= aux_tmp.c1;
	aux.c2 -= aux_tmp.c2;

	return aux;
}


#endif
