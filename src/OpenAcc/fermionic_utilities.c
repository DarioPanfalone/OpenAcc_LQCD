//Funzioni definite qui dentro:
/*
- per il prodotto scalare
    + scal_prod_loc_2double(in1,in2,idx,d_re,d_im)         --> compute the scalar product locally
    + scal_prod_global(in1,in2,result)                     --> compute the global scalar product
- per l'inverter (e per altre cose)
    + combine_in1xm2_minus_in2_minus_in3(in1,in2,in3,out)  --> compute out = in1*m2-in2-in3
    + combine_in1xm2_minus_in2(in1,in2,out)                --> compute out = in1*m2-in2
    + combine_in1_minus_in2(in1,in2,out)                   --> compute out = in1-in2
    + assign_in_to_out(in,out)                             --> set     out = in
    + combine_in1xfactor_plus_in2(in1,fattore,in2,out)     --> compute out = in1*fattore+in2
 */


#ifndef FERMIONIC_UTILITIES_C_
#define FERMIONIC_UTILITIES_C_

#ifdef __GNUC__
 #define __restrict
 #define _POSIX_C_SOURCE 200809L // not to have warning on posix memalign
#endif

#include "./geometry.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "./fermionic_utilities.h"

#include <stdlib.h>

#ifdef MULTIDEVICE
#include <mpi.h>
#include "../Mpi/multidev.h"

d_complex scal_prod_loc(const __restrict vec3_soa * const in_vect1,
        const __restrict vec3_soa * const in_vect2	)
{
  int t;
  d_complex res = 0.0 + 0.0I;
  double * res_RI_p;
  int check = posix_memalign((void **)&res_RI_p, 64, 2*sizeof(double));

  res_RI_p[0]=0.0;
  res_RI_p[1]=0.0;
  double resR = 0.0;
  double resI = 0.0;

  int loc_h0, loc_1, loc_2, loc_3;
#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop independent gang(LOC_N3) // REDUCTION PRAGMAS ?
  for(loc_3=0; loc_3<LOC_N3; loc_3++) {
#pragma acc loop independent gang(LOC_N2/DIM_BLOCK_2) vector(DIM_BLOCK_2)
    for(loc_2=0; loc_2<LOC_N2; loc_2++) {
#pragma acc loop independent gang(LOC_N1/DIM_BLOCK_1) vector(DIM_BLOCK_1)
      for(loc_1=0; loc_1<LOC_N1; loc_1++) {
#pragma acc loop independent vector(DIM_BLOCK_0)
        for(loc_h0=0; loc_h0 <LOC_N0H; loc_h0++) {
            int loc_0, 
                lnh_0 , lnh_1, lnh_2, lnh_3;

            loc_0 = 2*loc_h0 + ((loc_1+loc_2+loc_3) & 0x1);
            lnh_0 = loc_0 + D0_HALO;
            lnh_1 = loc_1 + D1_HALO;
            lnh_2 = loc_2 + D2_HALO;
            lnh_3 = loc_3 + D3_HALO;
            int lnh_idxh = lnh_to_lnh_snum(lnh_0,lnh_1,lnh_2,lnh_3);

            d_complex color_sum  =  conj(in_vect1->c0[t]) *  in_vect2->c0[t] ;
            color_sum +=  conj(in_vect1->c1[t]) *  in_vect2->c1[t] ;
            color_sum +=  conj(in_vect1->c2[t]) *  in_vect2->c2[t] ;

            resR+=creal(color_sum);
            resI+=cimag(color_sum);



        }
      }
    }
  }
  res = resR+resI*1.0I;
  return res;
}
double real_scal_prod_loc(  const __restrict vec3_soa * const in_vect1,
        const __restrict vec3_soa * const in_vect2 )
{
  int t;
  double res_R_p;
  double resR = 0.0;

  int loc_h0, loc_1, loc_2, loc_3;
#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop independent gang(LOC_N3)  // REDUCTION PRAGMAS ?
  for(loc_3=0; loc_3<LOC_N3; loc_3++) {
#pragma acc loop independent gang(LOC_N2/DIM_BLOCK_2) vector(DIM_BLOCK_2)
    for(loc_2=0; loc_2<LOC_N2; loc_2++) {
#pragma acc loop independent gang(LOC_N1/DIM_BLOCK_1) vector(DIM_BLOCK_1)
      for(loc_1=0; loc_1<LOC_N1; loc_1++) {
#pragma acc loop independent vector(DIM_BLOCK_0)
        for(loc_h0=0; loc_h0 <LOC_N0H; loc_h0++) {
            int loc_0, 
                lnh_0 , lnh_1, lnh_2, lnh_3;

            loc_0 = 2*loc_h0 + ((loc_1+loc_2+loc_3) & 0x1);
            lnh_0 = loc_0 + D0_HALO;
            lnh_1 = loc_1 + D1_HALO;
            lnh_2 = loc_2 + D2_HALO;
            lnh_3 = loc_3 + D3_HALO;
            int lnh_idxh = lnh_to_lnh_snum(lnh_0,lnh_1,lnh_2,lnh_3);

            res_R_p=scal_prod_loc_1double(in_vect1,in_vect2,
                    (const int)lnh_idxh);
            resR+=res_R_p;
        }
      }
    }
  }
  return resR;
}
double l2norm2_loc_vol( const __restrict vec3_soa * const in_vect1)
{
int t;
  double res_R_p;
  double resR = 0.0;

  int loc_h0, loc_1, loc_2, loc_3;
#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop independent gang(LOC_N3)  // REDUCTION PRAGMAS ?
  for(loc_3=0; loc_3<LOC_N3; loc_3++) {
#pragma acc loop independent gang(LOC_N2/DIM_BLOCK_2) vector(DIM_BLOCK_2)
    for(loc_2=0; loc_2<LOC_N2; loc_2++) {
#pragma acc loop independent gang(LOC_N1/DIM_BLOCK_1) vector(DIM_BLOCK_1)
      for(loc_1=0; loc_1<LOC_N1; loc_1++) {
#pragma acc loop independent vector(DIM_BLOCK_0)
        for(loc_h0=0; loc_h0 <LOC_N0H; loc_h0++) {
            int loc_0, 
                lnh_0 , lnh_1, lnh_2, lnh_3;

            loc_0 = 2*loc_h0 + ((loc_1+loc_2+loc_3) & 0x1);
            lnh_0 = loc_0 + D0_HALO;
            lnh_1 = loc_1 + D1_HALO;
            lnh_2 = loc_2 + D2_HALO;
            lnh_3 = loc_3 + D3_HALO;
            int lnh_idxh = lnh_to_lnh_snum(lnh_0,lnh_1,lnh_2,lnh_3);

            res_R_p=l2norm2_loc(in_vect1,lnh_idxh);
            resR+=res_R_p;
        }
      }
    }
  }
  return resR;
}


d_complex scal_prod_loc_1Dcut(  
        const __restrict vec3_soa * const in_vect1,
        const __restrict vec3_soa * const in_vect2 )
{
    int t;
    double resR = 0.0;
    double resI = 0.0;

#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop reduction(+:resR) // DEBUG
    for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
        d_complex color_sum  =  conj(in_vect1->c0[t]) *  in_vect2->c0[t] ;
        color_sum +=  conj(in_vect1->c1[t]) *  in_vect2->c1[t] ;
        color_sum +=  conj(in_vect1->c2[t]) *  in_vect2->c2[t] ;

        resR+=creal(color_sum);
        resI+=cimag(color_sum);

    }
    return resR+resI*I;

}
double real_scal_prod_loc_1Dcut(  const __restrict vec3_soa * in_vect1,
        const __restrict vec3_soa * in_vect2			       )
{
  int t;
  double res_R_p;
  double resR = 0.0;

#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop reduction(+:resR) // DEBUG
  for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
    res_R_p=scal_prod_loc_1double(in_vect1,in_vect2,t);
    resR+=res_R_p;
  }
  return resR;
}
double l2norm2_loc_1Dcut( 
        const __restrict vec3_soa * const in_vect1		       )
{
  int t;
  double res_R_p;
  double resR = 0.0;

#pragma acc kernels present(in_vect1)
#pragma acc loop reduction(+:resR) // DEBUG
  for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
    res_R_p=l2norm2_loc(in_vect1,t);
    resR+=res_R_p;
  }
  return resR;
}

// MPI REQUIRING FUNCTIONS
d_complex scal_prod_global(  const __restrict vec3_soa * const in_vect1,
        const __restrict vec3_soa * const in_vect2 )
{
    // RETURNS THE TOTAL RESULT USING MPI

     d_complex total_res=0;
     d_complex local_res;

     local_res = scal_prod_loc_1Dcut(in_vect1,in_vect2); 
     
     if(mdevinfo.myrank ==0){
         int rank;
         d_complex temp_loc_res = local_res;
         total_res += temp_loc_res; 
         for(rank = 1; rank < NRANKS ; rank++ ){
             MPI_Recv(&temp_loc_res,2,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
             total_res += temp_loc_res;

         }
         for(rank = 1; rank < NRANKS ;rank++ )
             MPI_Send(&total_res,2,MPI_DOUBLE,rank,0,MPI_COMM_WORLD);

     }else{
          MPI_Send(&local_res,2,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
          MPI_Recv(&total_res,2,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
     }
     return total_res;
}
inline double real_scal_prod_global(  
        const __restrict vec3_soa * const in_vect1,
        const __restrict vec3_soa * const in_vect2 )
{
    // RETURNS THE TOTAL RESULT USING MPI
     double total_res,local_res;
     //local_res = real_scal_prod_loc(in_vect1,in_vect2); 
     local_res = real_scal_prod_loc_1Dcut(in_vect1,in_vect2); 
     MPI_Allreduce((void*)&local_res,(void*)&total_res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
     return total_res;
}
inline double l2norm2_global( 
        const __restrict vec3_soa * const in_vect1  )
{

     double total_res,local_res;
     //local_res = l2norm2_loc(in_vect1); 
     local_res = l2norm2_loc_1Dcut(in_vect1); 
     MPI_Allreduce((void*)&local_res,(void*)&total_res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
     return total_res;
}


#else // NO MULTIDEV


d_complex scal_prod_global(  const __restrict vec3_soa * in_vect1,
			     const __restrict vec3_soa * in_vect2 )
{
  int t;

  double resR = 0.0;
  double resI = 0.0;

#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop reduction(+:resR) reduction(+:resI)
  for(t=0; t<sizeh; t++) {
  d_complex color_sum  =  conj(in_vect1->c0[t]) *  in_vect2->c0[t] ;
  color_sum +=  conj(in_vect1->c1[t]) *  in_vect2->c1[t] ;
  color_sum +=  conj(in_vect1->c2[t]) *  in_vect2->c2[t] ;

  resR+=creal(color_sum);
  resI+=cimag(color_sum);

  }
  return resR+resI*I;
}
double real_scal_prod_global(  __restrict const  vec3_soa * in_vect1,
			       __restrict const vec3_soa * in_vect2 )
{
  int t;
  double res_R_p;
  double resR = 0.0;

#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop reduction(+:resR)
  for(t=0; t<sizeh; t++) {
    res_R_p=scal_prod_loc_1double(in_vect1,in_vect2,t);
    resR+=res_R_p;
  }
  return resR;
}
double l2norm2_global( __restrict  const vec3_soa * in_vect1 )
{
  int t;
  double res_R_p;
  double resR = 0.0;

#pragma acc kernels present(in_vect1)
#pragma acc loop reduction(+:resR)
  for(t=0; t<sizeh; t++) {
    res_R_p=l2norm2_loc(in_vect1,t);
    resR+=res_R_p;
  }
  return resR;
}
#endif



// NEXT, only "embarassingly parallel" functions, not requiring any reduction.

void combine_in1xfactor_plus_in2(
        __restrict const vec3_soa * in_vect1,
        const double factor,
        __restrict const vec3_soa * in_vect2,
        __restrict vec3_soa * const out)
{

#pragma acc kernels present(in_vect1) present(in_vect2) present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih]*factor)+in_vect2->c0[ih];
    out->c1[ih]=(in_vect1->c1[ih]*factor)+in_vect2->c1[ih];
    out->c2[ih]=(in_vect1->c2[ih]*factor)+in_vect2->c2[ih];
  }
}


void multiply_fermion_x_doublefactor( 
        __restrict vec3_soa * const in1,
        const double factor )
{

#pragma acc kernels present(in1)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in1->c0[ih] = factor * (in1->c0[ih]);
    in1->c1[ih] = factor * (in1->c1[ih]);
    in1->c2[ih] = factor * (in1->c2[ih]);
  }
}

void combine_add_factor_x_in2_to_in1( 
        __restrict vec3_soa * const in1,
        __restrict const vec3_soa * const in2, 
        double factor)
{
#pragma acc kernels present(in1) present(in2)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in1->c0[ih] += (factor) *  (in2->c0[ih]);
    in1->c1[ih] += (factor) *  (in2->c1[ih]);
    in1->c2[ih] += (factor) *  (in2->c2[ih]);
  }
}




void combine_in1xferm_mass2_minus_in2_minus_in3(   
        __restrict const vec3_soa * const in_vect1,
        double ferm_mass,
        __restrict const vec3_soa * const in_vect2,
        __restrict const vec3_soa * const in_vect3,
        __restrict vec3_soa * const out)
{ // l'ultimo argomento e' l'unico che viene modificato
#pragma acc kernels present(in_vect1) present(in_vect2) present(in_vect3) present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih]*ferm_mass)-in_vect2->c0[ih]-in_vect3->c0[ih];
    out->c1[ih]=(in_vect1->c1[ih]*ferm_mass)-in_vect2->c1[ih]-in_vect3->c1[ih];
    out->c2[ih]=(in_vect1->c2[ih]*ferm_mass)-in_vect2->c2[ih]-in_vect3->c2[ih];
  }
}


void combine_inside_loop( __restrict vec3_soa * const vect_out,
			  __restrict vec3_soa * const vect_r,
			  const __restrict vec3_soa * const vect_s,
			  const __restrict vec3_soa * const vect_p,
			  const double omega)
{
#pragma acc kernels present(vect_out) present(vect_r) present(vect_s) present(vect_p)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    //out+=omega*p
    vect_out->c0[ih] += (vect_p->c0[ih]*omega);
    vect_out->c1[ih] += (vect_p->c1[ih]*omega);
    vect_out->c2[ih] += (vect_p->c2[ih]*omega);
    //r-=omega*s
    vect_r->c0[ih]   -= (vect_s->c0[ih]*omega);
    vect_r->c1[ih]   -= (vect_s->c1[ih]*omega);
    vect_r->c2[ih]   -= (vect_s->c2[ih]*omega);
  }
}

void combine_in1xferm_mass_minus_in2(
        __restrict const vec3_soa * const in_vect1,
        double ferm_mass2,
        __restrict vec3_soa * const in_vect2)
{
#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in_vect2->c0[ih]=(in_vect1->c0[ih]*ferm_mass2)-in_vect2->c0[ih];
    in_vect2->c1[ih]=(in_vect1->c1[ih]*ferm_mass2)-in_vect2->c1[ih];
    in_vect2->c2[ih]=(in_vect1->c2[ih]*ferm_mass2)-in_vect2->c2[ih];
  }
}

void combine_in1_minus_in2(  __restrict const vec3_soa * in_vect1,
                             __restrict const vec3_soa * in_vect2,
                            __restrict vec3_soa * out)
{
#pragma acc kernels present(in_vect1) present(in_vect2) present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih])-in_vect2->c0[ih];
    out->c1[ih]=(in_vect1->c1[ih])-in_vect2->c1[ih];
    out->c2[ih]=(in_vect1->c2[ih])-in_vect2->c2[ih];
  }
}

void assign_in_to_out(
        __restrict const vec3_soa * in_vect1,
        __restrict vec3_soa * out)
{
  SETINUSE(out);
#pragma acc kernels present(in_vect1)  present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih]);
    out->c1[ih]=(in_vect1->c1[ih]);
    out->c2[ih]=(in_vect1->c2[ih]);
  }

}
// Altre funzioni aggiunte per l'algebra lineare multishift "versatilizzato"
void set_vec3_soa_to_zero( __restrict vec3_soa* const fermion)
{
  //  printf(" puntatore incriminato (%p) \n",fermion);
  SETINUSE((fermion));
#pragma acc kernels present(fermion)
#pragma acc loop independent 
   // printf("Setting to zero memory from %p to %p ", fermion->c0 , &(fermion->c0[sizeh-1])+1);
   // printf(", %p to %p  ", fermion->c1 , &(fermion->c1[sizeh-1])+1);
   // printf(", %p to %p\n", fermion->c2 , &(fermion->c2[sizeh-1])+1);
    for(int i= 0; i < sizeh ; i++){
        fermion->c0[i]=0;
        fermion->c1[i]=0;
        fermion->c2[i]=0;
    }
}
void multiple_combine_in1_minus_in2x_factor_back_into_in1( 
        __restrict vec3_soa * const out, 
        __restrict const vec3_soa * const in, 
        const int maxiter, 
        __restrict const int * const flag, 
        __restrict const double * const omegas)
{
  int ia, ih;
  
#pragma acc kernels present(in) present(out) copyin(omegas[0:maxiter]) copyin(flag[0:maxiter])
  {
#pragma acc cache(omegas[0:maxiter])
#pragma acc cache(flag[0:maxiter])
#pragma acc loop independent gang
    for (ia=0; ia<maxiter; ia++) {
#pragma acc loop independent gang vector(512)
      for (ih=0; ih<sizeh; ih++) {
	if (flag[ia] == 1) {
	  double factor = omegas[ia];
	  out[ia].c0[ih] -= (factor)*(in[ia].c0[ih]);
	  out[ia].c1[ih] -= (factor)*(in[ia].c1[ih]);
	  out[ia].c2[ih] -= (factor)*(in[ia].c2[ih]);
	} //if
      } //ih
    } //ia
  } //acc kernels
}
void multiple1_combine_in1_x_fact1_plus_in2_x_fact2_back_into_in1( 
        __restrict vec3_soa * const in1,
        int maxiter, 
        __restrict const int * const flag, 
        __restrict const double * const gammas, 
        __restrict const vec3_soa * const in2, 
        __restrict const double * const zeta_iii ) 
{

  int ia, ih;

  #pragma acc kernels present(in1) present(in2) copyin(gammas[0:maxiter]) copyin(flag[0:maxiter]) copyin(zeta_iii[0:maxiter])
  {

  #pragma acc cache(gammas[0:maxiter]) 
  #pragma acc cache(zeta_iii[0:maxiter])
  #pragma acc cache(flag[0:maxiter]) 

  #pragma acc loop independent gang
  for (ia=0; ia<maxiter; ia++) {
    #pragma acc loop independent gang vector(512)
    for(ih=0; ih<sizeh; ih++) {

      if (flag[ia] == 1) {

        double fact1 = gammas[ia];
        double fact2 = zeta_iii[ia];

        in1[ia].c0[ih] = fact1 * (in1[ia].c0[ih]) + fact2 * (in2->c0[ih]);
        in1[ia].c1[ih] = fact1 * (in1[ia].c1[ih]) + fact2 * (in2->c1[ih]);
        in1[ia].c2[ih] = fact1 * (in1[ia].c2[ih]) + fact2 * (in2->c2[ih]);

      } //if flag

    } //ih

  } //ia

  } //acc data

}
void combine_in1_x_fact1_minus_in2_back_into_in2( 
        __restrict const vec3_soa * const in1, 
        double fact1, 
        __restrict vec3_soa * const in2 )
{

  int ih;


  #pragma acc kernels present(in1) present(in2) 
  {


#pragma acc loop independent gang vector(512)
    for(ih=0; ih<sizeh; ih++) {


        in2->c0[ih] = fact1 * (in1->c0[ih]) - (in2->c0[ih]);
        in2->c1[ih] = fact1 * (in1->c1[ih]) - (in2->c1[ih]);
        in2->c2[ih] = fact1 * (in1->c2[ih]) - (in2->c2[ih]);

    } //ih
  } //acc data

}
void combine_in1_minus_in2_allxfact( 
        __restrict const vec3_soa * const in1, 
        __restrict const vec3_soa * const in2, 
        double fact,
        __restrict vec3_soa * const out )
{
  int ih;
#pragma acc kernels present(in1) present(in2) present(out) 
  {
#pragma acc loop independent gang vector(512)
    for(ih=0; ih<sizeh; ih++) {
        out->c0[ih] = fact * (in1->c0[ih] - in2->c0[ih]);
        out->c1[ih] = fact * (in1->c1[ih] - in2->c1[ih]);
        out->c2[ih] = fact * (in1->c2[ih] - in2->c2[ih]);
    } //ih
  } //acc data
}

#endif


