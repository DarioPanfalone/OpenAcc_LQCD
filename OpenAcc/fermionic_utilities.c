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


static inline void scal_prod_loc_2double( const __restrict vec3_soa * const in_vect1,
                                          const __restrict vec3_soa * const in_vect2,
                                          const int idx_vect,
                                          double *res_ri
                                          ) {

  d_complex sum  =  conj(in_vect1->c0[idx_vect]) *  in_vect2->c0[idx_vect] ;
  sum +=  conj(in_vect1->c1[idx_vect]) *  in_vect2->c1[idx_vect] ;
  sum +=  conj(in_vect1->c2[idx_vect]) *  in_vect2->c2[idx_vect] ;

  res_ri[0] = creal(sum);
  res_ri[1] = cimag(sum);
}
static inline double scal_prod_loc_1double( const __restrict vec3_soa * const in_vect1,
                                          const __restrict vec3_soa * const in_vect2,
                                          const int idx_vect
                                          ) {
  // questa routine potrebbe essere migliorata perche per adesso
  // nel calcolare la parte reale del prodotto scalare si calcola anche
  // la parte immaginaria che poi viene buttata via.
  // forse si puo riscrivere il tutto in modo piu esplicito facendogli fare
  // solo i conti che contribuiscono al calcolo della parte reale.
  /*
  d_complex sum  =  conj(in_vect1->c0[idx_vect]) *  in_vect2->c0[idx_vect] ;
  sum +=  conj(in_vect1->c1[idx_vect]) *  in_vect2->c1[idx_vect] ;
  sum +=  conj(in_vect1->c2[idx_vect]) *  in_vect2->c2[idx_vect] ;
  */

  double sum  = creal(in_vect1->c0[idx_vect]) * creal(in_vect2->c0[idx_vect]) + cimag(in_vect1->c0[idx_vect]) * cimag(in_vect2->c0[idx_vect]);
         sum += creal(in_vect1->c1[idx_vect]) * creal(in_vect2->c1[idx_vect]) + cimag(in_vect1->c1[idx_vect]) * cimag(in_vect2->c1[idx_vect]);
         sum += creal(in_vect1->c2[idx_vect]) * creal(in_vect2->c2[idx_vect]) + cimag(in_vect1->c2[idx_vect]) * cimag(in_vect2->c2[idx_vect]);
  return sum;
}

static inline double l2norm2_loc( const __restrict vec3_soa * const in_vect1,
				  const int idx_vect
				  ) {
  double sum = creal(in_vect1->c0[idx_vect])*creal(in_vect1->c0[idx_vect])+cimag(in_vect1->c0[idx_vect])*cimag(in_vect1->c0[idx_vect]);
  sum += creal(in_vect1->c1[idx_vect])*creal(in_vect1->c1[idx_vect])+cimag(in_vect1->c1[idx_vect])*cimag(in_vect1->c1[idx_vect]);
  sum += creal(in_vect1->c2[idx_vect])*creal(in_vect1->c2[idx_vect])+cimag(in_vect1->c2[idx_vect])*cimag(in_vect1->c2[idx_vect]);
  return sum;
}



d_complex scal_prod_global(  const __restrict vec3_soa * const in_vect1,
			     const __restrict vec3_soa * const in_vect2
			     ){
  int t;
  d_complex res = 0.0 + 0.0I;
  double * res_RI_p;
  posix_memalign((void **)&res_RI_p, 64, 2*sizeof(double));

  res_RI_p[0]=0.0;
  res_RI_p[1]=0.0;
  double resR = 0.0;
  double resI = 0.0;

#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop reduction(+:resR) reduction(+:resI)
  for(t=0; t<sizeh; t++) {
    scal_prod_loc_2double(in_vect1,in_vect2,t,res_RI_p);
    resR+=res_RI_p[0];
    resI+=res_RI_p[1];
  }
  res = resR+resI*1.0I;
  return res;
}

double real_scal_prod_global(  const __restrict vec3_soa * const in_vect1,
			       const __restrict vec3_soa * const in_vect2
			       ){
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

double l2norm2_global( const __restrict vec3_soa * const in_vect1
		       ){
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


void  scal_prod_openacc(const vec3COM_soa *in1,const vec3COM_soa *in2, double *pre, double *pim){
  vec3_soa * ferm1_acc;
  vec3_soa * ferm2_acc;
  posix_memalign((void **)&ferm1_acc, ALIGN, sizeof(vec3_soa));
  posix_memalign((void **)&ferm2_acc, ALIGN, sizeof(vec3_soa));

  convert_vec3COM_soa_to_vec3_soa(in1,ferm1_acc);
  convert_vec3COM_soa_to_vec3_soa(in2,ferm2_acc);
  d_complex result;

#pragma acc data copyin(ferm1_acc[0:1]) copyin(ferm2_acc[0:1])
  {
    result = scal_prod_global(ferm1_acc,ferm2_acc);
  }

  *pre = creal(result);
  *pim = cimag(result);

  free(ferm1_acc);
  free(ferm2_acc);
}


void combine_in1xfactor_plus_in2(  const __restrict vec3_soa * const in_vect1,
				   const double factor,
                                   const __restrict vec3_soa * const in_vect2,
                                   __restrict vec3_soa * const out){

#pragma acc kernels present(in_vect1) present(in_vect2) present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih]*factor)+in_vect2->c0[ih];
    out->c1[ih]=(in_vect1->c1[ih]*factor)+in_vect2->c1[ih];
    out->c2[ih]=(in_vect1->c2[ih]*factor)+in_vect2->c2[ih];
  }
}


void multiply_fermion_x_doublefactor(  __restrict vec3_soa * const in1,
                                       const double factor
				       ){

#pragma acc kernels present(in1)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in1->c0[ih] = factor * (in1->c0[ih]);
    in1->c1[ih] = factor * (in1->c1[ih]);
    in1->c2[ih] = factor * (in1->c2[ih]);
  }
}

void combine_add_factor_x_in2_to_in1( __restrict vec3_soa * const in1,const __restrict vec3_soa * const in2, double factor){
#pragma acc kernels present(in1) present(in2)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in1->c0[ih] += (factor) *  (in2->c0[ih]);
    in1->c1[ih] += (factor) *  (in2->c1[ih]);
    in1->c2[ih] += (factor) *  (in2->c2[ih]);
  }
}



void combine_in1xm2_minus_in2_minus_in3(  const __restrict vec3_soa * const in_vect1,
                                          const __restrict vec3_soa * const in_vect2,
                                          const __restrict vec3_soa * const in_vect3,
                                          __restrict vec3_soa * const out){
#pragma acc kernels present(in_vect1) present(in_vect2) present(in_vect3) present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih]*mass2)-in_vect2->c0[ih]-in_vect3->c0[ih];
    out->c1[ih]=(in_vect1->c1[ih]*mass2)-in_vect2->c1[ih]-in_vect3->c1[ih];
    out->c2[ih]=(in_vect1->c2[ih]*mass2)-in_vect2->c2[ih]-in_vect3->c2[ih];
  }
}

void combine_before_loop(  const __restrict vec3_soa * const vect_in,
			   const __restrict vec3_soa * const vect_out,
			   __restrict vec3_soa * const vect_r,
			   __restrict vec3_soa * const vect_s,
			   __restrict vec3_soa * const vect_p){
#pragma acc kernels present(vect_in) present(vect_out) present(vect_r) present(vect_s) present(vect_p)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    //s=m^2*out-s
    vect_s->c0[ih] = (vect_out->c0[ih]*mass2)-vect_s->c0[ih];
    vect_s->c1[ih] = (vect_out->c1[ih]*mass2)-vect_s->c1[ih];
    vect_s->c2[ih] = (vect_out->c2[ih]*mass2)-vect_s->c2[ih];
    //r=in-s
    vect_r->c0[ih] =  vect_in->c0[ih]-vect_s->c0[ih];
    vect_r->c1[ih] =  vect_in->c1[ih]-vect_s->c1[ih];
    vect_r->c2[ih] =  vect_in->c2[ih]-vect_s->c2[ih];
    //p=r
    vect_p->c0[ih] =  vect_r->c0[ih];
    vect_p->c1[ih] =  vect_r->c1[ih];
    vect_p->c2[ih] =  vect_r->c2[ih];
  }
}

void combine_inside_loop( __restrict vec3_soa * const vect_out,
			  __restrict vec3_soa * const vect_r,
			  const __restrict vec3_soa * const vect_s,
			  const __restrict vec3_soa * const vect_p,
			  const double omega){
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
void combine_in1xm2_minus_in2(const __restrict vec3_soa * const in_vect1,
                              __restrict vec3_soa * const in_vect2){
#pragma acc kernels present(in_vect1) present(in_vect2)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in_vect2->c0[ih]=(in_vect1->c0[ih]*mass2)-in_vect2->c0[ih];
    in_vect2->c1[ih]=(in_vect1->c1[ih]*mass2)-in_vect2->c1[ih];
    in_vect2->c2[ih]=(in_vect1->c2[ih]*mass2)-in_vect2->c2[ih];
  }
}

/*
//OLD VERSION WITH 3 ARGUMENTS
void combine_in1xm2_minus_in2(  const __restrict vec3_soa * const in_vect1,
                                const __restrict vec3_soa * const in_vect2,
                                __restrict vec3_soa * const out){
#pragma acc kernels present(in_vect1) present(in_vect2) present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih]*mass2)-in_vect2->c0[ih];
    out->c1[ih]=(in_vect1->c1[ih]*mass2)-in_vect2->c1[ih];
    out->c2[ih]=(in_vect1->c2[ih]*mass2)-in_vect2->c2[ih];
  }
}
*/

void combine_in1_minus_in2(  const __restrict vec3_soa * const in_vect1,
                             const __restrict vec3_soa * const in_vect2,
                             __restrict vec3_soa * const out){
#pragma acc kernels present(in_vect1) present(in_vect2) present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih])-in_vect2->c0[ih];
    out->c1[ih]=(in_vect1->c1[ih])-in_vect2->c1[ih];
    out->c2[ih]=(in_vect1->c2[ih])-in_vect2->c2[ih];
  }
}

void assign_in_to_out(  const __restrict vec3_soa * const in_vect1,
                        __restrict vec3_soa * const out){
#pragma acc kernels present(in_vect1)  present(out)
#pragma acc loop independent 
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=(in_vect1->c0[ih]);
    out->c1[ih]=(in_vect1->c1[ih]);
    out->c2[ih]=(in_vect1->c2[ih]);
  }

}


void set_shiftmulti_to_zero( const int iorder, const int ips, __restrict ACC_ShiftMultiFermion * const out){

#pragma acc kernels present(out)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    out->shiftmulti[iorder][ips].c0[ih] = 0.0 + 0.0*I;
    out->shiftmulti[iorder][ips].c1[ih] = 0.0 + 0.0*I;
    out->shiftmulti[iorder][ips].c2[ih] = 0.0 + 0.0*I;
  }

}

void extract_pseudo_and_assign_to_fermion( const  __restrict ACC_MultiFermion * const in, const int ips,  __restrict vec3_soa * const out){
#pragma acc kernels present(in) present(out)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=in->multi[ips].c0[ih];
    out->c1[ih]=in->multi[ips].c1[ih];
    out->c2[ih]=in->multi[ips].c2[ih];
  }
}

void extract_pseudo_and_assign_to_shiftfermion(const __restrict ACC_MultiFermion * const in, int ips, __restrict ACC_ShiftFermion * const out, int ia){
#pragma acc kernels present(in) present(out)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    out->shift[ia].c0[ih]=in->multi[ips].c0[ih];
    out->shift[ia].c1[ih]=in->multi[ips].c1[ih];
    out->shift[ia].c2[ih]=in->multi[ips].c2[ih];
  }
}


void combine_shiftmulti_minus_shiftfermion_x_factor_back_into_shiftmulti(__restrict ACC_ShiftMultiFermion * const out,const  __restrict ACC_ShiftFermion * const in,int ips, int ia, double factor){
#pragma acc kernels present(in) present(out)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    out->shiftmulti[ia][ips].c0[ih] -= (factor)*(in->shift[ia].c0[ih]);
    out->shiftmulti[ia][ips].c1[ih] -= (factor)*(in->shift[ia].c1[ih]);
    out->shiftmulti[ia][ips].c2[ih] -= (factor)*(in->shift[ia].c2[ih]);
  }
}

void combine_shiftferm_x_fact1_plus_ferm_x_fact2_back_into_shiftferm( __restrict ACC_ShiftFermion * const in1, int ia, double fact1, const __restrict vec3_soa * const in2, double fact2){
#pragma acc kernels present(in1) present(in2)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in1->shift[ia].c0[ih] = (fact1) * (in1->shift[ia].c0[ih]) + (fact2) * (in2->c0[ih]);
    in1->shift[ia].c1[ih] = (fact1) * (in1->shift[ia].c1[ih]) + (fact2) * (in2->c1[ih]);
    in1->shift[ia].c2[ih] = (fact1) * (in1->shift[ia].c2[ih]) + (fact2) * (in2->c2[ih]);
  }
}

void extract_from_shiftmulti_and_assign_to_fermion( const  __restrict ACC_ShiftMultiFermion * const in, const int ia, const int ips,  __restrict vec3_soa * const out){
#pragma acc kernels present(in) present(out)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    out->c0[ih]=in->shiftmulti[ia][ips].c0[ih];
    out->c1[ih]=in->shiftmulti[ia][ips].c1[ih];
    out->c2[ih]=in->shiftmulti[ia][ips].c2[ih];
  }
}

void combine_in1_x_fact_minus_in2_minus_multiin3_back_into_in1( __restrict vec3_soa * const in1,double fact, const __restrict vec3_soa * const in2, const __restrict ACC_MultiFermion * const in3,int ips){
#pragma acc kernels present(in1) present(in2) present(in3)
#pragma acc loop independent
  for(int ih=0; ih<sizeh; ih++) {
    in1->c0[ih] = (in1->c0[ih]) * fact - (in2->c0[ih]) - (in3->multi[ips].c0[ih]);
    in1->c1[ih] = (in1->c1[ih]) * fact - (in2->c1[ih]) - (in3->multi[ips].c1[ih]);
    in1->c2[ih] = (in1->c2[ih]) * fact - (in2->c2[ih]) - (in3->multi[ips].c2[ih]);
  }
}


#endif


