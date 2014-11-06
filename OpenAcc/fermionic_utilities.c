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
  d_complex sum  =  conj(in_vect1->c0[idx_vect]) *  in_vect2->c0[idx_vect] ;
  sum +=  conj(in_vect1->c1[idx_vect]) *  in_vect2->c1[idx_vect] ;
  sum +=  conj(in_vect1->c2[idx_vect]) *  in_vect2->c2[idx_vect] ;

  return creal(sum);
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


void  scal_prod_openacc(vec3COM_soa *in1,vec3COM_soa *in2, double *pre, double *pim){
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
