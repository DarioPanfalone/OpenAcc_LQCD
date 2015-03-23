extern "C" {
  void  apply_Deo_openacc(const su3COM_soa *conf,const vec3COM_soa *in,vec3COM_soa *out);
}
extern "C" {
  void  apply_Doe_openacc(const su3COM_soa *conf,const vec3COM_soa *in,vec3COM_soa *out);
}

extern "C" {
  void  scal_prod_openacc(const vec3COM_soa *in1,const vec3COM_soa *in2, double *pre, double *pim);
}

extern "C" {
  //  void invert_openacc_full(su3COM_soa  *conf,vec3COM_soa *out,vec3COM_soa *in,double res,vec3COM_soa *trialSol = NULL);
  void invert_openacc_full(const su3COM_soa  *conf,vec3COM_soa *out,const vec3COM_soa *in,double res,const vec3COM_soa *trialSol);
}

extern "C" {
  void find_min_max_openacc(const  su3COM_soa  *conf, const vec3COM_soa *gaussian1,  const vec3COM_soa *gaussian2,double *minmax);
}

extern "C" {
  //  void multips_invert_openacc_full(const su3COM_soa  *conf,COM_ShiftMultiFermion *out, const COM_MultiFermion *in, REAL res, COM_RationalApprox approx);
  void multips_invert_openacc_full(const su3COM_soa  *conf,COM_ShiftMultiFermion *out, COM_MultiFermion *in, REAL res, COM_RationalApprox *approx);
}
extern "C" {
  void first_inv_approx_calc_openacc(const su3COM_soa  *conf,COM_MultiFermion *out, const COM_MultiFermion *in, double res,const COM_RationalApprox *approx);
}

extern "C" {
  void  calc_plaquette_openacc(const su3COM_soa *conf);
}
extern "C" {
  //  void  calc_staples_openacc(const su3COM_soa *conf);
  void  calc_staples_openacc(const su3COM_soa *conf,su3COM_soa *COM_staples);
}

extern "C" {
  void fermion_force_openacc(const su3COM_soa  *conf, tamatCOM_soa *out, const COM_MultiFermion *in, double res, const COM_RationalApprox *approx);
}

extern "C" {
  void  calc_ipdot_gauge_openacc(const su3COM_soa *conf,tamatCOM_soa * com_ipdot);
}


extern "C" {
  void mom_exp_times_conf_openacc(su3COM_soa *conf,const thmatCOM_soa * com_mom);
}
