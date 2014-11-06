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
