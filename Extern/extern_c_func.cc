extern "C" {
  void  apply_Deo_openacc(su3COM_soa *conf,vec3COM_soa *in,vec3COM_soa *out);
}
extern "C" {
  void  apply_Doe_openacc(su3COM_soa *conf,vec3COM_soa *in,vec3COM_soa *out);
}

extern "C" {
  void  scal_prod_openacc(vec3COM_soa *in1,vec3COM_soa *in2, double *pre, double *pim);
}

