extern "C" {
  void UPDATE_ACC(su3COM_soa *conf,double residue_metro,double residue_md,const COM_RationalApprox *approx1,const COM_RationalApprox *approx2,const COM_RationalApprox *approx3,const COM_MultiFermion *in,thmatCOM_soa * com_mom);
}

extern "C" {
  void UPDATE_ACC_UNOSTEP(su3COM_soa *conf,double res_metro, double res_md,COM_RationalApprox *approx1,COM_RationalApprox *approx2,COM_RationalApprox *approx3,int id_iter);
}

extern "C" {
  void THERM_UPDATE_ACC_UNOSTEP_NOMETRO(su3COM_soa *conf,double res_metro, double res_md,COM_RationalApprox *approx1,COM_RationalApprox *approx2,COM_RationalApprox *approx3,int id_iter);
}

extern "C" {
  int UPDATE_ACC_UNOSTEP_METRO(su3COM_soa *conf,double res_metro, double res_md,COM_RationalApprox *approx1,COM_RationalApprox *approx2,COM_RationalApprox *approx3,int id_iter,int acc);
}

