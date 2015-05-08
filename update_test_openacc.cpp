#include <iostream>
#include <ctime>
//constants, global variables, geometric parameters et cetera 
#include "./Include/global_const.cc"
#include "./Include/global_macro.cc"
#include "./Include/global_var.cc"

//#define DEBUG_MODE  //useful in some cases
#include "./Init/init.cc"
#include "./Meas/gaugemeas.cc"

using namespace std;



int main(int argc,char **argv){

  cout.precision(18);
  cerr.precision(18);
  rationalapprox_calc();
  
  init(0);
  
  cout << "Initialized Random Gauge Matrix.\n";
  
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////// CONFRONTO DINAMICA MOLECOLARE    /////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //    void multistep_2MN_ACC(su3COM_soa *conf,double res,const COM_RationalApprox *approx,const COM_MultiFermion *in,thmatCOM_soa * com_mom);
  
  //////  TRADUCO LA CONFIGURAZIONE    //////////////////////////////////////////////////////
  su3COM_soa conf_soaCOM[8];
  for(int index=0;index<8;index++)   gauge_conf->conf_aos_to_soaCOM(&conf_soaCOM[index],index);

  //////  ESTRAGGO I MOMENTI E LI TRADUCO     //////////////////////////////////////////////////////
  create_momenta();
  thmatCOM_soa momenta_soaCOM[8];
  for(int index=0;index<8;index++)   gauge_momenta->mom_aos_to_soaCOM(&momenta_soaCOM[index],index);
  cout << "Created momenta\n";
  
  //////  ESTRAGGO IL MULTIFERMIONE PHI E LO TRADUCO    ////////////////////////////////////////////
  create_phi();
  COM_MultiFermion    COMMON_phi;
  fermion_phi->ferm_Multi_to_MultiCOM(&COMMON_phi);
  cout << "Created phi\n";
  
  RationalApprox approx;
  COM_RationalApprox *COM_approx;
  COM_approx = new COM_RationalApprox[1];
  
  cout << "HERE 0 \n";
  //  use_stored=1;//  --> then in the following line the eigenvalues will be computed
  use_stored=1;//  DEBUG CHOICE
  cout << "HERE 1 \n";
  approx.first_inv_approx_coeff();
  cout << "HERE 2 \n";
  use_stored=0;//  --> then in the following the stored eigenvals  will be kept
  cout << "Before Printing Approx \n";
  //  cerr << approx << endl;
  cout << "Before Approx conversion \n";
  convert_RationalApprox_to_COM_RationalApprox(&COM_approx[0],approx);
  cout << "After Approx conversion \n";
  cout << "Rat approx done\n";
  
  //////  CALCOLO IL MULTIFERMIONE CHI CON OPENACC E LO TRADUCO PER DARLO (ANCHE) IN PASTO ALLA PARTE SERIALE   //////////////////
  COM_MultiFermion    COMMON_chi;

  first_inv_approx_calc_openacc(conf_soaCOM,&COMMON_chi,&COMMON_phi,residue_metro,COM_approx); // RESIDUO_METRO
  fermion_chi->ferm_MultiCOM_to_Multi(&COMMON_chi);

  //  first_inv_approx_calc(residue_metro); // RESIDUO_METRO
  //  fermion_chi->ferm_Multi_to_MultiCOM(&COMMON_chi);
  cout << "Created chi\n";
  
  
  //////  PRENDO L'APPROX DA USARE NELLA MD    ///////////////////////////////////////////////
  approx.md_inv_approx_coeff();
  convert_RationalApprox_to_COM_RationalApprox(&COM_approx[0],approx);
  clock_t time_start, time_finish;
  

  calc_plaquette_openacc(conf_soaCOM);
  double ps,pt;
  gauge_conf->calc_plaq(ps,pt);
  cout << "Plaquette CPU BEFORE    =  " << (ps+pt)*0.5 << endl;
  cout << "LINK 0 BEFORE  " << endl << gauge_conf->u_work[0] << endl;
  cout << "MOME 0 BEFORE  " << endl << gauge_momenta->momenta[0] << endl;
  cout << "CHIFERM 0 0 BEFORE   " << fermion_chi->fermion[0][0] << endl;


  //////////////////////////////// MD CPU //////////////////////////////////////////////////////

  time_start=clock();
  multistep_2MN();
  time_finish=clock();
  cout << "Time for Update with CPU = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";

  gauge_conf->calc_plaq_uwork(ps,pt);
  cout << "Plaquette CPU AFTER CPU UPDATE   =  " << (ps+pt)*0.5 << endl;
  cout << "LINK 0 AFTER CPU UPDATE  " << endl << gauge_conf->u_work[0] << endl;
  cout << "MOME 0 AFTER CPU UPDATE  " << endl << gauge_momenta->momenta[0] << endl;
  cout << "CHIFERM 0 0 BEFORE   " << fermion_chi->fermion[0][0] << endl;

  //////////////////////////////// MD ACC //////////////////////////////////////////////////////

  time_start=clock();
  multistep_2MN_ACC(conf_soaCOM,residue_md,COM_approx,&COMMON_chi,momenta_soaCOM);
  time_finish=clock();
  cout << "Time for Update with OPENACC = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";

  Conf *gauge_conf_bis;
  gauge_conf_bis=new Conf(0);

  for(int index=0;index<8;index++)   gauge_conf_bis->conf_soaCOM_to_aos(&conf_soaCOM[index],index);

  /*
  cout << gauge_conf->u_work[0] << endl;  
  gauge_conf->calc_plaq(ps,pt);
  cout << "Plaquette CPU AFTER ACC UPDATE   =  " << (ps+pt)*0.5 << endl;
  cout << "LINK 0 AFTER ACC UPDATE  " << endl << gauge_conf->u_work[0] << endl;
  */

  Su3 su3_check;
  int partial_sum = size;
  int indice;
  double normetta;
  for(int i=0;i<partial_sum;i++){
    normetta = 0.0;
    for(int mu=0;mu<4;mu++){
      indice = i + mu*size;
      su3_check  = gauge_conf_bis->u_work[indice];
      su3_check -= gauge_conf->u_work[indice];
      normetta += su3_check.l2norm2();
    }
    d_vector1[i] = normetta;
  }

  double differ;
  global_sum(d_vector1,partial_sum);
  differ = sqrt(d_vector1[0])*sqrt(1.0/4.0/8.0/((double)(partial_sum)));  // 4 directions, 9 components

  cout << "Delta Updated Confs / d.o.f. = " << differ <<"\n";

  calc_plaquette_openacc(conf_soaCOM);

  delete gauge_conf_bis;
  end();
  return 0;
}

