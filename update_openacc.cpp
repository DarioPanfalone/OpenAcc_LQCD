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
  
  //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////////////////////

  // NVIDIA GPUs
  acc_device_t my_device_type = acc_device_nvidia;

  // AMD GPUs
  // acc_device_t my_device_type = acc_device_radeon;

  // Intel XeonPhi
  //acc_device_t my_device_type = acc_device_xeonphi;

  // Select device ID
  int dev_index = 0;

  SELECT_INIT_ACC_DEVICE(my_device_type, dev_index);
  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////// CONFRONTO DINAMICA MOLECOLARE    /////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  su3COM_soa conf_soaCOM[8];

  RationalApprox approx1;
  COM_RationalApprox *COM_approx1;
  COM_approx1 = new COM_RationalApprox[1];
  RationalApprox approx2;
  COM_RationalApprox *COM_approx2;
  COM_approx2 = new COM_RationalApprox[1];
  RationalApprox approx3;
  COM_RationalApprox *COM_approx3;
  COM_approx3 = new COM_RationalApprox[1];
  COM_RationalApprox *COM_approx_mother1;
  COM_approx_mother1 = new COM_RationalApprox[1];
  COM_RationalApprox *COM_approx_mother2;
  COM_approx_mother2 = new COM_RationalApprox[1];
  COM_RationalApprox *COM_approx_mother3;
  COM_approx_mother3 = new COM_RationalApprox[1];
  clock_t time_start, time_finish;
  double ps,pt;
  int r,i;
  Su3 auxm,auxm_bis;
  double A,B,RC,IC,RD,ID,RE,IE;

  Su3 su3_check;
  int partial_sum = size;
  int indice;
  double normetta;
  double differ;

  convert_RationalApprox_to_COM_RationalApprox(&COM_approx_mother1[0],*first_inv_approx_norm_coeff);
  convert_RationalApprox_to_COM_RationalApprox(&COM_approx_mother2[0],*md_inv_approx_norm_coeff);
  convert_RationalApprox_to_COM_RationalApprox(&COM_approx_mother3[0],*last_inv_approx_norm_coeff);

  //////  TRADUCO LA CONFIGURAZIONE    //////////////////////////////////////////////////////
  for(int index=0;index<8;index++)   gauge_conf->conf_aos_to_soaCOM(&conf_soaCOM[index],index);

  ////////////////   THERMALIZATION   ////////////////////////////////////////////////////////////
  for(int id_iter=0;id_iter<1;id_iter++){
  time_start=clock();
  THERM_UPDATE_ACC_UNOSTEP_NOMETRO(conf_soaCOM,residue_metro,residue_md,COM_approx_mother1,COM_approx_mother2,COM_approx_mother3,id_iter);
  //  UPDATE_ACC_UNOSTEP(conf_soaCOM,residue_metro,residue_md,COM_approx_mother1,COM_approx_mother2,COM_approx_mother3,id_iter);
  time_finish=clock();
  cout << "Id_iter_therm (NEW ROUTINE) = " << id_iter << "   " << "Time for Update nometro with OPENACC = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
  }

  ////////////////   METROTEST   ////////////////////////////////////////////////////////////
  int accettate=0;
  for(int id_iter=0;id_iter<0;id_iter++){
    time_start=clock();
    accettate = UPDATE_ACC_UNOSTEP_METRO(conf_soaCOM,residue_metro,residue_md,COM_approx_mother1,COM_approx_mother2,COM_approx_mother3,id_iter,accettate);
    time_finish=clock();
    cout << "Id_iter_metro = " << id_iter << "   " << "Time for Update simetro with OPENACC = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
  }

  end();

  //////  OPENACC CONTEXT SHUTDOWN   //////////////////////////////////////////////////////

  SHUTDOWN_ACC_DEVICE(my_device_type);  

  
  return 0;

}

