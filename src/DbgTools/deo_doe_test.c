#include <stdio.h>
#include <stdlib.h>

#include "openacc.h"
#include "../Include/common_defines.h"
#include "../OpenAcc/fermion_matrix.c"
#include "../OpenAcc/alloc_vars.c"
#include "../OpenAcc/random_assignement.c"
#include "../OpenAcc/dbgtools.c"

//double casuale(void);


int main(){


  initrand(111);
  fflush(stdout);
  printf("INIZIO DEL PROGRAMMA \n");
  su3_soa  * conf_acc;
  int  allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di conf_acc \n");
  conf_acc->status = IN_USE;
  printf("Allocazione della configurazione : OK \n");

  // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
  init_ferm_params();



  mem_alloc();
  printf("Allocazione della memoria : OK \n");
  initialize_global_variables();
  printf("init vars : OK \n");
  compute_nnp_and_nnm_openacc();
  printf("nn computation : OK \n");
#ifdef BACKFIELD
  init_backfield();
  print_double_soa(u1_back_field_phases,"backfield");
  printf("u1_backfield initialization : OK \n");
#endif

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
  printf("Device Selected : OK \n");


  // init conf
  //generate_Conf_cold(conf_acc,0.05);
  //print_su3_soa(conf_acc, "conf_acc");
  //printf("Cold Gauge Conf Generated : OK \n");
  read_su3_soa(conf_acc,"indaddr_conf_acc");
  printf("Cold Gauge Conf READ : OK \n");

  conf_id_iter=0;


  // init fermion
  //generate_vec3_soa_gauss(ferm_chi_acc);
  //print_vec3_soa(ferm_chi_acc,"ferm_chi_acc");
  read_vec3_soa(ferm_chi_acc,"indaddr_ferm_chi_acc" );

#pragma acc data   copy(conf_acc[0:8]) copy(ferm_chi_acc) copyout(ferm_phi_acc)
  {
  acc_Doe(conf_acc, ferm_phi_acc, ferm_chi_acc, fermions_parameters,u1_back_field_phases ) ;


  print_vec3_soa(ferm_phi_acc,"ferm_phi_acc");
  }

  SHUTDOWN_ACC_DEVICE(my_device_type);
  mem_free();

  return 0; 




}
