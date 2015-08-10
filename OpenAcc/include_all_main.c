double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);


#include "openacc.h"
#include "../RationalApprox/rationalapprox.c"
#include "./struct_c_def.c"
#include "./alloc_vars.c"
#include "./fermionic_utilities.c"
#include "./su3_utilities.c"
#include "./random_assignement.c"
#include "./fermion_matrix.c"
#include "./inverter_full.c"
#include "./find_min_max.c"
#include "./inverter_multishift_full.c"
#include "./md_integrator.c"
#include "./md_integrator_soloopenacc.c"

int main(){

  printf("INIZIO DEL PROGRAMMA \n");

  su3_soa  * conf_acc;
  int  allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_conf_acc \n");
  printf("Allocazione della configurazione : OK \n");

  mem_alloc();
  printf("Allocazione della memoria : OK \n");
  initialize_global_variables();
  printf("init vars : OK \n");
  compute_nnp_and_nnm_openacc();
  printf("nn computation : OK \n");
  init_backfield();
  printf("u1_backfield initialization : OK \n");

  // RATIONAL APPROX COEFFS CALCLUATION ///////////////////////////////////////////////////////////
  init_ferm_params();
  printf("ratapprox computation : OK \n");
  /////////////////////////////////////////////////////////////////////////////////////////////////


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
  /////////////////////////////////////////////////////////////////////////////////////////////////


  // INIZIALIZZAZIONE DELLA CONFIGURAZIONE
  generate_Conf_cold(conf_acc);
  /////////////////////////////////////////////////////////////////////////////////////////////////

#pragma acc data   copy(conf_acc[0:8]) copyin(u1_back_field_phases[0:8]) create(ipdot_acc[0:8]) create(aux_conf_acc[0:8]) create(ferm_chi_acc[0:NPS_tot]) create(ferm_phi_acc[0:NPS_tot])  create(ferm_out_acc[0:NPS_tot]) create(ferm_shiftmulti_acc[0:max_ps*max_approx_order]) create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1])  create(k_p_shiftferm[0:max_approx_order]) create(momenta[0:8]) copyin(nnp_openacc) copyin(nnm_openacc) create(local_sums[0:2]) create(d_local_sums[0:1])  copyin(fermions_parameters[0:NDiffFlavs])
    {

    int accettate=0;
    ////////////////   THERMALIZATION   /////////////////////////////////////////////////////////////
    for(int id_iter=0;id_iter<10;id_iter++){
      accettate = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,residue_metro,residue_md,id_iter,accettate,0);
#pragma acc update host(conf_acc[0:8])
    }
    ////////////////   METROTEST   //////////////////////////////////////////////////////////////////
    accettate=0;
    for(int id_iter=0;id_iter<10;id_iter++){
      accettate = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,residue_metro,residue_md,id_iter,accettate,1);
#pragma acc update host(conf_acc[0:8])
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    }// end pragma acc data


  //////  OPENACC CONTEXT CLOSING    //////////////////////////////////////////////////////////////
  SHUTDOWN_ACC_DEVICE(my_device_type);
  /////////////////////////////////////////////////////////////////////////////////////////////////


  free(conf_acc);
  mem_free();


  return 0;
}

