double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);


#include "openacc.h"
#include "./struct_c_def.c"
#include "./alloc_vars.c"
#include "../RationalApprox/rationalapprox.c"
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

  su3_soa  * conf_acc;
  int  allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_conf_acc \n");
  mem_alloc();
  printf("Allocazione della memoria : OK \n");
  initialize_global_variables();
  printf("init vars : OK \n");
  compute_nnp_and_nnm_openacc();
  printf("nn computation : OK \n");
  init_backfield();
  printf("u1_backfield initialization : OK \n");

  // RATIONAL APPROX COEFFS CALCLUATION ///////////////////////////////////////////////////////////
  // la 1 e' la "first_inv_approx_norm_coeff"
  // la 2 e' la "md_inv_approx_norm_coeff"
  // la 3 e' la "last_inv_approx_norm_coeff"
  rationalapprox_read_soloopenacc(&approx_mother1[0],&approx_mother2[0],&approx_mother3[0]);
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

#pragma acc data copy(conf_acc[0:8]) create(momenta[0:8]) create(aux_conf_acc[0:8]) create(ferm_chi_acc[0:1]) create(ferm_phi_acc[0:1]) create(approx1[0:1]) create(approx2[0:1])  create(approx3[0:1])  create(ferm_out_acc[0:1])  create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1]) create(k_p_shiftferm[0:1]) create(ferm_shiftmulti_acc[0:1]) create(ipdot_acc[0:8])  copyin(nnp_openacc) copyin(nnm_openacc) create(local_sums[0:2]) create(d_local_sums[0:1]) copy(approx_mother1[0:1])  copy(approx_mother2[0:1])  copy(approx_mother3[0:1])  create(minmaxeig[0:2]) copyin(u1_back_field_phases[0:8])
    {

    ////////////////   THERMALIZATION   /////////////////////////////////////////////////////////////
    for(int id_iter=0;id_iter<10;id_iter++){
      THERM_UPDATE_SOLOACC_NOMETRO(conf_acc,residue_metro,residue_md,id_iter);
#pragma acc update host(conf_acc[0:8])
    }
    ////////////////   METROTEST   //////////////////////////////////////////////////////////////////
    int accettate=0;
    for(int id_iter=0;id_iter<10;id_iter++){
      accettate = UPDATE_SOLOACC_UNOSTEP_METRO(conf_acc,residue_metro,residue_md,id_iter,accettate);
#pragma acc update host(conf_acc[0:8])
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    }


  //////  OPENACC CONTEXT CLOSING    //////////////////////////////////////////////////////////////
  SHUTDOWN_ACC_DEVICE(my_device_type);
  /////////////////////////////////////////////////////////////////////////////////////////////////


  free(conf_acc);
  mem_free();


  return 0;
}

