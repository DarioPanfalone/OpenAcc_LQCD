double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);


#include "openacc.h"
#include "./struct_c_def.c"
#include "./rat_approx_rescale.c"
#include "./fermionic_utilities.c"
#include "./su3_utilities.c"
#include "./random_assignement.c"
#include "./fermion_matrix.c"
#include "./inverter_full.c"
#include "./find_min_max.c"
#include "./inverter_multishift_full.c"
#include "./md_integrator.c"
/*
#include "./md_integrator_soloopenacc.c"


int main(){

  COM_RationalApprox *approx_mother1;
  COM_RationalApprox *approx_mother2;
  COM_RationalApprox *approx_mother3;
  su3_soa  * conf_acc;

  initialize_global_variables();
  compute_nnp_and_nnm_openacc();
  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  posix_memalign((void **)&approx1, ALIGN, sizeof(COM_RationalApprox));
  posix_memalign((void **)&approx2, ALIGN, sizeof(COM_RationalApprox));
  posix_memalign((void **)&approx3, ALIGN, sizeof(COM_RationalApprox));





  return 0;
}
*/
