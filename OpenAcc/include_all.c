double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);


#include "openacc.h"
#include "./struct_c_def.c"
#include "./rat_approx_rescale.c"
#include "./fermionic_utilities.c"
//#include "../Rand/random.c"
#include "./su3_utilities.c"
#include "./random_assignement.c"
#include "./fermion_matrix.c"
#include "./inverter_full.c"
#include "./find_min_max.c"
#include "./inverter_multishift_full.c"
#include "./md_integrator.c"


