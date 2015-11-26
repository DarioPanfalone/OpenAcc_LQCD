#ifndef FERMION_PARAMETERS_H
#define FERMION_PARAMETERS_H

#include "../RationalApprox/rationalapprox.h"
#ifdef __GNUC__
 #define __USE_XOPEN2K
 #include <stdlib.h>
#endif

#ifndef FERMION_PARAMETERS_C_
#define EXT_TO_FERMION_PARAMETERS extern
#else
#define EXT_TO_FERMION_PARAMETERS 
#endif

typedef struct ferm_param_t{
  double ferm_charge;
  double ferm_mass;
  double ferm_im_chem_pot;
  int degeneracy;
  int number_of_ps;
  int index_of_the_first_ps;
  RationalApprox approx_fi_mother; // first inv   -> mother
  RationalApprox approx_md_mother; // md approx   -> mother
  RationalApprox approx_li_mother; // last inv    -> mother
  RationalApprox approx_fi;        // first inv
  RationalApprox approx_md;        // md approx
  RationalApprox approx_li;        // last inv
} ferm_param;

EXT_TO_FERMION_PARAMETERS int NDiffFlavs;
EXT_TO_FERMION_PARAMETERS int NPS_tot;
EXT_TO_FERMION_PARAMETERS int max_ps;
EXT_TO_FERMION_PARAMETERS ferm_param *fermions_parameters;




//FERMION PARAMETERS

#define  APPROX_METRO 19
#define  APPROX_MD 9
#define LAMBDA_MIN_METRO 4.0e-7  // rational approx valid on [lambda_min_metro, 1.0]
#define  LAMBDA_MIN_MD 4.0e-7  // rational approx valid on [lambda_min_metro, 1.0]
#define  GMP_REMEZ_PRECISION 100  // The precision that gmp uses

void init_ferm_params();


#endif

