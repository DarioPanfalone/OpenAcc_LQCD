#ifndef FERMION_PARAMETERS_H
#define FERMION_PARAMETERS_H

#include "../RationalApprox/rationalapprox.h"
#ifdef __GNUC__
 #define __USE_XOPEN2K
 #include <stdlib.h>
#endif

#define ALIGN 128

//FERMION PARAMETERS
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


int NDiffFlavs;
int NPS_tot;
int max_ps;
ferm_param *fermions_parameters;

#define  APPROX_METRO 19
#define  APPROX_MD 9
#define LAMBDA_MIN_METRO 4.0e-7  // rational approx valid on [lambda_min_metro, 1.0]
#define  LAMBDA_MIN_MD 4.0e-7  // rational approx valid on [lambda_min_metro, 1.0]
#define  GMP_REMEZ_PRECISION 100  // The precision that gmp uses

void init_ferm_params();


#endif

