#ifndef FERMION_PARAMETERS_C_
#define FERMION_PARAMETERS_C_

#include "./fermion_parameters.h"

#define ALIGN 128

int NDiffFlavs;
int NPS_tot;
int max_ps;
ferm_param *fermions_parameters;

void init_ferm_params(){



  NDiffFlavs = 3;  // the number of different quark flavours

  int allocation_check; 

  // al posto di 128 c'era ALIGN, solo che qui questa variabile non è ancora definita (viene fatto in struct_c_def)
  allocation_check =  posix_memalign((void **)&fermions_parameters, ALIGN, NDiffFlavs*sizeof(ferm_param));   //  -->  4*size phases (as many as links)
  if(allocation_check != 0)  printf("Errore nella allocazione di fermions_parameters \n");

  ferm_param *up,*down,*strange;
  up = &fermions_parameters[0];
  down = &fermions_parameters[1];
  strange = &fermions_parameters[2];

  up->ferm_charge       = -1.0;   // up    charge
  up->ferm_mass         = 0.00362345;    //0.075;  // up    mass
  up->ferm_im_chem_pot  = 0.0;    // up    chem pot
  up->degeneracy        = 1;      // up    degeneracy
  up->number_of_ps      = 1;      // up    number of pseudo fermions

  down->ferm_charge       = 2.0;    // down  charge
  down->ferm_mass         = 0.00362345;    //0.075;  // down  mass
  down->ferm_im_chem_pot  = 0.0;    // down  chem pot
  down->degeneracy        = 1;      // down  degeneracy
  down->number_of_ps      = 1;      // down  number of pseudo fermions

  strange->ferm_charge       = -1.0;    // strange  charge
  strange->ferm_mass         = 0.102;    //0.075;  // strange  mass
  strange->ferm_im_chem_pot  = 0.0;    // strange  chem pot
  strange->degeneracy        = 1;      // strange  degeneracy
  strange->number_of_ps      = 1;      // strange  number of pseudo fermions

  NPS_tot = 0;
  max_ps = up->number_of_ps;
  for(int i=0;i<NDiffFlavs;i++){
    // compute the total number of ps
    NPS_tot += fermions_parameters[i].number_of_ps;
    // compute the max number of ps among the various flavs
    if(fermions_parameters[i].number_of_ps>=max_ps) max_ps = fermions_parameters[i].number_of_ps;
    // deterime the offset (where does the ps of the flavour i starts?)
    if(i==0){
      fermions_parameters[i].index_of_the_first_ps=0;
    }else{
      fermions_parameters[i].index_of_the_first_ps = fermions_parameters[i-1].index_of_the_first_ps + fermions_parameters[i-1].number_of_ps;
    }
  }

  printf("NPS_tot = %d \n",NPS_tot);
  printf("max_ps = %d \n",max_ps);


  
  for(int i=0;i<NDiffFlavs;i++){
    ferm_param *quark = &fermions_parameters[i];
    quark->approx_fi_mother.exponent_num =  +quark->degeneracy;
    quark->approx_md_mother.exponent_num =  -quark->degeneracy;
    quark->approx_li_mother.exponent_num =  -quark->degeneracy;

    quark->approx_fi_mother.exponent_den =   quark->number_of_ps*8;
    quark->approx_md_mother.exponent_den =   quark->number_of_ps*4;
    quark->approx_li_mother.exponent_den =   quark->number_of_ps*4;

    quark->approx_fi_mother.approx_order =  APPROX_METRO;
    quark->approx_md_mother.approx_order =  APPROX_MD;
    quark->approx_li_mother.approx_order =  APPROX_METRO;

    quark->approx_fi_mother.lambda_min =  LAMBDA_MIN_METRO;
    quark->approx_md_mother.lambda_min =  LAMBDA_MIN_MD;
    quark->approx_li_mother.lambda_min =  LAMBDA_MIN_METRO;

    quark->approx_fi_mother.lambda_max =  1.0;
    quark->approx_md_mother.lambda_max =  1.0;
    quark->approx_li_mother.lambda_max =  1.0;

    quark->approx_fi_mother.lambda_max =  1.0;
    quark->approx_md_mother.lambda_max =  1.0;
    quark->approx_li_mother.lambda_max =  1.0;

    quark->approx_fi_mother.gmp_remez_precision = GMP_REMEZ_PRECISION;
    quark->approx_md_mother.gmp_remez_precision = GMP_REMEZ_PRECISION;
    quark->approx_li_mother.gmp_remez_precision = GMP_REMEZ_PRECISION;

    // copy everything also in the daughter approxs
    quark->approx_fi.exponent_num =   quark->approx_fi_mother.exponent_num;
    quark->approx_md.exponent_num =   quark->approx_md_mother.exponent_num;
    quark->approx_li.exponent_num =   quark->approx_li_mother.exponent_num;
    quark->approx_fi.exponent_den =   quark->approx_fi_mother.exponent_den;
    quark->approx_md.exponent_den =   quark->approx_md_mother.exponent_den;
    quark->approx_li.exponent_den =   quark->approx_li_mother.exponent_den;
    quark->approx_fi.approx_order =   quark->approx_fi_mother.approx_order;
    quark->approx_md.approx_order =   quark->approx_md_mother.approx_order;
    quark->approx_li.approx_order =   quark->approx_li_mother.approx_order;
    quark->approx_fi.gmp_remez_precision =   quark->approx_fi_mother.gmp_remez_precision;
    quark->approx_md.gmp_remez_precision =   quark->approx_md_mother.gmp_remez_precision;
    quark->approx_li.gmp_remez_precision =   quark->approx_li_mother.gmp_remez_precision;

    // READ THE RAT APPROXS FROM THE FILES
    rationalapprox_read(&(quark->approx_fi_mother));
    rationalapprox_read(&(quark->approx_md_mother));
    rationalapprox_read(&(quark->approx_li_mother));


  }

  

}

#endif
