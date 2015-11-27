#ifndef FERMION_PARAMETERS_C_
#define FERMION_PARAMETERS_C_

#include "./fermion_parameters.h"

#define ALIGN 128

void init_ferm_params(){



  NDiffFlavs = 3;  // the number of different quark flavours

  int allocation_check; 

  // al posto di 128 c'era ALIGN, solo che qui questa variabile non è ancora definita (viene fatto in struct_c_def)
  allocation_check =  posix_memalign((void **)&fermions_parameters, ALIGN, NDiffFlavs*sizeof(ferm_param));   //  -->  4*size phases (as many as links)
  if(allocation_check != 0)  printf("Errore nella allocazione di fermions_parameters \n");

  fermions_parameters[0].ferm_charge       = -1.0;   // up    charge
  fermions_parameters[0].ferm_mass         = 0.0011723;    //0.075;  // up    mass
  fermions_parameters[0].ferm_im_chem_pot  = 0.0;    // up    chem pot
  fermions_parameters[0].degeneracy        = 1;      // up    degeneracy
  fermions_parameters[0].number_of_ps      = 1;      // up    number of pseudo fermions

  fermions_parameters[1].ferm_charge       = 2.0;    // down  charge
  fermions_parameters[1].ferm_mass         = 0.0011723;    //0.075;  // down  mass
  fermions_parameters[1].ferm_im_chem_pot  = 0.0;    // down  chem pot
  fermions_parameters[1].degeneracy        = 1;      // down  degeneracy
  fermions_parameters[1].number_of_ps      = 1;      // down  number of pseudo fermions

  fermions_parameters[2].ferm_charge       = -1.0;    // strange  charge
  fermions_parameters[2].ferm_mass         = 0.033;    //0.075;  // strange  mass
  fermions_parameters[2].ferm_im_chem_pot  = 0.0;    // strange  chem pot
  fermions_parameters[2].degeneracy        = 1;      // strange  degeneracy
  fermions_parameters[2].number_of_ps      = 1;      // strange  number of pseudo fermions

  NPS_tot = 0;
  max_ps = fermions_parameters[0].number_of_ps;
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
    fermions_parameters[i].approx_fi_mother.exponent_num =  +fermions_parameters[i].degeneracy;
    fermions_parameters[i].approx_md_mother.exponent_num =  -fermions_parameters[i].degeneracy;
    fermions_parameters[i].approx_li_mother.exponent_num =  -fermions_parameters[i].degeneracy;

    fermions_parameters[i].approx_fi_mother.exponent_den =   fermions_parameters[i].number_of_ps*8;
    fermions_parameters[i].approx_md_mother.exponent_den =   fermions_parameters[i].number_of_ps*4;
    fermions_parameters[i].approx_li_mother.exponent_den =   fermions_parameters[i].number_of_ps*4;

    fermions_parameters[i].approx_fi_mother.approx_order =  APPROX_METRO;
    fermions_parameters[i].approx_md_mother.approx_order =  APPROX_MD;
    fermions_parameters[i].approx_li_mother.approx_order =  APPROX_METRO;

    fermions_parameters[i].approx_fi_mother.lambda_min =  LAMBDA_MIN_METRO;
    fermions_parameters[i].approx_md_mother.lambda_min =  LAMBDA_MIN_MD;
    fermions_parameters[i].approx_li_mother.lambda_min =  LAMBDA_MIN_METRO;

    fermions_parameters[i].approx_fi_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_md_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_li_mother.lambda_max =  1.0;

    fermions_parameters[i].approx_fi_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_md_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_li_mother.lambda_max =  1.0;

    fermions_parameters[i].approx_fi_mother.gmp_remez_precision = GMP_REMEZ_PRECISION;
    fermions_parameters[i].approx_md_mother.gmp_remez_precision = GMP_REMEZ_PRECISION;
    fermions_parameters[i].approx_li_mother.gmp_remez_precision = GMP_REMEZ_PRECISION;

    // copy everything also in the daughter approxs
    fermions_parameters[i].approx_fi.exponent_num =   fermions_parameters[i].approx_fi_mother.exponent_num;
    fermions_parameters[i].approx_md.exponent_num =   fermions_parameters[i].approx_md_mother.exponent_num;
    fermions_parameters[i].approx_li.exponent_num =   fermions_parameters[i].approx_li_mother.exponent_num;
    fermions_parameters[i].approx_fi.exponent_den =   fermions_parameters[i].approx_fi_mother.exponent_den;
    fermions_parameters[i].approx_md.exponent_den =   fermions_parameters[i].approx_md_mother.exponent_den;
    fermions_parameters[i].approx_li.exponent_den =   fermions_parameters[i].approx_li_mother.exponent_den;
    fermions_parameters[i].approx_fi.approx_order =   fermions_parameters[i].approx_fi_mother.approx_order;
    fermions_parameters[i].approx_md.approx_order =   fermions_parameters[i].approx_md_mother.approx_order;
    fermions_parameters[i].approx_li.approx_order =   fermions_parameters[i].approx_li_mother.approx_order;
    fermions_parameters[i].approx_fi.gmp_remez_precision =   fermions_parameters[i].approx_fi_mother.gmp_remez_precision;
    fermions_parameters[i].approx_md.gmp_remez_precision =   fermions_parameters[i].approx_md_mother.gmp_remez_precision;
    fermions_parameters[i].approx_li.gmp_remez_precision =   fermions_parameters[i].approx_li_mother.gmp_remez_precision;

    // READ THE RAT APPROXS FROM THE FILES
    rationalapprox_read(&(fermions_parameters[i].approx_fi_mother));
    rationalapprox_read(&(fermions_parameters[i].approx_md_mother));
    rationalapprox_read(&(fermions_parameters[i].approx_li_mother));


  }
}

#endif
