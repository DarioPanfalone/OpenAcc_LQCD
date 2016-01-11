#ifndef FERMION_PARAMETERS_C_
#define FERMION_PARAMETERS_C_

#include "./fermion_parameters.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/alloc_vars.h"
#include "./markowchain.h"
#include <string.h>
#include <math.h>

#define ALIGN 128


int NDiffFlavs;// set in init.c, from input file
int NPS_tot;
int max_ps;
ferm_param *fermions_parameters;// set in init.c, from input file

int init_ferm_params(ferm_param *fermion_settings){

    int errorstatus = 0;
    
    
  printf("Initializing fermions...\n");
    
  NPS_tot = 0;
  max_ps = 0;

  // calculation of NPS_tot, max_ps,index_of_the_first_ps; 
  for(int i=0;i<NDiffFlavs;i++){
    // compute the total number of ps
    NPS_tot += fermion_settings[i].number_of_ps;
    // compute the max number of ps among the various flavs
    if(fermion_settings[i].number_of_ps>=max_ps) max_ps = fermion_settings[i].number_of_ps;
    // deterime the offset (where does the ps of the flavour i starts?)
    if(i==0){
      fermion_settings[i].index_of_the_first_ps=0;
    }else{
      fermion_settings[i].index_of_the_first_ps = fermion_settings[i-1].index_of_the_first_ps + fermion_settings[i-1].number_of_ps;
    }
  }

  printf("NPS_tot = %d \n",NPS_tot);
  printf("max_ps = %d \n",max_ps);

  // Rational Approximation related stuff
  for(int i=0;i<NDiffFlavs;i++){
    ferm_param *quark = &fermion_settings[i];
    quark->approx_fi_mother.exponent_num =  +quark->degeneracy;
    quark->approx_md_mother.exponent_num =  -quark->degeneracy;
    quark->approx_li_mother.exponent_num =  -quark->degeneracy;

    quark->approx_fi_mother.exponent_den =   quark->number_of_ps*8;
    quark->approx_md_mother.exponent_den =   quark->number_of_ps*4;
    quark->approx_li_mother.exponent_den =   quark->number_of_ps*4;

    quark->approx_fi_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/mkwch_pars.expected_max_eigenvalue;
    quark->approx_md_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/mkwch_pars.expected_max_eigenvalue;
    quark->approx_li_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/mkwch_pars.expected_max_eigenvalue;

    quark->approx_fi_mother.lambda_max =  1.0;
    quark->approx_md_mother.lambda_max =  1.0;
    quark->approx_li_mother.lambda_max =  1.0;

    quark->approx_fi_mother.error =  mkwch_pars.residue_metro/
        pow(mkwch_pars.expected_max_eigenvalue,
                (double) quark->approx_fi_mother.exponent_num/
                quark->approx_fi_mother.exponent_den );
    quark->approx_md_mother.error =  mkwch_pars.residue_md/
        pow(mkwch_pars.expected_max_eigenvalue,
                (double) quark->approx_md_mother.exponent_num/
                quark->approx_md_mother.exponent_den );
        
    quark->approx_li_mother.error =  mkwch_pars.residue_metro/
        pow(mkwch_pars.expected_max_eigenvalue,
                (double) quark->approx_li_mother.exponent_num/
                quark->approx_li_mother.exponent_den );

    quark->approx_fi_mother.gmp_remez_precision = 100;
    quark->approx_md_mother.gmp_remez_precision = 100;
    quark->approx_li_mother.gmp_remez_precision = 100;

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
    quark->approx_fi.gmp_remez_precision =
        quark->approx_fi_mother.gmp_remez_precision;
    quark->approx_md.gmp_remez_precision =   
        quark->approx_md_mother.gmp_remez_precision;
    quark->approx_li.gmp_remez_precision =   
        quark->approx_li_mother.gmp_remez_precision;

    // READ THE RAT APPROXS FROM THE FILES
    errorstatus += rationalapprox_read(&(quark->approx_fi_mother));
    errorstatus += rationalapprox_read(&(quark->approx_md_mother));
    errorstatus += rationalapprox_read(&(quark->approx_li_mother));

  }

  return errorstatus;

}


void init_all_u1_phases(bf_param bfpars, ferm_param *fpar  )
{
   

  for(int i=0;i<NDiffFlavs;i++){

      fpar[i].phases = &u1_back_phases[i*8];





  
  }
 





}


#endif
