// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h

#ifndef FERM_MEAS_C 
#define FERM_MEAS_C

#define ALIGN 128

#include <stdio.h>
#include <stdlib.h>
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/inverter_full.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/random_assignement.h"
#include "./ferm_meas.h"
#include "../Include/fermion_parameters.h"
#include "../OpenAcc/fermion_matrix.h"
#include "../OpenAcc/fermionic_utilities.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "../Include/common_defines.h"
#ifdef STOUT_FERMIONS
#include "../OpenAcc/stouting.h"
#endif
#include "../OpenAcc/action.h"

char fermionic_outfilename[50];
char fermionic_outfile_header[100];

// vedi tesi LS F.Negro per ragguagli (Appendici)
void eo_inversion(su3_soa *tconf_acc,
		  ferm_param * tfermions_parameters,
                  double res,
		  vec3_soa * in_e,     // z2 noise
		  vec3_soa * in_o,     // z2 noise
		  vec3_soa * out_e,
		  vec3_soa * out_o,
		  vec3_soa * phi_e,    // parking variable
		  vec3_soa * phi_o,    // parking variable
                  vec3_soa * trialSolution,       // initial vector for the inversion 
		  vec3_soa * tloc_r,    // parking variable for the inverter      
		  vec3_soa * tloc_h,    // parking variable for the inverter
		  vec3_soa * tloc_s,    // parking variable for the inverter
                  vec3_soa * tloc_p){   // parking variable for the inverter




	    acc_Deo(tconf_acc, phi_e, in_o,tfermions_parameters->phases);
	    combine_in1_x_fact1_minus_in2_back_into_in2(in_e, tfermions_parameters->ferm_mass , phi_e);
	    ker_invert_openacc(tconf_acc,tbackfield,tfermions_parameters,
				      out_e,phi_e,res,trialSolution,
				      tloc_r,tloc_h,tloc_s,tloc_p);
	    acc_Doe(tconf_acc, phi_o, out_e,tfermions_parameters->phases);
	    combine_in1_minus_in2_allxfact(in_o,phi_o,(double)1/tfermions_parameters->ferm_mass,out_o);

    
}// end eo_inversion


d_complex chiral_condensate(vec3_soa * rnd_e,
	       vec3_soa * rnd_o,
	       vec3_soa * chi_e,
	       vec3_soa * chi_o){

  return scal_prod_global(rnd_o,chi_o) + scal_prod_global(rnd_e,chi_e);
     
}






void perform_chiral_measures( su3_soa * tconf_acc,
			      double_soa * tbackfield,
			      ferm_param * tfermions_parameters,
			      double res,
			      FILE *out_file){
  vec3_soa * rnd_e,* rnd_o;
  vec3_soa * chi_e,* chi_o;
  vec3_soa * phi_e,* phi_o;
  vec3_soa * trial_sol;

  su3_soa * conf_to_use;

#ifdef STOUT_FERMIONS
    SETREQUESTED(gstout_conf_acc_arr);
    stout_wrapper(tconf_acc ,gstout_conf_acc_arr);
    conf_to_use = &gstout_conf_acc_arr[8*(act_params.stout_steps-1)];
#else
    conf_to_use = tconf_acc;
#endif

  int allocation_check;
  allocation_check =  posix_memalign((void **)&rnd_e, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di rnd_e \n");
  allocation_check =  posix_memalign((void **)&rnd_o, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di rnd_o \n");
  allocation_check =  posix_memalign((void **)&phi_e, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di phi_e \n");
  allocation_check =  posix_memalign((void **)&phi_o, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di phi_o \n");
  allocation_check =  posix_memalign((void **)&chi_e, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di chi_e \n");
  allocation_check =  posix_memalign((void **)&chi_o, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di chi_o \n");
  allocation_check =  posix_memalign((void **)&trial_sol, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di trial_sol \n");


  generate_vec3_soa_z2noise(rnd_e);
  generate_vec3_soa_z2noise(rnd_o);
  generate_vec3_soa_gauss(trial_sol);
  d_complex chircond = 0.0 + 0.0*I;

#pragma acc data create(phi_e[0:1]) create(phi_o[0:1]) create(chi_e[0:1]) create(chi_o[0:1]) copyin(rnd_e[0:1]) copyin(rnd_o[0:1]) copyin(trial_sol[0:1])
  {
    // i fermioni ausiliari kloc_* sono quelli GLOBALI !!!
    eo_inversion(conf_to_use,tbackfield,tfermions_parameters,res,rnd_e,rnd_o,chi_e,chi_o,phi_e,phi_o,trial_sol,kloc_r,kloc_h,kloc_s,kloc_p);
    chircond = chiral_condensate(rnd_e,rnd_o,chi_e,chi_o);
    double factor = tfermions_parameters->degeneracy*0.25/size;
    fprintf(out_file,"%.16lf\t%.16lf\t",creal(chircond)*factor,cimag(chircond)*factor);
  }


  free(rnd_e);
  free(rnd_o);
  free(phi_e);
  free(phi_o);
  free(chi_e);
  free(chi_o);
  free(trial_sol);


}


#endif
