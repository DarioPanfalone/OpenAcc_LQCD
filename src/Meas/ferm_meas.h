// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h
#ifndef FERM_MEAS_H 
#define FERM_MEAS_H

#include <stdio.h>
#include "../OpenAcc/struct_c_def.h"

extern char fermionic_outfilename[50];
extern char fermionic_outfile_header[100];

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
                  vec3_soa * tloc_p);

d_complex chiral_condensate(vec3_soa * rnd_e, vec3_soa * rnd_o,
	       vec3_soa * chi_e, vec3_soa * chi_o);


void perform_chiral_measures( su3_soa * tconf_acc,
			      ferm_param * tfermions_parameters, double res,
			      FILE *out_file);


#endif
