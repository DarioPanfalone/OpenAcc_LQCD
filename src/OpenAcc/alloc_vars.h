#ifndef ALLOC_VARS_H_
#define ALLOC_VARS_H_

#include "../Include/common_defines.h"
#include "./struct_c_def.h"

// used in the dynamical allocation of structures

#if !defined(ALLOC_VARS_C_) && !defined(ONE_FILE_COMPILATION)
#define EXT_TO_ALLOC_VARS extern
#else
#define EXT_TO_ALLOC_VARS 
#endif

EXT_TO_ALLOC_VARS su3_soa  * conf_acc_bkp; // the old stored conf that will be recovered if the metro test fails.
EXT_TO_ALLOC_VARS su3_soa  * aux_conf_acc; // auxiliary 
EXT_TO_ALLOC_VARS su3_soa  * auxbis_conf_acc; // auxiliary 
EXT_TO_ALLOC_VARS double_soa * u1_back_field_phases; // BACKRGROUND EM FIELD
EXT_TO_ALLOC_VARS thmat_soa * momenta;// GAUGE FIELD EVOLUTION
EXT_TO_ALLOC_VARS tamat_soa * ipdot_acc;// GAUGE FIELD EVOLUTION
EXT_TO_ALLOC_VARS su3_soa * gconf_as_fermionmatrix; // conf to use in either cases in fermion related computation (with or without stouting)

// STOUTING 
#ifdef STOUT_FERMIONS
EXT_TO_ALLOC_VARS su3_soa * gstout_conf_acc; // max stouted conf, just pointer
EXT_TO_ALLOC_VARS su3_soa * gstout_conf_acc_arr; // all stouting steps except the zeroth
EXT_TO_ALLOC_VARS su3_soa * glocal_staples;
EXT_TO_ALLOC_VARS tamat_soa * gipdot;
EXT_TO_ALLOC_VARS tamat_soa * aux_ta; // aggiunta per il calcolo della forza stoutata
EXT_TO_ALLOC_VARS thmat_soa * aux_th; // aggiunta per il calcolo della forza stoutata
#endif

// FERMIONS

EXT_TO_ALLOC_VARS vec3_soa * ferm_chi_acc; // questo e' il chi [NPS_tot]
EXT_TO_ALLOC_VARS vec3_soa * ferm_phi_acc; // questo e' il phi [NPS_tot]
EXT_TO_ALLOC_VARS vec3_soa * ferm_out_acc; // questo e' uno ausiliario [NPS_tot]
EXT_TO_ALLOC_VARS vec3_soa * ferm_shiftmulti_acc; // ausiliario per l'invertitore multishift [max_ps*MAX_APPROX_ORDER]
EXT_TO_ALLOC_VARS vec3_soa * kloc_r;  // vettore ausiliario
EXT_TO_ALLOC_VARS vec3_soa * kloc_h;  // vettore ausiliario
EXT_TO_ALLOC_VARS vec3_soa * kloc_s;  // vettore ausiliario
EXT_TO_ALLOC_VARS vec3_soa * kloc_p;  // vettore ausiliario
EXT_TO_ALLOC_VARS vec3_soa * k_p_shiftferm; // ausiliario [max_nshift=MAX_APPROX_ORDER]


// LOCAL SUMS
EXT_TO_ALLOC_VARS dcomplex_soa * local_sums;
EXT_TO_ALLOC_VARS double_soa * d_local_sums;

void mem_alloc();

void mem_free();

#endif
