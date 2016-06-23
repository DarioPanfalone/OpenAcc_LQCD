#ifndef MAGNETIC_SUSCEPTIBILITY_UTILITIES_C_
#define MAGNETIC_SUSCEPTIBILITY_UTILITIES_C_

#include "./magnetic_susceptibility_utilities.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/backfield.h"

void idphase_dbz(double_soa * idphi_dbz_re, double_soa * idphi_dbz_im, // dphi_dbz_re = 0 
        ferm_param *tfermion_parameters){

    bf_param zero_point = (bf_param){0,0,0,0,0,0};
    bf_param bz1_point  = zero_point;
    bz1_point.bz += 1;
   
                                      
    calc_u1_phases(idphi_dbz_im,bz1_point ,0,tfermion_parameters->ferm_charge);// mu = 0 
    // using idphi_dbz_re as temp
    calc_u1_phases(idphi_dbz_re,zero_point,0,tfermion_parameters->ferm_charge);// mu = 0 
   
    phase_diff_in_place(idphi_dbz_im,idphi_dbz_re); // in place 
    set_double_soa_to_zero(idphi_dbz_re); 

}


// FINITE differences in e^(i phase_1) - e^(i phase_0)
void delta_u_delta_bz(double_soa * delta_u_re,double_soa * delta_u_im,
        double_soa *tp1, double_soa *tp2,
        ferm_param *tfermion_parameters, bf_param ext_field){

    bf_param bz1_point  = ext_field;
    bz1_point.bz += 1;
    
    calc_u1_phases(tp1,bz1_point,0,tfermion_parameters->ferm_charge);// mu = 0 
    calc_u1_phases(tp2,ext_field,0,tfermion_parameters->ferm_charge);// mu = 0 

    u1_diff( delta_u_re, delta_u_im, tp1, tp2);
    
}




#endif
