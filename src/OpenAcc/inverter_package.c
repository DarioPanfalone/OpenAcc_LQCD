#ifndef INVERTER_PACKAGE_C_
#define INVERTER_PACKAGE_C_

#include "./sp_struct_c_def.h"
#include "./struct_c_def.h"
#include "./inverter_package.h"



void setup_inverter_package_dp(inverter_package * ip, 
        su3_soa * u, 
        vec3_soa * ferm_shift_temp,
        int nshift,
        vec3_soa * loc_r,
        vec3_soa * loc_h,
        vec3_soa * loc_s,
        vec3_soa * loc_p){

    ip->u = (__restrict) u;
    ip->ferm_shift_temp = ferm_shift_temp;
    ip->nshifts = nshifts;
    ip->loc_r = (__restrict) loc_r;
    ip->loc_h = (__restrict) loc_h;
    ip->loc_s = (__restrict) loc_s;
    ip->loc_p = (__restrict) loc_p;
    
}

void setup_inverter_package_sp(inverter_package * ip, 
        su3_soa_f * u_f, 
        vec3_soa_f * ferm_shift_temp_f,
        int nshift,
        vec3_soa_f * loc_r_f,
        vec3_soa_f * loc_h_f,
        vec3_soa_f * loc_s_f,
        vec3_soa_f * loc_p_f,
        vec3_soa_f * out_f)
{

    ip->u_f = (__restrict) u_f;
    ip->ferm_shift_temp_f = ferm_shift_temp_f;
    ip->nshifts = nshifts;
    ip->loc_r_f = (__restrict) loc_r_f;
    ip->loc_h_f = (__restrict) loc_h_f;
    ip->loc_s_f = (__restrict) loc_s_f;
    ip->loc_p_f = (__restrict) loc_p_f;
    ip->out_f   = (__restrict)   out_f;
    
}






#endif
