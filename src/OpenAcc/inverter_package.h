#ifndef INVERTER_PACKAGE_H_
#define INVERTER_PACKAGE_H_

#include "./struct_c_def.h"
#include  "./sp_struct_c_def.h"


typedef struct inverter_package_t{

    const __restrict su3_soa * u;
    const __restrict su3_soa_f * u_f;
    __restrict vec3_soa_f * ferm_shift_temp_f;  //parking
    int nshifts;
    __restrict vec3_soa * loc_r;            //parking
    __restrict vec3_soa * loc_h;            //parking
    __restrict vec3_soa * loc_s;            //parking
    __restrict vec3_soa * loc_p;            //parking
    __restrict vec3_soa_f * loc_r_f;            //parking
    __restrict vec3_soa_f * loc_h_f;            //parking
    __restrict vec3_soa_f * loc_s_f;            //parking
    __restrict vec3_soa_f * loc_p_f;            //parking

} inverter_package;

void setup_inverter_package_dp(inverter_package * ip, 
        su3_soa * u, 
        vec3_soa * ferm_shift_temp,
        int nshift,
        vec3_soa * loc_r,
        vec3_soa * loc_h,
        vec3_soa * loc_s,
        vec3_soa * loc_p);


void setup_inverter_package_sp(inverter_package * ip, 
        su3_soa_f * u_f, 
        vec3_soa_f * ferm_shift_temp_f,
        int nshift,
        vec3_soa_f * loc_r_f,
        vec3_soa_f * loc_h_f,
        vec3_soa_f * loc_s_f,
        vec3_soa_f * loc_p_f);





#endif
