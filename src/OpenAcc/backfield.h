#ifndef BACKFIELD_H_
#define BACKFIELD_H_

#include "./struct_c_def.h"
#include "./backfield_parameters.h"

// phases, unbounded, without multiplication by 2pi 
void calc_u1_phases_unb_no2pi(double_soa * phases,bf_param bf_pars,
        double im_chem_pot, double ferm_charge);

// reduces phases to the -1 / +1 range
void rebound_u1_phases(double_soa * phases);

void mult_u1_phases(double_soa * phases,double factor);

//just wrapper of the three preceeding functions
void calc_u1_phases(double_soa * phases, bf_param bfpars, 
        double chpot,double charge);

void phase_diff_in_place(double_soa * inout, double_soa * sottraendo);

void set_double_soa_to_zero(double_soa * p);

void u1_diff(double_soa * out_re, double_soa * out_im,
        double_soa * phase_p,double_soa * phase_m );

static inline int KSphaseX(int x,int y,int z,int t){return 1;};
static inline int KSphaseY(int x,int y,int z,int t){
    return 1 - ( 2*(x & 0x1) );

}
static inline int KSphaseZ(int x,int y,int z,int t){
    return 1 - ( 2*((x+y) & 0x1) );
}
static inline int KSphaseT(int x,int y,int z,int t){ 
    return 1 - ( 2*((x+y+z) & 0x1) );
}








#endif
