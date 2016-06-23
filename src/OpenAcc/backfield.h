#ifndef BACKFIELD_H_
#define BACKFIELD_H_

#include "./struct_c_def.h"

// quanti di campo esterno

typedef struct bf_param_t{

 double ex,ey,ez,bx,by,bz;
 // maybe you want to add some other strange things, setting and so on

} bf_param;

extern bf_param backfield_parameters;


void calc_u1_phases(double_soa * phases, bf_param bfpars, 
        double chpot,double charge);

void phase_diff_in_place(double_soa * inout, double_soa * sottraendo);

void set_double_soa_to_zero(double_soa * p);

void u1_diff(double_soa * out_re, double_soa * out_im,
        double_soa * phase_p,double_soa * phase_m );

inline int KSphaseX(int x,int y,int z,int t){return 1;};
inline int KSphaseY(int x,int y,int z,int t){
    return 1 - ( 2*(x & 0x1) );

}
inline int KSphaseZ(int x,int y,int z,int t){
    return 1 - ( 2*((x+y) & 0x1) );
}
inline int KSphaseT(int x,int y,int z,int t){ 
    return 1 - ( 2*((x+y+z) & 0x1) );
}








#endif
