#ifndef BACKFIELD_H_
#define BACKFIELD_H_

#include "./struct_c_def.h"

// quanti di campo esterno

typedef struct bf_param_t{

 double ex,ey,ez,bx,by,bz;
 // maybe you want to add some other strange things, setting and so on

} bf_param;

extern bf_param backfield_parameters;


void init_backfield(double_soa * tu1_back_field_phases, bf_param bfpars);


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


/********
 * The next function produces a background field that includes all
 * phases
 * 1. The EM field phases (e^iQA, already containing the charge) .
 * 2. The staggered phases
 * 3. The antiperiodic boundary condition in one of the four direction
 * 4. The chemical potential related phases that must be applied in the 
 *    "T" direction 
 ***********************************************************************/
void init_fermion_backfield(bf_param bf_pars, ferm_param *fermion_parameters);


void mult_backfield_by_stag_phases(double_soa * phases);
#endif
