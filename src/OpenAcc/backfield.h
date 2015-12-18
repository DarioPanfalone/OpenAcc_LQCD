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


#endif
