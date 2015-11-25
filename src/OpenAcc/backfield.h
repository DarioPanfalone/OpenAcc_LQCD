#ifndef BACKFIELD_H_
#define BACKFIELD_H_

#include "./struct_c_def.h"


// quanti di campo esterno
#define  bx_quantum 0.0
#define  ex_quantum 0.0
#define  by_quantum 0.0
#define  ey_quantum 0.0
#define  bz_quantum 0.0
#define  ez_quantum 0.0

void init_backfield(double_soa * tu1_back_field_phases);


#endif
