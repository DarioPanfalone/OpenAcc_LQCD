#ifndef BACKFIELD_H_
#define BACKFIELD_H_

#include "./struct_c_def.h"

// quanti di campo esterno
extern double  bx_quantum;
extern double  ex_quantum;
extern double  by_quantum;
extern double  ey_quantum;
extern double  bz_quantum;
extern double  ez_quantum;

void set_field_quanta();

void init_backfield(double_soa * tu1_back_field_phases);


#endif
