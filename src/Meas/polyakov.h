#ifndef POLYAKOV_H_
#define POLYAKOV_H_

#include "../OpenAcc/struct_c_def.h"

// LOOPS in the various directions
d_complex polyakov_loop0(__restrict const su3_soa * const u);
d_complex polyakov_loop1(__restrict const su3_soa * const u);
d_complex polyakov_loop2(__restrict const su3_soa * const u);
d_complex polyakov_loop3(__restrict const su3_soa * const u);

extern d_complex (*polyakov_loop[4])(__restrict const su3_soa * const u);



#endif
