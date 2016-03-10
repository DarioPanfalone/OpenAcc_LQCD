#ifndef POLYAKOV_H_
#define POLYAKOV_H_

#include "../OpenAcc/struct_c_def.h"

// LOOPS in the various directions
extern d_complex (*polyakov_loop[4])(__restrict const su3_soa * const u);

#endif
