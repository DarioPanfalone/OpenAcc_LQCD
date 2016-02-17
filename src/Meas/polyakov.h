#ifndef POLYAKOV_H_
#define POLYAKOV_H_

#include "../OpenAcc/struct_c_def.h"
double polyakov_loop(	__restrict su3_soa * const u,
			__restrict su3_soa * const loopplk	);

#endif
