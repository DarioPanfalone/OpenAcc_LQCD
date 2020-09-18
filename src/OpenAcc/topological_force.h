#ifndef TOPO_FORCE_H
#define TOPO_FORCE_H

#include <stdlib.h>
#include "./alloc_vars.h"
#include "./action.h"



double compute_topodynamical_potential_der(__restrict const su3_soa * u, double Q);

double compute_topological_force_internal(__restrict const su3_soa * const u,__restrict su3_soa * const F);

double compute_topo_force(su3_soa * const u, const double Q);


#endif
