#ifndef TOPO_FORCE_H
#define TOPO_FORCE_H

#include <stdlib.h>
#include <stdio.h>
#include "./struct_c_def.h"


double compute_topodynamical_potential_der(double Q);

void four_leaves(su3_soa * const leaves, su3_soa * const u);

double compute_topological_force_internal(__restrict const su3_soa * const u,__restrict su3_soa * const F);

void calc_loc_topo_staples(__restrict su3_soa * const u, __restrict su3_soa * const staples);




#endif
