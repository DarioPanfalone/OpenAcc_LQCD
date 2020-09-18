#ifndef TOPO_ACTION_H
#define TOPO_ACTION_H

#include <stdlib.h>
#include "../Meas/gauge_meas.h"
#include "./alloc_vars.h"
#include "./action.h"


#ifndef STOUTING_H
 #include "./stouting.h"
#endif

int load(const char *path, double barr, double width, double* grid, int ngrid);


double topodynamical_pot(double Q);


double topodynamical_pot_der(double Q);


double compute_topo_action(su3_soa * const u);




#endif
