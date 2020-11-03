#ifndef TOPOLOGICAL_ACTION_H
#define TOPOLOGICAL_ACTION_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "./struct_c_def.h"

//read COPIATO DA NISSA													\
SI DEVE CREARE UN FILE CHE CONTENGA UN INSIEME DI VALORI CHE DEFINISCONO IL POTENZIALE


void load_topo(const char *path,const double barr,const double width, double *grid,const int ngrid);


double topodynamical_pot(double Q);


double compute_topo_action(__restrict su3_soa * u
#ifdef STOUT_TOPO
			   ,__restrict su3_soa * const tstout_conf_acc_arr
#endif
			   );




#endif
