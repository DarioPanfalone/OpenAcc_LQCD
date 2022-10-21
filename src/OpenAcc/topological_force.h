#ifndef TOPOLOGICAL_FORCE_H
#define TOPOLOGICAL_FORCE_H

#include <stdlib.h>
#include <stdio.h>
#include "./struct_c_def.h"

double compute_topodynamical_potential_der(double Q);

void four_leaves(su3_soa * const leaves,__restrict const su3_soa * const u);

double compute_topological_force_internal(__restrict const su3_soa * const u,__restrict su3_soa * const F);

void calc_loc_topo_staples(__restrict su3_soa * u,
													 __restrict su3_soa * const staples);

void topo_staples(__restrict su3_soa * u,__restrict su3_soa * const staples, double norm);

void calc_ipdot_topo(__restrict const su3_soa * const tconf_acc,  
#ifdef STOUT_TOPO
										 __restrict su3_soa * const tstout_conf_acc_arr,
#endif
										 __restrict su3_soa  * taux_conf_acc,
										 __restrict su3_soa * const local_staples,
										 __restrict tamat_soa * const tipdot);
  


#endif
