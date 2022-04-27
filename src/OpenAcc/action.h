#ifndef ACTION_H_
#define ACTION_H_

#include "../Include/common_defines.h"

typedef struct action_param_t{

    double beta;
    int stout_steps;
    double stout_rho; // AT PRESENT, ONLY FOR CHECKING PURPOSES

	int topo_action;
	double barrier;
	double width; // MAYBE IT GIVE REDUNDANT INFO: width = 2*barrier/(tot lines in topo_file)
	char topo_file_path[20];
	int topo_stout_steps;
	double topo_rho;
} action_param;

extern double grid[];

extern action_param act_params;
double gl_stout_rho, gl_topo_rho;

#define BETA_BY_THREE  (act_params.beta*ONE_BY_THREE)

#endif
