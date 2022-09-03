#ifndef ACTION_H_
#define ACTION_H_

#include "../Include/common_defines.h"

typedef struct action_param_t{

	double beta;
	int stout_steps;
	double stout_rho;

	int topo_action;
	double barrier;
	double width;
	char topo_file_path[20];
	int topo_stout_steps;
	double topo_rho;
} action_param;

extern double grid[];

extern action_param act_params;
double gl_stout_rho, gl_topo_rho;

#define BETA_BY_THREE  (act_params.beta*ONE_BY_THREE)

#endif
