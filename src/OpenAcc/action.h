#ifndef ACTION_H_
#define ACTION_H_

#include "../Include/common_defines.h"

typedef struct action_param_t{

    double beta;
    int stout_steps;
    double stout_rho; // AT PRESENT, ONLY FOR CHECKING PURPOSES

} action_param;



extern action_param act_params;

#define BETA_BY_THREE  (act_params.beta*ONE_BY_THREE)

#endif
