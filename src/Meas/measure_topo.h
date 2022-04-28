#ifndef MEASURE_TOPO_H_
#define MEASURE_TOPO_H_

#include "../Include/common_defines.h"

typedef struct meastopo_param_t{

    int meascool;
    char pathcool[20];
    int coolmeasstep;
    int cool_measinterval;
    int cooleach;
    
    int measstout;
    char pathstout[20];
    double measrhostout; //ONLY FOR CHECKING
    int stoutmeasstep;
    int stout_measinterval;
    int stouteach;
    
} meastopo_param;

extern meastopo_param meastopo_params;

#endif
