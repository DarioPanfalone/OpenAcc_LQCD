#ifndef BACKFIELD_PARAMETERS_H
#define BACKFIELD_PARAMETERS_H

// external field quanta
typedef struct bf_param_t{

	double ex,ey,ez,bx,by,bz;

} bf_param;

extern bf_param backfield_parameters;

#endif

