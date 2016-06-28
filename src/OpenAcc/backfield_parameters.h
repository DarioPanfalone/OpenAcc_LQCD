#ifndef BACKFIELD_PARAMETERS_H
#define BACKFIELD_PARAMETERS_H

// quanti di campo esterno

typedef struct bf_param_t{

 double ex,ey,ez,bx,by,bz;
 // maybe you want to add some other strange things, setting and so on

} bf_param;

extern bf_param backfield_parameters;

#endif

