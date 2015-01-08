
#ifndef GLOBAL_VAR_CC_
#define GLOBAL_VAR_CC_

#include "./global_macro.cc"

int update_iteration;

// eigenvalues
REAL min_stored=0.0;
REAL max_stored=0.0;
int use_stored=1;  // when =0 stored eigenvalues are used

// staggered phases , allocated and set in Geometry/geometry.cc :: init_geo()
int *eta;

// auxiliary global fermions to be used in inverters
class Fermion;
Fermion *loc_r;
Fermion *loc_h;
Fermion *loc_s;
Fermion *loc_p;


// next neighbours, allocated and set in Geometry/geometry.cc :: init_geo()
long int **nnp;    
long int **nnm;   

// configuration
class Conf;
Conf *gauge_conf;

// vectors for global sum (also used in inverters)
double *d_vector1;


//used in the CUDA version, to pack the configuration on the device
float *gauge_field_packed;
int *shift_table;

float *chi_packed;//not actually used but necessary in one of the cudainit functions 


#endif
