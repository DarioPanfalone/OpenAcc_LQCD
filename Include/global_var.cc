
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

// momenta derivative times complex_I
class Ipdot;
Ipdot *gauge_ipdot;

// momenta 
class Momenta;
Momenta *gauge_momenta;


// normalized coefficients for rational approximations
class RationalApprox;
RationalApprox *first_inv_approx_norm_coeff;
RationalApprox *md_inv_approx_norm_coeff;
RationalApprox *last_inv_approx_norm_coeff;
RationalApprox *meas_inv_coeff;

// auxiliary global fermions to be used in inverters
class Fermion;
Fermion *loc_r;
Fermion *loc_h;
Fermion *loc_s;
Fermion *loc_p;

class MultiFermion;
MultiFermion *fermion_phi;
MultiFermion *fermion_chi;

// to be used in shifted inverter and fermion force calculation
class ShiftFermion;
ShiftFermion *p_shiftferm;
class ShiftMultiFermion;
ShiftMultiFermion *fermion_shiftmulti;


// next neighbours, allocated and set in Geometry/geometry.cc :: init_geo()
long int **nnp;    
long int **nnm;   

// configuration
class Conf;
Conf *gauge_conf;

// staples
class Staples;
Staples *gauge_staples;


// vectors for global sum (also used in inverters)
double *d_vector1;
double *d_vector2;


//used in the CUDA version, to pack the configuration on the device
float *gauge_field_packed;
int *shift_table;

float *chi_packed;//not actually used but necessary in one of the cudainit functions 


#endif
