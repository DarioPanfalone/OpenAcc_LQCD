#ifndef PARAMETERS_CC_
#define PARAMETERS_CC_

#include "./global_macro.cc"

#define DIM_BLOCK_X 8 // This should divide (nx/2)
#define DIM_BLOCK_Y 8 // This should divide ny
#define DIM_BLOCK_Z 8  // This should divide nz*nt

// lattice dimensions
const int nx=16;
const int ny=16;
const int nz=16;
const int nt=16;



// run parameters
const int no_flavours=2; // number of quark species
const REAL beta=5.55;
const REAL mass=0.02;//0.01335; //0.04303;



// random seed (if 0 time machine is used)
const int rand_seed=1;

// fermion temporal bounday condition   =0 --> antiperiodic, else periodic
const int ferm_temp_bc=1;               


// start parameter
const int start=2;  //=0 ordered start, =1 random start, =2 start from saved conf

// RHMC parameters
const int max_cg=100000;             // maximum number of iteration in CG inverter

const REAL inv_single_double_prec=1.0e-05;//1.0e-08;  // if stopping residue <inv_single_double_prec inverter use double prec, else single
// const REAL inv_single_double_prec=10.0;  // if stopping residue <inv_single_double_prec inverter use double prec, else single

  //Rational approximations for Metropolis
  const int gpu_device_to_use=0;

#endif
