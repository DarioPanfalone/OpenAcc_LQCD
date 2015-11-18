#ifndef COMMON_DEFINES_H_
#define COMMON_DEFINES_H_

#include "../Include/fermion_parameters.h"

//#define BACKFIELD
//#define IMCHEMPOT

// se BACKFIELD o IMCHEMPOT sono definiti allora nell'applicazione della matrice
// di dirac usa la routine che moltiplica anche per la fase opportuna, altrimenti
// in modo hard coded usa l'altra routine moltiplica link per fermione e basta.
// Per farlo viene definita la variabile PHASE_MAT_VEC_MULT o meno.

#ifdef BACKFIELD
  #define PHASE_MAT_VEC_MULT 
#endif

#ifndef PHASE_MAT_VEC_MULT
  #ifdef IMCHEMPOT
    #define PHASE_MAT_VEC_MULT 
  #endif
#endif

#define DIM_BLOCK_X 8 // This should divide (nx/2)
#define DIM_BLOCK_Y 8 // This should divide ny
#define DIM_BLOCK_Z 8  // This should divide nz*nt

// lattice dimensions
#define nx 4
#define ny 4
#define nz 4
#define nt 4


#define sizehh nx*ny*nz*nt/2 

#define ANTIPERIODIC_T_BC  // else periodic time bc are taken

//#define TIMING_ALL // if defined many computation times are printed in the output

#define GAUGE_ACT_TLSM
//#define GAUGE_ACT_WILSON

#define STOUT_FERMIONS
#ifdef STOUT_FERMIONS
#define STOUT_STEPS 2
#endif
#define RHO 0.15

#define ONE_BY_THREE   0.33333333333333333333333
#define ONE_BY_SIX     0.16666666666666666666666
#define beta_by_three  (beta*ONE_BY_THREE)


#ifdef GAUGE_ACT_TLSM
#define GAUGE_ACTION   1     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         (5.0*ONE_BY_THREE)
#define C_ONE          (-0.25*ONE_BY_THREE)
#endif

#ifdef GAUGE_ACT_WILSON
#define GAUGE_ACTION   0     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         1.0
#define C_ONE          0.0
#endif

#define beta 3.55
  // 3.7  //5.35


const int start_opt=0;// 0 --> COLD START; 1 --> START FROM SAVED CONF
int conf_id_iter;
int ITERATIONS=1; // the code will generate new <ITERATIONS> confs, from <conf_id_iter+1> to <conf_id_iter+ITERATIONS>
int therm_ITERATIONS = 30; // the first <therm_ITERATIONS> of the history will be thermalization updates

int save_conf_every=10;

// quanti di campo esterno
const double bx_quantum=0.0;
const double ex_quantum=0.0;
const double by_quantum=0.0;
const double ey_quantum=0.0;
const double bz_quantum=0.0;
const double ez_quantum=1.0;


const char  *nome_file_gauge_output = "gauge_meas.dat";
const char  *nome_file_ferm_output  = "ferm_meas.dat";


#define max_cg 10000
#define no_md 8 // number of MD steps
#define gauge_scale 4  // Update fermions every gauge_scale gauge updates
const double residue_metro=1.0e-8;//-8    // stopping residual for CG
const double residue_md=1.0e-6;//-5    // stopping residual for CG

#endif

