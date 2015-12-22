#ifndef COMMON_DEFINES_H_
#define COMMON_DEFINES_H_


//#define BACKFIELD
//#define IMCHEMPOT

// se BACKFIELD o IMCHEMPOT sono definiti allora nell'applicazione della matrice
// di dirac usa la routine che moltiplica anche per la fase opportuna, altrimenti
// in modo hard coded usa l'altra routine moltiplica link per fermione e basta.
// Per farlo viene definita la variabile PHASE_MAT_VEC_MULT o meno.

#if defined(BACKFIELD) || defined (IMCHEMPOT)
  #define PHASE_MAT_VEC_MULT 
#endif

#define DIM_BLOCK_X 8 // This should divide (nx/2)
#define DIM_BLOCK_Y 8 // This should divide ny
#define DIM_BLOCK_Z 8  // This should divide nz*nt

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


extern int verbosity_lv;


#endif

