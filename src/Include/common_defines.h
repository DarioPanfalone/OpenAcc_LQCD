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

//keep these commented if you are using M100
/*
#define conj __builtin_conj
#define creal __builtin_creal
#define cimag __builtin_cimag
#define conjf __builtin_conjf
#define crealf __builtin_crealf
#define cimagf __builtin_cimagf
*/

#define DIM_BLOCK_X 8 // This should divide (nx/2)
#define DIM_BLOCK_Y 8 // This should divide ny
#define DIM_BLOCK_Z 8  // This should divide nz*nt

//#define TIMING_ALL // if defined many computation times are printed in the output

//**************************************//
//ACHTUNG !!!!!! comment the one you don't wanna use and decoment the other one

//#define GAUGE_ACT_TLSM
#define GAUGE_ACT_WILSON

//**************************************//

#define STOUT_FERMIONS
#ifdef STOUT_FERMIONS
#define STOUT_STEPS 2
#endif

#define STOUT_TOPO

#define TOPO_RHO 0.1
#define FERMION_RHO 0.15

#define TOPO_RHOF 0.1f
#define FERMION_RHOF 0.15f

extern int TOPO_GLOBAL_DONT_TOUCH;
#define RHO (TOPO_GLOBAL_DONT_TOUCH == 1? TOPO_RHO : FERMION_RHO)
#define RHOF (TOPO_GLOBAL_DONT_TOUCH == 1? TOPO_RHOF : FERMION_RHOF)

#define ONE_BY_THREE   0.33333333333333333333333
#define ONE_BY_THREEF   0.33333333333333333333333f
#define ONE_BY_SIX     0.16666666666666666666666
#define ONE_BY_SIXF     0.16666666666666666666666f


#ifdef GAUGE_ACT_TLSM
#define GAUGE_ACTION   1     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         (5.0*ONE_BY_THREE)
#define C_ZEROF         (5.0f*ONE_BY_THREEF)
#define C_ONE          (-0.25*ONE_BY_THREE)
#define C_ONEF          (-0.25f*ONE_BY_THREEF)
#endif

#ifdef GAUGE_ACT_WILSON
#define GAUGE_ACTION   0     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         1.0
#define C_ZEROF         1.0f
#define C_ONE          0.0
#define C_ONEF          0.0f
#endif


extern int verbosity_lv;


// macro to check the outcome of 'fscanf', reporting the file and line where the problem
// happened.
 
#define CHECKREAD(expr,should_read) \
{int read = expr ;if(read != should_read) { \
        printf("%s:%d, Error, not read expected number of entries : %d read vs %d should_read\n.", __FILE__, __LINE__ , read,should_read);\
        exit(1);}}\


#endif

