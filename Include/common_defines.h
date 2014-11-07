#ifndef COMMON_DEFINES_H_
#define COMMON_DEFINES_H_

/*****************************************
 * This file is included both in the c++ *
 * and the Openacc Version               *
 * ***************************************/

#define DIM_BLOCK_X 8 // This should divide (nx/2)
#define DIM_BLOCK_Y 8 // This should divide ny
#define DIM_BLOCK_Z 8  // This should divide nz*nt

// lattice dimensions
#define nx 16
#define ny 8
#define nz 4
#define nt 2

#define sizehh nx*ny*nz*nt/2 

#define mass 0.2 

typedef struct COM_t{
  double Re;
  double Im;
} COM;


typedef struct vec3COM_soa_t {
  COM c0[sizehh];
  COM c1[sizehh];
  COM c2[sizehh];
} vec3COM_soa;

typedef struct vec3COM_t {
  COM c0;
  COM c1;
  COM c2;
} vec3COM;

typedef struct su3COM_soa_t {
  vec3COM_soa r0;
  vec3COM_soa r1;
  vec3COM_soa r2;
} su3COM_soa;

#endif


