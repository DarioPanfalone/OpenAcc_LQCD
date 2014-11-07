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

#define mass 0.2 

#endif


