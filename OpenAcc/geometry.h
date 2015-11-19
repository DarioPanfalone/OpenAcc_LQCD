#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "../Include/common_defines.h"

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

#define vol1 nx
#define vol2 (ny * vol1)
#define vol3 (nz * vol2)
#define vol4 (nt * vol3)
#define nxh (nx >> 1) // nx/2
#define nyh (ny >> 1)
#define nzh (nz >> 1)
#define nth (nt >> 1)

#define size vol4
#define size2 (2*size)
#define size3 (3*size)
#define sizeh (size / 2)
#define no_links (4 * vol4)


#pragma acc routine seq
static inline int snum_acc(int x, int y, int z, int t) {
  int ris;
  ris = x + (y*vol1) + (z*vol2) + (t*vol3);
  return ris/2;   // <---  /2 Pay attention to even/odd  (see init_geo)
}

// TABLES FOR THE NEAREST NEIGHBOURS       
// nnp[site_half][dir][par] = nearest neighbour in the positive direction "dir"            
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"         
// nnm[site_half][dir][par] = nearest neighbour in the negative direction "dir"            
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"         
// temporarily defined and computed here (should be moved elsewhere!)      
int nnp_openacc[sizeh][4][2];
int nnm_openacc[sizeh][4][2];

void compute_nnp_and_nnm_openacc(void);

#endif
