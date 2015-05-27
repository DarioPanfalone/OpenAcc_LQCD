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
#define ny 16
#define nz 16
#define nt 16
#define sizehh nx*ny*nz*nt/2 

#define ANTIPERIODIC_T_BC  // else periodic time bc are taken

#define mass 0.32
#define beta 6.0

#define therm_updates 0
#define max_cg 100000
#define max_approx_order 19
#define no_ps 2

#define no_md 15 // number of MD steps
#define use_multistep 1 // =0 does not use multistep,   =1 2MN_multistep,   =2 4MN_multistep
#define gauge_scale 10  // Update fermions every gauge_scale gauge updates




typedef struct COM_t{
  double Re;
  double Im;
} COM;


typedef struct vec3COM_soa_t {
  COM c0[sizehh];
  COM c1[sizehh];
  COM c2[sizehh];
} vec3COM_soa;

typedef struct tamatCOM_soa_t {
  COM c01[sizehh];
  COM c02[sizehh];
  COM c12[sizehh];
  double rc00[sizehh];
  double rc11[sizehh];
} tamatCOM_soa;

typedef struct thmatCOM_soa_t {
  COM c01[sizehh];
  COM c02[sizehh];
  COM c12[sizehh];
  double rc00[sizehh];
  double rc11[sizehh];
} thmatCOM_soa;

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


//SHIFT FERMIONS   
typedef struct COM_ShiftFermion_t{
  vec3COM_soa shift[max_approx_order];
} COM_ShiftFermion;

//MULTI FERMIONS   
typedef struct COM_MultiFermion_t{
  vec3COM_soa multi[no_ps];
} COM_MultiFermion;


//SHIFT MULTI FERMIONS     
typedef struct COM_ShiftMultiFermion_t{
  vec3COM_soa shiftmulti[max_approx_order][no_ps];
} COM_ShiftMultiFermion;

typedef struct COM_RationalApprox_t{
  double COM_min_epsilon;
  int COM_approx_order;
  double COM_RA_a0;
  double COM_RA_a[max_approx_order];
  double COM_RA_b[max_approx_order];
}COM_RationalApprox;


#endif


