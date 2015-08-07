#ifndef COMMON_DEFINES_H_
#define COMMON_DEFINES_H_

/*****************************************
 * This file is included both in the c++ *
 * and the Openacc Version               *
 * ***************************************/

#define BACKFIELD
#define IMCHEMPOT

// se BACKFIELD o IMCHEMPOT sono definiti allora nell'applicazione della matrice
// di dirac usa la routine che moltiplica anche per la fase opportuna, altrimenti
// in modo hard coded usa l'altra routine moltiplica link per fermione e basta.
// Per farlo viene definita la variabile PHASE_MAT_VEC_MULT o meno.
#ifdef BACKFIELD
  #define PHASE_MAT_VEC_MULT 
#else
  #ifdef IMCHEMPOT
    #define PHASE_MAT_VEC_MULT 
  #endif
#endif

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

#define mass 0.03
#define beta 6.0

const int no_flavours=2; // number of quark species                                                                                                          


#define max_approx_order 19
#define max_ps 2
double approx_metro=19;
double approx_md=9;
const double lambda_min_metro=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double lambda_min_md=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double residue_metro=1.0e-8;    // stopping residual for CG
const double residue_md=1.0e-5;    // stopping residual for CG

// quanti di campo esterno
const double bx_quantum=0.0;
const double ex_quantum=0.0;
const double by_quantum=0.0;
const double ey_quantum=0.0;
const double bz_quantum=0.0;
const double ez_quantum=0.0;

#define no_ps 2
#define therm_updates 0 // not used 
#define max_cg 100000


#define no_md 13 // number of MD steps
#define use_multistep 1 // =0 does not use multistep,   =1 2MN_multistep,   =2 4MN_multistep
#define gauge_scale 5  // Update fermions every gauge_scale gauge updates


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


typedef struct ferm_param_t{
  double ferm_charge;
  double ferm_mass;
  int degeneracy;
  int number_of_ps;
  COM_RationalApprox approx1; // first inv  // prima c'era *approx1, ora ho messo approx1  e basta
  COM_RationalApprox approx2; // md approx
  COM_RationalApprox approx3; // last inv
} ferm_param;


int NDiffFlavs;
int NPS_tot;
ferm_param *fermions_parameters;

void init_ferm_params(){

  NDiffFlavs = 2;  // the number of different quark flavours

  int allocation_check;
  allocation_check =  posix_memalign((void **)&fermions_parameters, ALIGN, NDiffFlavs*sizeof(ferm_param));   //  -->  4*size phases (as many as links)
  if(allocation_check != 0)  printf("Errore nella allocazione di fermions_parameters \n");

  fermions_parameters[0].ferm_charge       = -1.0;   // up    charge
  fermions_parameters[0].ferm_mass         = mass;   // up    mass
  fermions_parameters[0].ferm_im_chem_pot  = 0.0;    // up    chem pot
  fermions_parameters[0].degeneracy        = 1;      // up    degeneracy
  fermions_parameters[0].number_of_ps      = 1;      // up    number of pseudo fermions

  fermions_parameters[1].ferm_charge       = 2.0;    // down  charge
  fermions_parameters[1].ferm_mass         = mass;   // down  mass
  fermions_parameters[1].ferm_im_chem_pot  = 0.0;    // down  chem pot
  fermions_parameters[1].degeneracy        = 1;      // down  degeneracy
  fermions_parameters[1].number_of_ps      = 1;      // down  number of pseudo fermions

  NPS_tot = 0;
  for(int i=0;i<NDiffFlavs;i++)
    NPS_tot += fermions_parameters[i].number_of_ps;
  

}


#endif


