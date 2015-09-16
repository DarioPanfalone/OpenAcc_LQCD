#ifndef COMMON_DEFINES_H_
#define COMMON_DEFINES_H_

#include "../RationalApprox/rationalapprox.c"

/*****************************************
 * This file is included both in the c++ *
 * and the Openacc Version               *
 * ***************************************/

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

#define beta 5.35 


const int no_flavours=2; // number of quark species
const int start_opt=0;// 0 --> COLD START; 1 --> START FROM SAVED CONF
int conf_id_iter;
int ITERATIONS=20; // the code will generate new <ITERATIONS> confs, from <conf_id_iter+1> to <conf_id_iter+ITERATIONS>
int therm_ITERATIONS = 10; // the first <therm_ITERATIONS> of the history will be thermalization updates

int save_conf_every=10000;

#define max_approx_order 19
int approx_metro=19;
int approx_md=9;
const double lambda_min_metro=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double lambda_min_md=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double residue_metro=1.0e-8;    // stopping residual for CG
const double residue_md=1.0e-5;    // stopping residual for CG
const int gmp_remez_precision=100; // The precision that gmp uses

// quanti di campo esterno
const double bx_quantum=0.0;
const double ex_quantum=0.0;
const double by_quantum=0.0;
const double ey_quantum=0.0;
const double bz_quantum=0.0;
const double ez_quantum=0.0;

#define no_ps 2
#define therm_updates 0 // not used 
#define max_cg 10000


#define no_md 8 // number of MD steps
#define use_multistep 1 // =0 does not use multistep,   =1 2MN_multistep,   =2 4MN_multistep
#define gauge_scale 4  // Update fermions every gauge_scale gauge updates



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


typedef struct ferm_param_t{
  double ferm_charge;
  double ferm_mass;
  double ferm_im_chem_pot;
  int degeneracy;
  int number_of_ps;
  int index_of_the_first_ps;
  RationalApprox approx_fi_mother; // first inv   -> mother
  RationalApprox approx_md_mother; // md approx   -> mother
  RationalApprox approx_li_mother; // last inv    -> mother
  RationalApprox approx_fi;        // first inv
  RationalApprox approx_md;        // md approx
  RationalApprox approx_li;        // last inv
} ferm_param;


int NDiffFlavs;
int NPS_tot;
int max_ps;
ferm_param *fermions_parameters;

char *nome_file_ferm_output;

void init_ferm_params(){

  nome_file_ferm_output = (char*)malloc(sizeof("ferm_meas.dat"));
  strcpy(nome_file_ferm_output,"ferm_meas.dat");


  NDiffFlavs = 2;  // the number of different quark flavours

  int allocation_check; 
  // al posto di 128 c'era ALIGN, solo che qui questa variabile non è ancora definita (viene fatto in struct_c_def)
  allocation_check =  posix_memalign((void **)&fermions_parameters, 128, NDiffFlavs*sizeof(ferm_param));   //  -->  4*size phases (as many as links)
  if(allocation_check != 0)  printf("Errore nella allocazione di fermions_parameters \n");

  fermions_parameters[0].ferm_charge       = -1.0;   // up    charge
  fermions_parameters[0].ferm_mass         = 0.075;  // up    mass
  fermions_parameters[0].ferm_im_chem_pot  = 0.0;    // up    chem pot
  fermions_parameters[0].degeneracy        = 1;      // up    degeneracy
  fermions_parameters[0].number_of_ps      = 1;      // up    number of pseudo fermions

  fermions_parameters[1].ferm_charge       = 2.0;    // down  charge
  fermions_parameters[1].ferm_mass         = 0.075;  // down  mass
  fermions_parameters[1].ferm_im_chem_pot  = 0.0;    // down  chem pot
  fermions_parameters[1].degeneracy        = 1;      // down  degeneracy
  fermions_parameters[1].number_of_ps      = 1;      // down  number of pseudo fermions

  NPS_tot = 0;
  max_ps = fermions_parameters[0].number_of_ps;
  for(int i=0;i<NDiffFlavs;i++){
    // compute the total number of ps
    NPS_tot += fermions_parameters[i].number_of_ps;
    // compute the max number of ps among the various flavs
    if(fermions_parameters[i].number_of_ps>=max_ps) max_ps = fermions_parameters[i].number_of_ps;
    // deterime the offset (where does the ps of the flavour i starts?)
    if(i==0){
      fermions_parameters[i].index_of_the_first_ps=0;
    }else{
      fermions_parameters[i].index_of_the_first_ps = fermions_parameters[i-1].index_of_the_first_ps + fermions_parameters[i-1].number_of_ps;
    }
  }

  printf("NPS_tot = %d \n",NPS_tot);
  printf("max_ps = %d \n",max_ps);


  
  for(int i=0;i<NDiffFlavs;i++){
    fermions_parameters[i].approx_fi_mother.exponent_num =  +fermions_parameters[i].degeneracy;
    fermions_parameters[i].approx_md_mother.exponent_num =  -fermions_parameters[i].degeneracy;
    fermions_parameters[i].approx_li_mother.exponent_num =  -fermions_parameters[i].degeneracy;

    fermions_parameters[i].approx_fi_mother.exponent_den =   fermions_parameters[i].number_of_ps*8;
    fermions_parameters[i].approx_md_mother.exponent_den =   fermions_parameters[i].number_of_ps*4;
    fermions_parameters[i].approx_li_mother.exponent_den =   fermions_parameters[i].number_of_ps*4;

    fermions_parameters[i].approx_fi_mother.approx_order =  approx_metro;
    fermions_parameters[i].approx_md_mother.approx_order =  approx_md;
    fermions_parameters[i].approx_li_mother.approx_order =  approx_metro;

    fermions_parameters[i].approx_fi_mother.lambda_min =  lambda_min_metro;
    fermions_parameters[i].approx_md_mother.lambda_min =  lambda_min_md;
    fermions_parameters[i].approx_li_mother.lambda_min =  lambda_min_metro;

    fermions_parameters[i].approx_fi_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_md_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_li_mother.lambda_max =  1.0;

    fermions_parameters[i].approx_fi_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_md_mother.lambda_max =  1.0;
    fermions_parameters[i].approx_li_mother.lambda_max =  1.0;

    fermions_parameters[i].approx_fi_mother.gmp_remez_precision = gmp_remez_precision;
    fermions_parameters[i].approx_md_mother.gmp_remez_precision = gmp_remez_precision;
    fermions_parameters[i].approx_li_mother.gmp_remez_precision = gmp_remez_precision;

    // copy everything also in the daughter approxs
    fermions_parameters[i].approx_fi.exponent_num =   fermions_parameters[i].approx_fi_mother.exponent_num;
    fermions_parameters[i].approx_md.exponent_num =   fermions_parameters[i].approx_md_mother.exponent_num;
    fermions_parameters[i].approx_li.exponent_num =   fermions_parameters[i].approx_li_mother.exponent_num;
    fermions_parameters[i].approx_fi.exponent_den =   fermions_parameters[i].approx_fi_mother.exponent_den;
    fermions_parameters[i].approx_md.exponent_den =   fermions_parameters[i].approx_md_mother.exponent_den;
    fermions_parameters[i].approx_li.exponent_den =   fermions_parameters[i].approx_li_mother.exponent_den;
    fermions_parameters[i].approx_fi.approx_order =   fermions_parameters[i].approx_fi_mother.approx_order;
    fermions_parameters[i].approx_md.approx_order =   fermions_parameters[i].approx_md_mother.approx_order;
    fermions_parameters[i].approx_li.approx_order =   fermions_parameters[i].approx_li_mother.approx_order;
    fermions_parameters[i].approx_fi.gmp_remez_precision =   fermions_parameters[i].approx_fi_mother.gmp_remez_precision;
    fermions_parameters[i].approx_md.gmp_remez_precision =   fermions_parameters[i].approx_md_mother.gmp_remez_precision;
    fermions_parameters[i].approx_li.gmp_remez_precision =   fermions_parameters[i].approx_li_mother.gmp_remez_precision;

    // READ THE RAT APPROXS FROM THE FILES
    rationalapprox_read(&(fermions_parameters[i].approx_fi_mother));
    rationalapprox_read(&(fermions_parameters[i].approx_md_mother));
    rationalapprox_read(&(fermions_parameters[i].approx_li_mother));


  }

  

}


#endif

