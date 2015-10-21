#ifndef STRUCT_C_DEF_
#define STRUCT_C_DEF_


#ifndef COMPLEX_C_
#define COMPLEX_C_
#include <complex.h>
#endif

#include <stdlib.h>
//#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include "../Include/common_defines.h"

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

//#define mass2 mass*mass
// used in the dynamical allocation of structures
#define ALIGN 128

#define ONE_BY_THREE   0.33333333333333333333333
#define ONE_BY_SIX     0.16666666666666666666666
#define beta_by_three  beta*ONE_BY_THREE


#ifdef GAUGE_ACT_TLSM
#define GAUGE_ACTION   1     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         5.0*ONE_BY_THREE
#define C_ONE          -0.25*ONE_BY_THREE
#endif

#ifdef GAUGE_ACT_WILSON
#define GAUGE_ACTION   0     // 0 --> Wilson; 1 --> tree level Symanzik
#define C_ZERO         1.0
#define C_ONE          0.0
#endif


#define acc_twopi 2.0 * M_PI 

// strutture native c con i complessi c --> quelle che verranno utilizzate nelle routine che verranno accelerate da openacc
typedef double complex  d_complex;

int no_md_acc,gauge_scale_acc;
double epsilon_acc;
d_complex ieps_acc,iepsh_acc;

void initialize_global_variables(){
  no_md_acc = no_md ;
  gauge_scale_acc = gauge_scale;
  epsilon_acc = 1.0/((double)(no_md_acc));
  ieps_acc  = 0.0 + (epsilon_acc) * 1.0I;
  iepsh_acc = 0.0 + (epsilon_acc) * 0.5 * 1.0I;
}


static inline int snum_acc(int x, int y, int z, int t) {
  int ris;
  ris = x + (y*vol1) + (z*vol2) + (t*vol3);
  return ris/2;   // <---  /2 Pay attention to even/odd  (see init_geo)
}


typedef struct vec3_soa_t {
    int flag;
    d_complex c0[sizeh];
    d_complex c1[sizeh];
    d_complex c2[sizeh];
} vec3_soa;

void checkfree_vec3soa(vec3_soa* fermion){
    if (fermion->flag !=0) 
        printf("CHECK: FERMION AT %p IS NOT FREE.\n", fermion);
}



typedef struct dcomplex_soa_t {
    int flag;
    d_complex c[sizeh];
} dcomplex_soa;

void checkfree_dcomplex_soa(dcomplex_soa* x){
    if (x->flag !=0) 
        printf("CHECK: DCOMPLEX_SOA AT %p IS NOT FREE.\n", x);
}
typedef struct double_soa_t {
    int flag;
    double d[sizeh];
} double_soa;

void checkfree_dcomplex_soa(double_soa* x){
    if (x->flag !=0) 
        printf("CHECK: DOUBLE_SOA AT %p IS NOT FREE.\n", x);
}
typedef struct single_su3_t {
  d_complex comp[3][3];
} single_su3;

typedef struct vec3_t {
  d_complex c0;
  d_complex c1;
  d_complex c2;
} vec3;

typedef struct su3_soa_t {
  int flag;  
  vec3_soa r0;
  vec3_soa r1;
  vec3_soa r2;
} su3_soa;
void checkfree_su3soa(su_3soa* M){
    if (M->flag !=0) 
        printf("CHECK: GAUGE MATRIX AT %p IS NOT FREE.\n", M);
}


typedef struct tamat_soa_t {
  int flag;
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double rc00[sizeh];   // Im(comp_00)
  double rc11[sizeh];   // Im(comp_11)
} tamat_soa;

void checkfree_tamat_soa(tamat_soa* M){
    if (M->flag !=0) 
        printf("CHECK: TAMAT_SOA AT %p IS NOT FREE.\n", M);
}

typedef struct single_tamat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double rc00;   // Im(comp_00)
  double rc11;   // Im(comp_11)
} single_tamat;

typedef struct thmat_soa_t {
  int flag;
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double rc00[sizeh];   // Re(comp_00)
  double rc11[sizeh];   // Re(comp_11)
} thmat_soa;

void checkfree_thmat_soa(thmat_soa* M){
    if (M->flag !=0) 
        printf("CHECK: THMAT_SOA AT %p IS NOT FREE.\n", M);
}




typedef struct single_thmat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double rc00;   // Re(comp_00)
  double rc11;   // Re(comp_11)
} single_thmat;









// TABLES FOR THE NEAREST NEIGHBOURS       
// nnp[site_half][dir][par] = nearest neighbour in the positive direction "dir"            
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"         
// nnm[site_half][dir][par] = nearest neighbour in the negative direction "dir"            
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"         
// temporarily defined and computed here (should be moved elsewhere!)      
int nnp_openacc[sizeh][4][2];
int nnm_openacc[sizeh][4][2];
void compute_nnp_and_nnm_openacc(){
  int x, y, z, t,parity;
  for(t=0; t<nt; t++) {
    for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
        for(x=0; x < nx; x++) {
          int  xm, ym, zm, tm, xp, yp, zp, tp, idxh;
          idxh = snum_acc(x,y,z,t);
          parity = (x+y+z+t) % 2;

          xm = x - 1;
          xm = xm + (((xm >> 31) & 0x1) * nx);
          ym = y - 1;
          ym = ym + (((ym >> 31) & 0x1) * ny);
          zm = z - 1;
          zm = zm + (((zm >> 31) & 0x1) * nz);
          tm = t - 1;
          tm = tm + (((tm >> 31) & 0x1) * nt);

          nnm_openacc[idxh][0][parity] = snum_acc(xm,y,z,t);
          nnm_openacc[idxh][1][parity] = snum_acc(x,ym,z,t);
          nnm_openacc[idxh][2][parity] = snum_acc(x,y,zm,t);
          nnm_openacc[idxh][3][parity] = snum_acc(x,y,z,tm);

          xp = x + 1;
          xp *= (((xp-nx) >> 31) & 0x1);
          yp = y + 1;
          yp *= (((yp-ny) >> 31) & 0x1);
          zp = z + 1;
          zp *= (((zp-nz) >> 31) & 0x1);
          tp = t + 1;
          tp *= (((tp-nt) >> 31) & 0x1);

          nnp_openacc[idxh][0][parity] = snum_acc(xp,y,z,t);
          nnp_openacc[idxh][1][parity] = snum_acc(x,yp,z,t);
          nnp_openacc[idxh][2][parity] = snum_acc(x,y,zp,t);
          nnp_openacc[idxh][3][parity] = snum_acc(x,y,z,tp);

        }
      }
    }
  }

}


#endif

