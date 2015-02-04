
#ifndef STRUCT_C_DEF_
#define STRUCT_C_DEF_


#ifndef COMPLEX_C_
#define COMPLEX_C_
#include <complex.h>
#endif

#include <math.h>
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

#define mass2 mass*mass


// used in the dynamical allocation of structures
#define ALIGN 128

#define ONE_BY_THREE 0.333333333333333333

static inline int snum_acc(int x, int y, int z, int t) {
  int ris;
  ris = x + (y*vol1) + (z*vol2) + (t*vol3);
  return ris/2;   // <---  /2 Pay attention to even/odd  (see init_geo)
}


// strutture native c con i complessi c --> quelle che verranno utilizzate nelle routine che verranno accelerate da openacc
typedef double complex  d_complex;

typedef struct vec3_soa_t {
  d_complex c0[sizeh];
  d_complex c1[sizeh];
  d_complex c2[sizeh];
} vec3_soa;

typedef struct dcomplex_soa_t {
  d_complex c[sizeh];
} dcomplex_soa;



typedef struct vec3_t {
  d_complex c0;
  d_complex c1;
  d_complex c2;
} vec3;

typedef struct su3_soa_t {
  vec3_soa r0;
  vec3_soa r1;
  vec3_soa r2;
} su3_soa;

//SHIFT FERMIONS
typedef struct ACC_ShiftFermion_t{
  vec3_soa shift[max_approx_order];
} ACC_ShiftFermion;

//MULTI FERMIONS
typedef struct ACC_MultiFermion_t{
  vec3_soa multi[no_ps];
} ACC_MultiFermion;

//SHIFT MULTI FERMIONS
typedef struct ACC_ShiftMultiFermion_t{
  vec3_soa shiftmulti[max_approx_order][no_ps];
} ACC_ShiftMultiFermion;


typedef struct tamat_soa_t {
  d_complex c01[sizeh];
  d_complex c02[sizeh];
  d_complex c12[sizeh];
  double rc00[sizeh];
  double rc11[sizeh];
} tamat_soa;




// funzioni di conversione:    su3_soa   ==>   su3COM_soa 
void convert_su3_soa_to_su3COM_soa(su3_soa *in, su3COM_soa *out){
  int i;
  for(i =0 ; i < sizeh ; i++){
    out->r0.c0[i].Re = creal(in->r0.c0[i]);
    out->r1.c0[i].Re = creal(in->r1.c0[i]);
    out->r2.c0[i].Re = creal(in->r2.c0[i]);
    out->r0.c1[i].Re = creal(in->r0.c1[i]);
    out->r1.c1[i].Re = creal(in->r1.c1[i]);
    out->r2.c1[i].Re = creal(in->r2.c1[i]);
    out->r0.c2[i].Re = creal(in->r0.c2[i]);
    out->r1.c2[i].Re = creal(in->r1.c2[i]);
    out->r2.c2[i].Re = creal(in->r2.c2[i]);

    out->r0.c0[i].Im = cimag(in->r0.c0[i]);
    out->r1.c0[i].Im = cimag(in->r1.c0[i]);
    out->r2.c0[i].Im = cimag(in->r2.c0[i]);
    out->r0.c1[i].Im = cimag(in->r0.c1[i]);
    out->r1.c1[i].Im = cimag(in->r1.c1[i]);
    out->r2.c1[i].Im = cimag(in->r2.c1[i]);
    out->r0.c2[i].Im = cimag(in->r0.c2[i]);
    out->r1.c2[i].Im = cimag(in->r1.c2[i]);
    out->r2.c2[i].Im = cimag(in->r2.c2[i]);
  }
}

// funzioni di conversione:    su3COM_soa   ==>   su3_soa 
void convert_su3COM_soa_to_su3_soa(su3COM_soa *in, su3_soa *out){
  int i;
  for( i =0 ; i < sizeh ; i++){
    out->r0.c0[i] = (in->r0.c0[i].Re) + (in->r0.c0[i].Im) * 1.0I;
    out->r1.c0[i] = (in->r1.c0[i].Re) + (in->r1.c0[i].Im) * 1.0I;
    out->r2.c0[i] = (in->r2.c0[i].Re) + (in->r2.c0[i].Im) * 1.0I;

    out->r0.c1[i] = (in->r0.c1[i].Re) + (in->r0.c1[i].Im) * 1.0I;
    out->r1.c1[i] = (in->r1.c1[i].Re) + (in->r1.c1[i].Im) * 1.0I;
    out->r2.c1[i] = (in->r2.c1[i].Re) + (in->r2.c1[i].Im) * 1.0I;

    out->r0.c2[i] = (in->r0.c2[i].Re) + (in->r0.c2[i].Im) * 1.0I;
    out->r1.c2[i] = (in->r1.c2[i].Re) + (in->r1.c2[i].Im) * 1.0I;
    out->r2.c2[i] = (in->r2.c2[i].Re) + (in->r2.c2[i].Im) * 1.0I;
  }
}

// funzioni di conversione:    vec3_soa  ==>  vec3COM_soa
void convert_vec3_soa_to_vec3COM_soa(vec3_soa *in, vec3COM_soa *out){
  int i;
  for( i =0 ; i < sizeh ; i++){
    out->c0[i].Re = creal(in->c0[i]);
    out->c1[i].Re = creal(in->c1[i]);
    out->c2[i].Re = creal(in->c2[i]);

    out->c0[i].Im = cimag(in->c0[i]);
    out->c1[i].Im = cimag(in->c1[i]);
    out->c2[i].Im = cimag(in->c2[i]);
  }
}

// funzioni di conversione:    vec3COM_soa  ==>  vec3_soa
void convert_vec3COM_soa_to_vec3_soa(vec3COM_soa *in, vec3_soa *out){
  int i;
  for( i =0 ; i < sizeh ; i++){
    out->c0[i] = (in->c0[i].Re) + (in->c0[i].Im) * 1.0I;
    out->c1[i] = (in->c1[i].Re) + (in->c1[i].Im) * 1.0I;
    out->c2[i] = (in->c2[i].Re) + (in->c2[i].Im) * 1.0I;
  }
}


// funzioni di conversione:    ACC_ShiftFermion  ==>  COM_ShiftFermion          
void convert_ACC_ShiftFermion_to_COM_ShiftFermion(ACC_ShiftFermion *in, COM_ShiftFermion *out){
  int i,ia;
  for(ia =0 ; ia < max_approx_order ; ia++){
    for( i =0 ; i < sizeh ; i++){
      out->shift[ia].c0[i].Re = creal(in->shift[ia].c0[i]);
      out->shift[ia].c1[i].Re = creal(in->shift[ia].c1[i]);
      out->shift[ia].c2[i].Re = creal(in->shift[ia].c2[i]);

      out->shift[ia].c0[i].Im = cimag(in->shift[ia].c0[i]);
      out->shift[ia].c1[i].Im = cimag(in->shift[ia].c1[i]);
      out->shift[ia].c2[i].Im = cimag(in->shift[ia].c2[i]);
    }
  }
}

// funzioni di conversione:    COM_ShiftFermion  ==>  ACC_ShiftFermion          
void convert_COM_ShiftFermion_to_ACC_ShiftFermion(COM_ShiftFermion *in, ACC_ShiftFermion *out){
  int i,ia;
  for(ia =0 ; ia < max_approx_order ; ia++){
    for( i =0 ; i < sizeh ; i++){
      out->shift[ia].c0[i] = (in->shift[ia].c0[i].Re) + (in->shift[ia].c0[i].Im) * 1.0I;
      out->shift[ia].c1[i] = (in->shift[ia].c1[i].Re) + (in->shift[ia].c1[i].Im) * 1.0I;
      out->shift[ia].c2[i] = (in->shift[ia].c2[i].Re) + (in->shift[ia].c2[i].Im) * 1.0I;
    }
  }
}

// funzioni di conversione:    ACC_MultiFermion  ==>  COM_MultiFermion          
void convert_ACC_MultiFermion_to_COM_MultiFermion(ACC_MultiFermion *in, COM_MultiFermion *out){
  int i,ips;
  for(ips =0 ; ips < no_ps ; ips++){
    for( i =0 ; i < sizeh ; i++){
      out->multi[ips].c0[i].Re = creal(in->multi[ips].c0[i]);
      out->multi[ips].c1[i].Re = creal(in->multi[ips].c1[i]);
      out->multi[ips].c2[i].Re = creal(in->multi[ips].c2[i]);

      out->multi[ips].c0[i].Im = cimag(in->multi[ips].c0[i]);
      out->multi[ips].c1[i].Im = cimag(in->multi[ips].c1[i]);
      out->multi[ips].c2[i].Im = cimag(in->multi[ips].c2[i]);
    }
  }
}

// funzioni di conversione:    COM_MultiFermion  ==>  ACC_MultiFermion          
void convert_COM_MultiFermion_to_ACC_MultiFermion(COM_MultiFermion *in, ACC_MultiFermion *out){
  int i,ips;
  for(ips =0 ; ips < no_ps ; ips++){
    for( i =0 ; i < sizeh ; i++){
      out->multi[ips].c0[i] = (in->multi[ips].c0[i].Re) + (in->multi[ips].c0[i].Im) * 1.0I;
      out->multi[ips].c1[i] = (in->multi[ips].c1[i].Re) + (in->multi[ips].c1[i].Im) * 1.0I;
      out->multi[ips].c2[i] = (in->multi[ips].c2[i].Re) + (in->multi[ips].c2[i].Im) * 1.0I;
    }
  }
}



// funzioni di conversione:    ACC_ShiftMultiFermion  ==>  COM_ShiftMultiFermion
void convert_ACC_ShiftMultiFermion_to_COM_ShiftMultiFermion(ACC_ShiftMultiFermion *in, COM_ShiftMultiFermion *out){
  int i,ips,ia;
  for(ips =0 ; ips < no_ps ; ips++){
    for(ia =0 ; ia < max_approx_order ; ia++){
      for( i =0 ; i < sizeh ; i++){
        out->shiftmulti[ia][ips].c0[i].Re = creal(in->shiftmulti[ia][ips].c0[i]);
        out->shiftmulti[ia][ips].c1[i].Re = creal(in->shiftmulti[ia][ips].c1[i]);
        out->shiftmulti[ia][ips].c2[i].Re = creal(in->shiftmulti[ia][ips].c2[i]);

        out->shiftmulti[ia][ips].c0[i].Im = cimag(in->shiftmulti[ia][ips].c0[i]);
        out->shiftmulti[ia][ips].c1[i].Im = cimag(in->shiftmulti[ia][ips].c1[i]);
        out->shiftmulti[ia][ips].c2[i].Im = cimag(in->shiftmulti[ia][ips].c2[i]);
      }
    }
  }
}

// funzioni di conversione:    COM_ShiftMultiFermion  ==>  ACC_ShiftMultiFermion
void convert_COM_ShiftMultiFermion_to_ACC_ShiftMultiFermion(COM_ShiftMultiFermion *in, ACC_ShiftMultiFermion *out){
  int i,ips,ia;
  for(ia =0 ; ia < max_approx_order ; ia++){
    for(ips =0 ; ips < no_ps ; ips++){
      for( i =0 ; i < sizeh ; i++){
        out->shiftmulti[ia][ips].c0[i] = (in->shiftmulti[ia][ips].c0[i].Re) + (in->shiftmulti[ia][ips].c0[i].Im) * 1.0I;
        out->shiftmulti[ia][ips].c1[i] = (in->shiftmulti[ia][ips].c1[i].Re) + (in->shiftmulti[ia][ips].c1[i].Im) * 1.0I;
        out->shiftmulti[ia][ips].c2[i] = (in->shiftmulti[ia][ips].c2[i].Re) + (in->shiftmulti[ia][ips].c2[i].Im) * 1.0I;
      }
    }
  }
}

#endif
