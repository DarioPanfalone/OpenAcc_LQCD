#ifndef STRUCT_C_DEF_
#define STRUCT_C_DEF_


#ifndef COMPLEX_C_
#define COMPLEX_C_
#include <complex.h>
#endif

#include <stdlib.h>
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

#define ONE_BY_THREE 0.33333333333333333333333
#define ONE_BY_SIX   0.16666666666666666666666
#define beta_by_three beta*ONE_BY_THREE

#define acc_pi2 2.0 * M_PI 

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
  d_complex c0[sizeh];
  d_complex c1[sizeh];
  d_complex c2[sizeh];
} vec3_soa;

typedef struct dcomplex_soa_t {
  d_complex c[sizeh];
} dcomplex_soa;

typedef struct double_soa_t {
  double d[sizeh];
} double_soa;

typedef struct single_su3_t {
  d_complex comp[3][3];
} single_su3;

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
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double rc00[sizeh];   // Im(comp_00)
  double rc11[sizeh];   // Im(comp_11)
} tamat_soa;
typedef struct single_tamat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double rc00;   // Im(comp_00)
  double rc11;   // Im(comp_11)
} single_tamat;

typedef struct thmat_soa_t {
  d_complex c01[sizeh]; // comp_01
  d_complex c02[sizeh]; // comp_02
  d_complex c12[sizeh]; // comp_12
  double rc00[sizeh];   // Re(comp_00)
  double rc11[sizeh];   // Re(comp_11)
} thmat_soa;
typedef struct single_thmat_t {
  d_complex c01; // comp_01
  d_complex c02; // comp_02
  d_complex c12; // comp_12
  double rc00;   // Re(comp_00)
  double rc11;   // Re(comp_11)
} single_thmat;






// funzioni di conversione:    thmat_soa   ==>   thmatCOM_soa 
void convert_thmat_soa_to_thmatCOM_soa(thmat_soa *in, thmatCOM_soa *out){
  int i;
  for(i =0 ; i < sizeh ; i++){
    out->c01[i].Re = creal(in->c01[i]);
    out->c02[i].Re = creal(in->c02[i]);
    out->c12[i].Re = creal(in->c12[i]);

    out->c01[i].Im = cimag(in->c01[i]);
    out->c02[i].Im = cimag(in->c02[i]);
    out->c12[i].Im = cimag(in->c12[i]);

    out->rc00[i] = in->rc00[i];
    out->rc11[i] = in->rc11[i];
  }
}

// funzioni di conversione:    thmatCOM_soa   ==>   thmat_soa 
void convert_thmatCOM_soa_to_thmat_soa(thmatCOM_soa *in, thmat_soa *out){
  int i;
  for( i =0 ; i < sizeh ; i++){
    out->c01[i] = in->c01[i].Re + in->c01[i].Im * 1.0I;
    out->c02[i] = in->c02[i].Re + in->c02[i].Im * 1.0I;
    out->c12[i] = in->c12[i].Re + in->c12[i].Im * 1.0I;
    out->rc00[i] = in->rc00[i];
    out->rc11[i] = in->rc11[i];
  }
}

// funzioni di conversione:    tamat_soa   ==>   tamatCOM_soa 
void convert_tamat_soa_to_tamatCOM_soa(tamat_soa *in, tamatCOM_soa *out){
  int i;
  for(i =0 ; i < sizeh ; i++){
    out->c01[i].Re = creal(in->c01[i]);
    out->c02[i].Re = creal(in->c02[i]);
    out->c12[i].Re = creal(in->c12[i]);

    out->c01[i].Im = cimag(in->c01[i]);
    out->c02[i].Im = cimag(in->c02[i]);
    out->c12[i].Im = cimag(in->c12[i]);

    out->rc00[i] = in->rc00[i];
    out->rc11[i] = in->rc11[i];
  }
}

// funzioni di conversione:    tamatCOM_soa   ==>   tamat_soa 
void convert_tamatCOM_soa_to_tamat_soa(tamatCOM_soa *in, tamat_soa *out){
  int i;
  for( i =0 ; i < sizeh ; i++){
    out->c01[i] = in->c01[i].Re + in->c01[i].Im * 1.0I;
    out->c02[i] = in->c02[i].Re + in->c02[i].Im * 1.0I;
    out->c12[i] = in->c12[i].Re + in->c12[i].Im * 1.0I;
    out->rc00[i] = in->rc00[i];
    out->rc11[i] = in->rc11[i];
  }
}
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
void convert_su3COM_soa_to_su3_soa(const su3COM_soa *in, su3_soa *out){
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
void convert_vec3COM_soa_to_vec3_soa(const vec3COM_soa *in, vec3_soa *out){
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
void convert_COM_MultiFermion_to_ACC_MultiFermion(const COM_MultiFermion *in, ACC_MultiFermion *out){
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

