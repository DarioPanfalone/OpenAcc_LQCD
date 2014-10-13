#ifndef COMPLEX_C_
#define COMPLEX_C_
#include <complex.h>
#endif


#define DIM_BLOCK_X 8 // This should divide (nx/2)
#define DIM_BLOCK_Y 8 // This should divide ny
#define DIM_BLOCK_Z 8  // This should divide nz*nt
#define nx 16
#define ny 16
#define nz 16
#define nt 16
#define vol1 nx
#define vol2 (ny * vol1)
#define vol3 (nz * vol2)
#define vol4 (nt * vol3)
//#define nxh (nx >> 1) // nx/2
//#define nyh (ny >> 1)
//#define nzh (nz >> 1)
//#define nth (nt >> 1)

#define nxh 8
#define nyh 8
#define nzh 8
#define nth 8

#define size vol4
#define size2 (2*size)
#define size3 (3*size)
//#define sizeh (size / 2)
#define sizeh 32768 // che sarebbe 16^4/2
#define no_links (4 * vol4)

#define ALIGN 128

static inline int snum_acc(int x, int y, int z, int t) {
  int ris;
  ris = x + (y*vol1) + (z*vol2) + (t*vol3);
  return ris/2;   // <---  /2 Pay attention to even/odd  (see init_geo)
}


// strutture intermedie che ci servono per traghettare le conf ed il resto da c a c++
typedef struct COM_t{
  double Re;
  double Im;
} COM;


typedef struct vec3COM_soa_t {
  COM c0[sizeh];
  COM c1[sizeh];
  COM c2[sizeh];
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



// strutture native c con i complessi c --> quelle che verranno utilizzate nelle routine che verranno accelerate da openacc
typedef double complex  d_complex;

typedef struct vec3_soa_t {
  d_complex c0[sizeh];
  d_complex c1[sizeh];
  d_complex c2[sizeh];
} vec3_soa;

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

