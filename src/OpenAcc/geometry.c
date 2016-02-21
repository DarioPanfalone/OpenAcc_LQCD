#ifndef GEOMETRY_C_
#define GEOMETRY_C_

#include "./geometry.h"


geom_parameters geom_par;

int nd[4]={nd0,nd1,nd2,nd3};
int vol3s[4]={vol4/nd0,vol4/nd1,vol4/nd2,vol4/nd3};

void compute_nnp_and_nnm_openacc(void){
  int d0, d1, d2, d3,parity;
  for(d3=0; d3<nd3; d3++) {
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
        for(d0=0; d0 < nd0; d0++) {
          int  d0m, d1m, d2m, d3m, d0p, d1p, d2p, d3p, idxh;
          idxh = snum_acc(d0,d1,d2,d3);
          parity = (d0+d1+d2+d3) % 2;

          d0m = d0 - 1;
          d0m = d0m + (((d0m >> 31) & 0x1) * nd0);
          d1m = d1 - 1;
          d1m = d1m + (((d1m >> 31) & 0x1) * nd1);
          d2m = d2 - 1;
          d2m = d2m + (((d2m >> 31) & 0x1) * nd2);
          d3m = d3 - 1;
          d3m = d3m + (((d3m >> 31) & 0x1) * nd3);

          nnm_openacc[idxh][0][parity] = snum_acc(d0m,d1,d2,d3);
          nnm_openacc[idxh][1][parity] = snum_acc(d0,d1m,d2,d3);
          nnm_openacc[idxh][2][parity] = snum_acc(d0,d1,d2m,d3);
          nnm_openacc[idxh][3][parity] = snum_acc(d0,d1,d2,d3m);

          d0p = d0 + 1;
          d0p *= (((d0p-nd0) >> 31) & 0x1);
          d1p = d1 + 1;
          d1p *= (((d1p-nd1) >> 31) & 0x1);
          d2p = d2 + 1;
          d2p *= (((d2p-nd2) >> 31) & 0x1);
          d3p = d3 + 1;
          d3p *= (((d3p-nd3) >> 31) & 0x1);

          nnp_openacc[idxh][0][parity] = snum_acc(d0p,d1,d2,d3);
          nnp_openacc[idxh][1][parity] = snum_acc(d0,d1p,d2,d3);
          nnp_openacc[idxh][2][parity] = snum_acc(d0,d1,d2p,d3);
          nnp_openacc[idxh][3][parity] = snum_acc(d0,d1,d2,d3p);

        }
      }
    }
  }
}



#endif
