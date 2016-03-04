#ifndef GEOMETRY_C_
#define GEOMETRY_C_

#include "./geometry.h"


geom_parameters geom_par;

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

void set_geom_glv(geom_parameters* gp){


    gp->nd[0]= nd0; gp->nd[1] = nd1; gp->nd[2] = nd2; gp->nd[3] = nd3;
    gp->nranks[0] = NRANKS_D0; gp->nranks[1] = NRANKS_D1; 
    gp->nranks[2] = NRANKS_D2; gp->nranks[3] = NRANKS_D3;

    gp->vol3s[0] = vol4/nd0; gp->vol3s[1] = vol4/nd1;
    gp->vol3s[2] = vol4/nd2; gp->vol3s[3] = vol4/nd3;

    gp->xyztmap[0] = gp->xmap;  gp->xyztmap[1] = gp->ymap;
    gp->xyztmap[2] = gp->zmap;  gp->xyztmap[3] = gp->tmap;

    gp->d0123map[gp->xmap] = 0;   gp->d0123map[gp->ymap] = 1;
    gp->d0123map[gp->zmap] = 2;   gp->d0123map[gp->tmap] = 3;


}


#endif
