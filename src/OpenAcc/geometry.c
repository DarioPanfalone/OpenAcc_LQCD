#ifndef GEOMETRY_C_
#define GEOMETRY_C_

#include "./geometry.h"


geom_parameters geom_par;


void compute_nnp_and_nnm_openacc(void){
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
