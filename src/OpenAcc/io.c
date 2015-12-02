#ifndef IO_C_
#define IO_C_

#include "./io.h"
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "./geometry.h"
#include "./struct_c_def.h"


void print_su3_soa_ASCII(su3_soa * const conf, const char* nomefile,int conf_id_iter){

    FILE *fp;
    fp = fopen(nomefile,"w");
    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",nx,ny,nz,nt,conf_id_iter);
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r0.c0[i]),cimag(conf[q].r0.c0[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r0.c1[i]),cimag(conf[q].r0.c1[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r0.c2[i]),cimag(conf[q].r0.c2[i]));//
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r1.c0[i]),cimag(conf[q].r1.c0[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r1.c1[i]),cimag(conf[q].r1.c1[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r1.c2[i]),cimag(conf[q].r1.c2[i]));//
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r2.c0[i]),cimag(conf[q].r2.c0[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r2.c1[i]),cimag(conf[q].r2.c1[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(conf[q].r2.c2[i]),cimag(conf[q].r2.c2[i]));
        }
    }
    fclose(fp);
}
void read_su3_soa_ASCII(su3_soa * conf, const char* nomefile,int * conf_id_iter ){

  FILE *fp;
  fp = fopen(nomefile,"r");

  int nxt,nyt,nzt,ntt;
  fscanf(fp,"%d\t%d\t%d\t%d\t%d\n",&nxt,&nyt,&nzt,&ntt,conf_id_iter);
  if((nx!=nxt)||(ny!=nyt)||(nz!=nzt)||(nz!=nzt)){
    printf(" Errore: DIMENSIONI DELLA CONFIGURAZIONE LETTA DIVERSE DA QUELLE ATTESE\n");
    abort();
  }

  for(int q = 0 ; q < 8 ; q++)
    for(int i = 0 ; i < sizeh ; i++){
      double re,im;
      //      fscanf(fp, "%.18lf\t%.18lf\n",&re,&im);
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r0.c0[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r0.c1[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r0.c2[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r1.c0[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r1.c1[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r1.c2[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r2.c0[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r2.c1[i] = re + im * I;
      fscanf(fp, "%lf\t%lf\n",&re,&im);conf[q].r2.c2[i] = re + im * I;
      
    }
  fclose(fp);
  
}

#endif
