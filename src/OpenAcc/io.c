#ifndef IO_C_
#define IO_C_

#include "./io.h"
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "./geometry.h"
#include "./struct_c_def.h"
#include "../Include/common_defines.h"
#include "./single_types.h"



void print_su3_soa_ASCII(su3_soa * const conf, const char* nomefile,int conf_id_iter){

    FILE *fp;
    fp = fopen(nomefile,"w");
    if(! fp ){
        printf("ERROR, %s unreadable.\n",nomefile);
        exit(1);
    }

    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",nx,ny,nz,nt,conf_id_iter);
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r0.c0[i]),cimag(conf[q].r0.c0[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r0.c1[i]),cimag(conf[q].r0.c1[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r0.c2[i]),cimag(conf[q].r0.c2[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r1.c0[i]),cimag(conf[q].r1.c0[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r1.c1[i]),cimag(conf[q].r1.c1[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r1.c2[i]),cimag(conf[q].r1.c2[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r2.c0[i]),cimag(conf[q].r2.c0[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r2.c1[i]),cimag(conf[q].r2.c1[i])); 
            fprintf(fp, "%.18lf %.18lf\n",creal(conf[q].r2.c2[i]),cimag(conf[q].r2.c2[i])); 
        }
    }
    fclose(fp);
}
int read_su3_soa_ASCII(su3_soa * conf, const char* nomefile,int * conf_id_iter ){


  FILE *fp;
  fp = fopen(nomefile,"r");
  if(!fp){
      printf("Gauge configuration file %s not readable.\n",nomefile);
      *conf_id_iter = -1;
      return 1;
  }
  else{



  int nxt,nyt,nzt,ntt;
  CHECKREAD(fscanf(fp,"%d\t%d\t%d\t%d\t%d\n",&nxt,&nyt,&nzt,&ntt,conf_id_iter),5);
  if((nx!=nxt)||(ny!=nyt)||(nz!=nzt)||(nz!=nzt)){
    printf(" Errore: DIMENSIONI DELLA CONFIGURAZIONE LETTA DIVERSE DA QUELLE ATTESE\n");
    abort();
  }

  for(int q = 0 ; q < 8 ; q++){
    for(int i = 0 ; i < sizeh ; i++){
      double re,im;
      single_su3 m;double det;
      //      fscanf(fp, "%.18lf\t%.18lf\n",&re,&im);
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][0]=conf[q].r0.c0[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][1]=conf[q].r0.c1[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[0][2]=conf[q].r0.c2[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][0]=conf[q].r1.c0[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][1]=conf[q].r1.c1[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[1][2]=conf[q].r1.c2[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][0]=conf[q].r2.c0[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][1]=conf[q].r2.c1[i] = re + im * I;
      CHECKREAD(fscanf(fp, "%lf %lf\n",&re,&im),2);m.comp[2][2]=conf[q].r2.c2[i] = re + im * I;
      det = detSu3(&m);
      if((float)det==(float)-1){
          printf("Warning in read_su3_soa_ASCII(), Det M = -1.\n");
        conf[q].r0.c0[i] = -conf[q].r0.c0[i];
        conf[q].r0.c1[i] = -conf[q].r0.c1[i];
        conf[q].r0.c2[i] = -conf[q].r0.c2[i];
        conf[q].r1.c0[i] = -conf[q].r1.c0[i];
        conf[q].r1.c1[i] = -conf[q].r1.c1[i];
        conf[q].r1.c2[i] = -conf[q].r1.c2[i];
        conf[q].r2.c0[i] = -conf[q].r2.c0[i];
        conf[q].r2.c1[i] = -conf[q].r2.c1[i];
        conf[q].r2.c2[i] = -conf[q].r2.c2[i];

      }else if ((float)det!=(float)1) printf("Matrice non unitaria, determinante %.18lf\n",det);
      
    }
  }
  fclose(fp);
  return 0;
  }
}

#endif
