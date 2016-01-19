#ifndef DBGTOOLS_OPENACC_
#define DBGTOOLS_OPENACC_

#include <complex.h>

#include "../OpenAcc/geometry.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/single_types.h"

void print_vec3_soa(vec3_soa * const fermion, const char* nomefile){

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int i = 0 ; i < sizeh ; i++){
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c0[i]),cimag(fermion->c0[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c1[i]),cimag(fermion->c1[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c2[i]),cimag(fermion->c2[i]));
    }
    fclose(fp);

}

int read_vec3_soa(vec3_soa * fermion, const char* nomefile){


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("vec3_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{

    for(int i = 0 ; i < sizeh ; i++){
        double re,im;
        fscanf(fp, "%lf\t%lf\n",&re,&im);fermion->c0[i] = re + im * I;
        fscanf(fp, "%lf\t%lf\n",&re,&im);fermion->c1[i] = re + im * I;
        fscanf(fp, "%lf\t%lf\n",&re,&im);fermion->c2[i] = re + im * I;
    }
    fclose(fp);
    return 0;
    }

}

void print_su3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter){

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
void print_1su3_soa(su3_soa * const conf, const char* nomefile){

    FILE *fp;
    fp = fopen(nomefile,"w");
    fprintf(fp,"%d\t%d\t%d\t%d\t\n",nx,ny,nz,nt);
    for(int q = 0 ; q < 1 ; q++){ // q = 0 !!!!!
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

void read_su3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter ){

  FILE *fp;
  fp = fopen(nomefile,"r");

  int nxt,nyt,nzt,ntt;
  fscanf(fp,"%d\t%d\t%d\t%d\t%d\n",&nxt,&nyt,&nzt,&ntt,conf_id_iter);
  if((nx!=nxt)||(ny!=nyt)||(nz!=nzt)||(nz!=nzt)){
    printf(" Errore: DIMENSIONI DELLA CONFIGURAZIONE LETTA DIVERSE DA QUELLE ATTESE\n");
    abort();
  }

  double max_error = 0;
  for(int q = 0 ; q < 8 ; q++)
    for(int i = 0 ; i < sizeh ; i++){
      double re,im;
      //      fscanf(fp, "%.18lf\t%.18lf\n",&re,&im);
      single_su3 m;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[0][0] =conf[q].r0.c0[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[0][1] =conf[q].r0.c1[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[0][2] =conf[q].r0.c2[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[1][0] =conf[q].r1.c0[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[1][1] =conf[q].r1.c1[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[1][2] =conf[q].r1.c2[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[2][0] =conf[q].r2.c0[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[2][1] =conf[q].r2.c1[i] =re+im*I;
      fscanf(fp, "%lf %lf\n",&re,&im); m.comp[2][2] =conf[q].r2.c2[i] =re+im*I;
      double error = 1-fabs(detSu3(&m));
      if (error>max_error) max_error = error;
      
    }
  printf("Max deviation from unitarity - 1-|detM|: %e (Wrong,but staggered phases are here).\n", max_error);
  fclose(fp);
      
}

void print_tamat_soa(tamat_soa * const ipdot, const char* nomefile){
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c01[i]),cimag(ipdot[q].c01[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c02[i]),cimag(ipdot[q].c02[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c12[i]),cimag(ipdot[q].c12[i]));//
            fprintf(fp, "%.18lf\n",ipdot[q].rc00[i]);
            fprintf(fp, "%.18lf\n",ipdot[q].rc11[i]);
        }
    }
    fclose(fp);
}

int read_tamat_soa(tamat_soa * ipdot, const char* nomefile){

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Tamat_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            double re,im;
            fscanf(fp, "%lf\t%lf\n",&re,&im);ipdot[q].c01[i] = re + im * I;
            fscanf(fp, "%lf\t%lf\n",&re,&im);ipdot[q].c02[i] = re + im * I;
            fscanf(fp, "%lf\t%lf\n",&re,&im);ipdot[q].c12[i] = re + im * I;//
            fscanf(fp, "%lf\n",&(ipdot[q].rc00[i]));
            fscanf(fp, "%lf\n",&(ipdot[q].rc11[i]));
        }
    }
    fclose(fp);
    return 0;
    }
}


void print_thmat_soa(thmat_soa * const ipdot, const char* nomefile){
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c01[i]),cimag(ipdot[q].c01[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c02[i]),cimag(ipdot[q].c02[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c12[i]),cimag(ipdot[q].c12[i]));//
            fprintf(fp, "%.18lf\n",ipdot[q].rc00[i]);
            fprintf(fp, "%.18lf\n",ipdot[q].rc11[i]);
        }
    }
    fclose(fp);
}
void print_1thmat_soa(thmat_soa * const ipdot, const char* nomefile){
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 1 ; q++){// q = 1 !!!!
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c01[i]),cimag(ipdot[q].c01[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c02[i]),cimag(ipdot[q].c02[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c12[i]),cimag(ipdot[q].c12[i]));//
            fprintf(fp, "%.18lf\n",ipdot[q].rc00[i]);
            fprintf(fp, "%.18lf\n",ipdot[q].rc11[i]);
        }
    }
    fclose(fp);
}

int read_thmat_soa(thmat_soa * ipdot, const char* nomefile){

    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Thmat_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            double re,im;
            fscanf(fp, "%lf\t%lf\n",&re,&im);ipdot[q].c01[i] = re + im * I;
            fscanf(fp, "%lf\t%lf\n",&re,&im);ipdot[q].c02[i] = re + im * I;
            fscanf(fp, "%lf\t%lf\n",&re,&im);ipdot[q].c12[i] = re + im * I;//
            fscanf(fp, "%lf\n",&(ipdot[q].rc00[i]));
	    fscanf(fp, "%lf\n",&(ipdot[q].rc11[i]));
        }
    }
    fclose(fp);
    return 0;
    }
}



void print_double_soa(double_soa * const backfield, const char* nomefile){

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\n",backfield[q].d[i]);
        }
    }
    fclose(fp);
}

void print_1double_soa(double_soa * const vettore, const char* nomefile){

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int i = 0 ; i < sizeh ; i++){
      fprintf(fp, "%.18lf\n",vettore->d[i]);
    }
    fclose(fp);
}


#endif
