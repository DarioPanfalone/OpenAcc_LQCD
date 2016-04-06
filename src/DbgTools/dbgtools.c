#ifndef DBGTOOLS_OPENACC_C_
#define DBGTOOLS_OPENACC_C_

#include <complex.h>

#include "../OpenAcc/geometry.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/single_types.h"
#include "../OpenAcc/io.h"

void save_gl_fermion(global_vec3_soa * const fermion, 
        const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int i = 0 ; i < GL_SIZEH ; i++){
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c0[i]),cimag(fermion->c0[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c1[i]),cimag(fermion->c1[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c2[i]),cimag(fermion->c2[i]));
    }

    fclose(fp);

}

int read_gl_fermion(global_vec3_soa * fermion, const char* nomefile)
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("vec3_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{
        if(verbosity_lv > 2) 
            printf("Reading vec3_soa %s\n", nomefile );

    for(int i = 0 ; i < GL_SIZEH ; i++){
        double re,im;
        CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c0[i] = re + im * I;
        CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c1[i] = re + im * I;
        CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c2[i] = re + im * I;
    }
    fclose(fp);
    return 0;
    }

}


void print_vec3_soa_wrapper(vec3_soa * const fermion, 
        const char* nomefile)
 {

    printf("MPI%02d - Saving whole fermion...\n", devinfo.myrank);

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subfermion_from_rank(ferm_rw,irank);
        recv_loc_subfermion_from_buffer(ferm_rw,fermion,0);
        save_gl_fermion(ferm_rw, nomefile);
    }
    else  send_lnh_subfermion_to_master(fermion,devinfo.myrank);

#else 
//    recv_loc_subfermion_from_buffer(ferm_rw,fermion,0);
    save_gl_fermion(ferm_rw, nomefile);
#endif


}


int read_vec3_soa_wrapper(vec3_soa * fermion, const char* nomefile)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        if(verbosity_lv > 2)
            printf("MPI%02d - reading global fermion \n",devinfo.myrank );
        error = read_gl_fermion(ferm_rw, nomefile);
        MPI_Bcast((void*) &error,1,MPI_INT,0,MPI_COMM_WORLD);

        if(!error) {
            send_lnh_subfermion_to_buffer(ferm_rw,fermion,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subfermion_to_rank(ferm_rw,irank);


        }
        else 
        if(verbosity_lv > 2)
            printf("MPI%02d - no fermion sent!\n",devinfo.myrank );

    }
    else{
        if(verbosity_lv > 2)
            printf("MPI%02d - receiving fermion \n",devinfo.myrank );
        
        MPI_Bcast((void*) &error,1,MPI_INT,0,MPI_COMM_WORLD);
        if(!error){ 
            receive_lnh_subfermion_from_master(fermion);
        }
        else 
        if(verbosity_lv > 2)
            printf("MPI%02d - no fermion received!\n",devinfo.myrank );
    }

#else 
    error = read_gl_fermion(ferm_rw, nomefile);
    if(!error)  send_lnh_subfermion_to_buffer(ferm_rw,fermion,0);
#endif

    return error;

}




void dbg_print_su3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter)
{

    FILE *fp;
    fp = fopen(nomefile,"w");


    int nx = geom_par.gnx;
    int ny = geom_par.gny;
    int nz = geom_par.gnz;
    int nt = geom_par.gnt;




    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",nx,ny,nz,nt,conf_id_iter);
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            // rebuilding the 3rd row
            single_su3 aux;
            single_su3_from_su3_soa(&conf[q],i,&aux);
            rebuild3row(&aux);


            for(int r=0; r < 3 ; r++)
                for(int c=0; c < 3 ; c++) 
                    fprintf(fp, "%.18lf %.18lf\n",
                            creal(aux.comp[r][c]),cimag(aux.comp[r][c])); 
        }
    }
    fclose(fp);
}
int dbgread_su3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter )
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Gauge configuration file %s not readable.\n",nomefile);
        *conf_id_iter = -1;
        return 1;
    }
    else{

        int nx = geom_par.gnx;
        int ny = geom_par.gny;
        int nz = geom_par.gnz;
        int nt = geom_par.gnt;


        int nxt,nyt,nzt,ntt;
        int minus_det_count = 0;
        CHECKREAD(fscanf(fp,"%d\t%d\t%d\t%d\t%d\n",&nxt,&nyt,&nzt,&ntt,conf_id_iter),5);
        if((nx!=nxt)||(ny!=nyt)||(nz!=nzt)||(nz!=nzt)){
            printf("Error, configuration dimensions not compatible with code.\n");
            exit(1);
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
                if(fabs(1+det) < 0.005 ){ // DEBUG, the limit should be FAR stricter.
                    if(verbosity_lv > 5)  printf("Warning in read_su3_soa_ASCII(), Det M = -1.\n");
                    minus_det_count ++;
                    conf[q].r0.c0[i] = -conf[q].r0.c0[i];
                    conf[q].r0.c1[i] = -conf[q].r0.c1[i];
                    conf[q].r0.c2[i] = -conf[q].r0.c2[i];
                    conf[q].r1.c0[i] = -conf[q].r1.c0[i];
                    conf[q].r1.c1[i] = -conf[q].r1.c1[i];
                    conf[q].r1.c2[i] = -conf[q].r1.c2[i];
                    conf[q].r2.c0[i] = -conf[q].r2.c0[i];
                    conf[q].r2.c1[i] = -conf[q].r2.c1[i];
                    conf[q].r2.c2[i] = -conf[q].r2.c2[i];

                }else if (fabs(fabs(det)-1.0)>0.005 ){  // DEBUG, SAME HERE
                    printf("Non unitary matrix read, Det M =  %.18lf\n",det);
                    double pippo =  fabs(det)-1.0;
                    printf("Diff %f,    ", pippo);
                    pippo = fabs(pippo);
                    printf("Diff %f\n", pippo);
                }

            }
        }

        if(minus_det_count !=0 ) {
            printf("Conf read has some matrices with determinant = -1.\n");
            printf("This might be due to the configuration being saved\n");
            printf("is multiplied by the staggered phase.\n");
            printf("In such case, we would expect %d matrices to have determinant = -1\n",
                    nx*ny*nz*nt*3/2 );
            printf("The count of matrices with determinant = -1 is %d. \n", minus_det_count);





        }
        fclose(fp);
        return 0;
    }
}

void dbgprint_gl3_soa(su3_soa * const conf, const char* nomefile,int conf_id_iter)
{

    FILE *fp;
    fp = fopen(nomefile,"w");


    int nx = geom_par.gnx;
    int ny = geom_par.gny;
    int nz = geom_par.gnz;
    int nt = geom_par.gnt;




    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",nx,ny,nz,nt,conf_id_iter);
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            // rebuilding the 3rd row
            single_su3 aux;
            single_gl3_from_su3_soa(&conf[q],i,&aux);


            for(int r=0; r < 3 ; r++)
                for(int c=0; c < 3 ; c++) 
                    fprintf(fp, "%.18lf %.18lf\n",
                            creal(aux.comp[r][c]),cimag(aux.comp[r][c])); 
        }
    }
    fclose(fp);
}
int dbgread_gl3_soa(su3_soa * conf, const char* nomefile,int * conf_id_iter )
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("Gauge configuration file %s not readable.\n",nomefile);
        *conf_id_iter = -1;
        return 1;
    }
    else{

        int nx = geom_par.gnx;
        int ny = geom_par.gny;
        int nz = geom_par.gnz;
        int nt = geom_par.gnt;


        int nxt,nyt,nzt,ntt;
        int minus_det_count = 0;
        CHECKREAD(fscanf(fp,"%d\t%d\t%d\t%d\t%d\n",&nxt,&nyt,&nzt,&ntt,conf_id_iter),5);
        if((nx!=nxt)||(ny!=nyt)||(nz!=nzt)||(nz!=nzt)){
            printf("Error, configuration dimensions not compatible with code.\n");
            exit(1);
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
                if(fabs(1+det) < 0.005 ){ // DEBUG, the limit should be FAR stricter.
                    if(verbosity_lv > 5)  printf("Warning in read_su3_soa_ASCII(), Det M = -1.\n");
                    minus_det_count ++;
                    conf[q].r0.c0[i] = -conf[q].r0.c0[i];
                    conf[q].r0.c1[i] = -conf[q].r0.c1[i];
                    conf[q].r0.c2[i] = -conf[q].r0.c2[i];
                    conf[q].r1.c0[i] = -conf[q].r1.c0[i];
                    conf[q].r1.c1[i] = -conf[q].r1.c1[i];
                    conf[q].r1.c2[i] = -conf[q].r1.c2[i];
                    conf[q].r2.c0[i] = -conf[q].r2.c0[i];
                    conf[q].r2.c1[i] = -conf[q].r2.c1[i];
                    conf[q].r2.c2[i] = -conf[q].r2.c2[i];

                }else if (fabs(fabs(det)-1.0)>0.005 ){  // DEBUG, SAME HERE
                    printf("Non unitary matrix read, Det M =  %.18lf\n",det);
                    double pippo =  fabs(det)-1.0;
                    printf("Diff %f,    ", pippo);
                    pippo = fabs(pippo);
                    printf("Diff %f\n", pippo);
                }

            }
        }

        if(minus_det_count !=0 ) {
            printf("Conf read has some matrices with determinant = -1.\n");
            printf("This might be due to the configuration being saved\n");
            printf("is multiplied by the staggered phase.\n");
            printf("In such case, we would expect %d matrices to have determinant = -1\n",
                    nx*ny*nz*nt*3/2 );
            printf("The count of matrices with determinant = -1 is %d. \n", minus_det_count);





        }
        fclose(fp);
        return 0;
    }
}


void print_vec3_soa(vec3_soa * const fermion, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int i = 0 ; i < sizeh ; i++){
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c0[i]),cimag(fermion->c0[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c1[i]),cimag(fermion->c1[i]));
        fprintf(fp, "%.18lf\t%.18lf\n",creal(fermion->c2[i]),cimag(fermion->c2[i]));
    }
    fclose(fp);

}

int read_vec3_soa(vec3_soa * fermion, const char* nomefile)
{


    FILE *fp;
    fp = fopen(nomefile,"r");
    if(!fp){
        printf("vec3_soa File %s not found.\n", nomefile );
        return 1;
    }
    else{
        if(verbosity_lv > 2) 
            printf("Reading vec3_soa %s\n", nomefile );

    for(int i = 0 ; i < sizeh ; i++){
        double re,im;
        CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c0[i] = re + im * I;
        CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c1[i] = re + im * I;
        CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);fermion->c2[i] = re + im * I;
    }
    fclose(fp);
    return 0;
    }

}


void print_1su3_soa(su3_soa * const conf, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");
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


void print_tamat_soa(tamat_soa * const ipdot, const char* nomefile)
{
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c01[i]),cimag(ipdot[q].c01[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c02[i]),cimag(ipdot[q].c02[i]));
            fprintf(fp, "%.18lf\t%.18lf\n",creal(ipdot[q].c12[i]),cimag(ipdot[q].c12[i]));//
            fprintf(fp, "%.18lf\n",ipdot[q].ic00[i]);
            fprintf(fp, "%.18lf\n",ipdot[q].ic11[i]);
        }
    }
    fclose(fp);
}

int read_tamat_soa(tamat_soa * ipdot, const char* nomefile)
{

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
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c01[i] = re + im * I;
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c02[i] = re + im * I;
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c12[i] = re + im * I;//
            CHECKREAD(fscanf(fp, "%lf\n",&(ipdot[q].ic00[i])),2);
            CHECKREAD(fscanf(fp, "%lf\n",&(ipdot[q].ic11[i])),2);
        }
    }
    fclose(fp);
    return 0;
    }
}


void print_thmat_soa(thmat_soa * const ipdot, const char* nomefile)
{
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
void print_1thmat_soa(thmat_soa * const ipdot, const char* nomefile)
{
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

int read_thmat_soa(thmat_soa * ipdot, const char* nomefile)
{

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
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c01[i] = re + im * I;
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c02[i] = re + im * I;
            CHECKREAD(fscanf(fp, "%lf\t%lf\n",&re,&im),2);ipdot[q].c12[i] = re + im * I;//
            CHECKREAD(fscanf(fp, "%lf\n",&(ipdot[q].rc00[i])),1);
     	    CHECKREAD(fscanf(fp, "%lf\n",&(ipdot[q].rc11[i])),1);
        }
    }
    fclose(fp);
    return 0;
    }
}



void print_double_soa(double_soa * const backfield, const char* nomefile)
{

    
    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int dir = 0 ; dir < 4 ; dir++)
        for(int d3 = 0 ; d3 < nd3; d3++)
        for(int d2 = 0 ; d2 < nd2; d2++)
        for(int d1 = 0 ; d1 < nd1; d1++)
        for(int d0 = 0 ; d0 < nd0; d0++){
            int parity = (d0+d1+d2+d3)%2;

            int idxh = snum_acc(d0,d1,d2,d3);

            fprintf(fp,
                    "dir: %d (d0,d1,d2,d3): %d %d %d %d value: %.18lf\n",
                    dir, d0,d1,d2,d3,
                    backfield[geom_par.xyztmap[dir]+parity].d[idxh]);
    }
    fclose(fp);
}

void print_1double_soa(double_soa * const vettore, const char* nomefile)
{

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int i = 0 ; i < sizeh ; i++){
      fprintf(fp, "%.18lf\n",vettore->d[i]);
    }
    fclose(fp);
}


void transpose_vec3_soa(vec3_soa * vecin, vec3_soa *vecout, int xmap, int ymap, int zmap, int tmap){

// WARNING : NOT TESTED
    int x,y,z,t;
    int d[4], idxh, newidxh;

    for(d[3]=0; d[3] < nd3; d[3]++) for(d[2]=0; d[2] < nd2; d[2]++)
        for(d[1]=0; d[1] < nd1; d[1]++) for(d[0]=0; d[0] < nd0; d[0]++){

                    idxh = snum_acc(d[0],d[1],d[2],d[3]);
                    int tnd[4] = {nd0,nd1,nd2,nd3};
                    x = d[xmap];int tnx = tnd[xmap] ;
                    y = d[ymap];int tny = tnd[ymap] ;
                    z = d[zmap];int tnz = tnd[zmap] ;
                    t = d[tmap];int tnt = tnd[tmap] ;
                    newidxh = (x+tnx*(y+tny*(z+tnz*t)))/2;
                    
                    vecout->c0[newidxh] = vecin->c0[idxh];
                    vecout->c1[newidxh] = vecin->c1[idxh];
                    vecout->c2[newidxh] = vecin->c2[idxh];

        }

}

void transpose_su3_soa(su3_soa * matin, su3_soa *matout, int xmap, int ymap, int zmap, int tmap){

// WARNING : NOT TESTED
    int dir,parity;
    for(dir = 1; dir < 4 ; dir++) for(parity = 1; parity < 2 ; parity++)
    {
        int dirmod = geom_par.xyztmap[dir];
        transpose_vec3_soa(&(matin[2*dirmod+parity].r0),
                &(matout[2*dir+parity].r0),xmap,ymap,zmap,tmap);
        transpose_vec3_soa(&(matin[2*dirmod+parity].r1),
                &(matout[2*dir+parity].r1),xmap,ymap,zmap,tmap);
        transpose_vec3_soa(&(matin[2*dirmod+parity].r2),
                &(matout[2*dir+parity].r2),xmap,ymap,zmap,tmap);

    }

}


#endif
