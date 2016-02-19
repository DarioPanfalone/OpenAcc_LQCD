#ifndef IO_C_
#define IO_C_

// FOR ILDG IO
#define __USE_XOPEN2K
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <sys/types.h>
#include "endian.h"
#define MAXILDGHEADERS 10


#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include "./io.h"
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
        int minus_det_count = 0;
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


// STUFF FOR READING AND WRITING ILDG FILES

// ENDIANNESS SWAP
void doubleswe(double  * doubletoswap){

    char *p= (char*) doubletoswap; char t[4];
    t[0] = p[7]; p[7] = p[0] ; p[0] = t[0];
    t[1] = p[6]; p[6] = p[1] ; p[1] = t[1];
    t[2] = p[5]; p[5] = p[2] ; p[2] = t[2];
    t[3] = p[4]; p[4] = p[3] ; p[3] = t[3];

}
void uint64swe(uint64_t  * inttoswap){

    char *p= (char*) inttoswap; char t[4];
    t[0] = p[7]; p[7] = p[0] ; p[0] = t[0];
    t[1] = p[6]; p[6] = p[1] ; p[1] = t[1];
    t[2] = p[5]; p[5] = p[2] ; p[2] = t[2];
    t[3] = p[4]; p[4] = p[3] ; p[3] = t[3];

}
void uint32swe(uint32_t  * inttoswap){

    char *p= (char*) inttoswap; char t[2];
    t[0] = p[3]; p[3] = p[0] ; p[0] = t[0];
    t[1] = p[2]; p[2] = p[1] ; p[1] = t[1];

}
void uint16swe(uint16_t  * inttoswap){

    char *p= (char*) inttoswap; char t;
    t = p[1]; p[1] = p[0] ; p[0] = t;

}
typedef struct ILDG_header_t
{
    uint32_t magic_no; // LIME magic number 0x456789ab
    uint16_t version;
    uint16_t mbme_flag;
    uint64_t data_length;
    char type[128];
}ILDG_header;

int read_su3_soa_ildg_binary(su3_soa * conf, const char* nomefile,int * conf_id_iter )
{
    FILE *fg;int reads;
    char header[1000];
    *conf_id_iter = 1000 ; // random number


    ILDG_header ildg_headers[MAXILDGHEADERS]; // all ildg headers
    // note _FILE_OFFSET_BITS 64 
    off_t ildg_header_ends_positions[MAXILDGHEADERS]; // point in files 
    // where the message starts

    int ildg_format_index =-1;   // where, in ildg_headers, is the ildg_format header
    int ildg_binary_data_index =-1;// where, in ildg_headers, is the ildg_binary_data header

    // ILDG confs are BIG ENDIAN
    int conf_machine_endianness_disagreement = (__BYTE_ORDER == __LITTLE_ENDIAN) ;

    fg = fopen(nomefile,"r+");
    if(!fg){
        printf("Gauge configuration ILDG file %s not readable.\n",nomefile);
        *conf_id_iter = -1;
        return 1;
    } 

    int alldimfound = 0;
    fseek(fg,0,SEEK_SET);
    // Find all headers
    int i=0;
    while((ildg_format_index==-1 ||  ildg_binary_data_index==-1) && fg){
        reads = fread(&ildg_headers[i],sizeof(ILDG_header),1,fg);
        if(reads!= 1 ){
            printf("Error in reading file: %s ,%d\n",__FILE__,__LINE__);
            exit(1);
        }




        ildg_header_ends_positions[i] = ftello(fg);
        if(strcmp(ildg_headers[i].type,"ildg-format")==0) ildg_format_index = i;
        if(strcmp(ildg_headers[i].type,"ildg-binary-data")==0) ildg_binary_data_index = i;

        if(conf_machine_endianness_disagreement){
            uint32swe(&ildg_headers[i].magic_no);
            uint16swe(&ildg_headers[i].version);
            uint16swe(&ildg_headers[i].mbme_flag);
            uint64swe(&ildg_headers[i].data_length);
        }

        if(verbosity_lv>3){
            printf("header type: %s, data length:%"PRIu64", magic number: %"PRIu32"\n",
                    ildg_headers[i].type, ildg_headers[i].data_length,
                    ildg_headers[i].magic_no);
            printf("version: %d, mbme_flag:%d\n",
                    ildg_headers[i].version, ildg_headers[i].mbme_flag);

        }





        // padding to 8 bytes
        off_t missing_bytes = 
            (off_t)(ildg_headers[i].data_length%8==0?0:8-ildg_headers[i].data_length%8);
        // note _FILE_OFFSET_BITS 64 
        off_t mod_data_length =  ildg_headers[i].data_length + missing_bytes ;
        printf("%"PRIu64" vs %"PRIu64"\n", (uint64_t) mod_data_length, ildg_headers[i].data_length);

        // going to the next record  
        fseeko(fg,mod_data_length,SEEK_CUR);

        i++;
    }


    // read ildg-format for 
    char ildg_format_str[1000];
    if(verbosity_lv>2)
        printf("Reading ildg-format...\n");


    fseeko(fg,ildg_header_ends_positions[ildg_format_index],SEEK_SET);
    reads = fread(ildg_format_str,1,ildg_headers[ildg_format_index].data_length,fg);
    if(reads!= ildg_headers[ildg_format_index].data_length){
        printf("Error in reading file: %s ,%d\n",__FILE__,__LINE__);
        exit(1);
    }
    char * strfnx = strstr(ildg_format_str,"<lx>");
    char * strfny = strstr(ildg_format_str,"<ly>");
    char * strfnz = strstr(ildg_format_str,"<lz>");
    char * strfnt = strstr(ildg_format_str,"<lt>");
    int nx_r,ny_r,nz_r,nt_r;
    int allfound = 1;
    if(strfnx!= NULL){sscanf(strfnx," <lx> %d </lx> ",&nx_r);printf("nx: %d\n",nx_r);}
    else{allfound = 0;}
    if(strfny!= NULL){sscanf(strfny," <ly> %d </ly> ",&ny_r);printf("ny: %d\n",ny_r);}
    else{allfound = 0;}
    if(strfnz!= NULL){sscanf(strfnz," <lz> %d </lz> ",&nz_r);printf("nz: %d\n",nz_r);}
    else{allfound = 0;}
    if(strfnt!= NULL){sscanf(strfnt," <lt> %d </lt> ",&nt_r);printf("nt: %d\n",nt_r);}
    else{allfound = 0;}
    if(!allfound){
        printf("Error, %s:%d : lx,ly,lz or lt not found in \"ildg-format\"\n",
                __FILE__,__LINE__);
        printf("Please check. Exiting now.\n");
        exit(1);
    }



    // read ildg-binary-data (su3 gauge conf)
    if(verbosity_lv>2)
        printf("Reading ildg-binary-data...\n");



    fseeko(fg,ildg_header_ends_positions[ildg_binary_data_index],SEEK_SET);
    int x,y,z,t,dir;
    for(t=0;t<nt_r;t++) for(z=0;z<nz_r;z++)
        for(y=0;y<ny_r;y++) for(x=0;x<nx_r;x++)
        {
            int idxh = snum_acc(x,y,z,t);
            int parity = (x+y+z+t)%2;

            for(dir=0;dir<4;dir++){

                single_su3 m;          
                reads = fread((void*)m.comp,sizeof(double),18,fg);
                if(reads!= 18){
                    printf("Error in reading file: %s ,%d\n",__FILE__,__LINE__);
                    exit(1);
                }


                // check enddiannes
                if(conf_machine_endianness_disagreement){
                    int irow,icol;
                    for(irow=0;irow<2;irow++) for(icol=0;icol<3;icol++){
                        double re = creal(m.comp[irow][icol]);
                        double im = cimag(m.comp[irow][icol]);
                        doubleswe(&re);
                        doubleswe(&im);
                        m.comp[irow][icol] = re + im*I;
                    }
                    single_su3_into_su3_soa(&conf[2*dir+parity],idxh,&m);
                }

            }

        }

    fclose(fg);

    return 0;

}

void print_su3_soa_ildg_binary(su3_soa * const conf, const char* nomefile,
        int conf_id_iter)
{
    printf("Writing ildg conf...\n");
    int conf_machine_endianness_disagreement = (__BYTE_ORDER == __LITTLE_ENDIAN) ;
    int writes; FILE *fp;
    fp = fopen(nomefile,"w");
    if(! fp ){
        printf("ERROR, %s unreadable.\n",nomefile);
        exit(1);
    }

    char ildg_format_str[1000];  

    sprintf(ildg_format_str,
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
            <ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n\
            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n\
            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n\
            <version>1.0</version>\n\
            <field>su3gauge</field>\n\
            <precision>64</precision>\n\
            <lx>%d</lx>\n\
            <ly>%d</ly>\n\
            <lz>%d</lz>\n\
            <lt>%d</lt>\n\
            </ildgFormat>",nx,ny,nz,nt);
    uint64_t len = strlen(ildg_format_str) ;
    uint32_t magic_number =  0x456789ab;
    uint16_t version = 1;
    uint16_t mbme_flag = 0;
    uint64_t len64 = len;

    if(conf_machine_endianness_disagreement){
        uint32swe(&magic_number);
        uint16swe(&version);
        uint16swe(&mbme_flag);
        uint64swe(&len64);
    }
    ILDG_header ildg_format_header = 
        (ILDG_header){magic_number,version,mbme_flag,len64,"ildg-format"};
    fwrite(&ildg_format_header,sizeof(ILDG_header),1,fp);
    fwrite(ildg_format_str,1,len,fp);

    char pad[8] = {0,0,0,0,0,0,0,0};

    // padding to 8 bytes
    off_t missing_bytes = (off_t)(len%8==0 ? 0:8-len%8);
    fwrite(pad,1,missing_bytes,fp);


    len = nd0*nd1*nd2*nd3*4*3*3*2*8;
    len64 = len;
    if(verbosity_lv > 3) printf("Writing ildg-binary-data, %" PRIu64 "\n", len64);
    if(conf_machine_endianness_disagreement) uint64swe(&len64);

    ILDG_header ildg_binary_data_header = 
        (ILDG_header){magic_number,version,mbme_flag,len64,"ildg-binary-data"};
    fwrite(&ildg_binary_data_header,sizeof(ILDG_header),1,fp);

    int x,y,z,t,dir;
    for(t=0;t<nt;t++) for(z=0;z<nz;z++)
        for(y=0;y<ny;y++) for(x=0;x<nx;x++)
        {
            int idxh = snum_acc(x,y,z,t);
            int parity = (x+y+z+t)%2;

            for(dir=0;dir<4;dir++){

                single_su3 m;          
                single_su3_from_su3_soa(&conf[2*dir+parity],idxh,&m);
                rebuild3row(&m);

                // fix enddiannes
                if(conf_machine_endianness_disagreement){
                    int irow,icol;
                    for(irow=0;irow<3;irow++) for(icol=0;icol<3;icol++){
                        double re = creal(m.comp[irow][icol]);
                        double im = cimag(m.comp[irow][icol]);
                        doubleswe(&re);
                        doubleswe(&im);

                        writes =  fwrite((void*)&re,sizeof(double),1,fp);
                        if(writes!= 1){
                           printf("Error in writing file: %s ,%d\n",__FILE__,__LINE__);
                        exit(1);
                        }
                        writes =  fwrite((void*)&im,sizeof(double),1,fp);
                        if(writes!= 1){
                           printf("Error in writing file: %s ,%d\n",__FILE__,__LINE__);
                        exit(1);
                        }

                    }
                }
            }
        }

    // padding to 8 bytes
    missing_bytes = (off_t)(len%8==0?0:8-len%8);
    fwrite(pad,1,missing_bytes,fp);

    ILDG_header ildg_data_lfn_header = 
        (ILDG_header){0x456789ab,1,0,0,"ildg-data-lfn"};
    fwrite(&ildg_binary_data_header,sizeof(ILDG_header),1,fp);

    fclose(fp);
    

    return ;
}





#endif
