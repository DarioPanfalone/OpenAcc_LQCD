//
//  Mod_functions.c
//  
//
//  Created by Luca Parente on 02/06/21.
//

#include "Mod_functions.h"
#include "./struct_c_def.h"
#include "./alloc_settings.h"
#include "./alloc_vars.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "./geometry.h"
#include "../Mpi/multidev.h"






/////funzione aggiunta//////////////////////
//per ora la metto qui poi la sposto in struct_def.h

//this function does actually nothing useful, I wrote it down just for studying how bc have been implementated.
void counter_size_function(int d0,int d1, int d2, int d3){
    int id=0;
    for(d3=0; d3<nd3; d3++) {
        for(d2=0; d2<nd2; d2++) {
            for(d1=0; d1<nd1; d1++) {
                for(d0=0; d0 < nd0; d0++) { //i cicli sui siti del reticolo.
                    
                    id=snum_acc(d0,d1,d2,d3);
                    printf("(%d,%d,%d,%d):      id=%d\n",d0,d1,d2,d3,id);
                }
            }
        }
    }
    
    return;
}

/*
void init_k_values(su3_soa * conf,int c_r,int * pos_def){


    for(mu=0;mu<8;mu++){
        for(i=0; i<sizeh; i++ ){
            
            
            conf[mu].K.d[i]=c_r;
            }
    
    
    return;
}


void init_k(su3_soa * conf,int c_r,int def_axis,int * def_vet){
    
    int a0,af,b0,bf,c0,cf; //the cuboid defect region's sides. a0 start, af end.
    int pos_def[6];
  
    
    switch (def_axis) {
        case 0:
            printf("defect on x's boundary");
            
            pos_def[1]=snum_acc(0,def_vet[1],);
            
            
            break;
            
        case 1:
            printf("defect on y's boundary");
            
            break;
            
        case 2:
            printf("defect on z's boundary");
            
            break;
            
        case 3:
            printf("defect on t's boundary");
            
            break;
            
            
            
        default:
            printf("ERROR WRONG AXIS CHOICE!\n");
            break;
    }
    
    
    
    
    
  
}


*/




//function that tested the initialization of k_mu  by setting every value to one.
int init_k_test(su3_soa *conf_acc,int c_r){
    int kk2=0;
    int mu1=0;
    for(mu1=0;mu1<8;mu1++){
        for(kk2=0;kk2<sizeh;kk2++){
            if(conf_acc[mu1].K.d[kk2]!=(1 && c_r)){
                return 1;
            }
        printf("%d:ku[%d]:%f\n",mu1,kk2, conf_acc[mu1].K.d[kk2]);
        }
    }
    return 0;
}







//Outdated and useless function. Once used before a better mod of the setting_parser_file.c
/*
int n_replicas_reader(const char* input_filename){
    int value_nr;
    int trovato=0;
    int i=0;
    char riga[20];
    char riga2[20]="Replicas number";
    FILE *input = fopen(input_filename,"r"); //questo ovviamente apre il file
    printf("LETTURA DEL NUMERO DI REPLICHE\n");
    while(trovato==0){
        fgets(riga,20,input);
      //  printf("ecc %d\n",i);
 
        //i=i+1;
        if(strncmp(riga,riga2,15)==0){
          fscanf(input,"%d",&value_nr );
            trovato=1;
        }
        
        
        
        
        
    }
    
    
    
    fclose(input);
    
    return value_nr;
}
*/

