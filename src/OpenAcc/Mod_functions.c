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
    
    conf[mu].K.d[snum_acc(2*i,j,k,t)]
    return;
}
*/


int init_k(su3_soa * conf,double c_r,int def_axis,int * def_vet){
    int mu;
    int i,j,z,t;
    int res=0;
    int defect_volume;
    int parity;
    
    
    int counter=0;
   //il case considera le 4 possibili direzioni in cui viene fatto il defect.
    
    switch (def_axis) {
        case 0:
            printf("defect on x's boundary\n");
        for(mu=0;mu<4;mu++){
           for(t=0;t<nd3;t++) {
                for (z=0; z<nd2; z++){
                     for(j=0;j<nd1;j++){
                       for(i=0;i<nd0;i++){
                            parity = (i+j+z+t) % 2;
                           
                            if(j>=0 && j<def_vet[0] && z>=0 && z<def_vet[1] && t>=0 && t<def_vet[2] && i==((nd0)-1) && mu==0 )
                            {
                                
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;
                                  /*  printf("(%d,%d,%d,%d):  k_mu[%d]=%f (%d)parity\n",i,j,z,t,snum_acc(i,j,z,t),conf[2*mu].K.d[snum_acc(i,j,z,t)],parity);*/
                                } //inizializza il vettore
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;
                                    /*printf("(%d,%d,%d,%d):  k_mu[%d]=%f (%d)parity\n",i,j,z,t,snum_acc(i,j,z,t),conf[2*mu+1].K.d[snum_acc(i,j,z,t)],parity);*/
                                } //inizializza il vettore
                                    
                                
                                counter=counter+1;
                                
                                ;
                     /* printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                                 }
                                
                           
                        
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1; }
                                
                                
                                
                            } //else inizializza a 1.
                        
                     /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                       }
                    }
                }
            }
            
        }
        
            
            
            
            break;
            
        case 1:
            printf("defect on y's boundary\n");
            for(mu=0;mu<4;mu++){
                for(t=0;t<nd3;t++) {
                    for (z=0; z<nd2; z++){
                        for(j=0;j<nd1;j++){
                            for(i=0;i<nd0;i++){
                            if(i>=0 && i<def_vet[0] && z>=0 && z<def_vet[1] && t>=0 && t<def_vet[2] && j==nd1-1 && mu==1 )
                            {
                                 parity = (i+j+z+t) % 2;
                                
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                
                                
                                counter=counter+1;
                                
                            }
                                
                                
                                
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                
                                
                            } //else inizializza a 1.
                                
                                /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                                
                    }
              }
            }
          }
        }
            
            
            
            
            
            break;
            
        case 2:
            printf("defect on z's boundary\n");
            for(mu=0;mu<4;mu++){
                for(t=0;t<nd3;t++) {
                    for (z=0; z<nd2; z++){
                        for(j=0;j<nd1;j++){
                            for(i=0;i<nd0;i++){
                            if(i>=0 && i<def_vet[0] && j>=0 && j<def_vet[1] && t>=0 && t<def_vet[2] && z==nd2-1 && mu==2 )
                            {
                                 parity = (i+j+z+t) % 2;
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                
                                
                                counter=counter+1;
                                /* printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                            }
                                
                                
                                
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                
                                
                            } //else inizializza a 1.
                                
                                /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                    }
                }
            }
            
           }
      }
            
            
            
            break;
            
        case 3:
            printf("defect on t's boundary\n");
            for(mu=0;mu<4;mu++){
                for(t=0;t<nd3;t++) {
                    for (z=0; z<nd2; z++){
                        for(j=0;j<nd1;j++){
                            for(i=0;i<nd0;i++){
                            if(i>=0 && i<def_vet[0] && j>=0 && j<def_vet[1] && z>=0 && z<def_vet[2] && t==nd3-1 && mu==3 )
                            {
                                 parity = (i+j+z+t) % 2;
                                
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                if(parity!=0){conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=c_r;} //inizializza il vettore}
                                
                                
                                counter=counter+1;
                                /* printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                            }
                                
                                
                                
                            else{
                                if (parity==0){conf[2*mu].K.d[snum_acc(i,j,z,t)]=1;}
                                else{conf[2*mu+1].K.d[snum_acc(i,j,z,t)]=1;}
                                
                                
                                
                            } //else inizializza a 1.
                                
                                /*printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);*/
                    }
                }
            }
            
         }
        }
            break;
            
            
            
        default:
            printf("ERROR WRONG AXIS CHOICE!\n");
            res=1;
            break;
    }
    
    defect_volume=(def_vet[0])*(def_vet[1])*(def_vet[2]);
    
   printf("counter %d\n",counter);
    
    if(counter!=defect_volume){printf("wrong defect initialization!\n"); res=1;}
    
    return res;
  
}






//function that tested the initialization of k_mu  by setting every value to one. //OUTADATED
int init_k_test(su3_soa *conf_acc,int c_r){
    int kk2=0;
    int mu1=0;
    for(mu1=0;mu1<8;mu1++){
        for(kk2=0;kk2<sizeh;kk2++){
            if(conf_acc[mu1].K.d[kk2]!=(1 && c_r)){
                return 1;
            }
        printf("%d:ku[%d]:%f\n",mu1,kk2,conf_acc[mu1].K.d[kk2]);
        }
    }
    return 0;
}
/* work in progress
int defect_shifter_function(su3_soa * conf,int def_axis,int * def_vet){
        int mu;
    double aux;
    
    switch (def_axis){
            
         case 0: for(mu=0;mu<8;mu++){
             for(t=def_vet[4];t<def_vet[5];t++) {
                 for (z=def_vet[2]; z<def_vet[3]; z++){
                     for(j=def_vet[0];j<def_vet[1];j++){
             
             
                         conf[mu].K.d[snum_acc((nd0/2)-1),j,z,t)];
                         
            
                    }
                 }
             }
         }
            break;
            
            
        case 1:
            
    }
    
    
    
    
    return 0;
}


*/




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


void printing_k_mu(su3_soa * conf){
    int mu,t,z,j,i;
    
    for(mu=0;mu<8;mu++){
        for(t=0;t<nd3;t++) {
            for (z=0; z<nd2; z++){
                for(j=0;j<nd1;j++){
                    for(i=0;i<(nd0/2);i++){
    
    printf("(%d,%d,%d,%d):     k_mu[%d]=%f\n",2*i,j,z,t,snum_acc(2*i,j,z,t),conf[mu].K.d[snum_acc(2*i,j,z,t)]);
                    }
                }
            }
            
        }
    }
    return;
}

int replicas_swap_1(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet ){
    double aux;
    int i,j,k,t,mu,parity;
    int res=0;
    //test variable
    int counter=0;

    switch (def_axis){
            
        case 0:
            i=nd0-1;
                mu=0;
                for(t=0;t<def_vet[2];t++){
                    for(k=0;k<def_vet[1];k++){
                        for(j=0;j<def_vet[0];j++){
                            
                            parity = (i+j+k+t) % 2;
                            
                            //test
                            
                            
                            //K_mu_values swap
                            if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                           /* printf("beforrre (%d) %f %f",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );*/
                                conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                                conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                         /*   printf("aftermath (%d) %f %f\n",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );*/
                            }
                            if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                         /*   printf("beforrre (%d) %f %f",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );*/
                            conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                        /*    printf("aftermath (%d) %f %f\n",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );*/
                            }
                            //test
                            
                           /* counter=counter+1;*/
                        }
                    }
                }
            break;
            
        case 1:
            j=nd0-1;
            mu=1;
            for(t=0;t<def_vet[2];t++){
                for(k=0;k<def_vet[1];k++){
                    for(i=0;i<def_vet[0];i++){
                         parity = (i+j+k+t) % 2;
                        
                        //K_mu_values swap
                        if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                            conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                            conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        
                        
                        
            
                   }
                }
            }
            break;
            
        case 2:
            k=nd0-1;
            mu=2;
                    for(t=0;t<def_vet[2];t++){
                        for(k=0;k<def_vet[1];k++){
                            for(j=0;j<def_vet[0];j++){
                                
                                parity = (i+j+k+t) % 2;
                                
                                //K_mu_values swap
                                if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                                    printf("beforrre (%d) %f %f",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );
                                    conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                                    conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                                    printf("aftermath (%d) %f %f\n",counter,conf1[mu].K.d[snum_acc(i,j,k,t)],conf2[mu].K.d[snum_acc(i,j,k,t)] );
                                }
                                if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                                    printf("beforrre (%d) %f %f",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );
                                    conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                                    conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                                     printf("aftermath (%d) %f %f\n",counter,conf1[mu+1].K.d[snum_acc(i,j,k,t)],conf2[mu+1].K.d[snum_acc(i,j,k,t)] );
                                }
                                
                                
                            }
                        }
                    }
                                
            break;
            
        case 3:
            t=nd0-1;
            mu=3;
            for(t=0;t<def_vet[2];t++){
                for(k=0;k<def_vet[1];k++){
                    for(j=0;j<def_vet[0];j++){
                        
                        parity = (i+j+k+t) % 2;
                        
                        //K_mu_values swap
                        if (parity==0){aux=conf1[mu].K.d[snum_acc(i,j,k,t)];
                            conf1[mu].K.d[snum_acc(i,j,k,t)]=conf2[mu].K.d[snum_acc(i,j,k,t)];
                            conf2[mu].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        if(parity!=0){aux=conf1[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf1[mu+1].K.d[snum_acc(i,j,k,t)]=conf2[mu+1].K.d[snum_acc(i,j,k,t)];
                            conf2[mu+1].K.d[snum_acc(i,j,k,t)]=aux;
                        }
                        
                    }
                }
            }
                        
            break;
            

        default:
            printf("ERROR WRONG AXIS CHOICE! (SWAP)\n");
            res=1;
            
    }
    
    
    //test
   /* printf("counter:%d\n%",counter);*/
    
 
    return res;
}
//replicas_swap function: 2 confs defects are exchanged.
int replicas_swap(su3_soa * conf1,su3_soa * conf2,int def_axis,int * def_vet ){
    vec3_soa  aux;
    int aux_label;
        int res=0;
    int mu=0;
    
    for(mu=0;mu<8;mu++){
        printf("beforrre (%d) %.18lf %.18lf",mu,creal(conf1[mu].r1.c1[snum_acc(31,6,6,6)]),creal(conf2[mu].r1.c1[snum_acc(31,6,6,6)]));
        
        //label swap.
        aux_label=conf1[mu].label;
        conf1[mu].label=conf2[mu].label;
        conf2[mu].label=aux_label;
        
        
        //swap r0
        aux=conf1[mu].r0;
        conf1[mu].r0=conf2[mu].r0;
        conf2[mu].r0=aux;
        
        
        
        //swap r1
        aux=conf1[mu].r1;
        conf1[mu].r1=conf2[mu].r1;
        conf2[mu].r1=aux;
        
        //swap r2
        aux=conf1[mu].r2;
        conf1[mu].r2=conf2[mu].r2;
        conf2[mu].r2=aux;
        
    printf("aftermath (%d) %.18lf %.18lf\n",mu,creal(conf1[mu].r1.c1[snum_acc(31,6,6,6)]),creal(conf2[mu].r1.c1[snum_acc(31,6,6,6)]));

    }
    
    return res;
}


//function which print the confs'labels.
int label_print(su3_soa ** conf_hasen, int replicas_number,FILE *file,int step_number){
    int res=0;
    int i;
    
    fprintf(file,"#step: &d\n",step_number);
    
    for(i=0;i<replicas_number;i++){
        fprintf(file,"%d   %d\n",i,conf_hasen[i][0].label);
        
    }
  
    
    
    return res;
}
