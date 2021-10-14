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
#include "./plaquettes.h"
#include "./su3_utilities.h"
#include "../Include/common_defines.h"
#include "./single_types.h"
#include "./rettangoli.h"
#include "./plaquettes.h"
#include "../Include/debug.h"
#include "../Rand/random.h"

#include <time.h>





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












int init_k(su3_soa * conf, double c_r, int def_axis, int * def_vet, defect_info * def,int defect_info_config){

  int mu;
  int i;
  int parity,condition,condition_2=1,condition_3=1;
  int def_vet_4d[4];
  int perp_dir[4][3] = {{ 1, 2, 3}, { 0, 2, 3}, { 0, 1, 3}, { 0, 1, 2}};

  int map[4] = {geom_par.xmap,geom_par.ymap,geom_par.zmap,geom_par.tmap};
  int def_axis_mapped=map[def_axis];
if(defect_info_config==0){  def->def_axis_mapped=def_axis_mapped;}


    //INITIALIZZATION
    int nd[4]={nd0-1,nd1-1,nd2-1,nd3-1};
    if(defect_info_config==0){   
    for(int nu=0; nu<4;nu++){
        for (i=0; i<3; i++)
        {
            def->defect_swap_max[nu][i]=0;
            def->defect_swap_min[nu][i]=nd[i];
	     #ifdef GAUGE_ACT_TLSM
	    for(int j=0;j<2;j++){
	     def->defect_swap_max_TLSM[nu][i][j]=0;
            def->defect_swap_min_TLSM[nu][i][j] = nd[i];
	    }
             #endif


        }
	  def->defect_swap_max[nu][3]=D3_HALO;
            def->defect_swap_min[nu][3]=nd[3]-D3_HALO;
             #ifdef GAUGE_ACT_TLSM
	     for(int j=0;j<2;j++){

             def->defect_swap_max_TLSM[nu][3][j]=D3_HALO;
            def->defect_swap_min_TLSM[nu][3][j] = nd[3]-D3_HALO;
	     }
	    #endif


    }
    }



  //NB: def_axis map: 0 -> x, 1 -> y, 2 -> z, 3 -> t
  //same for def_vet_4d.
  //it should be kept in mind: when you choose def_axis, the def_vet mapping follow that order.
  //So, if def_axis = 0, def_vet[0] is defect extent along y, def_vet[1(2)] along z(t).
  //if def_axis = 1, def_vet[0] is along x, def_vet[1(2)] along z(t), and so on...
  //This is independent on the physical axis map on the lattice.
  //If, in input, defect_boundary = 0(1-2-3), def_axis will be x(y-z-t) EVER.


  def_vet_4d[def_axis]=1;
  for(i=0;i<3;i++){
    def_vet_4d[perp_dir[def_axis][i]]=def_vet[i];
      if(defect_info_config==0){ def->def_mapped_perp_dir[i]=perp_dir[def_axis_mapped][i];}
  }

  int x_mind = (geom_par.xmap==def_axis)? geom_par.gnx-1 : 0;
  int y_mind = (geom_par.ymap==def_axis)? geom_par.gny-1 : 0;
  int z_mind = (geom_par.zmap==def_axis)? geom_par.gnz-1 : 0;
  int t_mind = (geom_par.tmap==def_axis)? geom_par.gnt-1 : 0;

  int x_maxd = x_mind + def_vet_4d[0];
  int y_maxd = y_mind + def_vet_4d[1];
  int z_maxd = z_mind + def_vet_4d[2];
  int t_maxd = t_mind + def_vet_4d[3];
  
  int tnx = geom_par.gnx;
  int tny = geom_par.gny;
  int tnz = geom_par.gnz;
  int tnt = geom_par.gnt;
  int count=0;
  
  int d[4], idxh, x, y, z, t;
  for(d[3]=0; d[3] < nd3; d[3]++)
    for(d[2]=0; d[2] < nd2; d[2]++)
      for(d[1]=0; d[1] < nd1; d[1]++)
	for(d[0]=0; d[0] < nd0; d[0]++){
	  
	  idxh = snum_acc(d[0],d[1],d[2],d[3]);
	  
	  x = d[geom_par.xmap];
	  y = d[geom_par.ymap];
	  z = d[geom_par.zmap];
	  t = d[geom_par.tmap];

#ifdef MULTIDEVICE
	  x+= devinfo.origin_0123[geom_par.xmap]        //x is now physical x-coordinate
	    - devinfo.halo_widths0123[geom_par.xmap];   //(for every MPI Rank).
	  y+= devinfo.origin_0123[geom_par.ymap]
	    - devinfo.halo_widths0123[geom_par.ymap];
	  z+= devinfo.origin_0123[geom_par.zmap]
	    - devinfo.halo_widths0123[geom_par.zmap];
	  t+= devinfo.origin_0123[geom_par.tmap]
	    - devinfo.halo_widths0123[geom_par.tmap];

	  if(x>tnx-1)  x-= tnx ; if(x<0)  x+= tnx ; //this guarantees boundary conditions.
	  if(y>tny-1)  y-= tny ; if(y<0)  y+= tny ; 
	  if(z>tnz-1)  z-= tnz ; if(z<0)  z+= tnz ; 
	  if(t>tnt-1)  t-= tnt ; if(t<0)  t+= tnt ; 

#endif                  

int	  parity_fis = (x+y+z+t)%2; 

int 	  parity_log=(d[0]+d[1]+d[2]+d[3])%2;
 parity= parity_log;
//  parity= parity_log; 

	  for(mu=0;mu<4;mu++)
	    conf[2*mu+parity].K.d[idxh] = 1;

	  condition = (x >= x_mind) && (y >= y_mind) && (z >= z_mind) && (t >= t_mind) &&
	    (x < x_maxd) && (y < y_maxd) && (z < z_maxd) && (t < t_maxd);
	  
#ifdef MULTIDEVICE
	  condition_2=(  d[0]<=nd[0]-devinfo.halo_widths0123[0] &&d[1]<=nd[1]-devinfo.halo_widths0123[1] &&  d[2]<=nd[2]-devinfo.halo_widths0123[2] &&  d[3]<=nd[3]-devinfo.halo_widths0123[3]);

#endif

	  if(condition){
printf("rank %d,parity %d parity_fis %d parity_log %d %d %d %d %d\n",devinfo.myrank,parity,parity_fis,parity_log,d[0],d[1],d[2],d[3]);
if(defect_info_config==0){
     for(int nu=0; nu<4;nu++){
          for (i=0; i<4; i++)
         {
		 #ifdef MULTIDEVICE
		 condition_3=(d[i]>=devinfo.halo_widths0123[i]);
                 #endif

              if ( d[i] <def->defect_swap_min[nu][i]  &&  condition_3  && condition_2 ){
                   def->defect_swap_min[nu][i] = d[i];
		   #ifdef GAUGE_ACT_TLSM
		    for(int j=0;j<2;j++){
                   def->defect_swap_min_TLSM[j][nu][i] = d[i];
		    }
                   #endif
                        count ++;
	      }
              
              if ( d[i] >  def->defect_swap_max[nu][i] && condition_3 && condition_2) {
                  def->defect_swap_max[nu][i] = d[i];
		  #ifdef GAUGE_ACT_TLSM
		   for(int j=0;j<2;j++){
                  def->defect_swap_max_TLSM[j][nu][i] = d[i];
		   }	
                    #endif
                         count ++;
	      }
	  } 

        }
    }
	    conf[2*def_axis_mapped+parity].K.d[idxh] = c_r;
	    printf("%d %d %d %d || %d %d %d %d\n",d[0],d[1],d[2],d[3],x,y,z,t);
	    

	  }


	

	  
 

	}//lattice loop

    if(defect_info_config==0){
     if(count==0){
              for(int nu=0; nu<4;nu++){
                      for (i=0; i<4; i++){
                  def->defect_swap_min[nu][i] = 0;
                 def->defect_swap_max[nu][i] = 0; 

                  #ifdef GAUGE_ACT_TLSM
                   for(int j=0;j<2;j++){
		  def->defect_swap_min_TLSM[j][nu][i] = 0;
                  def->defect_swap_max_TLSM[j][nu][i] = 0;
		   }
                  #endif
              }
          }

     }	

     else{
              for(int nu=0; nu<4;nu++){
                      for (i=0; i<4; i++){
                 def->defect_swap_max[nu][i] += 1;

                  #ifdef GAUGE_ACT_TLSM
 for(int j=0;j<2;j++){
                  def->defect_swap_max_TLSM[j][nu][i] += 1;
 }
                  #endif
              }
          }

     }
    


if(count!=0){

     def->defect_swap_min[0][0] += -1;
     def->defect_swap_min[1][1] += -1;
     def->defect_swap_min[2][2] += -1;
if(0==devinfo.myrank){ def->defect_swap_min[3][3] += -1;}
                   #ifdef GAUGE_ACT_TLSM
     for(i=0;i<4;i++){
  if(0==devinfo.myrank || def_axis_mapped!=3 ){ def->defect_swap_min_TLSM[1][i][def_axis_mapped] +=-1;//plaq 2x1
//this set mu one step back

     }
     }

 for(i=0;i<3;i++){
 if(0==devinfo.myrank || perp_dir[def_axis_mapped][i]!=3 ){ def->defect_swap_min_TLSM[1][perp_dir[def_axis_mapped][i]][perp_dir[def_axis_mapped][i]] +=-1 ;}
if(0==devinfo.myrank || perp_dir[def_axis_mapped][i]!=3 ){ def->defect_swap_min_TLSM[0][perp_dir[def_axis_mapped][i]][perp_dir[def_axis_mapped][i]] +=-2;}
}	





                   #endif

     


}

   printf("rank %d count: %d\n",devinfo.myrank,count);

 printf("defect_swap_max: rank %d\n",devinfo.myrank);

    for(int nu=0;nu<4;nu++){
	     if(nu!=def->def_axis_mapped){
            printf("nu:%d||",nu);
        for(i=0;i<4;i++){
    printf("%d||",def->defect_swap_max[nu][i]);}
    printf("\n");
        }
    }

    printf("defect_swap_min: rank %d\n",devinfo.myrank);

    for(int nu=0;nu<4;nu++){
	     if(nu!=def->def_axis_mapped){
        printf("nu:%d||",nu);
        for(i=0;i<4;i++){
            printf("%d||",def->defect_swap_min[nu][i]);}
    printf("\n");
        }
    }
#ifdef GAUGE_ACT_TLSM
  
 printf("defect_swap_max_tlsm 1x2: rank %d\n",devinfo.myrank);

    for(int nu=0;nu<4;nu++){
            printf("nu:%d||",nu);
        for(i=0;i<4;i++){
    printf("%d||",def->defect_swap_max_TLSM[0][nu][i]);

        }
        printf("\n");
    }
                     

    printf("defect_swap_min_tlsm 1x2: rank %d\n",devinfo.myrank);

    for(int nu=0;nu<4;nu++){
        printf("nu:%d||",nu);
        for(i=0;i<4;i++){
            printf("%d||",def->defect_swap_min_TLSM[0][nu][i]);

        }
        printf("\n");
    }
    

 printf("defect_swap_max_tlsm 2x1: rank %d\n",devinfo.myrank);

    for(int nu=0;nu<4;nu++){
            printf("nu:%d||",nu);
        for(i=0;i<4;i++){
    printf("%d||",def->defect_swap_max_TLSM[1][nu][i]);

        }
        printf("\n");
    }


    printf("defect_swap_min_tlsm 2x1: rank %d\n",devinfo.myrank);

    for(int nu=0;nu<4;nu++){
        printf("nu:%d||",nu);
        for(i=0;i<4;i++){
            printf("%d||",def->defect_swap_min_TLSM[1][nu][i]);

        }
        printf("\n");
    }



#endif

    }
  //printf("MPI%d: condition satisfied %d times.\n",devinfo.myrank,count);
}

/*
//function that tested the initialization of k_mu  by setting every value to one. //OUTADATED
int init_k_test(su3_soa *conf_acc,double c_r){
    int kk2=0;
    int mu1=0;
    

    for(mu1=0;mu1<8;mu1++){
        for(kk2=0;kk2<sizeh;kk2++){
        
            if(conf_acc[mu1].K.d[kk2]!=1 && conf_acc[mu1].K.d[kk2]!=c_r){
                return 1;
            }
            
            if(conf_acc[mu1].K.d[kk2]==c_r){
                printf("trovato\n");
                printf("%d:ku[%d]:%f\n",mu1,kk2,conf_acc[mu1].K.d[kk2]);}
        }
         printf("%d:ku[%d]:%f\n",mu1,kk2,conf_acc[mu1].K.d[kk2]);
        printf("%d\n",kk2);
        
    }
    return 0;
}
*/
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

int replicas_swap(su3_soa * conf1,su3_soa * conf2){
    vec3_soa  aux;
    int aux_label;
        int res=0;
    int mu=0;
    
    for(mu=0;mu<8;mu++){
       // printf("beforrre (%d) %.18lf %.18lf",mu,creal(conf1[mu].r1.c1[snum_acc(31,6,6,6)]),creal(conf2[mu].r1.c1[snum_acc(31,6,6,6)]));
        
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
        
   // printf("aftermath (%d) %.18lf %.18lf\n",mu,creal(conf1[mu].r1.c1[snum_acc(31,6,6,6)]),creal(conf2[mu].r1.c1[snum_acc(31,6,6,6)]));

    }
    
    return res;
}


//function which print the confs'labels.
int label_print(su3_soa ** conf_hasen, int replicas_number,FILE *file,int step_number){
    int res=0;
    int i;
    
    fprintf(file,"%d    ",step_number);
    
    for(i=0;i<replicas_number;i++){
        fprintf(file,"%d    ",conf_hasen[i][0].label);
        
    }
  
    fprintf(file,"\n");
    
    return res;
}



double  calc_Delta_S_soloopenacc_SWAP(
                                      __restrict  su3_soa * const tconf_acc,
                                      __restrict  su3_soa * const tconf_acc2,
                                      __restrict su3_soa * const local_plaqs,
                                      dcomplex_soa * const tr_local_plaqs,defect_info * def, int improved )
{
    double result=0.0;
    double total_result=0.0;
    int mu,nu;
    mu=def->def_axis_mapped;
    int counter=0;


    for(counter=0;counter<3;counter++){
        // sommo i 6 risultati in tempo
        nu=def->def_mapped_perp_dir[counter];
                printf("(%d,%d)\n",mu,nu);
            
                result  += calc_Delta_S_Wilson_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def); //here ol the plaquettes of a specific plane's choice are computed.
            
            printf("wilson's part done\n");

            if(improved==1){
		    result*=C_ZERO;
		       #ifdef GAUGE_ACT_TLSM
                result  +=C_ONE*calc_Delta_S_Symanzik_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def);
                          #endif
		  printf("symanzik's part done\n");

            }
            
 
 	
  }	
#ifdef MULTIDEVICE
        MPI_Allreduce((void*)&result,(void*)&total_result,
                      1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
        total_result = result;
#endif
        return total_result;
    
    
}


//Function which chooses the plane and iterate only on useful them.

/*
double  calc_Delta_S_soloopenacc_SWAP(
                                   __restrict  su3_soa * const tconf_acc,
                                        __restrict  su3_soa * const tconf_acc2,
                                   __restrict su3_soa * const local_plaqs,
                                   dcomplex_soa * const tr_local_plaqs,int def_axis, int * def_vet, int improved )
{
    
    
    double result=0.0;
    double total_result=0.0;
    int mu;
    //int def_vet_2[3];
    
    
    // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)//(the couple has to be chosen excluding same direction ones.(4 2)binomial coefficient.
    
    switch (def_axis) {
        case 0:
            mu=0;  //TEST MOD
           
            for(int nu=mu+1;nu<4;nu++){
                // sommo i 6 risultati in tempo
               // printf("(%d,%d)\n",mu,nu);
                 if(improved==0){
                     
                   //  if(nu==1){def_vet_2[0]=def_vet[0]+1; def_vet_2[1]=def_vet[1];def_vet_2[2]=def_vet[2];}
                   //  if(nu==2){def_vet_2[0]=def_vet[0]; def_vet_2[1]=def_vet[1]+1;def_vet_2[2]=def_vet[2];}
                   //  if(nu==3){def_vet_2[0]=def_vet[0]; def_vet_2[1]=def_vet[1];def_vet_2[2]=def_vet[2]+1;}
                      
                result  += calc_Delta_S_Wilson_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet); //here ol the plaquettes of a specific plane's choice are computed.
                 }
                
                
                if(improved==1){
                    result  += calc_Delta_S_Symanzik_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet);
                
                    
                }
                
                
            }
            
            
            break;
            
        case 1:
            mu=1;
            for(int nu=0;nu<4;nu++){
                // sommo i 6 risultati in tempo
                if(nu!=mu){
                 //     printf("(%d,%d)\n",mu,nu);
                     if(improved==0){
                result  += calc_Delta_S_Wilson_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet); //here ol the plaquettes of a specific plane's choice are computed.
                     }
                    
                    if(improved==1){
                        result  += calc_Delta_S_Symanzik_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet);
                        
                        
                    }
                
                }
                
                
            }
            
            
            break;
      
        case 2:
            mu=2;
            for(int nu=0;nu<4;nu++){
                // sommo i 6 risultati in tempo
                if(nu!=mu){
                   //   printf("(%d,%d)\n",mu,nu);
                    if(improved==0){
                    result  += calc_Delta_S_Wilson_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet); //here all the plaquettes of a specific plane's choice are computed.
                    }
                    if(improved==1){
                        result  += calc_Delta_S_Symanzik_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet);
                        
                        
                    }
                }
                
                
            }
            
            
            break;
            
        case 3:
            mu=3;
            for(int nu=0;nu<3;nu++){
             //     printf("(%d,%d)\n",mu,nu);
                // sommo i 6 risultati in tempo
                 if(improved==0){
                    result  += calc_Delta_S_Wilson_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet); //here all the plaquettes of a specific plane's choice are computed.
                 }
                if(improved==1){
                    result  += calc_Delta_S_Symanzik_SWAP(tconf_acc,tconf_acc2,local_plaqs,tr_local_plaqs,mu,nu,def_axis,def_vet);
                    
                }
                
            }
            
            
            break;
            
            
        default:
            printf("DELTA_S_SWAP ERROR!\n");
            break;
    }
 
    
    
    
#ifdef MULTIDEVICE
    MPI_Allreduce((void*)&result,(void*)&total_result,
                  1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
    total_result = result;
#endif
    return total_result;
    
}
*/


double calc_Delta_S_Wilson_SWAP(
                                __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                __restrict const su3_soa * const w,
                                __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                
                                const int mu, const int nu, defect_info * def)
{



    double K_mu_nu; //MOD.
    double K_mu_nu2; //MOD.
    int d0,d1,d2,d3;
  //  int D_max[4];
   //int D_min[4];
    int i;
  printf("nu: %d\n",nu);
/*
 #pragma acc  enter data create(D_max[0:4])
#pragma acc  enter data create(D_min[0:4])
#pragma acc update device(D_max[0:4])
#pragma acc update device(D_max[0:4])

for(i=0;i<4;i++){

    D_min[i]=def->defect_swap_min[nu][i];
    D_max[i]=def->defect_swap_max[nu][i];
    printf("%d | %d | %d \n",i,D_min[i],D_max[i]);
}

*/
int D0_min,D1_min,D2_min,D3_min;
int D0_max,D1_max,D2_max,D3_max;




D0_min=def->defect_swap_min[nu][0];
D1_min=def->defect_swap_min[nu][1];
D2_min=def->defect_swap_min[nu][2];
D3_min=def->defect_swap_min[nu][3];

D0_max=def->defect_swap_max[nu][0];
D1_max=def->defect_swap_max[nu][1];
D2_max=def->defect_swap_max[nu][2];
D3_max=def->defect_swap_max[nu][3];



    for(i=0; i<sizeh;i++){
        tr_local_plaqs[0].c[i]=0;
        tr_local_plaqs[1].c[i]=0;
        
    }
#pragma acc update device(tr_local_plaqs[0:2])
 





printf("|%d %d||%d %d| |%d %d| |%d %d| \n",D0_min, D0_max, D1_min, D1_max,D2_min, D2_max,D3_min, D3_max);


#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)/* present(D_max) present(D_min)*/
#pragma acc loop independent gang(STAPGANG3)
  

                   for(d3=D3_min; d3< D3_max; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
            for(d2=D2_min; d2< D2_max; d2++) {
                for(d1=D1_min; d1< D1_max; d1++) {
                    for(d0=D0_min; d0< D0_max; d0++) {
           
/*  
                   for(d3=D_min[3]; d3< D_max[3]; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
            for(d2=D_min[2]; d2< D_max[2]; d2++) {
                for(d1=D_min[1]; d1< D_max[1]; d1++) {
                    for(d0=D_min[0]; d0< D_max[0]; d0++) {
  */            
  

			    
                int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                int parity; //parity
                int dir_muA,dir_nuB; //mu and nu directions.
                int dir_muC,dir_nuD;
                
                idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                
                
			   
                if(d1==-1){parity = (d0+d2+d3) % 2; idxh=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];parity=!parity;}
                if(d2==-1){parity = (d0+d1+d3) % 2; idxh=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];parity=!parity;}
                if(d0==-1){parity = (d3+d1+d2) % 2; idxh=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];parity=!parity;}
                if(d3==-1){parity = (d0+d1+d2) % 2; idxh=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];parity=!parity;}
              
                
       
                
                dir_muA = 2*mu +  parity;
                dir_muC = 2*mu + !parity;
                idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                
                dir_nuB = 2*nu + !parity;
                dir_nuD = 2*nu +  parity;
                idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                //       r+nu (C)  r+mu+nu
                //          +<---+
                // nu       |    ^
                // ^    (D) V    | (B)
                // |        +--->+
                // |       r  (A)  r+mu
                // +---> mu
                
                //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                
                //plaquette u
                mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                
                
                d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                
                
                
                
                // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                //MOD
                
                //K_mu_nu computation;
                K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                
                
                
                //MOD_END
                
                //SECOND MOD
                
                
                //plaquette w
                mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                
                
                //K_mu_nu computation;
                K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                
                
                
                d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                
                
                tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                
                tr_local_plaqs[parity].c[idxh]=(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                
                
                
 /* 																																																              
    if(K_mu_nu2 !=1){                printf("( %d %d %d  %d parity %d idxh %d K_mu_nu %lf K_mu_nu %lf mu %d nu %d k_mu(conf1) %f k_mu(conf2) %f\n ",d0,d1,d2,d3,parity,idxh,K_mu_nu,K_mu_nu2,mu,nu,u[dir_muA].K.d[idxh],w[dir_muA].K.d[idxh]);           
    
    printf("tracce: conf0 %f conf 2%f\n",creal(tr_local_plaqs[parity].c[idxh]));
    }
   */             
                
                
	       } //d0  
                
            }  // d1
        }  // d2
        
    }  // d3
double res_R_p = 0.0;
double res_I_p = 0.0;
double resR = 0.0;
int t;



#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
for(t=0; t  <sizeh; t++) {
    res_R_p += creal(tr_local_plaqs[0].c[t]); //even sites plaquettes
    
    res_R_p += creal(tr_local_plaqs[1].c[t]); //odd sites plaquettes
}

									


res_R_p=BETA_BY_THREE *res_R_p;

return res_R_p;
    
    
}
#ifdef GAUGE_ACT_TLSM

double calc_Delta_S_Symanzik_SWAP(
                                  __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
                                  __restrict const su3_soa * const w,
                                  __restrict su3_soa * const loc_plaq, //la placchetta locale.
                                  dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                  const int mu, const int nu, defect_info * def)
{

    double K_mu_nu; //MOD.
    double K_mu_nu2; //MOD.
    int d0, d1, d2, d3;
    int D0_min,D1_min,D2_min,D3_min;
    int D0_max,D1_max,D2_max,D3_max;



    //min e max init
    D0_max=def->defect_swap_max_TLSM[1][nu][def->def_axis_mapped];
    D1_max=def->defect_swap_max_TLSM[1][nu][def->def_mapped_perp_dir[0]];
    D2_max=def->defect_swap_max_TLSM[1][nu][def->def_mapped_perp_dir[1]];
    D3_max=def->defect_swap_max_TLSM[1][nu][def->def_mapped_perp_dir[2]];


    D0_min=def->defect_swap_min_TLSM[1][nu][def->def_axis_mapped];
    D1_min=def->defect_swap_min_TLSM[1][nu][def->def_mapped_perp_dir[0]];
    D2_min=def->defect_swap_min_TLSM[1][nu][def->def_mapped_perp_dir[1]];
    D3_min=def->defect_swap_min_TLSM[1][nu][def->def_mapped_perp_dir[2]];



    int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
    int parity; //parity
    int dir_muA,dir_nuB; //mu and nu directions.
    int dir_muC,dir_nuD;

    //rectangular adjoints
    int dir_muB,dir_muD,dir_nuC;


    int dir_muE,dir_nuF,dir_nuE; //rectangular adjoints
    int idxpmupmu,idxpmupnu;//2x1
    int idxpnupnu; //1x2








    int is;
    for(is=0; is<sizeh;is++){
        tr_local_plaqs[0].c[is]=0;
        tr_local_plaqs[1].c[is]=0;

    }
#pragma acc update device(tr_local_plaqs[0:2])
    printf("2x1 |%d %d||%d %d| |%d %d| |%d %d| \n",D0_min, D0_max, D1_min, D1_max,D2_min, D2_max,D3_min, D3_max);

#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)


                   ////-----------2x1---------------------////


    for(d3=D3_min; d3< D3_max; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
        for(d2=D2_min; d2< D2_max; d2++) {
            for(d1=D1_min; d1< D1_max; d1++) {
                for(d0=D0_min; d0< D0_max; d0++) {




                    idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                    parity = (d0+d1+d2+d3) % 2; //obviously the parity_term



                    if(d1==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d2+d3) % 2; idxh=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];parity=!parity;}
                    if(d2==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d1+d3) % 2; idxh=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];parity=!parity;}
                    if(d0==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d3+d1+d2) % 2; idxh=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];parity=!parity;}
                    if(d3==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d1+d2) % 2; idxh=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];parity=!parity;}



                    ////-----------2x1---------------------////

                        //THIRD MOD
                        dir_muA = 2*mu +  parity;
                        dir_muB = 2*mu + !parity;
                        dir_nuC = 2*nu +  parity;
                        dir_muD = 2*mu +  parity;
                        dir_muE = 2*mu + !parity;
                        dir_nuF = 2*nu +  parity;
                        idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                        idxpmupmu = nnp_openacc[idxpmu][mu][!parity];// r+2mu
                        idxpmupnu = nnp_openacc[idxpmu][nu][!parity];// r+mu+nu
                        idxpnu = nnp_openacc[idxh][nu][parity];// r+nu

                        //       r+nu r+mu+nu r+2mu+nu
                        //          +<---+<---+
                        // nu       | (E) (D) ^
                        // ^    (F) V (A) (B) |  (C)
                        // |        +--->+--->+
                        // |       r   r+mu r+2mu
                        // +---> mu

                        //ret u
                        mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                        mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                        d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);


                        //K_mu_nu computation;
                        double K_mu_nu_RET;
                        K_mu_nu_RET=(u[dir_muA].K.d[idxh])*(u[dir_muB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupmu])*(u[dir_muD].K.d[idxpmupnu])*(u[dir_muE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);



                        //ret w

                        mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                        mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                        d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);


                        //K_mu_nu2 computation;
                        double K_mu_nu_RET2;
                        K_mu_nu_RET2=(w[dir_muA].K.d[idxh])*(w[dir_muB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupmu])*(w[dir_muD].K.d[idxpmupnu])*(w[dir_muE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);

                        //FINAL SUM

                        tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET-K_mu_nu_RET2)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);












                        //THIRD MOD END.



                }//d0
            }  // d1
        }  // d2
    }  // d3



    double res_R_p = 0.0;
    double res_I_p = 0.0;
    double resR = 0.0;
    int t;


#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
    for(t=0; t  < sizeh; t++) {
        res_R_p += creal(tr_local_plaqs[0].c[t]); //even sites plaquettes

        res_R_p += creal(tr_local_plaqs[1].c[t]); //odd sites plaquettes
    }






        ////-----------1x2---------------------////
    //min e max init
    
    int idxh1;
    D0_max=def->defect_swap_max_TLSM[0][nu][def->def_axis_mapped];
    D1_max=def->defect_swap_max_TLSM[0][nu][def->def_mapped_perp_dir[0]];
    D2_max=def->defect_swap_max_TLSM[0][nu][def->def_mapped_perp_dir[1]];
    D3_max=def->defect_swap_max_TLSM[0][nu][def->def_mapped_perp_dir[2]];


    D0_min=def->defect_swap_min_TLSM[0][nu][def->def_axis_mapped];
    D1_min=def->defect_swap_min_TLSM[0][nu][def->def_mapped_perp_dir[0]];
    D2_min=def->defect_swap_min_TLSM[0][nu][def->def_mapped_perp_dir[1]];
    D3_min=def->defect_swap_min_TLSM[0][nu][def->def_mapped_perp_dir[2]];





    for(is=0; is<sizeh;is++){
        tr_local_plaqs[0].c[is]=0;
        tr_local_plaqs[1].c[is]=0;

    }
#pragma acc update device(tr_local_plaqs[0:2])

 printf("1x2|%d %d||%d %d| |%d %d| |%d %d| \n",D0_min, D0_max, D1_min, D1_max,D2_min, D2_max,D3_min, D3_max);
   
 for(d3=D3_min; d3< D3_max; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
        for(d2=D2_min; d2< D2_max; d2++) {
            for(d1=D1_min; d1< D1_max; d1++) {
                for(d0=D0_min; d0< D0_max; d0++) {

   idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                    parity = (d0+d1+d2+d3) % 2; //obviously the parity_term



    if(d1==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d2+d3) % 2; idxh=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];parity=!parity;}
    if(d2==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d1+d3) % 2; idxh=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];parity=!parity;}
    if(d0==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d3+d1+d2) % 2; idxh=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];parity=!parity;}
    if(d3==-1){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d1+d2) % 2; idxh=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];parity=!parity;}

    int idxh1;

    if(d1==-2){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d2+d3) % 2; idxh1=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];idxh=nnm_openacc[idxh1][1][!parity]; }
    if(d2==-2){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d1+d3) % 2; idxh1=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];idxh=nnm_openacc[idxh1][2][!parity];}
    if(d0==-2){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d3+d1+d2) % 2; idxh1=nnm_openacc[snum_acc(0,d1,d2,d3)][0][parity];idxh=nnm_openacc[idxh1][0][!parity];}
    if(d3==-2){printf("%d %d %d %d\n",d0,d1,d2,d3);parity = (d0+d1+d2) % 2; idxh1=nnm_openacc[snum_acc(d0,d1,d2,0)][3][parity];idxh=nnm_openacc[idxh1][3][!parity]; }



                     ////-----------1x2---------------------////



                        //FOURTH MOD
                        dir_muA = 2*mu +  parity;
                        dir_nuB = 2*nu + !parity;
                        dir_nuC = 2*nu +  parity;
                        dir_muD = 2*mu +  parity;
                        dir_nuE = 2*nu + !parity;
                        dir_nuF = 2*nu +  parity;

                        idxpmu = nnp_openacc[idxh][mu][parity];      //r+mu
                        idxpmupnu = nnp_openacc[idxpmu][nu][!parity];//r+mu+nu
                        idxpnu = nnp_openacc[idxh][nu][parity];      //r+nu
                        idxpnupnu = nnp_openacc[idxpnu][nu][!parity];//r+nu+nu
                        //            (D)
                        //    r+2nu +<---+ r+mu+2nu
                        //          |    ^
                        //      (E) V    | (C)
                        //     r+nu +    + r+mu+nu
                        // nu       |    ^
                        // ^    (F) V    | (B)
                        // |        +--->+
                        // |       r  (A)  r+mu
                        // +---> mu

                        //ret u
                        mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                        mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                        d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);

                        //K_mu_nu computation;
                        double K_mu_nu_RET3;
                        K_mu_nu_RET3=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupnu])*(u[dir_muD].K.d[idxpnupnu])*(u[dir_nuE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);

                        //ret w
                        mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                        mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                        d_complex  ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);


                        //K_mu_nu computation;
                        double K_mu_nu_RET4;
                        K_mu_nu_RET4=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupnu])*(w[dir_muD].K.d[idxpnupnu])*(w[dir_nuE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);




                        //FINAL SUM

                        tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET3-K_mu_nu_RET4)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);



                        //FOURTH MOD END.



                }//d0
            }  // d1
        }  // d2
    }  // d3




#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
    for(t=0; t  < sizeh; t++) {
        res_R_p += creal(tr_local_plaqs[0].c[t]); //even sites plaquettes

        res_R_p += creal(tr_local_plaqs[1].c[t]); //odd sites plaquettes
    }


    res_R_p=BETA_BY_THREE *res_R_p;

    return res_R_p;





    return 0;
}
#endif


/*
double calc_Delta_S_Wilson_SWAP(
    __restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
    __restrict const su3_soa * const w,
__restrict su3_soa * const loc_plaq, //la placchetta locale.
    dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
                                         
    const int mu, const int nu, const int def_axis, const int * const def_vet)
    {
        double K_mu_nu; //MOD.
        double K_mu_nu2; //MOD.
        int d0, d1, d2, d3;
        int D1,D2,D3;
        int D0s,D1s,D2s,D3s;
        
        
        
        
        D1=def_vet[0];
        D2=def_vet[1];
        D3=def_vet[2];
        
        
        
        D0s=0;
        D1s=0;
        D2s=0;
        D3s=0;
        
        if(nu==0){D0s=-1;}
        
        
        if(nu==1){D1s=-1;}
        if(nu==2){D2s=-1;}
        if(nu==3){D3s=-PLAQ_EXTENT;}
    
        
    
        
       
       // printf("%d %d %d\n",def_vet[0],def_vet[1],def_vet[2]);
        //printf("%d %d %d %d\n",nd0,nd1,nd2,nd3);

      //  printf("%d %d %d %d\n",D0s,D1s,D2s,D3s);
        int is;

#pragma acc update self(tr_local_plaqs[0:2])
	 if(0==devinfo.myrank){
        for(is=0; is<sizeh;is++){
            tr_local_plaqs[0].c[is]=0;
            tr_local_plaqs[1].c[is]=0;
            
        }}
#pragma acc update device(tr_local_plaqs[0:2])
        
 switch (def_axis){
            case 0:
          d0=nd0-1;
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
         
                for(d3=D3s+D3_HALO; d3<D3-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d2=D2s; d2<D2; d2++) {
                        for(d1=D1s; d1<D1; d1++) {
                          //  for(d0=0;d0<nd0;d0++){  //TEST MOD
                        
                            
                        
                        int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                        int parity; //parity
                        int dir_muA,dir_nuB; //mu and nu directions.
                        int dir_muC,dir_nuD;

                        idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                        parity = (d0+d1+d2+(d3-D3_HALO)) % 2; //obviously the parity_term

                            
                            if(d1==-1){idxh=snum_acc(d0,nd1-1,d2,d3);parity=(d0+(nd1-1)+d2+d3)%2;}
                            if(d2==-1){idxh=snum_acc(d0,d1,nd2-1,d3);parity=(d0+d1+(nd2-1)+d3)%2;}
                            if(d3==-1){idxh=snum_acc(d0,d1,d2,nd3-1);parity=(d0+d1+d2+(nd3-1))%2;}
                            
                            
                            if(d1==-1){ parity = (d0+d2+d3) % 2; idxh=nnm_openacc[snum_acc(d0,0,d2,d3)][1][parity];parity=!parity;}
                            if(d2==-1){parity = (d0+d1+d3) % 2; idxh=nnm_openacc[snum_acc(d0,d1,0,d3)][2][parity];parity=!parity;}
                            if(d3==-1+D3_HALO){parity = (d0+d1+d2+D3_HALO) % 2; idxh=nnm_openacc[snum_acc(d0,d1,d2,D3_HALO)][3][parity];parity=!parity;}
                            
                       // idxh=nnm_openacc[idxh][nu][parity]; // the previous one. //MOD
                        
                       // parity = 1-parity;
                            
                        dir_muA = 2*mu +  parity;
                        dir_muC = 2*mu + !parity;
                        idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                        
                        dir_nuB = 2*nu + !parity;
                        dir_nuD = 2*nu +  parity;
                        idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                        //       r+nu (C)  r+mu+nu
                        //          +<---+
                        // nu       |    ^
                        // ^    (D) V    | (B)
                        // |        +--->+
                        // |       r  (A)  r+mu
                        // +---> mu
                        
                        //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                        
                            //plaquette u
                        mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                        mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                          
                        
                        d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                        d_complex aux = creal(ciao)+cimag(ciao)*I;
                    
                            
                            
                        
                        // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                        //MOD
                        
                        //K_mu_nu computation;
                        K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                        
                        
                       
                        //MOD_END
                            
                        //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D

                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                           

                             d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            aux+=-creal(ciao2)-cimag(ciao2)*I;
                            
                           tr_local_plaqs[parity].c[idxh]+=K_mu_nu-K_mu_nu2)aux;
                            
  	  
                   
                         
                          
                
                           
                            //SECOND_MOD_END
                            
                            
                          
                            
                    //        }//d0 TEST
                       
                        
                }  // d1
            }  // d2
                   
        }  // d3
        
                break;
               
            case 1:
          d1=nd1-1;
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
         
                for(d3=D3s+D3_HALO; d3<D3-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d2=D2s; d2<D2; d2++) {
                        for(d0=D0s; d0<D1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                           
                            
                            
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            
                            //third mod
                            if(d0==-1){idxh=snum_acc(nd0-1,d1,d2,d3);parity=((nd0-1)+d1+d2+d3)%2;}
                            if(d2==-1){idxh=snum_acc(d0,d1,nd2-1,d3);parity=(d0+d1+(nd2-1)+d3)%2;}
                            if(d3==-1){idxh=snum_acc(d0,d1,d2,nd3-1);parity=(d0+d1+d2+(nd3-1))%2;}
                            
                          //  idxh=nnm_openacc[idxh][nu][parity]; // the previous one. //MOD
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                            //MOD
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            
                            //MOD_END
                            
                            //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                            
                            
                            d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                            
                            tr_local_plaqs[parity].c[idxh]=(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                            
                            //  printf("DELTA_K_MU_NU= %f-%f= %f  (%d,%d,%d,%d) \n",K_mu_nu,K_mu_nu2,(K_mu_nu-K_mu_nu2),d0,d1,d2,d3);
                            
                            
                            
                            
                            
                             
                             //SECOND_MOD_END
                            
                            
                            
                        }  // d1
                    }  // d2
                }  // d3
                break;
                
            case 2:
            d2=nd2-1;
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
         
                for(d3=D3s+D3_HALO; d3<D3-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d1=D1s; d1<D2; d1++) {
                        for(d0=D0s; d0<D1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                           
                            
                            
                            //third mod
                            if(d0==-1){idxh=snum_acc(nd0-1,d1,d2,d3);parity=((nd0-1)+d1+d2+d3)%2;}
                            if(d1==-1){idxh=snum_acc(d0,nd1-1,d2,d3);parity=(d0+(nd1-1)+d2+d3)%2;}
                            if(d3==-1){idxh=snum_acc(d0,d1,d2,nd3-1);parity=(d0+d1+d2+(nd3-1))%2;}
                            
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                            //MOD
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            
 
                            
                            
                            //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                            
                            
                            d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                            
                            tr_local_plaqs[parity].c[idxh]=(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                            
                            //  printf("DELTA_K_MU_NU= %f-%f= %f  (%d,%d,%d,%d) \n",K_mu_nu,K_mu_nu2,(K_mu_nu-K_mu_nu2),d0,d1,d2,d3);
                            
                            
                            
                            
                            
                            //SECOND_MOD_END
                            
                            
                            
                            
                            
                            
                        }  // d1
                    }  // d2
                }  // d3
                break;
                
            case 3:
          d3=nd3-1-D3_HALO;
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
         
                for(d2=D2s; d2<D3; d2++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d1=D1s; d1<D2; d1++) {
                        for(d0=D0s; d0<D1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            
                            //third mod
                            if(d0==-1){idxh=snum_acc(nd0-1,d1,d2,d3);parity=((nd0-1)+d1+d2+d3)%2;}
                            if(d1==-1){idxh=snum_acc(d0,nd1-1,d2,d3);parity=(d0+(nd1-1)+d2+d3)%2;}
                            if(d2==-1){idxh=snum_acc(d0,d1,nd2-1,d3);parity=(d0+d1+(nd2-1)+nd3)%2;}
                            
                            
                            
                            
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                            //MOD
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            
                            //Mod_end
                            
                            
                            
                            //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                            
                            
                            d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                            
                            tr_local_plaqs[parity].c[idxh]=(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                            
                            //  printf("DELTA_K_MU_NU= %f-%f= %f  (%d,%d,%d,%d) \n",K_mu_nu,K_mu_nu2,(K_mu_nu-K_mu_nu2),d0,d1,d2,d3);
                            
                            
                            
                            
                            
                            //SECOND_MOD_END
                            
                            
                            
                            
                            
                        }  // d1
                    }  // d2
                }  // d3
         
         
                break;


                
                
                
            default:
                printf("DELTA_S_SWAP 1x1 plaquette ERROR!\n");
                break;

                
                
    }
        double res_R_p = 0.0;
        double res_I_p = 0.0;
        double resR = 0.0;
        int t;
        

        
#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
        for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
            res_R_p += creal(tr_local_plaqs[0].c[t]); //even sites plaquettes
            
            res_R_p += creal(tr_local_plaqs[1].c[t]); //odd sites plaquettes
        }
        
        
        res_R_p=BETA_BY_THREE *res_R_p;
        for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
           // printf("rank %d  site %d value %f  \n",devinfo.myrank,t,creal(tr_local_plaqs[0].c[t]));
     
        }




        return res_R_p;
    }// closes routine

*/



/*
double calc_Delta_S_Symanzik_SWAP(
__restrict const su3_soa * const u,//for an unknown reason the vet conf is called u. this is a vector odf su3_soa.
__restrict const su3_soa * const w,
__restrict su3_soa * const loc_plaq, //la placchetta locale.
dcomplex_soa * const tr_local_plaqs, //complex number that states the value of the trace. Of course is a vector of the struct dcomplex_soa.
const int mu, const int nu, defect_info *def)
{
    

    
        double K_mu_nu; //MOD.
        double K_mu_nu2; //MOD.
        int d0, d1, d2, d3;
        int D1,D2,D3;
        int D0s,D1s,D2s,D3s;
        
        
        
        
        D1=def_vet[0];
        D2=def_vet[1];
        D3=def_vet[2];
        
        
    
        D0s=0;
        D1s=0;
        D2s=0;
        D3s=0;
    
    
        if(nu==0){D0s=-2;}
        
        
        if(nu==1){D1s=-2;}
        if(nu==2){D2s=-2;}
        if(nu==3){D3s=-2;}
    
        
    
        
        
        // printf("%d %d %d\n",def_vet[0],def_vet[1],def_vet[2]);
        //printf("%d %d %d %d\n",nd0,nd1,nd2,nd3);
        
      //  printf("%d %d %d %d\n",D0s,D1s,D2s,D3s);
    //printf("%d %d %d \n",D1,D2,D3);
        int is;
        for(is=0; is<sizeh;is++){
            tr_local_plaqs[0].c[is]=0;
            tr_local_plaqs[1].c[is]=0;
            
        }
#pragma acc update device(tr_local_plaqs[0:2])
        
        switch (def_axis){
            case 0:
              //mod_test  d0=nd0-1;
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
                
                for(d3=D3s+D3_HALO; d3<D3-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d2=D2s; d2<D2; d2++) {
                        for(d1=D1s; d1<D1; d1++) {
                             for(d0=nd0-2;d0<nd0;d0++){
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            //rectangular adjoints
                            int dir_muB,dir_muD,dir_nuC;
                         
                            
                            int dir_muE,dir_nuF,dir_nuE; //rectangular adjoints
                            int idxpmupmu,idxpmupnu;//2x1
                            int idxpnupnu; //1x2
                            

                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            
                            if(d1==-1){idxh=snum_acc(d0,nd1-1,d2,d3);parity=(d0+(nd1-1)+d2+d3)%2;}
                            if(d2==-1){idxh=snum_acc(d0,d1,nd2-1,d3);parity=(d0+d1+(nd2-1)+d3)%2;}
                            if(d3==-1){idxh=snum_acc(d0,d1,d2,nd3-1);parity=(d0+d1+d2+(nd3-1))%2;}
                                 
                                 
                            if(d1==-2){idxh=snum_acc(d0,nd1-2,d2,d3);parity=(d0+(nd1-2)+d2+d3)%2;}
                            if(d2==-2){idxh=snum_acc(d0,d1,nd2-2,d3);parity=(d0+d1+(nd2-2)+d3)%2;}
                            if(d3==-2){idxh=snum_acc(d0,d1,d2,nd3-2);parity=(d0+d1+d2+(nd3-2))%2;}
                                 
                            
                            
                       
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            //plaquette u
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            
                            
                            
                            // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                            //MOD
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            
                            //MOD_END
                            
                            //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                            
                            
                            d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                            
                            tr_local_plaqs[parity].c[idxh]=C_ZERO*(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                            
                            
                            
                            
                            
                            
                            
                            //SECOND_MOD_END
                            
                            
                            
                                 if(d1>-2 && d2>-2 && d3>-2){
                             //THIRD MOD
                            dir_muA = 2*mu +  parity;
                            dir_muB = 2*mu + !parity;
                            dir_nuC = 2*nu +  parity;
                            dir_muD = 2*mu +  parity;
                            dir_muE = 2*mu + !parity;
                            dir_nuF = 2*nu +  parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            idxpmupmu = nnp_openacc[idxpmu][mu][!parity];// r+2mu
                            idxpmupnu = nnp_openacc[idxpmu][nu][!parity];// r+mu+nu
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
                            
                            //       r+nu r+mu+nu r+2mu+nu
                            //          +<---+<---+
                            // nu       | (E) (D) ^
                            // ^    (F) V (A) (B) |  (C)
                            // |        +--->+--->+
                            // |       r   r+mu r+2mu
                            // +---> mu
                            
                            //ret u
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                            mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                            ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            //K_mu_nu computation;
                            double K_mu_nu_RET;
                            K_mu_nu_RET=(u[dir_muA].K.d[idxh])*(u[dir_muB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupmu])*(u[dir_muD].K.d[idxpmupnu])*(u[dir_muE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                            
                            
                            
                            //ret w
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                            mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                            ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            //K_mu_nu2 computation;
                            double K_mu_nu_RET2;
                            K_mu_nu_RET2=(w[dir_muA].K.d[idxh])*(w[dir_muB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupmu])*(w[dir_muD].K.d[idxpmupnu])*(w[dir_muE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                            
                            //FINAL SUM
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET-K_mu_nu_RET2)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                            
                            
                            
                            
                           
                            
                            
                            
                            
                            
                            
                            //THIRD MOD END.
                                 }
                            
                                 if(d0!=nd0-2){
                            //FOURTH MOD
                            dir_muA = 2*mu +  parity;
                            dir_nuB = 2*nu + !parity;
                            dir_nuC = 2*nu +  parity;
                            dir_muD = 2*mu +  parity;
                            dir_nuE = 2*nu + !parity;
                            dir_nuF = 2*nu +  parity;
                            
                            idxpmu = nnp_openacc[idxh][mu][parity];      //r+mu
                            idxpmupnu = nnp_openacc[idxpmu][nu][!parity];//r+mu+nu
                            idxpnu = nnp_openacc[idxh][nu][parity];      //r+nu
                            idxpnupnu = nnp_openacc[idxpnu][nu][!parity];//r+nu+nu
                            //            (D)
                            //    r+2nu +<---+ r+mu+2nu
                            //          |    ^
                            //      (E) V    | (C)
                            //     r+nu +    + r+mu+nu
                            // nu       |    ^
                            // ^    (F) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //ret u
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                            mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                            ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            //K_mu_nu computation;
                            double K_mu_nu_RET3;
                            K_mu_nu_RET3=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupnu])*(u[dir_muD].K.d[idxpnupnu])*(u[dir_nuE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                            
                            //ret w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                            mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                            ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            //K_mu_nu computation;
                            double K_mu_nu_RET4;
                            K_mu_nu_RET4=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupnu])*(w[dir_muD].K.d[idxpnupnu])*(w[dir_nuE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                            
                            
                            
                            
                            //FINAL SUM
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET3-K_mu_nu_RET4)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                            
                            
                            
                            
                            
                            
                            
                            //FOURTH MOD END.
                                 }
                             }//test_mod
                            
                        }  // d1
                    }  // d2
                }  // d3
                
                break;
                
            case 1:
                
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
                
                for(d3=D3s+D3_HALO; d3<D3-D3_HALO; d3++) {//what?
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d2=D2s; d2<D2; d2++) {
                        for(d1=nd1-2;d1<nd1;d1++){
                            
                          for(d0=D0s; d0<D1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                              //rectangular adjoints
                              int dir_muB,dir_muD,dir_nuC;
                              
                              
                              int dir_muE,dir_nuF,dir_nuE; //rectangular adjoints
                              int idxpmupmu,idxpmupnu;//2x1
                              int idxpnupnu; //1x2
                              
                              
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            
                            
                            
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            
                            //third mod
                            if(d0==-1){idxh=snum_acc(nd0-1,d1,d2,d3);parity=((nd0-1)+d1+d2+d3)%2;}
                            if(d2==-1){idxh=snum_acc(d0,d1,nd2-1,d3);parity=(d0+d1+(nd2-1)+d3)%2;}
                            if(d3==-1){idxh=snum_acc(d0,d1,d2,nd3-1);parity=(d0+d1+d2+(nd3-1))%2;}
                              
                            if(d0==-2){idxh=snum_acc(nd0-2,d1,d2,d3);parity=((nd0-2)+d1+d2+d3)%2;}
                            if(d2==-2){idxh=snum_acc(d0,d1,nd2-2,d3);parity=(d0+d1+(nd2-2)+d3)%2;}
                            if(d3==-2){idxh=snum_acc(d0,d1,d2,nd3-2);parity=(d0+d1+d2+(nd3-2))%2;}
                              
                            
                            
                              
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                            //MOD
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            
                            //MOD_END
                            
                            //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                            
                            
                            d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                            
                            tr_local_plaqs[parity].c[idxh]=C_ZERO*(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                            
                            //  printf("DELTA_K_MU_NU= %f-%f= %f  (%d,%d,%d,%d) \n",K_mu_nu,K_mu_nu2,(K_mu_nu-K_mu_nu2),d0,d1,d2,d3);
                            
                            
                              
                            
                            
                            //SECOND_MOD_END
                            
                              if( d2>-2 && d3>-2 && d0>-2){
                                  //THIRD MOD
                                  dir_muA = 2*mu +  parity;
                                  dir_muB = 2*mu + !parity;
                                  dir_nuC = 2*nu +  parity;
                                  dir_muD = 2*mu +  parity;
                                  dir_muE = 2*mu + !parity;
                                  dir_nuF = 2*nu +  parity;
                                  idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                                  idxpmupmu = nnp_openacc[idxpmu][mu][!parity];// r+2mu
                                  idxpmupnu = nnp_openacc[idxpmu][nu][!parity];// r+mu+nu
                                  idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
                                  
                                  //       r+nu r+mu+nu r+2mu+nu
                                  //          +<---+<---+
                                  // nu       | (E) (D) ^
                                  // ^    (F) V (A) (B) |  (C)
                                  // |        +--->+--->+
                                  // |       r   r+mu r+2mu
                                  // +---> mu
                                  
                                  //ret u
                                  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                  mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                  ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                  
                                  
                                  //K_mu_nu computation;
                                  double K_mu_nu_RET;
                                  K_mu_nu_RET=(u[dir_muA].K.d[idxh])*(u[dir_muB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupmu])*(u[dir_muD].K.d[idxpmupnu])*(u[dir_muE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                                  
                                  
                                  
                                  //ret w
                                  
                                  mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                  mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                  ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                  
                                  
                                  //K_mu_nu2 computation;
                                  double K_mu_nu_RET2;
                                  K_mu_nu_RET2=(w[dir_muA].K.d[idxh])*(w[dir_muB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupmu])*(w[dir_muD].K.d[idxpmupnu])*(w[dir_muE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                                  
                                  //FINAL SUM
                                  
                                  tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET-K_mu_nu_RET2)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  //THIRD MOD END.
                              }
                              
                              if(d1!=nd1-2){
                                  //FOURTH MOD
                                  dir_muA = 2*mu +  parity;
                                  dir_nuB = 2*nu + !parity;
                                  dir_nuC = 2*nu +  parity;
                                  dir_muD = 2*mu +  parity;
                                  dir_nuE = 2*nu + !parity;
                                  dir_nuF = 2*nu +  parity;
                                  
                                  idxpmu = nnp_openacc[idxh][mu][parity];      //r+mu
                                  idxpmupnu = nnp_openacc[idxpmu][nu][!parity];//r+mu+nu
                                  idxpnu = nnp_openacc[idxh][nu][parity];      //r+nu
                                  idxpnupnu = nnp_openacc[idxpnu][nu][!parity];//r+nu+nu
                                  //            (D)
                                  //    r+2nu +<---+ r+mu+2nu
                                  //          |    ^
                                  //      (E) V    | (C)
                                  //     r+nu +    + r+mu+nu
                                  // nu       |    ^
                                  // ^    (F) V    | (B)
                                  // |        +--->+
                                  // |       r  (A)  r+mu
                                  // +---> mu
                                  
                                  //ret u
                                  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                  mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                  ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                  
                                  //K_mu_nu computation;
                                  double K_mu_nu_RET3;
                                  K_mu_nu_RET3=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupnu])*(u[dir_muD].K.d[idxpnupnu])*(u[dir_nuE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                                  
                                  //ret w
                                  mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                  mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                  ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                  
                                  
                                  //K_mu_nu computation;
                                  double K_mu_nu_RET4;
                                  K_mu_nu_RET4=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupnu])*(w[dir_muD].K.d[idxpnupnu])*(w[dir_nuE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                                  
                                  
                                  
                                  
                                  //FINAL SUM
                                  
                                  tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET3-K_mu_nu_RET4)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  //FOURTH MOD END.
                            
                            
                         
                            
                              }
                            
                            
                             }//d0
                        }  // d1
                    }  // d2
                }  // d3
                break;
                
            case 2:
                
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
                
                for(d3=D3s+D3_HALO; d3<D3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d2=nd2-2;d2<nd2;d2++){
                        for(d1=D1s; d1<D2; d1++) {
                            for(d0=D0s; d0<D1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                                //rectangular adjoints
                                int dir_muB,dir_muD,dir_nuC;
                                
                                
                                int dir_muE,dir_nuF,dir_nuE; //rectangular adjoints
                                int idxpmupmu,idxpmupnu;//2x1
                                int idxpnupnu; //1x2
                                
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            
                            
                            
                            //third mod
                            if(d0==-1){idxh=snum_acc(nd0-1,d1,d2,d3);parity=((nd0-1)+d1+d2+d3)%2;}
                            if(d1==-1){idxh=snum_acc(d0,nd1-1,d2,d3);parity=(d0+(nd1-1)+d2+d3)%2;}
                            if(d3==-1){idxh=snum_acc(d0,d1,d2,nd3-1);parity=(d0+d1+d2+(nd3-1))%2;}
                                
                            if(d0==-2){idxh=snum_acc(nd0-2,d1,d2,d3);parity=((nd0-2)+d1+d2+d3)%2;}
                            if(d1==-2){idxh=snum_acc(d0,nd1-2,d2,d3);parity=(d0+(nd1-2)+d2+d3)%2;}
                            if(d3==-2){idxh=snum_acc(d0,d1,d2,nd3-2);parity=(d0+d1+d2+(nd3-2))%2;}
                            
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                            //MOD
                            
                                
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            
                            
                            //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                            
                            
                            d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                            
                            tr_local_plaqs[parity].c[idxh]=C_ZERO*(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                            
                            //  printf("DELTA_K_MU_NU= %f-%f= %f  (%d,%d,%d,%d) \n",K_mu_nu,K_mu_nu2,(K_mu_nu-K_mu_nu2),d0,d1,d2,d3);
                            
                            
                            
                            
                            
                            //SECOND_MOD_END
                            
                                if(d1>-2 && d0>-2 && d3>-2){
                                    //THIRD MOD
                                    dir_muA = 2*mu +  parity;
                                    dir_muB = 2*mu + !parity;
                                    dir_nuC = 2*nu +  parity;
                                    dir_muD = 2*mu +  parity;
                                    dir_muE = 2*mu + !parity;
                                    dir_nuF = 2*nu +  parity;
                                    idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                                    idxpmupmu = nnp_openacc[idxpmu][mu][!parity];// r+2mu
                                    idxpmupnu = nnp_openacc[idxpmu][nu][!parity];// r+mu+nu
                                    idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
                                    
                                    //       r+nu r+mu+nu r+2mu+nu
                                    //          +<---+<---+
                                    // nu       | (E) (D) ^
                                    // ^    (F) V (A) (B) |  (C)
                                    // |        +--->+--->+
                                    // |       r   r+mu r+2mu
                                    // +---> mu
                                    
                                    //ret u
                                    mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                    mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                    ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                    
                                    
                                    //K_mu_nu computation;
                                    double K_mu_nu_RET;
                                    K_mu_nu_RET=(u[dir_muA].K.d[idxh])*(u[dir_muB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupmu])*(u[dir_muD].K.d[idxpmupnu])*(u[dir_muE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                                    
                                    
                                    
                                    //ret w
                                    
                                    mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                    mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                    ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                    
                                    
                                    //K_mu_nu2 computation;
                                    double K_mu_nu_RET2;
                                    K_mu_nu_RET2=(w[dir_muA].K.d[idxh])*(w[dir_muB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupmu])*(w[dir_muD].K.d[idxpmupnu])*(w[dir_muE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                                    
                                    //FINAL SUM
                                    
                                    tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET-K_mu_nu_RET2)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    //THIRD MOD END.
                                }
                                
                                if(d2!=nd2-2){
                                    //FOURTH MOD
                                    dir_muA = 2*mu +  parity;
                                    dir_nuB = 2*nu + !parity;
                                    dir_nuC = 2*nu +  parity;
                                    dir_muD = 2*mu +  parity;
                                    dir_nuE = 2*nu + !parity;
                                    dir_nuF = 2*nu +  parity;
                                    
                                    idxpmu = nnp_openacc[idxh][mu][parity];      //r+mu
                                    idxpmupnu = nnp_openacc[idxpmu][nu][!parity];//r+mu+nu
                                    idxpnu = nnp_openacc[idxh][nu][parity];      //r+nu
                                    idxpnupnu = nnp_openacc[idxpnu][nu][!parity];//r+nu+nu
                                    //            (D)
                                    //    r+2nu +<---+ r+mu+2nu
                                    //          |    ^
                                    //      (E) V    | (C)
                                    //     r+nu +    + r+mu+nu
                                    // nu       |    ^
                                    // ^    (F) V    | (B)
                                    // |        +--->+
                                    // |       r  (A)  r+mu
                                    // +---> mu
                                    
                                    //ret u
                                    mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                    mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                    ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                    
                                    //K_mu_nu computation;
                                    double K_mu_nu_RET3;
                                    K_mu_nu_RET3=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupnu])*(u[dir_muD].K.d[idxpnupnu])*(u[dir_nuE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                                    
                                    //ret w
                                    mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                    mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                    mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                    ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                    
                                    
                                    //K_mu_nu computation;
                                    double K_mu_nu_RET4;
                                    K_mu_nu_RET4=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupnu])*(w[dir_muD].K.d[idxpnupnu])*(w[dir_nuE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                                    
                                    
                                    
                                    
                                    //FINAL SUM
                                    
                                    tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET3-K_mu_nu_RET4)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    //FOURTH MOD END.
                            
                                }
                            }//d0
                            
                        }  // d1
                    }  // d2
                }  // d3
                break;
                
            case 3:
             
#pragma acc kernels present(u) present(w) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
        for(d3=nd3-2-D3_HALO;d3<nd3-D3_HALO;d3++){
            for(d2=D2s; d2<D3; d2++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
                    for(d1=D1s; d1<D2; d1++) {
                        for(d0=D0s; d0<D1; d0++) {
                            
                            
                            
                            int idxh,idxpmu,idxpnu; //idxh is the half-lattice position, idxpmu and idxpnu the nearest neighbours.
                            int parity; //parity
                            int dir_muA,dir_nuB; //mu and nu directions.
                            int dir_muC,dir_nuD;
                            
                            //rectangular adjoints
                            int dir_muB,dir_muD,dir_nuC;
                            
                            
                            int dir_muE,dir_nuF,dir_nuE; //rectangular adjoints
                            int idxpmupmu,idxpmupnu;//2x1
                            int idxpnupnu; //1x2
                            
                            
                            idxh = snum_acc(d0,d1,d2,d3);// the site on the  half-lattice.
                            parity = (d0+d1+d2+d3) % 2; //obviously the parity_term
                            
                            //third mod
                            if(d0==-1){idxh=snum_acc(nd0-1,d1,d2,d3);parity=((nd0-1)+d1+d2+d3)%2;}
                            if(d1==-1){idxh=snum_acc(d0,nd1-1,d2,d3);parity=(d0+(nd1-1)+d2+d3)%2;}
                            if(d2==-1){idxh=snum_acc(d0,d1,nd2-1,d3);parity=(d0+d1+(nd2-1)+nd3)%2;}
                            
                            if(d0==-2){idxh=snum_acc(nd0-2,d1,d2,d3);parity=((nd0-2)+d1+d2+d3)%2;}
                            if(d1==-2){idxh=snum_acc(d0,nd1-2,d2,d3);parity=(d0+(nd1-2)+d2+d3)%2;}
                            if(d2==-2){idxh=snum_acc(d0,d1,nd2-2,d3);parity=(d0+d1+(nd2-2)+nd3)%2;}
                            
                            
                            
                            
                            dir_muA = 2*mu +  parity;
                            dir_muC = 2*mu + !parity;
                            idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                            
                            dir_nuB = 2*nu + !parity;
                            dir_nuD = 2*nu +  parity;
                            idxpnu = nnp_openacc[idxh][nu][parity];// r+nu //the table that states which is the nearest neighbour.
                            //       r+nu (C)  r+mu+nu
                            //          +<---+
                            // nu       |    ^
                            // ^    (D) V    | (B)
                            // |        +--->+
                            // |       r  (A)  r+mu
                            // +---> mu
                            
                            //(&u[dir_muA] & &u[dir_nuB] States which part of the the conf will be used. It is important to pass them as pointer, cause loc_plaq has to be modified.
                            
                            mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            d_complex ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            tr_local_plaqs[parity].c[idxh] = creal(ciao)+cimag(ciao)*I;
                            
                            // printf("%f +i%f ||",creal(tr_local_plaqs[parity].c[idxh]),cimag(tr_local_plaqs[parity].c[idxh])*I);
                            //MOD
                            
                            //K_mu_nu computation;
                            K_mu_nu=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_muC].K.d[idxpnu])*(u[dir_nuD].K.d[idxh]);
                            
                            
                            
                            //Mod_end
                            
                            
                            
                            //SECOND MOD
                            
                            
                            //plaquette w
                            mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_PLAQ = A * B
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muC],idxpnu);              // LOC_PLAQ = LOC_PLAQ * C
                            mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
                            
                            
                            //K_mu_nu computation;
                            K_mu_nu2=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_muC].K.d[idxpnu])*(w[dir_nuD].K.d[idxh]);
                            
                            
                            
                            d_complex ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                            
                            
                            tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]-creal(ciao2)-cimag(ciao2)*I;
                            
                            tr_local_plaqs[parity].c[idxh]=C_ZERO*(K_mu_nu-K_mu_nu2)*tr_local_plaqs[parity].c[idxh];
                            
                            //  printf("DELTA_K_MU_NU= %f-%f= %f  (%d,%d,%d,%d) \n",K_mu_nu,K_mu_nu2,(K_mu_nu-K_mu_nu2),d0,d1,d2,d3);
                            
                            
                            
                            
                            
                            //SECOND_MOD_END
                            
                            if(d1>-2 && d2>-2 && d0>-2){
                                //THIRD MOD
                                dir_muA = 2*mu +  parity;
                                dir_muB = 2*mu + !parity;
                                dir_nuC = 2*nu +  parity;
                                dir_muD = 2*mu +  parity;
                                dir_muE = 2*mu + !parity;
                                dir_nuF = 2*nu +  parity;
                                idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
                                idxpmupmu = nnp_openacc[idxpmu][mu][!parity];// r+2mu
                                idxpmupnu = nnp_openacc[idxpmu][nu][!parity];// r+mu+nu
                                idxpnu = nnp_openacc[idxh][nu][parity];// r+nu
                                
                                //       r+nu r+mu+nu r+2mu+nu
                                //          +<---+<---+
                                // nu       | (E) (D) ^
                                // ^    (F) V (A) (B) |  (C)
                                // |        +--->+--->+
                                // |       r   r+mu r+2mu
                                // +---> mu
                                
                                //ret u
                                mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                
                                
                                //K_mu_nu computation;
                                double K_mu_nu_RET;
                                K_mu_nu_RET=(u[dir_muA].K.d[idxh])*(u[dir_muB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupmu])*(u[dir_muD].K.d[idxpmupnu])*(u[dir_muE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                                
                                
                                
                                //ret w
                                
                                mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                
                                
                                //K_mu_nu2 computation;
                                double K_mu_nu_RET2;
                                K_mu_nu_RET2=(w[dir_muA].K.d[idxh])*(w[dir_muB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupmu])*(w[dir_muD].K.d[idxpmupnu])*(w[dir_muE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                                
                                //FINAL SUM
                                
                                tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET-K_mu_nu_RET2)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                //THIRD MOD END.
                            }
                            
                            if(d3!=nd3-2){
                                //FOURTH MOD
                                dir_muA = 2*mu +  parity;
                                dir_nuB = 2*nu + !parity;
                                dir_nuC = 2*nu +  parity;
                                dir_muD = 2*mu +  parity;
                                dir_nuE = 2*nu + !parity;
                                dir_nuF = 2*nu +  parity;
                                
                                idxpmu = nnp_openacc[idxh][mu][parity];      //r+mu
                                idxpmupnu = nnp_openacc[idxpmu][nu][!parity];//r+mu+nu
                                idxpnu = nnp_openacc[idxh][nu][parity];      //r+nu
                                idxpnupnu = nnp_openacc[idxpnu][nu][!parity];//r+nu+nu
                                //            (D)
                                //    r+2nu +<---+ r+mu+2nu
                                //          |    ^
                                //      (E) V    | (C)
                                //     r+nu +    + r+mu+nu
                                // nu       |    ^
                                // ^    (F) V    | (B)
                                // |        +--->+
                                // |       r  (A)  r+mu
                                // +---> mu
                                
                                //ret u
                                mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                ciao = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                
                                //K_mu_nu computation;
                                double K_mu_nu_RET3;
                                K_mu_nu_RET3=(u[dir_muA].K.d[idxh])*(u[dir_nuB].K.d[idxpmu])*(u[dir_nuC].K.d[idxpmupnu])*(u[dir_muD].K.d[idxpnupnu])*(u[dir_nuE].K.d[idxpnu])*(u[dir_nuF].K.d[idxh]);
                                
                                //ret w
                                mat1_times_mat2_into_mat3_absent_stag_phases(&w[dir_muA],idxh,&w[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
                                mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E
                                mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&w[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
                                ciao2 = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);
                                
                                
                                //K_mu_nu computation;
                                double K_mu_nu_RET4;
                                K_mu_nu_RET4=(w[dir_muA].K.d[idxh])*(w[dir_nuB].K.d[idxpmu])*(w[dir_nuC].K.d[idxpmupnu])*(w[dir_muD].K.d[idxpnupnu])*(w[dir_nuE].K.d[idxpnu])*(w[dir_nuF].K.d[idxh]);
                                
                                
                                
                                
                                //FINAL SUM
                                
                                tr_local_plaqs[parity].c[idxh]=tr_local_plaqs[parity].c[idxh]+C_ONE*(K_mu_nu_RET3-K_mu_nu_RET4)*( creal(ciao)+cimag(ciao)*I-creal(ciao2)-cimag(ciao2)*I);
                                
                                
                                
                                
                                
                                
                                
                                //FOURTH MOD END.
                            }
                             }//d0
                            
                        }  // d1
                    }  // d2
                }  // d3
                
                
                break;
                
                
                
                
                
            default:
                printf("DELTA_S_SWAP 2x1/1x2 plaquette ERROR!\n");
                break;
                
                
                
        }
        double res_R_p = 0.0;
        double res_I_p = 0.0;
        double resR = 0.0;
        int t;
        
        
        
#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
        for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
            res_R_p += creal(tr_local_plaqs[0].c[t]); //even sites plaquettes
            
            res_R_p += creal(tr_local_plaqs[1].c[t]); //odd sites plaquettes
        }
        
        
        res_R_p=BETA_BY_THREE *res_R_p;
        
        return res_R_p;
    
    

    
    
    printf("d\n");
    return 0;
}
*/























int metro_SWAP(su3_soa ** conf_hasenbusch,
               __restrict su3_soa * const loc_plaq, //la placchetta locale.
               dcomplex_soa * const tr_local_plaqs,
            int rep_indx1, int rep_indx2,defect_info * def )
{
    int gauge_param;
    
    gauge_param=GAUGE_ACTION;
    
    double p1,p2;
    double Delta_S_SWAP;
    int accettata=0;
    Delta_S_SWAP=calc_Delta_S_soloopenacc_SWAP(conf_hasenbusch[rep_indx1],conf_hasenbusch[rep_indx2],loc_plaq,tr_local_plaqs,def,gauge_param);
    printf("DELTA_SWAP:%f\n",Delta_S_SWAP);
    if(Delta_S_SWAP<0){
        accettata=1;
        printf(" p1 p2 :%f %f\n",p1,p2);

    }
    else
    {  p1=exp(-Delta_S_SWAP);
        if(debug_settings.do_norandom_test) p2=0; // NORANDOM
        else{   // NORMAL, RANDOM
            if(0==devinfo.myrank)p2=casuale();
       
#ifdef MULTIDEVICE
                MPI_Bcast((void*) &p2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                printf("MPI%02d p2 : %f, p1 %f \n",devinfo.myrank, p2,p1);
#endif       

	}
        if(p2<p1)
        {
        accettata=1;
            
        }
        
        else
        {
        accettata=0;
        // configuration reject
            
        }
    printf(" p1 p2 :%f %f\n",p1,p2);
    }
    
    
    
    if (accettata==1){
     replicas_swap(conf_hasenbusch[rep_indx1],conf_hasenbusch[rep_indx2]);
    }
   
    

    
    
    
    
    return accettata;
}


void All_Conf_SWAP(su3_soa ** conf_hasenbusch,
                   __restrict su3_soa * const loc_plaq, //la placchetta locale.
                   dcomplex_soa * const tr_local_plaqs,
                   
                   int replicas_number,  defect_info * def, int* swap_num,int * all_swap_vet,int * acceptance_vet ){
    double swap_order;

   if(0==devinfo.myrank){swap_order=casuale();}

    #ifdef MULTIDEVICE
                MPI_Bcast((void*) &swap_order,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                printf("MPI%02d swap_order : %f \n",devinfo.myrank,swap_order);
#endif


    int accettata=0;
    int i_counter, j_counter;
 
   
    if(swap_order<=0.5){
        
        
        for(i_counter=0;i_counter<replicas_number-1;i_counter++){
            printf("%d %d\n",i_counter,i_counter+1);
            
            accettata=metro_SWAP( conf_hasenbusch,loc_plaq,tr_local_plaqs, i_counter, i_counter+1, def);
            #pragma acc update device(conf_hasenbusch[0:replicas_number][0:8])
            
          
          
                
                printf("proposed: all_swap_vet %d\n",all_swap_vet[i_counter]);
                
            *swap_num=*swap_num+1;
           
            all_swap_vet[i_counter]++;
            //all_swap_vet[i_counter+1]++;
            printf("proposed: all_swap_vet %d\n",all_swap_vet[i_counter]);
            
            if (accettata==1){
                acceptance_vet[i_counter]++;
             //   acceptance_vet[i_counter+1]++;
            }
         
            
        }
        
        

    }
    else{
        
        for(i_counter=0;i_counter<replicas_number-1;i_counter++){
        printf("%d %d\n",replicas_number-1-i_counter,replicas_number-i_counter-2);
        
        accettata=metro_SWAP( conf_hasenbusch,loc_plaq,tr_local_plaqs, replicas_number-i_counter-1, replicas_number-i_counter-2, def);
#pragma acc update device(conf_hasenbusch[0:replicas_number][0:8])
            
            *swap_num=*swap_num+1;
        
            //all_swap_vet[replicas_number-i_counter-1]++;
            all_swap_vet[replicas_number-i_counter-2]++;
            
            if(accettata==1){
               // acceptance_vet[replicas_number-i_counter-1]++;
               acceptance_vet[replicas_number-i_counter-2]++;
            }
            
      
     }
        
        
    }
    
    return;
}

/*
void trasl_conf( __restrict const su3_soa *  const tconf_acc,
                 __restrict const su3_soa *  const taux_conf){
        printf("CONF 0 TRASL\n");
    
   set_su3_soa_to_su3_soa(tconf_acc,taux_conf);// conf_aux=conf_acc
    

    
    int dir;
    
    dir=rand()%4;
    

    
    printf("Mu is  %d\n",dir);
  
    
    
    int i,j,z,t; int idxh=0;
    int parity=0;
    int mu=0;
    int idxpnu=0;
    
    printf("conf e conf aux :%f || %f\n", creal( tconf_acc[2*mu].r0.c0[snum_acc(1,1,1,1)]),creal(taux_conf[2*mu].r0.c0[snum_acc(1,1,1,1)]));
  
    for(mu=0;mu<4;mu++){

        for(t=0;t<nd3;t++) {
            for (z=0; z<nd2;z++){
                for(j=0;j<nd1;j++){
                    for(i=0;i<nd0;i++){
        
                        idxh=snum_acc(i,j,z,t);
                        parity= (i+j+z+t) % 2;
            
            idxpnu = nnm_openacc[idxh][dir][parity];// r-nu //the table that states which is the nearest neighbour.
            
              
                        
            if(parity==0){
                tconf_acc[2*mu].r0.c0[idxh]=taux_conf[2*mu].r0.c0[idxpnu];
                tconf_acc[2*mu].r0.c1[idxh]=taux_conf[2*mu].r0.c1[idxpnu];
                tconf_acc[2*mu].r0.c2[idxh]=taux_conf[2*mu].r0.c2[idxpnu];
                
                tconf_acc[2*mu].r1.c0[idxh]=taux_conf[2*mu].r1.c0[idxpnu];
                tconf_acc[2*mu].r1.c1[idxh]=taux_conf[2*mu].r1.c1[idxpnu];
                tconf_acc[2*mu].r1.c2[idxh]=taux_conf[2*mu].r1.c2[idxpnu];
                
                
                tconf_acc[2*mu].r2.c0[idxh]=taux_conf[2*mu].r2.c0[idxpnu];
                tconf_acc[2*mu].r2.c1[idxh]=taux_conf[2*mu].r2.c1[idxpnu];
                tconf_acc[2*mu].r2.c2[idxh]=taux_conf[2*mu].r2.c2[idxpnu];
                
                        }
             if(parity==1){
                 
                 tconf_acc[2*mu+1].r0.c0[idxh]=taux_conf[2*mu+1].r0.c0[idxpnu];
                 tconf_acc[2*mu+1].r0.c1[idxh]=taux_conf[2*mu+1].r0.c1[idxpnu];
                 tconf_acc[2*mu+1].r0.c2[idxh]=taux_conf[2*mu+1].r0.c2[idxpnu];
                 
                 tconf_acc[2*mu+1].r1.c0[idxh]=taux_conf[2*mu+1].r1.c0[idxpnu];
                 tconf_acc[2*mu+1].r1.c1[idxh]=taux_conf[2*mu+1].r1.c1[idxpnu];
                 tconf_acc[2*mu+1].r1.c2[idxh]=taux_conf[2*mu+1].r1.c2[idxpnu];
                 
                 
                 tconf_acc[2*mu+1].r2.c0[idxh]=taux_conf[2*mu+1].r2.c0[idxpnu];
                 tconf_acc[2*mu+1].r2.c1[idxh]=taux_conf[2*mu+1].r2.c1[idxpnu];
                 tconf_acc[2*mu+1].r2.c2[idxh]=taux_conf[2*mu+1].r2.c2[idxpnu];
                 
                 
                 
             }
                        
                        
                        
                    }//close for x
                }//y
            }//z
        }//t
    }//mu
    
     printf("conf e conf aux :%f || %f\n", creal( tconf_acc[2*mu].r0.c0[snum_acc(1,1,1,1)]),creal(taux_conf[2*mu].r0.c0[snum_acc(1,1,1,1)]));
    
   // #pragma acc update device(tconf_acc[0:8])

    
    return;
}
*/
    
void trasl_conf( __restrict const su3_soa *  const tconf_acc,
                __restrict const su3_soa *  const taux_conf){
    printf("CONF 0 TRASL\n");
    
#ifdef MULTIDEVICE
    communicate_su3_borders(tconf_acc, GAUGE_HALO);
    #pragma acc update self(tconf_acc[0:8])
#endif
    
    set_su3_soa_to_su3_soa(tconf_acc,taux_conf);// conf_aux=conf_acc
    
    
    
 double dir0;
   int dir=0;
    
    
    
    if(0==devinfo.myrank){dir0=casuale();}
    
   #ifdef MULTIDEVICE
                MPI_Bcast((void*) &dir0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                printf("MPI%02d dir0 : %f \n",devinfo.myrank,dir0);
#endif


    

    if(dir0<=0.25){dir=0;}
    if(dir0>0.25 && dir0<=0.5){dir=1;}
    if(dir0>0.5 && dir0<=0.75){dir=2;}
    if(dir0>0.75){dir=3;}
    
    if(0==devinfo.myrank){
     printf("dir0 %f",dir0);
    printf("Mu is  %d\n",dir);
    
    
    
    printf("conf e conf aux :%f || %f\n", creal( tconf_acc[0].r0.c0[snum_acc(1,1,1,1)]),creal(taux_conf[0].r0.c0[snum_acc(1,1,1,1)]));
    }
    

         printf("dir0 %f,rank %d",dir0,devinfo.myrank);
    printf("Mu is  %d\n",dir);

printf("%d GAUGE_HALO %d D3_HALO \n",GAUGE_HALO,D3_HALO);


    set_su3_soa_to_su3_soa_trasl( taux_conf,tconf_acc, dir);
#pragma acc update device(tconf_acc[0:8])  
    
#ifdef MULTIDEVICE
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
    #pragma acc update self(tconf_acc[0:8])
#endif


    if(0==devinfo.myrank){
    printf("conf e conf aux :%f || %f\n", creal( tconf_acc[0].r0.c0[snum_acc(1,1,1,1)]),creal(taux_conf[0].r0.c0[snum_acc(1,1,1,1)]));
    
        printf("conf %f e nnp %f\n",creal( tconf_acc[0].r0.c0[snum_acc(1,1,1,1)]),creal(tconf_acc[0].r0.c0[nnp_openacc[snum_acc(1,1,1,1)][dir][1]]));
    }
    // #pragma acc update device(tconf_acc[0:8])
    
    
    return;
}


