#ifndef MULTIDEV_C_
#define MULTIDEV_C_

#include <stdio.h>
#include <stdlib.h>
#include "./geometry_multidev.h"
#include "../OpenAcc/geometry.h"

#include "./multidev.h"

dev_info devinfo;

#ifdef MULTIDEVICE
#include <mpi.h>

extern int verbosity_lv;


void pre_init_multidev1D(dev_info * mdi){
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &(mdi->myrank));
    MPI_Comm_size(MPI_COMM_WORLD, &(mdi->nranks));
    MPI_Get_processor_name(mdi->processor_name,&(mdi->namelen));

    if(mdi->nranks != NRANKS_D3){
        printf("#MPI%02d: NRANKS_D3 is different from nranks: no salamino? Exiting now\n",mdi->myrank);
        MPI_Finalize();
        exit(1);
        return;
    }
    
    if(verbosity_lv > 2){
        printf("MPI%02d - Called MPI_Init\n",mdi->myrank);
    }



}

void init_multidev1D(dev_info * mdi){
    mdi->myrank_L = (mdi->myrank + (mdi->nranks-1))%mdi->nranks;//SALAMINO
    mdi->myrank_R = (mdi->myrank + 1 ) % mdi->nranks;       //SALAMINO
    mdi->node_subrank =  mdi->myrank % mdi->proc_per_node;   //SALAMINO

    sprintf(mdi->myrankstr,"MPI%02d",mdi->myrank);

    printf("#MPI%02d: of \"%02d\" tasks running on host \"%s\", local rank: %d, rankL: %d, rankR: %d\n",
    mdi->myrank, mdi->nranks, mdi->processor_name, 
    mdi->node_subrank, mdi->myrank_L, mdi->myrank_R); //SALAMINO

    mdi->myrank4int = xyzt_rank(mdi->myrank); // ALL RANKS


    mdi->myrank4int = xyzt_rank(mdi->myrank);
    int myrank_vect[4];
    myrank_vect[0] = mdi->myrank4int.d0 ;
    myrank_vect[1] = mdi->myrank4int.d1 ; 
    myrank_vect[2] = mdi->myrank4int.d2 ; 
    myrank_vect[3] = mdi->myrank4int.d3 ; 

    int dir;
    for(dir=0;dir<4;dir++){
        int lrank_dir = (myrank_vect[dir]-1);
        if(lrank_dir == -1) lrank_dir = geom_par.nranks[dir]-1;
        int rrank_dir = (myrank_vect[dir]+1)% geom_par.nranks[dir];
        mdi->nnranks[dir][0] = lrank_dir;
        mdi->nnranks[dir][1] = rrank_dir;
    }

    mdi->gl_loc_origin4int = gl_loc_origin_from_rank(mdi->myrank);

    mdi->halo_widths0123[0] = D0_HALO ;  
    mdi->halo_widths0123[1] = D1_HALO ;  
    mdi->halo_widths0123[2] = D2_HALO ;  
    mdi->halo_widths0123[3] = D3_HALO ;  
    mdi->origin_0123[0]     = mdi->gl_loc_origin4int.d0;
    mdi->origin_0123[1]     = mdi->gl_loc_origin4int.d1;
    mdi->origin_0123[2]     = mdi->gl_loc_origin4int.d2;
    mdi->origin_0123[3]     = mdi->gl_loc_origin4int.d3;
    if(verbosity_lv > 2){
        printf("MPI%02d - Finished init_multidev1D\n", mdi->myrank);   
   
    }



}



void shutdown_multidev(){

    MPI_Finalize();     

}





#endif
#endif
