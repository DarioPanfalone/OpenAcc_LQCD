#ifndef MULTIDEV_C_
#define MULTIDEV_C_

#include <mpi.h>
#include <stdio.h>
#include "./geometry_multidev.h"

#ifdef MULTIDEVICE
#include "./multidev.h"

multidev_info mdevinfo;

void init_multidev1D(multidev_info * mdi){

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &(mdi->myrank));
    MPI_Comm_size(MPI_COMM_WORLD, &(mdi->nranks));
    MPI_Get_processor_name(mdi->processor_name,&(mdi->namelen));

    if(mdi->nranks != NRANKS_D3){
        printf("#MPI%02d: NRANKS_D3 is different from nranks: no salamino? Exiting now\n",mdi->myrank);
        MPI_Finalize();
        return;
    }
    
    mdi->myrank_L = (mdi->myrank + (mdi->nranks-1))%mdi->nranks;//SALAMINO
    mdi->myrank_R = (mdi->myrank + 1 ) % mdi->nranks;       //SALAMINO
    mdi->node_subrank =  mdi->myrank % mdi->proc_per_node;   //SALAMINO

    printf("#MPI%02d: of \"%02d\" tasks running on host \"%s\", local rank: %d, rankL: %d, rankR: %d\n",
    mdi->myrank, mdi->nranks, mdi->processor_name, 
    mdi->node_subrank, mdi->myrank_L, mdi->myrank_R); //SALAMINO

    mdi->myrank4int = xyzt_rank(mdi->myrank); // ALL RANKS
    setup_nnranks(mdi->nnranks,mdi->myrank);
    mdi->gl_loc_origin4int = gl_loc_origin_from_rank(mdi->myrank);

}



void shutdown_multidev(){

    MPI_Finalize();     

}





#endif
#endif
