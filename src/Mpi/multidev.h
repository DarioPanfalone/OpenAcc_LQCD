#ifndef MULTIDEV_H_
#define MULTIDEV_H_
#ifdef MULTIDEVICE

#include <mpi.h>
#include "./geometry_multidev.h"

typedef struct multidev_info_t{

    // FROM MPI INIT AND SIMILAR
    int myrank, nranks;
    int myrank_L,myrank_R; // nearest neighbours, salamino
    int node_subrank; // which card it will use?
    char processor_name[MPI_MAX_PROCESSOR_NAME]; int namelen;

    // GEOMETRIC QUANTITIES
    vec4int myrank4int;
    int nnranks[4][2];// ranks of nearest neighbour sublattices
    vec4int gl_loc_origin4int; = gl_loc_origin_from_rank(myrank);

    // FROM INPUT FILE
    int proc_per_node;
    int nranks_read;

} multidev_info ; 


extern multidev_info mdevinfo;

void init_multidev1D(multidev_info * mdi);

void shutdown_multidev();


#endif
#endif
