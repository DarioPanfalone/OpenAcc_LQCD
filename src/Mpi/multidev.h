#ifndef MULTIDEV_H_
#define MULTIDEV_H_

#include "./geometry_multidev.h"

#ifdef MULTIDEVICE
#include <mpi.h>
#endif

typedef struct dev_info_t{

    // FROM MPI INIT AND SIMILAR
    
    int single_dev_choice; // from input file
    int myrank, nranks;
    char myrankstr[50];
    int nranks_read;

#ifdef MULTIDEVICE
    int async_comm_fermion;
    int async_comm_gauge;

    int myrank_L,myrank_R; // nearest neighbours, salamino
    int node_subrank; // which card it will use?
    char processor_name[MPI_MAX_PROCESSOR_NAME]; int namelen;

    // GEOMETRIC QUANTITIES
    vec4int myrank4int;
    int nnranks[4][2];// ranks of nearest neighbour sublattices
    vec4int gl_loc_origin4int;
    // FROM INPUT FILE
    int proc_per_node;
    int halo_widths0123[4];
    int origin_0123[4];
#endif


} dev_info ; 


extern dev_info devinfo;


#ifdef MULTIDEVICE

void pre_init_multidev1D(dev_info * mdi);
void init_multidev1D(dev_info * mdi);

void shutdown_multidev();

#endif

#endif
