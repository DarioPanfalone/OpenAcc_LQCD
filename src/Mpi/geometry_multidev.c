#ifndef GEOMETRY_MULTIDEV_H_
#define GEOMETRY_MULTIDEV_H_

#ifdef MULTIDEVICE
#include "mpi.h"

void set_geom_glv_multidev(geom_par_multidev* gp){

    MPI_Comm_rank( MPI_COMM_WORLD, &(gp->rank));
    MPI_Comm_size( MPI_COMM_WORLD, &(gp->n_ranks));
    myrank4int;
    gl_loc_origin4int;

    nnranks[4][2];// ranks of nearest neighbour sublattices

}

#endif

#endif
