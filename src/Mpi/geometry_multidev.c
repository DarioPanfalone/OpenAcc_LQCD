#ifndef GEOMETRY_MULTIDEV_H_
#define GEOMETRY_MULTIDEV_H_


#ifdef MULTIDEVICE
#include "./geometry_multidev.h"

void setup_nnranks(int **nnranks, int myrank){

    int nranks_vect[4] = {NRANKS_D0, NRANKS_D1, NRANKS_D2, NRANKS_D3};

    myrank4int = xyzt_rank(myrank);
    int myrank_vect[4];
    myrank_vect[0] = myrank4int.d0 ;
    myrank_vect[1] = myrank4int.d1 ; 
    myrank_vect[2] = myrank4int.d2 ; 
    myrank_vect[3] = myrank4int.d3 ; 
    int dir;
    for(dir=0;dir<4;dir++){
        int lrank_dir = (myrank_vect[dir]-1);
        if(lrank_dir == -1) lrank_dir = nranks_vect[dir]-1;
        int rrank_dir = (myrank_vect[dir]+1)%nranks_vect[dir];
        nnranks[dir][0] = lrank_dir;
        nnranks[dir][1] = rrank_dir;
    }
}

#endif

#endif
