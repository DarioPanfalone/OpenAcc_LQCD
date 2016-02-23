#ifndef GEOMETRY_MULTIDEV_H_
#define GEOMETRY_MULTIDEV_H_

    // GEOMETRIC DEFINES
#define HALO_WIDTH 2  // for rectangles
//LOCAL lattice dimensions

#include "../../build/lattice_dimensions.h"

///HALO STRUCTURE (border widths)
// AUTOMATIC DEFINITION OF HALOS
#if NRANKS_D0 == 1 
#define D0_HALO 0
#else 
#define D0_HALO HALO_WIDTH
#define MULTIDEVICE
#endif

#if NRANKS_D1 == 1 
#define D1_HALO 0
#else 
#define D1_HALO HALO_WIDTH
#define MULTIDEVICE
#endif

#if NRANKS_D2 == 1 
#define D2_HALO 0
#else 
#define D2_HALO HALO_WIDTH
#define MULTIDEVICE
#endif

#if NRANKS_D3 == 1 
#define D3_HALO 0
#else 
#define D3_HALO HALO_WIDTH
#define MULTIDEVICE
#endif

// Local aNd Halo lattice dimensions 
// Hopefully automatically calculated by the compiler during optimization
//AUTOMATIC DEFINITIONS OF LNH STUFF
#define LNH_N0 (LOC_N0+2*D0_HALO)
#define LNH_N1 (LOC_N1+2*D1_HALO)
#define LNH_N2 (LOC_N2+2*D2_HALO)
#define LNH_N3 (LOC_N3+2*D3_HALO)

#define nd0 LNH_N0
#define nd1 LNH_N1
#define nd2 LNH_N2
#define nd3 LNH_N3

#ifdef MULTIDEVICE

#define NRANKS NRANKS_D0*NRANKS_D1*NRANKS_D2*NRANKS_D3

// LOCAL STUFF
#define LOC_VOL1  LOC_N0
#define LOC_VOL2 (LOC_N1 * LOC_VOL1)
#define LOC_VOL3 (LOC_N2 * LOC_VOL2)
#define LOC_VOL4 (LOC_N3 * LOC_VOL3)
#define LOC_N0H ( LOC_N0 >> 1) // nx/2
#define LOC_N1H ( LOC_N1 >> 1)
#define LOC_N2H ( LOC_N2 >> 1)
#define LOC_N3H ( LOC_N3 >> 1)

#define LOC_SIZE LOC_VOL4
#define LOC_SIZE2 (2*LOC_SIZE)
#define LOC_SIZE3 (3*LOC_SIZE)
#define LOC_SIZEH ((int)(LOC_SIZE / 2) + LOC_SIZE % 2)
#define LOC_SIZEH_A ((int)(LOC_SIZE / 2) + LOC_SIZE % 2) //no of sites having the same parity as the local origin
#define LOC_SIZEH_B ((int)(LOC_SIZE / 2)) // no of sites having opposite parity to the local origin
#define LOC_NO_LINKS (4 * LOC_VOL4)


#define LNH_VOL1  LNH_N0
#define LNH_VOL2 (LNH_N1 * LNH_VOL1)
#define LNH_VOL3 (LNH_N2 * LNH_VOL2)
#define LNH_VOL4 (LNH_N3 * LNH_VOL3)
#define LNH_N0H ( LNH_N0 >> 1) // nx/2
#define LNH_N1H ( LNH_N1 >> 1)
#define LNH_N2H ( LNH_N2 >> 1)
#define LNH_N3H ( LNH_N3 >> 1)

#define LNH_SIZE LNH_VOL4
#define LNH_SIZE2 (2*LNH_SIZE)
#define LNH_SIZE3 (3*LNH_SIZE)
#define LNH_SIZEH ((int)(LNH_SIZE / 2)+ LNH_SIZE % 2 )
#define LNH_SIZEH_A ((int)(LNH_SIZE / 2)+ LNH_SIZE % 2 ) //no of sites having the same parity as the local origin
#define LNH_SIZEH_B ((int)(LNH_SIZE / 2)) // no of sites having opposite parity to the local origin
#define LNH_NO_LINKS (4 * LNH_VOL4)

//SURFACE STRUCTURE (surface widths)
//AUTOMATIC DEFINITION OF BULK STUFF
#define D0_SURF D0_HALO 
#define D1_SURF D1_HALO
#define D2_SURF D2_HALO
#define D3_SURF D3_HALO

//Local Bulk lattice dimensions

#define LB_N0 (LOC_N0-2*D0_SURF)
#define LB_N1 (LOC_N1-2*D1_SURF)
#define LB_N2 (LOC_N2-2*D2_SURF)
#define LB_N3 (LOC_N3-2*D3_SURF)

#define LB_VOL1  LB_N0
#define LB_VOL2 (LB_N1 * LB_VOL1)
#define LB_VOL3 (LB_N2 * LB_VOL2)
#define LB_VOL4 (LB_N3 * LB_VOL3)
#define LB_N0H ( LB_N0 >> 1) // nx/2
#define LB_N1H ( LB_N1 >> 1)
#define LB_N2H ( LB_N2 >> 1)
#define LB_N3H ( LB_N3 >> 1)

#define LB_SIZE LB_VOL4
#define LB_SIZE2 (2*LB_SIZE)
#define LB_SIZE3 (3*LB_SIZE)
#define LB_SIZEH (LB_SIZE / 2)
#define LB_NO_LINKS (4 * LB_VOL4)

// GLOBAL quantities
//AUTOMATIC DEFINITIONS OF GLOBAL STUFF
#define GL_N0 (LOC_N0*NRANKS_D0)
#define GL_N1 (LOC_N1*NRANKS_D1)
#define GL_N2 (LOC_N2*NRANKS_D2)
#define GL_N3 (LOC_N3*NRANKS_D3)

#define GL_VOL1  GL_N0
#define GL_VOL2 (GL_N1 * GL_VOL1)
#define GL_VOL3 (GL_N2 * GL_VOL2)
#define GL_VOL4 (GL_N3 * GL_VOL3)
#define GL_N0H ( GL_N0 >> 1) // nx/2
#define GL_N1H ( GL_N1 >> 1)
#define GL_N2H ( GL_N2 >> 1)
#define GL_N3H ( GL_N3 >> 1)

#define GL_SIZE GL_VOL4
#define GL_SIZE2 (2*GL_SIZE)
#define GL_SIZE3 (3*GL_SIZE)
#define GL_SIZEH (GL_SIZE / 2)
#define GL_NO_LINKS (4 * GL_VOL4)

typedef struct vec4int_t{
int x,y,z,t;

} vec4int;


typedef struct geom_par_multidev_t{
    int n_ranks;//  = -1; //this should be correctly set by a call to MPI_init()
    int myrank;// = -1 ; //this should be correctly set by a call to MPI_init()
// following stuff set in init_LNH_geo() using myrank and NRANKS_X etc
    vec4int myrank4int;
    vec4int gl_loc_origin4int;
    int gl_orig_d0123[4];
    // "global x coordinate of the local origin"

} geom_par_multidev;

extern geom_par_multidev;
int nnranks[4][2];// ranks of nearest neighbour sublattices

void set_geom_glv_multidev(geom_par_multidev* gp);


static inline int gl_to_gl_snum(int gl_0, int gl_1, int gl_2, int gl_3){
//global coordinates to global 'snum' index

    int ris = gl_x + gl_y * GL_VOL1 + gl_z * GL_VOL2 + gl_t * GL_VOL3 ;
    return ris/2;// <---  /2 Pay attention to even/odd  (see init_geo) 

}
static inline int lnh_to_lnh_snum(int lnh_0, int lnh_1, int lnh_2, int lnh_3){
// local'n'halo coordinates to lnh 'snum' index

    //NOTE: Since lnh_0 goes from 0 to LNH_N0-1,
    //      here loc_0 goes from -D0_HALO to LOC_N0+D0_HALO-1

    int ris = lnh_0 + lnh_1 * LNH_VOL1 + lnh_2 * LNH_VOL2 + lnh_3 * LNH_VOL3;

    return ris/2;// <---  /2 Pay attention to even/odd  (see init_geo) 

}
static inline int loc_to_lnh_snum(int loc_0, int loc_1, int loc_2, int loc_3){
//local coordinates to loc'n'halo 'snum'index
    // Actually, the real memory layout is of the LNH type.
    // So, at the moment I don't see any real use for this function.
    
    loc_0+=D0_HALO;
    loc_1+=D1_HALO;
    loc_2+=D2_HALO;
    loc_3+=D3_HALO;

    return lnh_to_lnh_snum(loc_0,loc_1,loc_2,loc_3);
    // ^^ /2 Pay attention to even/odd  (see init_geo) 
}
static inline int lnh_to_gl_snum(int lnh_0, int lnh_1, int lnh_2, int lnh_3){

    lnh_0 += erte - D0_HALO;
    lnh_1 +=  - D1_HALO;
    lnh_2 +=  - D2_HALO;
    lnh_3 +=  - D3_HALO;

    // in case lnh_* is out of bonds
    lnh_0 = lnh_0 % GL_N0; lnh_0-= GL_N0 * ( lnh_0 >> 31 );
    lnh_1 = lnh_1 % GL_N1; lnh_1-= GL_N1 * ( lnh_1 >> 31 );
    lnh_2 = lnh_2 % GL_N2; lnh_2-= GL_N2 * ( lnh_2 >> 31 );
    lnh_3 = lnh_3 % GL_N3; lnh_3-= GL_N3 * ( lnh_3 >> 31 );

    return gl_to_gl_snum(lnh_0,lnh_1,lnh_2,lnh_3);
    // pay attention to even/odd

}
static inline int target_lnh_to_gl_snum(int lnh_x, int lnh_y, int lnh_z, int lnh_t,                          vec4int target_gl_loc_origin4int ){

    lnh_x += target_gl_loc_origin4int.x - X_HALO;
    lnh_y += target_gl_loc_origin4int.y - Y_HALO;
    lnh_z += target_gl_loc_origin4int.z - Z_HALO;
    lnh_t += target_gl_loc_origin4int.t - T_HALO;

    lnh_x = lnh_x % GL_NX - GL_NX * ( lnh_x % GL_NX >> 31 );
    lnh_y = lnh_y % GL_NY - GL_NY * ( lnh_y % GL_NY >> 31 );
    lnh_z = lnh_z % GL_NZ - GL_NZ * ( lnh_z % GL_NZ >> 31 );
    lnh_t = lnh_t % GL_NT - GL_NT * ( lnh_t % GL_NT >> 31 );

    return gl_to_gl_snum(lnh_x,lnh_y,lnh_z,lnh_t);
    // pay attention to even/odd

}
/* // ETA FUNCTIONS WHICH SHOULD NOT BE NECESSARY

static inline int loc_eta(int loc_x, int loc_y, int loc_z, int loc_t, int dir){

    loc_x +=  gl_loc_origin4int.x;
    loc_y +=  gl_loc_origin4int.y;
    loc_z +=  gl_loc_origin4int.z;
    loc_t +=  gl_loc_origin4int.t;
    return gl_eta(loc_x,loc_y,loc_z,loc_t,dir);

}
static inline int lnh_eta(int lnh_x, int lnh_y, int lnh_z, int lnh_t, int dir){
    lnh_x -= X_HALO; 
    lnh_y -= Y_HALO;
    lnh_z -= Z_HALO;
    lnh_t -= T_HALO;
    return loc_eta(lnh_x,lnh_y,lnh_z,lnh_t,dir);
}
*/
/***************************************************
links used according to this scheme

0           LNH_SIZE       2*LNH_SIZE      3*LNH_SIZE       LNH_NO_LINKS
|-------|-------|-------|-------|-------|-------|-------|-------|
    e       o        e      o       e       o       e       o
      x-dir           y-dir           z-dir           t-dir
NOTE:the number of even and odd sites on the local lattice can be different for the local sublattice if the local sublattice dimensions are all odd.
initialize geometry
periodic spatial bc are always assumed (this is relevant only if
the considered direction is not paralelized))
*/

vec4int xyzt_rank(int rank){

    vec4int rank4int;

    rank4int.x = rank % NRANKS_X;
    rank4int.y =( rank / NRANKS_X) % NRANKS_Y;
    rank4int.z =( rank / (NRANKS_X*NRANKS_Y)) % NRANKS_Z;
    rank4int.t =( rank / (NRANKS_X*NRANKS_Y*NRANKS_Z)) % NRANKS_T;

    return rank4int;

}
int rank_from_xyzt_rank(int rank_x, int rank_y, int rank_z, int rank_t){

    return rank_x + ( NRANKS_X * (rank_y + NRANKS_Y * (rank_z  + NRANKS_Z * rank_t)));

}
static inline vec4int gl_loc_origin_from_rank(int rank){
    vec4int res, rank4int = xyzt_rank(rank);
    res.x = rank4int.x * LOC_NX;
    res.y = rank4int.y * LOC_NY;
    res.z = rank4int.z * LOC_NZ;
    res.t = rank4int.t * LOC_NT;
    return res;

}
void init_LNH_geo(void){

    // NOTE: loc_x = lnh_x - X_HALO; 
    //       loc_xm = lnh_xm - X_HALO; 
    //       loc_xp = lnh_xp - X_HALO; 

    int lnh_num;
    int sum, rest;


    // local myrank stuff
    // NOTICE : These global variables are declared in ../Include/common_defines.h


    myrank4int = xyzt_rank(myrank);

    gl_loc_origin4int = gl_loc_origin_from_rank(myrank);

    int gl_coord_loc_origin_sum = gl_loc_origin4int.x + gl_loc_origin4int.y + 
        gl_loc_origin4int.z + gl_loc_origin4int.t;


    loc_origin_parity = gl_coord_loc_origin_sum % 2 ; 
    even_LNH_size = (int)( LNH_NX*LNH_NY*LNH_NZ*LNH_NT /2);
    if (loc_origin_parity == 0 ) even_LNH_size++;

    // The number of sites in the LNH lattice having the same parity
    // as the local origin is (int)( LNH_NX*LNH_NY*LNH_NZ*LNH_NT /2).
    // If the total number of sites in the LNH lattice is odd, the number of 
    // sites having opposite parity to the local origin is 
    // (int)( LNH_NX*LNH_NY*LNH_NZ*LNH_NT /2) + 1.
    // So, if the local origin is odd ( then loc_origin_parity==1) 
    // the number of even sites in the LNH lattice is 
    // (int)( LNH_NX*LNH_NY*LNH_NZ*LNH_NT /2) + 1.

}

void setup_nnranks(){

    int nranks_vect[4] = {NRANKS_X, NRANKS_Y, NRANKS_Z, NRANKS_T};

    myrank4int = xyzt_rank(myrank);
    int myrank_vect[4];
    myrank_vect[0] = myrank4int.x ;
    myrank_vect[1] = myrank4int.y ; 
    myrank_vect[2] = myrank4int.z ; 
    myrank_vect[3] = myrank4int.t ; 
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
