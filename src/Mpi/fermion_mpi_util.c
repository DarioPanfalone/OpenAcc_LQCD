#ifndef FERMION_MPI_UTIL_C
#define FERMION_MPI_UTIL_C


#include <string.h>
#include <stdlib.h>

#ifdef __PGI
#include "openacc.h"
#endif

#ifdef MPI
#include <mpi.h>
#endif

#include "./struct_c_def.c"
#include "../Geometry/geometry.cc"

#include "../DbgTools/debug_tools.c"

/*********************
 * INITIALIZATION    *
 * ******************/


    void load_global_vec3_soa_from_file(global_vec3_soa * fermion, const char * filename){

        FILE *fp;
        double ar, ai, br, bi, cr, ci;
        int i = 0;
        int error = 0;

        fp = fopen(filename, "rt");

        if (fp == NULL) {
            printf("Could not open file %s \n", filename);
            exit(-1);
        }

        while ((i < GL_SIZEH) && (!error)) {

            if (fscanf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf) \n", &ar, &ai, &br, &bi, &cr, &ci) == 6) {
                fermion->c0[i] = (ar + ai * I);
                fermion->c1[i] = (br + bi * I);
                fermion->c2[i] = (cr + ci * I);
            } else {
                printf("Read error... ");
                error = 1;
            }

            i++;

        }

        printf("Read %d vectors from file %s \n", i, filename);

        fclose(fp);

    }
void write_global_vec3_soa_to_file(const global_vec3_soa * fermion, const char * filename){

 FILE *fp;
  double ar, ai, br, bi, cr, ci;
  int idx;
  int i = 0;
  int error = 0;

  fp = fopen(filename, "wt");

  if (fp == NULL) {
    printf("Could not open file %s \n", filename);
    exit(-1);
  }


    while ( (i < GL_SIZEH) && (!error) ) {
    
      ar  = creal(fermion->c0[i]);
      br  = creal(fermion->c1[i]);
      cr  = creal(fermion->c2[i]);
      ai  = cimag(fermion->c0[i]);
      bi  = cimag(fermion->c1[i]);
      ci  = cimag(fermion->c2[i]);

      fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf)\n", ar, ai, br, bi, cr, ci);

     i++;

    }

  printf("Written %d vectors to file %s \n", i, filename);

  fclose(fp);

}
void read_lnh_vec3_soa_from_file(lnh_vec3_soa * u, const char *lnh_filename){

    // TO WRITE

}
void write_lnh_vec3_soa_to_file(const lnh_vec3_soa * lnh_fermion, const char * global_filename){

    FILE *fp;
    double ar, ai, br, bi, cr, ci;
    int i = 0;
    int error = 0;

    // defining a filename for the sublattice file 
    char filename[50] = "rank";
    char rank_str[5];
    snprintf(rank_str, 5*sizeof(char),"%d",myrank);
    char n_ranks_str[5];
    snprintf(n_ranks_str, 5*sizeof(char),"%d",n_ranks);

    strcat(filename, rank_str);
    strcat(filename, "of");
    strcat(filename, n_ranks_str);    
    strcat(filename, "_");
    strcat(filename,global_filename);//global filename appended

    fp = fopen(filename, "wt");
    
    if (fp == NULL) {
        printf("Could not open file %s \n", filename);
        exit(-1);
    }

    vec4int sl_glo = gl_loc_origin_from_rank(myrank);
  //Sub Lattice Global (coordinate of) Local Origin
    printf("#MPI%02d: Writing file %s, \n",myrank, filename);
    printf("%d %d %d %d\n", sl_glo.x, sl_glo.y, sl_glo.z, sl_glo.t); 
    
    while ( (i < LNH_SIZEH) && (!error) ) {

        ar  = creal(lnh_fermion->c0[i]);
        br  = creal(lnh_fermion->c1[i]);
        cr  = creal(lnh_fermion->c2[i]);
        ai  = cimag(lnh_fermion->c0[i]);
        bi  = cimag(lnh_fermion->c1[i]);
        ci  = cimag(lnh_fermion->c2[i]);

       fprintf(fp, "(%lf,%lf) (%lf,%lf) (%lf,%lf)\n", ar, ai, br, bi, cr, ci);
       i++;

    }

    printf("#MPI%02d: Written %d vectors to file %s \n",myrank, i, filename);

    fclose(fp);

}
void global_vec3_soa_init(int initmode, global_vec3_soa * fermion){
    if( initmode == 0) // Read from file
        load_global_vec3_soa_from_file(fermion, "fermion");
    else{
       
        int i = 0;
        while( i < GL_SIZEH ) {
            
            vec3 v;

            if (initmode == 2 ) v = zero_vec3;
            else if (initmode == 1 ) v = vec3_gauss();
            else if (initmode == 3 ) v = generator_vec3_debug();

            fermion->c0[i] = v.c0;
            fermion->c1[i] = v.c1;
            fermion->c2[i] = v.c2;

            i++;
        }
    }
}
void lnh_vec3_soa_init(int initmode, lnh_vec3_soa * fermion){
       
        int i = 0;
        while( i < LNH_SIZEH ) {
            
            vec3 v;

            if (initmode == 0 ) v = zero_vec3;
            else if (initmode == 1 ) v = vec3_gauss();
            else if (initmode == 3 ) v = generator_vec3_debug();

            fermion->c0[i] = v.c0;
            fermion->c1[i] = v.c1;
            fermion->c2[i] = v.c2;

            i++;
        }
}

// "BUFFER FUNCTIONS"
void send_lnh_subfermion_to_buffer(const global_vec3_soa * gl_soa_fermion, lnh_vec3_soa* lnh_fermion, int target_rank){
    // USE ONLY FROM MASTER RANK 
    // Notice that this function cannot tell the difference between
    // even and odd fermions
    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);

    // building sublattice duplicate, target_conf
    lnh_vec3_soa* target_vec3_soa = lnh_fermion; 

    int tg_lnh_xh,tg_lnh_y,tg_lnh_z,tg_lnh_t; //target-lnh coordinates

    
    // Copying all relevant vec3s into the sublattice
    for(tg_lnh_t=0;tg_lnh_t<LNH_NT; tg_lnh_t++)
    for(tg_lnh_z=0;tg_lnh_z<LNH_NZ; tg_lnh_z++)
    for(tg_lnh_y=0;tg_lnh_y<LNH_NY; tg_lnh_y++)
    for(tg_lnh_xh=0;tg_lnh_xh<LNH_NXH; tg_lnh_xh++){

        int tg_lnh_x = 2* tg_lnh_xh; // CHECK


        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
        int target_lnh_snum = lnh_to_lnh_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

        vec3 aux = vec3_from_global_vec3_soa(gl_soa_fermion,target_gl_snum);
        vec3_into_vec3_soa(aux,target_vec3_soa,target_lnh_snum);

    }

    //sending the subfermion
    //MPI_Send(target_su3_soa, 3*2*LNH_SIZEH , MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD); 
    // ^^ CHECK

}
void recv_loc_subfermion_from_buffer(global_vec3_soa * gl_soa_fermion, lnh_vec3_soa* lnh_fermion, int target_rank){
    // USE ONLY FROM MASTER RANK 
    // Notice that this function cannot tell the difference between
    // even and odd fermions

    //target sublattice information
    vec4int origin_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);
    // building sublattice duplicate, target_conf
    lnh_vec3_soa* origin_vec3_soa = lnh_fermion; 
    //receive the subfermion
//    MPI_Recv(origin_su3_soa, 3*2*LNH_SIZEH , MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // ^^ CHECK

    int or_loc_xh,or_loc_y,or_loc_z,or_loc_t; //target-lnh coordinates

    // NOTICE: only LOC sites are important, thus cycling only on terget local sites

    // Copying all relevant vec3s into the sublattice
    for(or_loc_t=0;or_loc_t<LOC_NT; or_loc_t++)
    for(or_loc_z=0;or_loc_z<LOC_NZ; or_loc_z++)
    for(or_loc_y=0;or_loc_y<LOC_NY; or_loc_y++)
    for(or_loc_xh=0;or_loc_xh<LOC_NXH; or_loc_xh++){


        int or_loc_x = 2* or_loc_xh; // CHECK
        int or_lnh_x,or_lnh_y,or_lnh_z,or_lnh_t; //target-lnh coordinates
        
        or_lnh_x = or_loc_x + X_HALO;
        or_lnh_y = or_loc_y + Y_HALO;
        or_lnh_z = or_loc_z + Z_HALO;
        or_lnh_t = or_loc_t + T_HALO;

        int origin_gl_snum = target_lnh_to_gl_snum(or_lnh_x, or_lnh_y, or_lnh_z, or_lnh_t, origin_gl_loc_origin_from_rank);
        int origin_lnh_snum = lnh_to_lnh_snum(or_lnh_x, or_lnh_y, or_lnh_z, or_lnh_t);

        vec3 aux = vec3_from_vec3_soa(origin_vec3_soa,origin_lnh_snum);
        vec3_into_global_vec3_soa(aux,gl_soa_fermion,origin_gl_snum);
    }
}

#ifdef MPI // TRUE MPI-USING FUNCTIONS
void send_lnh_subfermion_to_rank(const global_vec3_soa * gl_soa_fermion, int target_rank, int tag){

    // USE ONLY FROM MASTER RANK 
    // Notice that this function cannot tell the difference between
    // even and odd fermions
    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);

    // building sublattice duplicate, target_conf
    lnh_vec3_soa* target_vec3_soa = malloc(sizeof(lnh_vec3_soa)); 

    int tg_lnh_xh,tg_lnh_y,tg_lnh_z,tg_lnh_t; //target-lnh coordinates

    // Copying all relevant vec3s into the sublattice
    for(tg_lnh_t=0;tg_lnh_t<LNH_NT; tg_lnh_t++)
    for(tg_lnh_z=0;tg_lnh_z<LNH_NZ; tg_lnh_z++)
    for(tg_lnh_y=0;tg_lnh_y<LNH_NY; tg_lnh_y++)
    for(tg_lnh_xh=0;tg_lnh_xh<LNH_NXH; tg_lnh_xh++){

        int tg_lnh_x = 2* tg_lnh_xh; // CHECK


        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
        int target_lnh_snum = lnh_to_lnh_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

        vec3 aux = vec3_from_global_vec3_soa(gl_soa_fermion,target_gl_snum);
        vec3_into_vec3_soa(aux,target_vec3_soa,target_lnh_snum);

    }

    //sending the subfermion
    MPI_Send(target_vec3_soa, 3*2*LNH_SIZEH , MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD);

}
inline void recv_lnh_subfermion_from_master(lnh_vec3_soa * lnh_fermion, int tag){

    MPI_Recv(lnh_fermion, 3*2*LNH_SIZEH,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE );
    // ^^ CHECK
}
inline void send_lnh_subfermion_to_master(lnh_vec3_soa * lnh_fermion, int tag){

    MPI_Send(lnh_fermion, 3*2*LNH_SIZEH,MPI_DOUBLE,0,tag,MPI_COMM_WORLD );
    // ^^ CHECK
}
void recv_loc_subfermion_from_rank(const global_vec3_soa * gl_soa_fermion, int rank, int tag){
    // USE ONLY FROM MASTER RANK 
    // Notice that this function cannot tell the difference between
    // even and odd fermions

    //target sublattice information
    vec4int origin_gl_loc_origin_from_rank = gl_loc_origin_from_rank(rank);
    // building sublattice duplicate, target_conf
    lnh_vec3_soa* origin_vec3_soa = malloc(sizeof(lnh_vec3_soa)); 
    //receive the subfermion
    MPI_Recv(origin_vec3_soa, 3*2*LNH_SIZEH , MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // ^^ CHECK

    int or_loc_xh,or_loc_y,or_loc_z,or_loc_t; //target-lnh coordinates

    // NOTICE: only LOC sites are important, thus cycling only on terget local sites

    // Copying all relevant vec3s into the sublattice
    for(or_loc_t=0;or_loc_t<LOC_NT; or_loc_t++)
    for(or_loc_z=0;or_loc_z<LOC_NZ; or_loc_z++)
    for(or_loc_y=0;or_loc_y<LOC_NY; or_loc_y++)
    for(or_loc_xh=0;or_loc_xh<LOC_NXH; or_loc_xh++){


        int or_loc_x = 2* or_loc_xh; // CHECK
        int or_lnh_x,or_lnh_y,or_lnh_z,or_lnh_t; //target-lnh coordinates
        
        or_lnh_x = or_loc_x + X_HALO;
        or_lnh_y = or_loc_y + Y_HALO;
        or_lnh_z = or_loc_z + Z_HALO;
        or_lnh_t = or_loc_t + T_HALO;

        int origin_gl_snum = target_lnh_to_gl_snum(or_lnh_x, or_lnh_y, or_lnh_z, or_lnh_t, origin_gl_loc_origin_from_rank);
        int origin_lnh_snum = lnh_to_lnh_snum(or_lnh_x, or_lnh_y, or_lnh_z, or_lnh_t);

        vec3 aux = vec3_from_vec3_soa(origin_vec3_soa,origin_lnh_snum);
        vec3_into_global_vec3_soa(aux,gl_soa_fermion,origin_gl_snum);
    }
}
 void sendrecv_fermion_borders_1Dcut(lnh_vec3_soa *lnh_fermion, int rankL, int rankR){
   // NOTICE YOU HAVE TO SET MYRANK CORRECTLY TO USE THIS FUNCTION
  
  if(NRANKS_X != 1 || NRANKS_Y != 1 || NRANKS_Z != 1)
     printf("THIS SETUP IS NOT SALAMINO-LIKE!!!\n communication of fermion borders will FAIL!!\n");

/*   // PREAMBLE
   // This function is written taking the following assumptions:
   // 1. The domain is divided only along one direction, which is 
   //  the 'slowest' (strong assumption);
   // 2. That direction is the T direction (a bit weaker assumption, 
   //  changing that requires a complete redefinition of the site 
   //  ordering, which can be, perhaps, easily done by redefining the
   //  various '*snum' functions in Geometry/geometry.cc, and by 
   //  writing a tool to 'transpose' the configurations which have
   //  been written in a 'standard' ordering.*/

   //must be done for the three components of the fermion.

  int slab_sizeh = (LNH_NXH * LNH_NY * LNH_NZ * T_HALO); // NOTICE THERE IS LNH_NXH
  // no. of 'fermion' point in each slab)
  MPI_Status status;
#ifdef USE_MPI_CUDA_AWARE
#pragma acc host_data use_device(lnh_fermion)
  {
#endif
      d_complex *c[3] ;
      c[0] = lnh_fermion->c0;
      c[1] = lnh_fermion->c1;
      c[2] = lnh_fermion->c2;

      int ii;
      for(int ii =0; ii<3; ii++){

          // ASK FOR BACKS FIRST, THEN FACES
          int sendtag = ii;
          int recvtag = ii;

          d_complex *tmpc = c[ii];
#ifndef USE_MPI_CUDA_AWARE
#pragma acc update host(tmpc[slab_sizeh:slab_sizeh])
#endif
          MPI_Sendrecv((void*) &(c[ii][slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankL,sendtag,
                  (void*) &(c[ii][LNH_SIZEH-slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankR,recvtag,
                  MPI_COMM_WORLD, &status);
#ifndef USE_MPI_CUDA_AWARE
#pragma acc update device(tmpc[(LNH_SIZEH-slab_sizeh):slab_sizeh])
#endif

          sendtag = ii+3;
          recvtag = ii+3;

#ifndef USE_MPI_CUDA_AWARE
#pragma acc update host(tmpc[(LNH_SIZEH-2*slab_sizeh):slab_sizeh])
#endif
          MPI_Sendrecv((void*) &(c[ii][LNH_SIZEH-2*slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankR,sendtag,
                  (void*) c[ii],2*slab_sizeh,MPI_DOUBLE,
                  rankL,recvtag,
                  MPI_COMM_WORLD, &status);
#ifndef USE_MPI_CUDA_AWARE
#pragma acc update device(tmpc[0:slab_sizeh])
#endif
      }
#ifdef USE_MPI_CUDA_AWARE
  }
#endif
 }
 void sendrecv_fermion_borders_1Dcut_hostonly(lnh_vec3_soa *lnh_fermion, int rankL, int rankR){
   // NOTICE YOU HAVE TO SET MYRANK CORRECTLY TO USE THIS FUNCTION
  
  if(NRANKS_X != 1 || NRANKS_Y != 1 || NRANKS_Z != 1)
     printf("THIS SETUP IS NOT SALAMINO-LIKE!!!\n communication of fermion borders will FAIL!!\n");

/*   // PREAMBLE
   // This function is written taking the following assumptions:
   // 1. The domain is divided only along one direction, which is 
   //  the 'slowest' (strong assumption);
   // 2. That direction is the T direction (a bit weaker assumption, 
   //  changing that requires a complete redefinition of the site 
   //  ordering, which can be, perhaps, easily done by redefining the
   //  various '*snum' functions in Geometry/geometry.cc, and by 
   //  writing a tool to 'transpose' the configurations which have
   //  been written in a 'standard' ordering.*/

   //must be done for the three components of the fermion.

  int slab_sizeh = (LNH_NXH * LNH_NY * LNH_NZ * T_HALO); // NOTICE THERE IS LNH_NXH
  // no. of 'fermion' point in each slab)
  MPI_Status status;
  {
      d_complex *c[3] ;
      c[0] = lnh_fermion->c0;
      c[1] = lnh_fermion->c1;
      c[2] = lnh_fermion->c2;

      int ii;
      for(int ii =0; ii<3; ii++){

          // ASK FOR BACKS FIRST, THEN FACES
          int sendtag = ii;
          int recvtag = ii;

          d_complex *tmpc = c[ii];
          MPI_Sendrecv((void*) &(c[ii][slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankL,sendtag,
                  (void*) &(c[ii][LNH_SIZEH-slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankR,recvtag,
                  MPI_COMM_WORLD, &status);

          sendtag = ii+3;
          recvtag = ii+3;

          MPI_Sendrecv((void*) &(c[ii][LNH_SIZEH-2*slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankR,sendtag,
                  (void*) c[ii],2*slab_sizeh,MPI_DOUBLE,
                  rankL,recvtag,
                  MPI_COMM_WORLD, &status);
      }
  }
}
inline void communicate_fermion_borders(lnh_vec3_soa *lnh_fermion){ //WRAPPER

    // NOTICE: GEOMETRY MUST BE SET UP BEFORE!!
    MPI_Barrier(MPI_COMM_WORLD);
    sendrecv_fermion_borders_1Dcut(lnh_fermion,nnranks[3][0],nnranks[3][1] );
    MPI_Barrier(MPI_COMM_WORLD);
}
inline void communicate_fermion_borders_hostonly(lnh_vec3_soa *lnh_fermion){ //WRAPPER

    // NOTICE: GEOMETRY MUST BE SET UP BEFORE!!
    MPI_Barrier(MPI_COMM_WORLD);
    sendrecv_fermion_borders_1Dcut_hostonly(lnh_fermion,nnranks[3][0],nnranks[3][1] );
    MPI_Barrier(MPI_COMM_WORLD);
}

#ifdef USE_MPI_CUDA_AWARE  // This can work only if cuda aware MPI is used.
                           // If instead you wanted to use the update
                           // directive, this cannot work with MPI_Irecv
                           // because the directive should wait for the 
                           // completion the recv, but Irecv is not bloking!!
                          
void sendrecv_fermion_borders_1Dcut_async(lnh_vec3_soa *lnh_fermion, int rankL, int rankR,MPI_Request* send_border_requests, MPI_Request* recv_border_requests){
   // NOTICE YOU HAVE TO SET MYRANK CORRECTLY TO USE THIS FUNCTION
   // NOTICE send_border_requests recv_border_requests are both 
   // 6-elements long
  
  if(NRANKS_X != 1 || NRANKS_Y != 1 || NRANKS_Z != 1)
     printf("THIS SETUP IS NOT SALAMINO-LIKE!!!\n communication of fermion borders will FAIL!!\n");

   //SEE PREAMBLE FOR sendrecv_fermion_borders_1Dcut()
   //must be done for the three components of the fermion.

  int slab_sizeh = (LNH_NXH * LNH_NY * LNH_NZ * T_HALO); // NOTICE THERE IS LNH_NXH
  // no. of 'fermion' point in each slab)

#pragma acc host_data use_device(lnh_fermion)
  {
      d_complex *c[3] ;
      c[0] = lnh_fermion->c0;
      c[1] = lnh_fermion->c1;
      c[2] = lnh_fermion->c2;

      int ii;
      for(int ii =0; ii<3; ii++){

          // ASK FOR BACKS FIRST, THEN FACES
          int sendtag = ii;
          int recvtag = ii;

          d_complex *tmpc = c[ii];
          MPI_Isend((void*) &(c[ii][slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankL,sendtag,MPI_COMM_WORLD,&(send_border_requests[ii]));
          MPI_Irecv((void*) &(c[ii][LNH_SIZEH-slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankR,recvtag,MPI_COMM_WORLD,&(recv_border_requests[ii]));


          sendtag = 3+ii;
          recvtag = 3+ii;

          MPI_Isend((void*) &(c[ii][LNH_SIZEH-2*slab_sizeh]),2*slab_sizeh,MPI_DOUBLE,
                  rankR,sendtag,MPI_COMM_WORLD,&(send_border_requests[3+ii]));
          MPI_Irecv( (void*) c[ii],2*slab_sizeh,MPI_DOUBLE,
                  rankL,recvtag,MPI_COMM_WORLD,&(recv_border_requests[3+ii]));

      }
  }
}
#endif

inline void communicate_fermion_borders_async(lnh_vec3_soa *lnh_fermion, MPI_Request* send_border_requests, MPI_Request* recv_border_requests){ //WRAPPER
   
    // NOTICE: GEOMETRY MUST BE SET UP BEFORE!!
    // NOTICE send_border_requests recv_border_requests are both 
    // 6-elements long
    sendrecv_fermion_borders_1Dcut_async(lnh_fermion,nnranks[3][0],nnranks[3][1], send_border_requests, recv_border_requests );
}
#endif

/******************
 * COMMUNICATIONS *
 *****************/

void fake_sendrecv_fermion_borders_1Dcut(lnh_vec3_soa ** lnh_fermions){
   // NOTICE YOU HAVE TO CHANGE MYRANK TO USE THIS FUNCTION
  
  if(NRANKS_X != 1 || NRANKS_Y != 1 || NRANKS_Z != 1)
     printf("THIS SETUP IS NOT SALAMINO-LIKE!!!\n communication of fermion borders will FAIL!!\n");

   // PREAMBLE
   // This function is written taking the following assumptions:
   // 1. The domain is divided only along one direction, which is 
   //  the 'slowest' (strong assumption);
   // 2. That direction is the T direction (a bit weaker assumption, 
   //  changing that requires a complete redefinition of the site 
   //  ordering, which can be, perhaps, easily done by redefining the
   //  various '*snum' functions in Geometry/geometry.cc, and by 
   //  writing a tool to 'transpose' the configurations which have
   //  been written in a 'standard' ordering.

   //must be done for the three components of the fermion.

  int mpi_rankR, mpi_rankL;
  mpi_rankL = myrank - 1;
  if ( mpi_rankL == -1 ) mpi_rankL = NRANKS-1;
  mpi_rankR = myrank + 1;
  if ( mpi_rankR == NRANKS ) mpi_rankR = 0;
  int slab_sizeh = (LNH_NXH * LNH_NY * LNH_NZ * T_HALO); // NOTICE THERE IS LNH_NXH
  // no. of 'fermion' point in each slab)

//  MPI_Status status;

  // CHECK if in
  // pragma acc update host(lnh_fermion->c0[slab_size:slab_size])   
  // et simila
  // the size 'slab_size' is correct or if it must be multiplied 
  // by some factor

// Sendrecv c0
  printf("c0, slab_sizeh %d, myrank %d , L %d R %d, NRANKS %d\n", slab_sizeh, myrank, mpi_rankL, mpi_rankR, NRANKS);//DEBUG

  fakeMPI_send((void*) &(lnh_fermions[myrank]->c0[slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankL]->c0[LNH_SIZEH-slab_sizeh])); 
  fakeMPI_recv((void*) &(lnh_fermions[myrank]->c0[LNH_SIZEH-slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankR]->c0[slab_sizeh])); 
  fakeMPI_send((void*) &(lnh_fermions[myrank]->c0[LNH_SIZEH-2*slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) lnh_fermions[mpi_rankR]->c0); 
  fakeMPI_recv((void*) lnh_fermions[myrank]->c0, 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankL]->c0[LNH_SIZEH-2*slab_sizeh])); 

// Sendrecv c1

  printf("c1\n");//DEBUG
  fakeMPI_send((void*) &(lnh_fermions[myrank]->c1[slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankL]->c1[LNH_SIZEH-slab_sizeh])); 
  fakeMPI_recv((void*) &(lnh_fermions[myrank]->c1[LNH_SIZEH-slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankR]->c1[slab_sizeh])); 
  fakeMPI_send((void*) &(lnh_fermions[myrank]->c1[LNH_SIZEH-2*slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) lnh_fermions[mpi_rankR]->c1); 
  fakeMPI_recv((void*) lnh_fermions[myrank]->c1, 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankL]->c1[LNH_SIZEH-2*slab_sizeh])); 

// Sendrecv c2
  printf("c2\n");//DEBUG
   
  fakeMPI_send((void*) &(lnh_fermions[myrank]->c2[slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankL]->c2[LNH_SIZEH-slab_sizeh])); 
  fakeMPI_recv((void*) &(lnh_fermions[myrank]->c2[LNH_SIZEH-slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankR]->c2[slab_sizeh])); 
  fakeMPI_send((void*) &(lnh_fermions[myrank]->c2[LNH_SIZEH-2*slab_sizeh]), 2*slab_sizeh, sizeof(double), (void*) lnh_fermions[mpi_rankR]->c2); 
  fakeMPI_recv((void*) lnh_fermions[myrank]->c2, 2*slab_sizeh, sizeof(double), (void*) &(lnh_fermions[mpi_rankL]->c2[LNH_SIZEH-2*slab_sizeh])); 

}


#endif

