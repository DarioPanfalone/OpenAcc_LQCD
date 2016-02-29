#ifndef COMMUNICATIONS_C_
#define COMMUNICATIONS_C_


#ifdef MULTIDEVICE

#include "./multidevice.h"
#include "../OpenAcc/single_types.h"

THE FOLLOWING MUST BE FIXED // DEBUG

// fermions

void sendrecv_fermion_borders_1Dcut(vec3_soa *lnh_fermion,
        int rankL, int rankR, int thickness, int offset){
   // NOTICE YOU HAVE TO SET MYRANK CORRECTLY TO USE THIS FUNCTION
  
  if(NRANKS_D0 != 1 || NRANKS_D1 != 1 || NRANKS_D2 != 1)
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

  // no. of 'fermion' point in each slab)
  int slab_sizeh = (LNH_N0H * LNH_N1 * LNH_N2)*thickness;
  int offset_size =  (LNH_N0H * LNH_N1 * LNH_N2) * offset;
  // NOTICE THERE IS LNH_NXH
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
#pragma acc update host(tmpc[offset_size:slab_sizeh])
#endif
          MPI_Sendrecv((void*) &(c[ii][offset_size]),2*slab_sizeh,MPI_DOUBLE,
                  rankL,sendtag,
                  (void*) &(c[ii][sizeh-offset_size]),2*slab_sizeh,MPI_DOUBLE,
                  rankR,recvtag,
                  MPI_COMM_WORLD, &status);
#ifndef USE_MPI_CUDA_AWARE
#pragma acc update device(tmpc[(sizeh-offset_size):slab_sizeh])
#endif

          sendtag = ii+3;
          recvtag = ii+3;

#ifndef USE_MPI_CUDA_AWARE
#pragma acc update host(tmpc[(sizeh-offset_size-slab_sizeh):slab_sizeh])
#endif
          MPI_Sendrecv((void*) &(c[ii][sizeh-offset_size-slab_sizeh]),
                  2*slab_sizeh,MPI_DOUBLE,
                  rankR,sendtag,
                  (void*) &(c[ii][offset_size-slab_sizeh]),
                  2*slab_sizeh,MPI_DOUBLE,
                  rankL,recvtag,
                  MPI_COMM_WORLD, &status);
#ifndef USE_MPI_CUDA_AWARE
#pragma acc update device(tmpc[offset_size-slab_sizeh:slab_sizeh])
#endif
      }
#ifdef USE_MPI_CUDA_AWARE
  }
#endif
 }
inline void communicate_fermion_borders(vec3_soa *lnh_fermion){ //WRAPPER

    // NOTICE: GEOMETRY MUST BE SET UP BEFORE!!
    MPI_Barrier(MPI_COMM_WORLD);
    sendrecv_fermion_borders_1Dcut(lnh_fermion,
            multidev_info.myrank_L, multidev_info.myrank_R, 
            FERMION_HALO, HALO );
    MPI_Barrier(MPI_COMM_WORLD);
}

#ifdef USE_MPI_CUDA_AWARE  // This can work only if cuda aware MPI is used.
                           // If instead you wanted to use the update
                           // directive, this cannot work with MPI_Irecv
                           // because the directive should wait for the 
                           // completion the recv, but Irecv is not bloking!!
                          
void sendrecv_fermion_borders_1Dcut_async(vec3_soa *lnh_fermion, 
        int rankL, int rankR, 
        int thickness, int offset,
        MPI_Request* send_border_requests, 
        MPI_Request* recv_border_requests){
   // NOTICE YOU HAVE TO SET MYRANK CORRECTLY TO USE THIS FUNCTION
   // NOTICE send_border_requests recv_border_requests are both 
   // 6-elements long
  
  if(NRANKS_D0 != 1 || NRANKS_D1 != 1 || NRANKS_D2 != 1)
     printf("THIS SETUP IS NOT SALAMINO-LIKE!!!\n communication of fermion borders will FAIL!!\n");

   //SEE PREAMBLE FOR sendrecv_fermion_borders_1Dcut()
   //must be done for the three components of the fermion.

  // no. of 'fermion' point in each slab)
  int slab_sizeh = (LNH_N0H * LNH_N1 * LNH_N2 )*thickness;
  int offset_size =  (LNH_N0H * LNH_N1 * LNH_N2) * offset;
  // NOTICE THERE IS LNH_NXH

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
          MPI_Isend((void*) &(c[ii][offset_size]),2*slab_sizeh,MPI_DOUBLE,
                  rankL,sendtag,MPI_COMM_WORLD,
                  &(send_border_requests[ii]));
          MPI_Irecv((void*) &(c[ii][sizeh-offset_size]),2*slab_sizeh,
                  MPI_DOUBLE,
                  rankR,recvtag,MPI_COMM_WORLD,
                  &(recv_border_requests[ii]));


          sendtag = 3+ii;
          recvtag = 3+ii;

          MPI_Isend((void*) &(c[ii][sizeh-offset_size-slab_sizeh]),
                  2*slab_sizeh,MPI_DOUBLE,
                  rankR,sendtag,MPI_COMM_WORLD,
                  &(send_border_requests[3+ii]));
          MPI_Irecv( (void*) &(c[ii][offset_size-slab_sizeh]),
                  2*slab_sizeh,MPI_DOUBLE,
                  rankL,recvtag,MPI_COMM_WORLD,
                  &(recv_border_requests[3+ii]));

      }
  }
}

inline void communicate_fermion_borders_async(vec3_soa *lnh_fermion, 
        MPI_Request* send_border_requests, 
        MPI_Request* recv_border_requests){ //WRAPPER
   
    // NOTICE: GEOMETRY MUST BE SET UP BEFORE!!
    // NOTICE send_border_requests recv_border_requests are both 
    // 6-elements long
    sendrecv_fermion_borders_1Dcut_async(lnh_fermion,
            multidev_info.myrank_L, multidev_info.myrank_R,
            FERMION_HALO, HALO,
            send_border_requests,
            recv_border_requests );
}

#endif



// GAUGE COMMS

// borders
inline void communicate_su3_borders(su3_soa* lnh_conf, 
        int thickness, int offset){

    // NOTICE: GEOMETRY MUST BE SET UP BEFORE!!
    // note for async:
    // NOTICE send_border_requests recv_border_requests are both 
    // 12*8-elements long

    for(int c = 0 ; c < 8 ; c++){ // Remember lnh_conf has 8 components
        sendrecv_fermion_borders_1Dcut(&(lnh_conf[c].r0),
            multidev_info.myrank_L, multidev_info.myrank_R,
            thickness,offset);
        sendrecv_fermion_borders_1Dcut(&(lnh_conf[c].r1),
            multidev_info.myrank_L, multidev_info.myrank_R,
            thickness,offset);
    }

}
inline void communicate_gl3_borders(su3_soa* lnh_conf){

    // NOTICE: GEOMETRY MUST BE SET UP BEFORE!!
    // note for async:
    // NOTICE send_border_requests recv_border_requests are both 
    // 18*8-elements long

    for(int c = 0 ; c < 8 ; c++){ // Remember lnh_conf has 8 components
        sendrecv_fermion_borders_1Dcut(&(lnh_conf[c].r0),
            multidev_info.myrank_L, multidev_info.myrank_R,1,HALO);
        sendrecv_fermion_borders_1Dcut(&(lnh_conf[c].r1),
            multidev_info.myrank_L, multidev_info.myrank_R,1,HALO);
        sendrecv_fermion_borders_1Dcut(&(lnh_conf[c].r2),
            multidev_info.myrank_L, multidev_info.myrank_R,1,HALO);
    }

}


// chunks
void send_lnh_subconf_to_rank(global_su3_soa *gl_soa_conf, int target_rank){
    // USE ONLY FROM MASTER RANK
    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);
    // building sublattice duplicate, target_conf
    su3_soa* target_su3_soa = (su3_soa*) malloc(8*sizeof(su3_soa)); 

    int tg_lnh_x,tg_lnh_y,tg_lnh_z,tg_lnh_t,dir; //target-lnh coordinates
    // and link direction
    // Copying all relevant links into the sublattice
    for(dir =0; dir < 4; dir++)
        for(tg_lnh_t=0;tg_lnh_t<LNH_NT; tg_lnh_t++)
            for(tg_lnh_z=0;tg_lnh_z<LNH_NZ; tg_lnh_z++)
                for(tg_lnh_y=0;tg_lnh_y<LNH_NY; tg_lnh_y++)
                    for(tg_lnh_x=0;tg_lnh_x<LNH_NX; tg_lnh_x++){

                        //        int gtsp; // global target site parity
                        int tsprlo ; // target site parity respect (to his) local origin;

                        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
                        int target_lnh_snum = lnh_to_lnh_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

                        tsprlo = (X_HALO+Y_HALO+Z_HALO+T_HALO+ tg_lnh_t+tg_lnh_z+tg_lnh_y+tg_lnh_x)%2;

                        single_su3 aux;
                        su3_from_global_su3_soa(gl_soa_conf,dir,tsprlo,target_gl_snum);

                        su3_into_su3_soa(aux,target_su3_soa,dir,tsprlo,target_lnh_snum);

                    }

    //sending the subconfiguration
    MPI_Send(target_su3_soa, 2*4*(6*3)*LNH_SIZEH,MPI_DOUBLE, target_rank , target_rank, MPI_COMM_WORLD); // tag = target_rank
    // In case we remove the third line, this is possibly the thing to do
    //MPI_Send(target_su3_soa, 2*4*(6*2)*LNH_SIZEH,MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
    //Or maybe do multiple sends in a more complicated way
    // ^^ CHECK
    free(target_su3_soa);
}
void recv_loc_subconf_from_rank(global_su3_soa *gl_soa_conf, int target_rank, int tag){
    // USE ONLY FROM MASTER RANK
    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);
    /*
       int target_loc_origin_parity = (target_gl_loc_origin_from_rank.x +
       target_gl_loc_origin_from_rank.y +
       target_gl_loc_origin_from_rank.z +
       target_gl_loc_origin_from_rank.t)%2;// with the current setup, this should always be 0

       if(target_loc_origin_parity) printf("Problems\n");
       */
    // building sublattice duplicate, target_conf
    lnh_su3_soa* target_su3_soa = (lnh_su3_soa*) malloc(8*sizeof(lnh_su3_soa)); 
    MPI_Recv(target_su3_soa, 2*4*(6*3)*LNH_SIZEH,MPI_DOUBLE,target_rank,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int tg_loc_x,tg_loc_y,tg_loc_z,tg_loc_t,dir; //target-loc coordinates
    // and link direction
    // Copying all relevant links from the sublattice to the global lattice
    for(dir =0; dir < 4; dir++)
        for(tg_loc_t=0;tg_loc_t<LOC_NT; tg_loc_t++)
            for(tg_loc_z=0;tg_loc_z<LOC_NZ; tg_loc_z++)
                for(tg_loc_y=0;tg_loc_y<LOC_NY; tg_loc_y++)
                    for(tg_loc_x=0;tg_loc_x<LOC_NX; tg_loc_x++){
    
                        int tg_lnh_x,tg_lnh_y,tg_lnh_z,tg_lnh_t; //target-lnh coordinates
                        tg_lnh_x = tg_loc_x + X_HALO;
                        tg_lnh_y = tg_loc_y + Y_HALO;
                        tg_lnh_z = tg_loc_z + Z_HALO;
                        tg_lnh_t = tg_loc_t + T_HALO;

                        //        int gtsp; // global target site parity
                        int tsprlo ; // target site parity respect (to his) local origin;

                        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
                        int target_lnh_snum = lnh_to_lnh_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

                        tsprlo = (X_HALO+Y_HALO+Z_HALO+T_HALO+ tg_lnh_t+tg_lnh_z+tg_lnh_y+tg_lnh_x)%2;

                        su3 aux = su3_from_su3_soa(target_su3_soa,dir,tsprlo,target_lnh_snum);

                        su3_into_global_su3_soa(aux,gl_soa_conf,dir,tsprlo,target_gl_snum);

                    }

    free(target_su3_soa);
}
void send_lnh_subconf_to_master(su3_soa *lnh_soa_conf, int tag){
   //sending the subconfiguration
    MPI_Send(lnh_soa_conf, 2*4*(6*3)*LNH_SIZEH,MPI_DOUBLE, 0, tag , MPI_COMM_WORLD);
}
inline void receive_lnh_subconf_from_master(su3_soa* lnh_su3_conf){

    MPI_Recv(lnh_su3_conf, 2*4*(6*3)*LNH_SIZEH,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD, MPI_STATUS_IGNORE); // tag = myrank
    // In case we remove the third line possibly we have to do something different
    // ^^ CHECK
}

// only for the master rank
void send_lnh_subconf_to_buffer(global_su3_soa *gl_soa_conf, su3_soa *lnh_conf, int target_rank){
    // USE ONLY FROM MASTER RANK


    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);
    // building sublattice duplicate, target_conf
    lnh_su3_soa* target_su3_soa = lnh_conf;

    int tg_lnh_x,tg_lnh_y,tg_lnh_z,tg_lnh_t,dir; //target-lnh coordinates
    // and link direction
    // Copying all relevant links into the sublattice
    for(dir =0; dir < 4; dir++)
        for(tg_lnh_t=0;tg_lnh_t<LNH_NT; tg_lnh_t++)
            for(tg_lnh_z=0;tg_lnh_z<LNH_NZ; tg_lnh_z++)
                for(tg_lnh_y=0;tg_lnh_y<LNH_NY; tg_lnh_y++)
                    for(tg_lnh_x=0;tg_lnh_x<LNH_NX; tg_lnh_x++){

                        //        int gtsp; // global target site parity
                        int tsprlo ; // target site parity respect (to his) local origin;

                        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
                        int target_lnh_snum = lnh_to_lnh_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

                        tsprlo = (X_HALO+Y_HALO+Z_HALO+T_HALO+ tg_lnh_t+tg_lnh_z+tg_lnh_y+tg_lnh_x)%2;
                        //      gtsp = (target_loc_origin_parity + tsprlo )%2;

                        //      printf("%d %d  %d  %d  %d  %d  %d ",dir, tg_lnh_t, tg_lnh_z, tg_lnh_y, tg_lnh_x, target_gl_snum, tsprlo);
                        su3 aux = su3_from_global_su3_soa(gl_soa_conf,dir,tsprlo,target_gl_snum);

                        //    printf("ciao \n");
                        //        printf("ciao ");
                        //      print_su3(aux);
                        su3_into_su3_soa(aux,target_su3_soa,dir,tsprlo,target_lnh_snum);
                        //        printf("ciao \n");

                    }

    /*
    //sending the subconfiguration
    MPI_Send(target_su3_soa, 2*4*(6*3)*LNH_SIZEH,MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
    // In case we remove the third line, this is possibly the thing to do
    //MPI_Send(target_su3_soa, 2*4*(6*2)*LNH_SIZEH,MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
    //Or maybe do multiple sends in a more complicated way
    // ^^ CHECK
    */ // COMMENTED FOR TEST
}
void recv_loc_subconf_from_buffer(global_su3_soa *gl_soa_conf, su3_soa* lnh_conf, int target_rank){
    // USE ONLY FROM MASTER RANK
    //target sublattice information
    vec4int target_gl_loc_origin_from_rank = gl_loc_origin_from_rank(target_rank);
    // building sublattice duplicate, target_conf
    su3_soa* target_su3_soa = lnh_conf; 

    int tg_loc_x,tg_loc_y,tg_loc_z,tg_loc_t,dir; //target-loc coordinates
    // and link direction
    // Copying all relevant links from the sublattice to the global lattice
    for(dir =0; dir < 4; dir++)
        for(tg_loc_t=0;tg_loc_t<LOC_NT; tg_loc_t++)
            for(tg_loc_z=0;tg_loc_z<LOC_NZ; tg_loc_z++)
                for(tg_loc_y=0;tg_loc_y<LOC_NY; tg_loc_y++)
                    for(tg_loc_x=0;tg_loc_x<LOC_NX; tg_loc_x++){
    
                        int tg_lnh_x,tg_lnh_y,tg_lnh_z,tg_lnh_t; //target-lnh coordinates
                        tg_lnh_x = tg_loc_x + X_HALO;
                        tg_lnh_y = tg_loc_y + Y_HALO;
                        tg_lnh_z = tg_loc_z + Z_HALO;
                        tg_lnh_t = tg_loc_t + T_HALO;

                        //        int gtsp; // global target site parity
                        int tsprlo ; // target site parity respect (to his) local origin;

                        int target_gl_snum = target_lnh_to_gl_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t, target_gl_loc_origin_from_rank);
                        int target_lnh_snum = lnh_to_lnh_snum(tg_lnh_x, tg_lnh_y, tg_lnh_z, tg_lnh_t);

                        tsprlo = (X_HALO+Y_HALO+Z_HALO+T_HALO+ tg_lnh_t+tg_lnh_z+tg_lnh_y+tg_lnh_x)%2;

                        su3 aux = su3_from_su3_soa(lnh_conf,dir,tsprlo,target_lnh_snum);

                        su3_into_global_su3_soa(aux,gl_soa_conf,dir,tsprlo,target_gl_snum);

                    }

}

#endif




#endif
