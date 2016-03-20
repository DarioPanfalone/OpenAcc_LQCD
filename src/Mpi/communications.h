#ifndef COMMUNICATIONS_H_
#define COMMUNICATIONS_H_

#include "./geometry_multidev.h"
#include "./geometry_multidev.h"
#include "../OpenAcc/struct_c_def.h"

#define GAUGE_HALO HALO_WIDTH
#define FERMION_HALO 1

#ifdef MULTIDEVICE
#include <mpi.h>
// fermion border communications (only FERMION HALO thick)
void communicate_fermion_borders(vec3_soa *lnh_fermion);
void communicate_fermion_borders_async(vec3_soa *lnh_fermion, 
        MPI_Request* send_border_requests, MPI_Request* recv_border_requests);


// gauge conf border communication (only GAUGE_HALO thick)
void communicate_su3_borders(su3_soa* lnh_conf);
// communications of gl3 quantities (e.g. Sigma for stouting)
void communicate_gl3_borders(su3_soa* lnh_conf);

// used just at the beginning of MD trajectory (only GAUGE_HALO thick)
void communicate_thmat_soa_borders(thmat_soa* lnh_momenta);
// force communication (only GAUGE_HALO thick)
// this is performance critical and an async version should be produced
void communicate_tamat_soa_borders(tamat_soa* lnh_ipdot);

// communication of lattice chunks, for file IO only on master
// chunks, conf
void send_lnh_subconf_to_rank(global_su3_soa *gl_soa_conf, 
        int target_rank);
void recv_loc_subconf_from_rank(global_su3_soa *gl_soa_conf, 
        int target_rank, int tag);
void send_lnh_subconf_to_master(su3_soa *lnh_soa_conf, int tag);
void receive_lnh_subconf_from_master(su3_soa* lnh_su3_conf);

// chunks, fermions
void send_lnh_subfermion_to_rank(global_vec3_soa *gl_soa_ferm, 
        int target_rank);
void recv_loc_subfermion_from_rank(global_vec3_soa *gl_soa_ferm,
        int target_rank);


void send_lnh_subfermion_to_master(vec3_soa *lnh_ferm, int tag);
void receive_lnh_subfermion_from_master(vec3_soa* lnh_ferm);




#endif



// internal communication on master rank
void send_lnh_subconf_to_buffer(global_su3_soa *gl_soa_conf, 
        su3_soa *lnh_conf, int target_rank);
void recv_loc_subconf_from_buffer(global_su3_soa *gl_soa_conf,
        su3_soa* lnh_conf, int target_rank);

void send_lnh_subfermion_to_buffer(global_vec3_soa *gl_soa_ferm,
        vec3_soa *lnh_ferm, int target_rank);
void recv_loc_subfermion_from_buffer(global_vec3_soa *gl_soa_ferm,
        vec3_soa* lnh_ferm, int target_rank);





#endif
