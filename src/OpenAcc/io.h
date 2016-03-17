#ifndef IO_H_
#define IO_H_

#include "./struct_c_def.h"
#include "./alloc_vars.h"

#ifdef MULTIDEVICE
#include "../Mpi/multidev.h"
#endif 

#include "../Mpi/communications.h"



void print_su3_soa_ASCII(global_su3_soa * conf, const char* nomefile,int conf_id_iter);
void print_su3_soa_ildg_binary(global_su3_soa * conf, const char* nomefile,
        int conf_id_iter);


//WRAPPER
inline void save_conf(global_su3_soa * const conf, const char* nomefile,int conf_id_iter,
        int use_ildg){

#ifdef MULTIDEVICE
    if(devinfo.myrank != 0){
        printf("MPI%02d: Rank is not allowed to use this function!\n",devinfo.myrank );
        printf("ERROR: %s:%d\n",__FILE__, __LINE__);
        exit(1);
    }
#endif 

    if(use_ildg){
        printf("Using ILDG format.\n");
        print_su3_soa_ildg_binary(conf,nomefile,conf_id_iter);
    }
    else 
        print_su3_soa_ASCII(conf,nomefile,conf_id_iter);

}


// WRAPPER OF WRAPPER
inline void save_conf_wrapper(su3_soa* conf, const char* nomefile,int conf_id_iter, 
        int use_ildg){

#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        int irank;
        for(irank = 1 ; irank < devinfo.nranks; irank++)
            recv_loc_subconf_from_rank(conf_rw,irank,irank);
        recv_loc_subconf_from_buffer(conf_rw,conf,0);
        save_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
    }
    else  send_lnh_subconf_to_master(conf,devinfo.myrank);

#else 
    recv_loc_subconf_from_buffer(conf_rw,conf,0);
    save_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
#endif


}


int read_su3_soa_ASCII(global_su3_soa * conf, const char* nomefile,int * conf_id_iter );
int read_su3_soa_ildg_binary(global_su3_soa * conf, const char* nomefile,
        int * conf_id_iter );


//WRAPPER
inline int read_conf(global_su3_soa * conf, const char* nomefile,int * conf_id_iter,
        int use_ildg ){

#ifdef MULTIDEVICE
    if(devinfo.myrank != 0){
        printf("MPI%02d: Rank is not allowed to use this function!\n",devinfo.myrank );
        printf("ERROR: %s:%d\n",__FILE__, __LINE__);
        exit(1);
    }
#endif 

    printf("Reading configuration. ");
    if(use_ildg){ 
        printf("Using ILDG format.\n");
        return read_su3_soa_ildg_binary(conf,nomefile,conf_id_iter);
    }
    else{
        printf("Using ascii format.\n");
        return read_su3_soa_ASCII(conf,nomefile,conf_id_iter );
    }


}


// WRAPPER OF WRAPPER
inline int read_conf_wrapper(su3_soa* conf, const char* nomefile,int * conf_id_iter, 
        int use_ildg)
{

    int error =  0;
#ifdef MULTIDEVICE

    if(devinfo.myrank == 0){
        error = read_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
        if(!error) {
            send_lnh_subconf_to_buffer(conf_rw,conf,0);
            int irank;
            for(irank = 1 ; irank < devinfo.nranks; irank++)
                send_lnh_subconf_to_rank(conf_rw,irank);
        }
    }
    else receive_lnh_subconf_from_master(conf);

#else 
    error = read_conf(conf_rw, nomefile,conf_id_iter, use_ildg);
    if(!error)  send_lnh_subconf_to_buffer(conf_rw,conf,0);
#endif

    return error;

}


// READ AND WRITE FOR THE POOR MAN
//void write_conf_binary_chunks(su3_soa * conf, const char *rootname,int nchunks);
//void read_conf_binary_chunks(su3_soa * conf, const char *rootname,int chunk);

#endif 
