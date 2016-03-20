#include <stdio.h>
#include <stdlib.h>

#ifndef __GNUC__
#include "openacc.h"
#endif

#include "../Include/common_defines.h"
#include "../Include/setting_file_parser.h"
#include "../Include/fermion_parameters.h"
#include "../OpenAcc/fermion_matrix.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/io.h"
#include "../OpenAcc/random_assignement.h"
#include "../Rand/random.h"
#include "./dbgtools.h"
#include "../OpenAcc/action.h"
#include "../Include/markowchain.h"

#include "../Mpi/multidev.h"
#ifdef MULTIDEVICE
#include <mpi.h>
#endif

//double casuale(void);


int id_iter;
int verbosity_lv;

int main(int argc, char* argv[]){

    const char confname[50] = "test_conf";
    const char fermionname[50] = "test_fermion";
    const char fermionname_doe[50] = "test_fermion_result_doe";
    const char fermionname_deo[50] = "test_fermion_result_deo";
    const char fermionname_fulldirac[50] = "test_fermion_result_fulldirac";
    char myconfname             [50];
    char myfermionname          [50];
    char myfermionname_doe      [50];
    char myfermionname_deo      [50];
    char myfermionname_fulldirac[50];


    act_params.stout_steps = 0;
    NPS_tot = 1;

#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
#endif
    fflush(stdout);
    printf("DEODOE test\n");
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    printf("WELCOME! \n");
    // READ input file.
    set_global_vars_and_fermions_from_input_file(argv[1]);

    if(NPS_tot == 0) 
    {
        printf("ERROR, at least a fermion flavour is needed.\n");
        exit(1);
    }
    //
    initrand((unsigned int) mkwch_pars.seed+devinfo.myrank);
    verbosity_lv = mkwch_pars.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters))
        printf("Ignoring issues in init_ferm_params,\
                this is a deo-doe test.\n");
#ifdef MULTIDEVICE
    init_multidev1D(&devinfo);
#else
    devinfo.myrank = 0;
    devinfo.nranks = 1;
#endif

    mem_alloc();
    printf("Allocazione della memoria : OK \n");
    compute_nnp_and_nnm_openacc();
    printf("nn computation : OK \n");
    init_all_u1_phases(backfield_parameters,fermions_parameters);
    printf("u1_backfield initialization : OK \n");
    sprintf(myconfname             ,"%s_MPI%02d",confname             ,devinfo.myrank);
    sprintf(myfermionname          ,"%s_MPI%02d",fermionname          ,devinfo.myrank);
    sprintf(myfermionname_doe      ,"%s_MPI%02d",fermionname_doe      ,devinfo.myrank);
    sprintf(myfermionname_deo      ,"%s_MPI%02d",fermionname_deo      ,devinfo.myrank);
    sprintf(myfermionname_fulldirac,"%s_MPI%02d",fermionname_fulldirac,devinfo.myrank);



#ifndef __GNUC__
    //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////////////////////
    // NVIDIA GPUs
    acc_device_t my_device_type = acc_device_nvidia;
    // AMD GPUs
    // acc_device_t my_device_type = acc_device_radeon;
    // Intel XeonPhi
    //acc_device_t my_device_type = acc_device_xeonphi;
    // Select device ID
    int dev_index = 0;
    SELECT_INIT_ACC_DEVICE(my_device_type, dev_index);
    printf("Device Selected : OK \n");
#endif
    int conf_id_iter = 0;
    if(!read_conf_wrapper(conf_acc,mkwch_pars.save_conf_name,
                &conf_id_iter,mkwch_pars.use_ildg)){
        // READS ALSO THE conf_id_iter
        printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
                devinfo.myrank, mkwch_pars.save_conf_name);

    }
    else{
        generate_Conf_cold(conf_acc,mkwch_pars.eps_gen);
        printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
                devinfo.myrank);
        save_conf_wrapper(conf_acc,mkwch_pars.save_conf_name,0,0);
        conf_id_iter=0;
    }


    // init fermion
    //generate_vec3_soa_gauss(ferm_chi_acc);
    //print_vec3_soa(ferm_chi_acc,"ferm_chi_acc");
    if(!read_vec3_soa_wrapper(ferm_chi_acc,fermionname )){
        printf("MPI%02d - Fermion READ : OK \n",devinfo.myrank);
    }else{
        printf("MPI%02d - Fermion READ : FAIL, generating\n",devinfo.myrank);
        generate_vec3_soa_gauss(ferm_chi_acc);
        printf("MPI%02d - Fermion GENERATED : OK\n",devinfo.myrank);
        printf("MPI%02d - Writing file %s.\n",devinfo.myrank, fermionname);
        print_vec3_soa_wrapper(ferm_chi_acc,fermionname);
    }

    communicate_fermion_borders(ferm_chi_acc);
    print_vec3_soa(ferm_chi_acc,myfermionname);
    //    ferm_phi_acc->c0[0] = 5;
    //#pragma acc data  copyin(conf_acc[0:8]) copyin(ferm_chi_acc[0:1])\
    create(ferm_phi_acc[0:1])  copyin(u1_back_phases[0:8*NDiffFlavs]) \
        create(kloc_s[0:1])
#pragma acc data  copy(conf_acc[0:8]) copy(ferm_chi_acc[0:1])\
        copy(ferm_phi_acc[0:1])  copy(u1_back_phases[0:8*NDiffFlavs]) \
        copy(kloc_s[0:1])
        {

            printf("Multiplication by Doe...");
            acc_Doe(conf_acc, ferm_phi_acc, ferm_chi_acc, fermions_parameters[0].phases);
            printf("Writing file %s.\n", fermionname_doe);

#pragma acc update host(ferm_phi_acc[0:1])
            print_vec3_soa_wrapper(ferm_phi_acc,fermionname_doe);
            print_vec3_soa(ferm_phi_acc,myfermionname_doe);

            printf("Multiplication by Deo...");
            acc_Deo(conf_acc, ferm_phi_acc, ferm_chi_acc, fermions_parameters[0].phases) ;
            printf("Writing file %s.\n", fermionname_deo);
#pragma acc update host(ferm_phi_acc[0:1])
            print_vec3_soa_wrapper(ferm_phi_acc,fermionname_deo);
            print_vec3_soa(ferm_phi_acc,myfermionname_deo);

            printf("Multiplication by M^\\dagM+m^2...");
            fermion_matrix_multiplication(conf_acc, ferm_phi_acc, ferm_chi_acc, kloc_s, &fermions_parameters[0]) ;
            printf("Writing file %s.\n", fermionname_fulldirac);
#pragma acc update host(ferm_phi_acc[0:1])
            print_vec3_soa_wrapper(ferm_phi_acc,fermionname_fulldirac);
            print_vec3_soa(ferm_phi_acc,myfermionname_fulldirac);
        }
#ifndef __GNUC__
    SHUTDOWN_ACC_DEVICE(my_device_type);
#endif
#ifdef MULTIDEVICE
    MPI_Finalize();
#endif 

    mem_free();

    return 0; 


}
