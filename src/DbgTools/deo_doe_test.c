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
#include "../Include/montecarlo_parameters.h"
#include "../Include/debug.h"
#include "../OpenAcc/deviceinit.h" 

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
    const char fermionname_doe[50] = "test_fermion_result_doe2";
    const char fermionname_deo[50] = "test_fermion_result_deo2";
    const char fermionname_fulldirac[50] = "test_fermion_result_fulldirac2";
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

    initrand((unsigned int) mc_params.seed+devinfo.myrank);
    verbosity_lv = debug_settings.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters))
        printf("Ignoring issues in init_ferm_params,\
                this is a deo-doe test.\n");

    if(NPS_tot == 0) 
    {
        printf("Scusami ma sono handicappato e ho bisogno che metti almeno un fermione nell'input file.\n"); // GOLIARDIA
        exit(1);
    }
    else {
        printf("buono capo grazie per il fermone! %d\n", NPS_tot); // GOLIARDIA
    }
    //

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
    printf("MPI%02d: Selecting device.\n");
#ifdef MULTIDEVICE
    select_init_acc_device(my_device_type, devinfo.myrank%devinfo.proc_per_node);
#else
    select_init_acc_device(my_device_type, devinfo.single_dev_choice);
#endif
    printf("Device Selected : OK \n");
#endif
    int conf_id_iter = 0;
    if(!read_conf_wrapper(conf_acc,mc_params.save_conf_name,
                &conf_id_iter,debug_settings.use_ildg)){
        // READS ALSO THE conf_id_iter
        printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
                devinfo.myrank, mc_params.save_conf_name);

    }
    else{
        generate_Conf_cold(conf_acc,mc_params.eps_gen);
        printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
                devinfo.myrank);
        save_conf_wrapper(conf_acc,mc_params.save_conf_name,0,debug_settings.use_ildg);
        if(debug_settings.use_ildg)
            printf("MPI%02d: You're using ILDG format.\n", devinfo.myrank);
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


#ifdef MULTIDEVICE
    communicate_fermion_borders_hostonly(ferm_chi_acc);
#endif

    print_vec3_soa(ferm_chi_acc,myfermionname);




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
            printf("MPI%02d: Ce l'hai quasi fatta, gino!\n", devinfo.myrank); // GOLIARDIA
        }
#ifndef __GNUC__
    shutdown_acc_device(my_device_type);
#endif
    printf("il tuo test del cazzo e' QUASI finito, %d %d\n", FERMION_HALO, HALO_WIDTH); // GOLIARDIA
#ifdef MULTIDEVICE
    MPI_Finalize();
#endif 

    mem_free();

    printf("Bene, il tuo test del cazzo e' finito, sei contento adesso?\n"); // GOLIARDIA

    return 0; 


}
