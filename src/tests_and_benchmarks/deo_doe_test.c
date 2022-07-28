#include <stdio.h>
#include <stdlib.h>

#ifndef __GNUC__
#include "openacc.h"
#endif
#include <sys/time.h>


#include "../Include/common_defines.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"
#include "../Include/setting_file_parser.h"
#include "../Include/tell_geom_defines.h"
#include "../Mpi/multidev.h"
#include "../OpenAcc/action.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/deviceinit.h" 
#include "../OpenAcc/fermion_matrix.h"
#include "../OpenAcc/float_double_conv.h"
#include "../OpenAcc/io.h"
#include "../OpenAcc/random_assignement.h"
#include "../OpenAcc/sp_alloc_vars.h"
#include "../OpenAcc/sp_fermion_matrix.h"
#include "../Rand/random.h"
#include "../DbgTools/dbgtools.h"
#include "../DbgTools/sp_dbgtools.h"
#include "./test_and_benchmarks.h"
#include "../OpenAcc/alloc_settings.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#include "../Mpi/sp_communications.h"
#include <mpi.h>
#endif


#include "../Include/stringify.h"

//double casuale(void);


int id_iter;
int verbosity_lv;

int main(int argc, char* argv[]){

    const char confname[50] = "test_conf";
    const char fermionname[50] = "test_fermion";
    const char fermionname_doe[50] = "test_fermion_result_doe2";
    const char fermionname_deo[50] = "test_fermion_result_deo2";
    const char fermionname_fulldirac[50] = "test_fermion_result_fulldirac2";
    const char fermionname_f[50] = "sp_test_fermion";
    const char fermionname_doe_f[50] = "sp_test_fermion_result_doe2";
    const char fermionname_deo_f[50] = "sp_test_fermion_result_deo2";
    const char fermionname_fulldirac_f[50] = "sp_test_fermion_result_fulldirac2";

    char myconfname             [50];
    char myfermionname          [50];
    char myfermionname_doe      [50];
    char myfermionname_deo      [50];
    char myfermionname_fulldirac[50];

    //char myconfname_f             [50];
    char myfermionname_f          [50];
    char myfermionname_doe_f      [50];
    char myfermionname_deo_f      [50];
    char myfermionname_fulldirac_f[50];





    act_params.stout_steps = 0;
    alloc_info.NPS_tot = 1;
    //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////////////////////
#ifndef __GNUC__
    acc_device_t my_device_type = acc_device_nvidia; // NVIDIA GPUs
    // acc_device_t my_device_type = acc_device_radeon; // AMD GPUs
    //acc_device_t my_device_type = acc_device_xeonphi; // Intel XeonPhi
#ifdef MULTIDEVICE
    {
        char* localRankStr = NULL;
        if ((localRankStr = getenv("OMPI_COMM_WORLD_RANK")) != NULL){
            printf("LocalRankStr: %s\n",localRankStr);
            devinfo.myrank = atoi(localRankStr);
        }
    }
    select_init_acc_device(my_device_type, devinfo.myrank%devinfo.proc_per_node);
    printf("MPI%02d: Device Selected : OK \n", devinfo.myrank );
#else
    select_init_acc_device(my_device_type, devinfo.single_dev_choice);
    printf("MPI%02d: Device Selected : OK \n", devinfo.single_dev_choice );
#endif
#endif


#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
#endif
    if(0 == devinfo.myrank){
        fflush(stdout);
        printf("DEODOE test\n");
        printf("commit: %s\n", xstr(COMMIT_HASH) );
        // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
        printf("WELCOME! \n");
    }
    // READ input file.
    //
    set_global_vars_and_fermions_from_input_file(argv[1]);


    if(! test_settings.parametersAreSet){
        if(0 == devinfo.myrank){
            printf("ERROR : Test parameters are not set in %s. Is a proper TestSetting section present?\n", argv[1]);
            printf("Exiting now!\n");
        }
        exit(1);
    }
    initrand((unsigned int) mc_params.seed+devinfo.myrank);
    verbosity_lv = debug_settings.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters) && (0==devinfo.myrank))
        printf("Ignoring issues in init_ferm_params,\
                this is a deo-doe test.\n");

    if(alloc_info.NPS_tot == 0) 
    {
        if(0 == devinfo.myrank){
            printf("ERROR: For technical reasons, at least a fermion in the input file is needed\n");
            printf("       (namely, the allocation of external U(1) phases.\n");
        }
        exit(1);
    }
    else {
        if(0 == devinfo.myrank)
            printf("N fermions found: %d, only one is necessary.\n", alloc_info.NPS_tot);
    }

    printf("Setting number of shifts to 1 - no shifts are needeed for this benchmark.\n");
    printf("                                but using 1 just to avoid errors.\n");
    alloc_info.maxNeededShifts = 1; // DEBUG
    printf("Setting stout level to zero - no stout is needed for this benchmark.\n");
    alloc_info.stoutAllocations = 0;
    act_params.stout_steps = 0;





#ifdef MULTIDEVICE
    init_multidev1D(&devinfo);
#else
    devinfo.myrank = 0;
    devinfo.nranks = 1;
#endif

    if(0==devinfo.myrank) print_geom_defines();
    mem_alloc_core();
    mem_alloc_extended();
    mem_alloc_core_f();
    mem_alloc_extended_f();

    printf("Allocazione della memoria : OK \n");
    compute_nnp_and_nnm_openacc();
#pragma acc enter data copyin(nnp_openacc)
#pragma acc enter data copyin(nnm_openacc)

    printf("nn computation : OK \n");
    init_all_u1_phases(backfield_parameters,fermions_parameters);
#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])
    printf("u1_backfield initialization : OK \n");
    sprintf(myconfname             ,"%s_MPI%02d",confname             ,devinfo.myrank);
    sprintf(myfermionname          ,"%s_MPI%02d",fermionname          ,devinfo.myrank);
    sprintf(myfermionname_doe      ,"%s_MPI%02d",fermionname_doe      ,devinfo.myrank);
    sprintf(myfermionname_deo      ,"%s_MPI%02d",fermionname_deo      ,devinfo.myrank);
    sprintf(myfermionname_fulldirac,"%s_MPI%02d",fermionname_fulldirac,devinfo.myrank);
    //sprintf(myconfname_f             ,"%s_MPI%02d",confname             ,devinfo.myrank);
    sprintf(myfermionname_f          ,"%s_MPI%02d",fermionname_f          ,devinfo.myrank);
    sprintf(myfermionname_doe_f      ,"%s_MPI%02d",fermionname_doe_f      ,devinfo.myrank);
    sprintf(myfermionname_deo_f      ,"%s_MPI%02d",fermionname_deo_f      ,devinfo.myrank);
    sprintf(myfermionname_fulldirac_f,"%s_MPI%02d",fermionname_fulldirac_f,devinfo.myrank);




    // INITIALIZING GAUGE CONFIGURATION
    int conf_id_iter = 0;
    if(!read_conf_wrapper(conf_acc[0],mc_params.save_conf_name,
                &conf_id_iter,debug_settings.use_ildg)){
        // READS ALSO THE conf_id_iter
        printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
                devinfo.myrank, mc_params.save_conf_name);

    }
    else{
        generate_Conf_cold(conf_acc[0],mc_params.eps_gen);
        printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
                devinfo.myrank);
        save_conf_wrapper(conf_acc[0],mc_params.save_conf_name,0,debug_settings.use_ildg);
        if(debug_settings.use_ildg)
            printf("MPI%02d: You're using ILDG format.\n", devinfo.myrank);
        conf_id_iter=0;
    }
#pragma acc update device(conf_acc[0:1][0:8])


    // init fermion
    if(!read_vec3_soa_wrapper(ferm_chi_acc,fermionname )){
        printf("MPI%02d - Fermion READ : OK \n",devinfo.myrank);
    }else{
        printf("MPI%02d - Fermion READ : FAIL, generating\n",devinfo.myrank);
        generate_vec3_soa_gauss(ferm_chi_acc);
        printf("MPI%02d - Fermion GENERATED : OK\n",devinfo.myrank);
        printf("MPI%02d - Writing file %s.\n",devinfo.myrank, fermionname);
        print_vec3_soa_wrapper(ferm_chi_acc,fermionname);
    }


#pragma acc update device(ferm_chi_acc[0:1])



    // double precision
    if(0 == devinfo.myrank){
        printf("####################\n");
        printf("# DOUBLE PRECISION #\n");
        printf("####################\n");
    }
    struct timeval t0,t1,t2,t3,t4,t5;
    int r;
    if(0 == devinfo.myrank)
        printf("Multiplication by Doe, %d times...\n", test_settings.deoDoeIterations);
    gettimeofday(&t0,NULL);
    for(r=0; r<test_settings.deoDoeIterations; r++)
        acc_Doe(conf_acc[0], ferm_phi_acc, ferm_chi_acc, fermions_parameters[0].phases);
    gettimeofday(&t1,NULL);

    if(test_settings.saveResults){
        if(0 == devinfo.myrank)
            printf("Writing file %s.\n", fermionname_doe);
#pragma acc update host(ferm_phi_acc[0:1])
        print_vec3_soa_wrapper(ferm_phi_acc,fermionname_doe);
        print_vec3_soa(ferm_phi_acc,myfermionname_doe);
    }

    if(0 == devinfo.myrank)
        printf("Multiplication by Deo, %d times...\n", test_settings.deoDoeIterations);
    gettimeofday(&t2,NULL);
    for(r=0; r<test_settings.deoDoeIterations; r++)
        acc_Deo(conf_acc[0], ferm_phi_acc, ferm_chi_acc,fermions_parameters[0].phases);
    gettimeofday(&t3,NULL);

    if(test_settings.saveResults){
        if(0 == devinfo.myrank)
            printf("Writing file %s.\n", fermionname_deo);
#pragma acc update host(ferm_phi_acc[0:1])
        print_vec3_soa_wrapper(ferm_phi_acc,fermionname_deo);
        print_vec3_soa(ferm_phi_acc,myfermionname_deo);
    }

    if(0 == devinfo.myrank)
        printf("Multiplication by M^\\dagM+m^2, %d times...\n", test_settings.deoDoeIterations);
    gettimeofday(&t4,NULL);
    for(r=0; r<test_settings.deoDoeIterations; r++)
        fermion_matrix_multiplication(conf_acc[0], ferm_phi_acc, 
                ferm_chi_acc, kloc_s, &fermions_parameters[0]) ;
    gettimeofday(&t5,NULL);

    if(test_settings.saveResults){
        if(0 == devinfo.myrank)
            printf("Writing file %s.\n", fermionname_fulldirac);
#pragma acc update host(ferm_phi_acc[0:1])
        print_vec3_soa_wrapper(ferm_phi_acc,fermionname_fulldirac);
        print_vec3_soa(ferm_phi_acc,myfermionname_fulldirac);
    }

    double dt_doe = (double)(t1.tv_sec - t0.tv_sec) + 
        ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);
    double dt_deo = (double)(t3.tv_sec - t2.tv_sec) + 
        ((double)(t3.tv_usec - t2.tv_usec)/1.0e6);
    double dt_dirac = (double)(t5.tv_sec - t4.tv_sec) + 
        ((double)(t5.tv_usec - t4.tv_usec)/1.0e6);

    if(0 == devinfo.myrank){
        printf("Time for 1 application of Doe           : %e\n",
                dt_doe/test_settings.deoDoeIterations);
        printf("Time for 1 application of Deo           : %e\n", 
                dt_deo/test_settings.deoDoeIterations);
        printf("Time for 1 application of Dirac Operator: %e\n", 
                dt_dirac/test_settings.deoDoeIterations);
    }
    // conversion to single precision
    convert_double_to_float_su3_soa(conf_acc[0],conf_acc_f[0]);
    convert_double_to_float_vec3_soa(ferm_chi_acc,ferm_chi_acc_f);

    // single precision
    if(0 == devinfo.myrank){
        printf("####################\n");
        printf("# SINGLE PRECISION #\n");
        printf("####################\n");
    }
    struct timeval t0_f,t1_f,t2_f,t3_f,t4_f,t5_f;
    //int r;
    if(0 == devinfo.myrank)
        printf("Multiplication by Doe, %d times...\n", test_settings.deoDoeIterations);
    gettimeofday(&t0_f,NULL);
    for(r=0; r<test_settings.deoDoeIterations; r++)
        acc_Doe_f(conf_acc_f[0], ferm_phi_acc_f, ferm_chi_acc_f, fermions_parameters[0].phases_f);
    gettimeofday(&t1_f,NULL);


    if(test_settings.saveResults){
        if(0 == devinfo.myrank)
            printf("Writing file %s.\n", fermionname_doe_f);
#pragma acc update host(ferm_phi_acc_f[0:1])
        print_vec3_soa_wrapper_f(ferm_phi_acc_f,fermionname_doe_f);
        print_vec3_soa_f(ferm_phi_acc_f,myfermionname_doe_f);
    }

    if(0 == devinfo.myrank)
        printf("Multiplication by Deo, %d times...\n", test_settings.deoDoeIterations);
    gettimeofday(&t2_f,NULL);
    for(r=0; r<test_settings.deoDoeIterations; r++)
        acc_Deo_f(conf_acc_f[0], ferm_phi_acc_f, ferm_chi_acc_f, fermions_parameters[0].phases_f) ;
    gettimeofday(&t3_f,NULL);

    if(test_settings.saveResults){
        if(0 == devinfo.myrank)
            printf("Writing file %s.\n", fermionname_deo_f);
#pragma acc update host(ferm_phi_acc_f[0:1])
        print_vec3_soa_wrapper_f(ferm_phi_acc_f,fermionname_deo_f);
        print_vec3_soa_f(ferm_phi_acc_f,myfermionname_deo_f);
    }

    if(0 == devinfo.myrank)
        printf("Multiplication by M^\\dagM+m^2, %d times...\n", test_settings.deoDoeIterations);
    gettimeofday(&t4_f,NULL);
    for(r=0; r<test_settings.deoDoeIterations; r++)
        fermion_matrix_multiplication_f(conf_acc_f[0], ferm_phi_acc_f, 
                ferm_chi_acc_f, kloc_s_f, &fermions_parameters[0]) ;
    gettimeofday(&t5_f,NULL);

    if(test_settings.saveResults){
        if(0 == devinfo.myrank)
            printf("Writing file %s.\n", fermionname_fulldirac_f);
#pragma acc update host(ferm_phi_acc_f[0:1])
        print_vec3_soa_wrapper_f(ferm_phi_acc_f,fermionname_fulldirac_f);
        print_vec3_soa_f(ferm_phi_acc_f,myfermionname_fulldirac_f);
    }

    double dt_doe_f = (double)(t1_f.tv_sec - t0_f.tv_sec) + 
        ((double)(t1_f.tv_usec - t0_f.tv_usec)/1.0e6);
    double dt_deo_f = (double)(t3_f.tv_sec - t2_f.tv_sec) + 
        ((double)(t3_f.tv_usec - t2_f.tv_usec)/1.0e6);
    double dt_dirac_f = (double)(t5_f.tv_sec - t4_f.tv_sec) + 
        ((double)(t5_f.tv_usec - t4_f.tv_usec)/1.0e6);

    if(0 == devinfo.myrank){
        printf("Time for 1 application of Doe           : %e\n",
                dt_doe_f/test_settings.deoDoeIterations);
        printf("Time for 1 application of Deo           : %e\n", 
                dt_deo_f/test_settings.deoDoeIterations);
        printf("Time for 1 application of Dirac Operator: %e\n", 
                dt_dirac_f/test_settings.deoDoeIterations);
    }

    printf("MPI%02d: Test completed.\n",devinfo.myrank);
    
    mem_free_extended_f();
    mem_free_core_f();
    mem_free_extended();
    mem_free_core();

#ifndef __GNUC__
    shutdown_acc_device(my_device_type);
#endif
    printf("MPI%02d: Before MPI_Finalize()...\n",devinfo.myrank);
#ifdef MULTIDEVICE
    MPI_Finalize();
#endif 


    return 0; 


}
