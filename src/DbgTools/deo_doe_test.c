#include <stdio.h>
#include <stdlib.h>

#ifndef __GNUC__
#include "openacc.h"
#endif

#include "../Include/common_defines.h"
#include "../Include/init.h"
#include "../Include/fermion_parameters.h"
#include "../OpenAcc/fermion_matrix.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/io.h"
#include "../OpenAcc/random_assignement.h"
#include "../Rand/random.h"
#include "./dbgtools.h"
#include "../OpenAcc/action.h"
#include "../Include/markowchain.h"

//double casuale(void);


int id_iter;
int verbosity_lv;

int main(int argc, char* argv[]){

    const char confname[] = "test_conf";
    const char fermionname[] = "test_fermion";
    const char fermionname_doe[] = "test_fermion_result_doe";
    const char fermionname_deo[] = "test_fermion_result_deo";

    act_params.stout_steps = 0;
    NPS_tot = 1;


    initrand(111);
    fflush(stdout);
    printf("DEODOE test\n");

    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    printf("WELCOME! \n");
    // READ input file.
    set_global_vars_and_fermions_from_input_file(argv[1]);
    //

    initrand((unsigned int) mkwch_pars.seed);
    verbosity_lv = mkwch_pars.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters))
        printf("Ignoring issues in init_ferm_params,\
this is a deo-doe test.\n");

    mem_alloc();
    printf("Allocazione della memoria : OK \n");
    compute_nnp_and_nnm_openacc();
    printf("nn computation : OK \n");
    init_all_u1_phases(backfield_parameters,fermions_parameters);
    printf("u1_backfield initialization : OK \n");
 

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

    // init conf
    //generate_Conf_cold(conf_acc,0.05);
    //print_su3_soa(conf_acc, "conf_acc");
    //printf("Cold Gauge Conf Generated : OK \n");
    if(!read_su3_soa_ASCII(conf_acc,confname, &id_iter)){
        printf("Cold Gauge Conf READ : OK \n");
    }else{
        id_iter = 0;
        printf("Cold Gauge Conf READ : FAIL, generating\n");
        generate_Conf_cold(conf_acc,mkwch_pars.eps_gen); 
        printf("Cold Gauge Conf GENERATED: OK \n");
        printf("Writing file %s.\n", confname);
        print_su3_soa_ASCII(conf_acc,confname, id_iter);
    }


    // init fermion
    //generate_vec3_soa_gauss(ferm_chi_acc);
    //print_vec3_soa(ferm_chi_acc,"ferm_chi_acc");
    if(!read_vec3_soa(ferm_chi_acc,fermionname )){
        printf("Fermion READ : OK \n");
    }else{
        printf("Fermion READ : FAIL, generating\n");
        generate_vec3_soa_gauss(ferm_chi_acc);
        printf("Fermion GENERATED : OK\n");
        printf("Writing file %s.\n", fermionname);
        print_vec3_soa(ferm_chi_acc,fermionname);
    }


#pragma acc data   copy(conf_acc[0:8]) copy(ferm_chi_acc) copyout(ferm_phi_acc) 
    {
        printf("Multiplication by Doe...");
        acc_Doe(conf_acc, ferm_phi_acc, ferm_chi_acc, fermions_parameters[0].phases) ;
        printf("Writing file %s.\n", fermionname_doe);
        print_vec3_soa(ferm_phi_acc,fermionname_doe);

        printf("Multiplication by Deo...");
        acc_Deo(conf_acc, ferm_phi_acc, ferm_chi_acc, fermions_parameters[0].phases) ;
        printf("Writing file %s.\n", fermionname_deo);
        print_vec3_soa(ferm_phi_acc,fermionname_deo);


    }

#ifndef __GNUC__
    SHUTDOWN_ACC_DEVICE(my_device_type);
#endif

    mem_free();

    return 0; 


}
