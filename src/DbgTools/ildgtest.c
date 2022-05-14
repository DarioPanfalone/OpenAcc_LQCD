#include "../Include/setting_file_parser.h"
#include "../OpenAcc/io.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"

int conf_id_iter;
int verbosity_lv;


su3_soa  * conf_hasenbusch[0]; // the gauge configuration.

int main(int argc, char* argv[]){

    printf("WELCOME! \n");
    // READ input file.
    set_global_vars_and_fermions_from_input_file(argv[1]);
    //
    verbosity_lv = mkwch_pars.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters)) exit(1);

#define ALLOCCHECK(control_int,var)  if(control_int != 0 ) \
    printf("\tError in  allocation of %s . \n", #var);\
    else if(verbosity_lv > 2) printf("\tAllocation of %s : OK , %p\n", #var, var );\

#define ALIGN 128
    int allocation_check =  posix_memalign((void **)&conf_hasenbusch[0], ALIGN, 8*sizeof(su3_soa));
    ALLOCCHECK(allocation_check, conf_hasenbusch[0]);



    printf("Reading ildg conf\n");
    if(read_su3_soa_ildg_binary(conf_hasenbusch[0],"confildgtest",&conf_id_iter))
        exit(1);
    printf("Writing ASCII conf, %d\n", conf_id_iter);
    // READS ALSO THE conf_id_iter
    print_su3_soa_ASCII(conf_hasenbusch[0],mkwch_pars.save_conf_name,conf_id_iter); 

    /*
    printf("Writing ildg conf\n");
    print_su3_soa_ildg_binary(conf_hasenbusch[0],"conf.out.ildg",conf_id_iter);
    */



    free(conf_hasenbusch[0]);


    return 0;

}
