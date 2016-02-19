#include "../Include/init.h"
#include "../OpenAcc/io.h"
#include "../Include/fermion_parameters.h"
#include "../OpenAcc/alloc_vars.h"
#include "../Include/markowchain.h"

int conf_id_iter;
int verbosity_lv;

int main(int argc, char* argv[]){

    printf("WELCOME! \n");
    // READ input file.
    set_global_vars_and_fermions_from_input_file(argv[1]);
    //
    verbosity_lv = mkwch_pars.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters)) exit(1);


    mem_alloc();

    printf("Reading ildg conf\n");
    if(read_su3_soa_ildg_binary(conf_acc,"confildgtest",&conf_id_iter))
        exit(1);
    printf("Writing ASCII conf\n");
    print_su3_soa_ASCII(conf_acc,mkwch_pars.save_conf_name,conf_id_iter); // READS ALSO THE conf_id_iter

    printf("Writing ildg conf\n");
    print_su3_soa_ildg_binary(conf_acc,"morite_tutti",conf_id_iter);

    mem_free();

    return 0;

}





