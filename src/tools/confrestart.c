
#include "../OpenAcc/io.h"
#include "./struct_c_def.h"
#include <stdlib.h>

int main(int argc, char ** argv){

    global_su3_soa * conf = (global_su3_soa * ) malloc(8*sizeof(global_su3_soa));
    int conf_id_iter;
    int new_conf_id_iter = 0;


    printf("Reading conf %s...\n", argv[1]);
    read_su3_soa_ildg_binary(conf,argv[1],&conf_id_iter);
           
    printf("conf_id_iter is %d.\n", conf_id_iter);

    printf("Writing conf in file %s with conf_id_iter = %d", argv[2], new_conf_id_iter);
    print_su3_soa_ildg_binary(conf,argv[2],new_conf_id_iter);

    free(conf);

}
