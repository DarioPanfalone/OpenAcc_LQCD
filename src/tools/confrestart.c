
#include "../OpenAcc/io.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/geometry.h"
#include "../Mpi/geometry_multidev.h"
#include <stdlib.h>
#include <stdio.h>

#include "../Include/stringify.h"

int verbosity_lv;

int main(int argc, char ** argv){

    global_su3_soa * conf = (global_su3_soa * ) malloc(8*sizeof(global_su3_soa));
    int conf_id_iter;
    int new_conf_id_iter = 0;

    geom_par.gnx = GL_N0;
    geom_par.gny = GL_N1;
    geom_par.gnz = GL_N2;
    geom_par.gnt = GL_N3;


    printf("Reading conf %s...\n", argv[1]);
    read_su3_soa_ildg_binary(conf,argv[1],&conf_id_iter);
           
    printf("conf_id_iter is %d.\n", conf_id_iter);

    printf("Writing conf in file %s with conf_id_iter = %d", argv[2], new_conf_id_iter);
    print_su3_soa_ildg_binary(conf,argv[2],new_conf_id_iter);

    free(conf);
    printf("(commit: %s)\n", xstr(COMMIT_HASH) );

}
