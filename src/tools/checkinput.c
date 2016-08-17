

#include "../Include/setting_file_parser.h"
#include "../Mpi/multidev.h"

#include "stdio.h"

int verbosity_lv;

int main(int argc, char* argv[]){


    devinfo.myrank = 0;
    printf("Just reading input file and calculating hash...\n") ;
    set_global_vars_and_fermions_from_input_file(argv[1]);

    return 0 ; 

}
