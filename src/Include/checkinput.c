

#include "./input_file_parser.h"

#include "stdio.h"

int verbosity_lv;

int main(int argc, char* argv[]){

    printf("Just reading input file and calculating hash...\n") ;
    set_global_vars_and_fermions_from_input_file(argv[1]);

    return 0 ; 

}
