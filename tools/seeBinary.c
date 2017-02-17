#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]){

    char *fileToRead = argv[1]; 
    unsigned int charsToRead = atoi(argv[2]);

    FILE * file = fopen(fileToRead,"r");

    unsigned char * array;
    array = (unsigned char*)    malloc(charsToRead);


    fread(array,1,charsToRead,file);


    int i;
    for(i=0;i<charsToRead;++i){

        if(i%32 == 0 ) printf("\n %03d",i/32);
        if((i%32)% 4 == 0 ) printf(" ");
        printf("%02x",array[i]);
    }

    printf("\n");
    free(array);

    return 0;
}
