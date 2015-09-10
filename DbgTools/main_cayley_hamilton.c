

#include "../OpenAcc/cayley_hamilton.c"
//#include "../Rand/random.c"
#include <stdlib.h>
#include <time.h>  


int main(){

    single_thmat P;
    single_su3 Q;

//    initrand(2);
    srand(time(NULL));

    P.rc00 =0 ; (double)rand()/((double)(RAND_MAX));
    P.rc11 =0 ; (double)rand()/((double)(RAND_MAX));
    P.c01 = 0 ; (double)rand()/((double)(RAND_MAX));
    P.c02 = 0 ; (double)rand()/((double)(RAND_MAX));
    P.c12 = 1 ; (double)rand()/((double)(RAND_MAX));
    printf("P.rc00  = %f \n", P.rc00);
    printf("P.rc11  = %f \n", P.rc11);
    printf("P.c01   = %f \n", P.c01 );
    printf("P.c02   = %f \n", P.c02 );
    printf("P.c12   = %f \n", P.c12 );




    printf("Det(Q)=%f\n",detSu3(&Q) );

    CH_exponential(&Q,&P);

    print_su3_stdout(&Q);

    printf("Det(Q)=%f\n",detSu3(&Q) );

    return 0;

}
