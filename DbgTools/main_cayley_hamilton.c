

#include "../OpenAcc/cayley_hamilton.c"
//#include "../Rand/random.c"
#include <stdlib.h>
#include <time.h>  

double crand_casuale(){
  return  1.0-2.0*(double)rand()/((double)(RAND_MAX));
}

int main(){

    single_thmat P;
    single_tamat PA;
    single_su3 Q,QA;

//    initrand(2);
    srand(time(NULL));

    P.rc00 =crand_casuale();
    P.rc11 =crand_casuale();
    P.c01  =crand_casuale() + I*crand_casuale();
    P.c02  =crand_casuale() + I*crand_casuale();
    P.c12  =crand_casuale() + I*crand_casuale();
    printf("P.rc00  = %f \n", P.rc00);
    printf("P.rc11  = %f \n", P.rc11);
    printf("P.c01   = %f \n", P.c01 );
    printf("P.c02   = %f \n", P.c02 );
    printf("P.c12   = %f \n", P.c12 );

    PA.rc00 = - P.rc00 ;// (double)rand()/((double)(RAND_MAX));
    PA.rc11 = - P.rc11 ;// (double)rand()/((double)(RAND_MAX));
    PA.c01  = (-1.0*I) * P.c01  ;// (double)rand()/((double)(RAND_MAX));
    PA.c02  = (-1.0*I) * P.c02  ;// (double)rand()/((double)(RAND_MAX));
    PA.c12  = (-1.0*I) * P.c12  ;// (double)rand()/((double)(RAND_MAX));

    d_complex determ = detSu3(&Q);
    printf("Det(Q) = %.18lf %.18lf \n",creal(determ), cimag(determ) );

    CH_exponential(&Q,&P);
    print_su3_stdout(&Q);

    determ = detSu3(&Q);
    printf("Det(Q) = %.18lf %.18lf \n",creal(determ), cimag(determ) );

    CH_exponential_antihermitian(&QA,&PA);
    print_su3_stdout(&QA);
    determ = detSu3(&QA);
    printf("Det(QA)= %.18lf %.18lf \n",creal(determ), cimag(determ) );

    return 0;

}
