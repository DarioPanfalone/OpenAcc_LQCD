
#include "../OpenAcc/cayley_hamilton.h"
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/single_types.h"
#include "../Rand/random.h"

#include <stdio.h>
#include <sys/time.h>

int main(){

    d_complex c01,c02,c12;
    double rc00,rc11;

    single_su3 RES1,RES2,MOM,AUX;
    single_tamat QA;

    initrand(2);

    rc00 =casuale();
    rc11 =casuale();
    c01  =casuale() + I*casuale();
    c02  =casuale() + I*casuale();
    c12  =casuale() + I*casuale();

    printf("rc00  =  %e \n", rc00);
    printf("rc11  =  %e \n", rc11);
    printf("c01   = (%e , %e )\n", creal(c01), cimag(c01));
    printf("c02   = (%e , %e )\n", creal(c02), cimag(c02));
    printf("c12   = (%e , %e )\n", creal(c12), cimag(c12));

    // Taylor test

    MOM.comp[0][0] = rc00 * ( 1.0I);
    MOM.comp[0][1] = c01 * ( 1.0I);
    MOM.comp[0][2] = c02 * ( 1.0I);
    MOM.comp[1][0] = - conj(MOM.comp[0][1]);
    MOM.comp[1][1] = rc11 * ( 1.0I);
    MOM.comp[1][2] = c12 * ( 1.0I);
    MOM.comp[2][0] = - conj(MOM.comp[0][2]);
    MOM.comp[2][1] = - conj(MOM.comp[1][2]);
    MOM.comp[2][2] = - MOM.comp[0][0] - MOM.comp[1][1];

    // CH test
    QA.ic00 = -rc00;
    QA.ic11 = -rc11;
    QA.c01  = -c01 *1.0I;
    QA.c02  = -c02 *1.0I;
    QA.c12  = -c12 *1.0I;

    struct timeval t0,t1,t2,t3;
    int iter;
    int maxiter = 1000000;

    gettimeofday(&t0,NULL);
    for(iter =0; iter < maxiter ; iter++){
        MOM.comp[0][0] += (iter%2-1)*0.05;
        taylor_exponential_su3(&RES1,&MOM,&AUX);
    }
    gettimeofday(&t1,NULL);

    double dt = (double)(t1.tv_sec - t0.tv_sec) + 
        ((double)(t1.tv_usec - t0.tv_usec)/1.0e6);

    printf("Taylor, time %e:\n", dt/maxiter);
    print_su3_stdout(&RES1);

    printf("determinant-1: %e  +  %e I\n",creal(detSu3(&RES1))-1, cimag(detSu3(&RES1)) );
    
    maxiter = 1000000;

    gettimeofday(&t2,NULL);
    for(iter =0; iter < maxiter ; iter++){
        QA.ic00 += (iter%2-1)*0.05;
        CH_exponential_antihermitian_nissalike(&RES2,&QA);
    }
    gettimeofday(&t3,NULL);
    dt = (double)(t3.tv_sec - t2.tv_sec) + 
        ((double)(t3.tv_usec -t2.tv_usec)/1.0e6);


    printf("Cayley Hamilton, time %e:\n", dt/maxiter);
    print_su3_stdout(&RES2);
    rebuild3row(&RES2);

    printf("determinant-1: %e  +  %e I\n",creal(detSu3(&RES2))-1, cimag(detSu3(&RES2)) );

    int i,j;
    double diff=0;
    d_complex Z;
    for(i=0; i< 2 ; i++)
        for(j=0; j< 3 ; j++){
            Z = RES1.comp[i][j] - RES2.comp[i][j];
            diff += creal(Z)*creal(Z) + cimag(Z)*cimag(Z);
        }

    printf("Difference: %e \n", diff);

    return 0;



}


