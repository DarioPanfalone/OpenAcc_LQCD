#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "./rationalapprox.h"
#include"./Remez/bigfloat.h"
#include"./Remez/alg_remez.h"


//using namespace std;

const int gmp_remez_precision=100; // The precision that gmp uses

//Rational approximations for Metropolis
const double lambda_min=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double residue=1.0e-8;    // stopping residual for CG

int main(int argc, char **argv){
    if(argc!=5  && argc != 6){
        printf("Use as arguments:    tolerance exponent_num exponent_den lambda_min (optional: min approx order)\n");
        return 0;
    }

    RationalApprox* approx = new RationalApprox;

    double goal_error = atof(argv[1]);
    approx->error = 1;

    if(argc == 6 ) approx->approx_order = atoi(argv[5]);
    else approx->approx_order = 7;

    while(approx->error > goal_error){
        approx->exponent_num = atoi(argv[2]);
        approx->exponent_den = atoi(argv[3]);
        approx->gmp_remez_precision = gmp_remez_precision; // The precision that gmp uses
        approx->lambda_min = atof(argv[4]);
        approx->lambda_max = 1.0;              // The bounds of the approximation

        // The error from the approximation (the relative error is minimised
        // - if another error minimisation is requried, then line 398 in
        // alg_remez.C is where to change it)

        // The partial fraction expansion takes the form 
        // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])


       fflush(stdout);
        // Instantiate the Remez class
        AlgRemez remez1(approx->lambda_min,approx->lambda_max,approx->gmp_remez_precision);
        // Generate the required approximation
        approx->error = remez1.generateApprox(approx->approx_order,approx->approx_order,abs(approx->exponent_num),abs(approx->exponent_den));
         printf("\nf(x) = (x)^(%d/%d), approximation on [%e, %e], order %d, error = %e, ",
                approx->exponent_num, approx->exponent_den,approx->lambda_min,
                approx->lambda_max, approx->approx_order, approx->error);

        // Find the partial fraction expansion of the approximation 
        // to the function x^{y/z} (this only works currently for 
        // the special case that n = d)

        if(approx->error > goal_error) {
            printf(" not precise enough, increasing order.\n");
            approx->approx_order++;
        }
        else {
            printf("OK!\n");
            if(approx->exponent_num*approx->exponent_den > 0)
                remez1.getPFE(approx->RA_a,approx->RA_b,&(approx->RA_a0));
            else
                remez1.getIPFE(approx->RA_a,approx->RA_b,&(approx->RA_a0));
        }
    }
    char * nomefile = rational_approx_filename(goal_error, approx->exponent_num, approx->exponent_den, approx->lambda_min);

    rationalapprox_save(nomefile, approx);

}
